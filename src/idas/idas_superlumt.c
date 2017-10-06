/*
 * -----------------------------------------------------------------
 * $Revision: 4938 $
 * $Date: 2016-09-21 14:33:08 -0700 (Wed, 21 Sep 2016) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the IDASSUPERLUMT linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "sundials/sundials_math.h"
#include "sundials/sundials_superlumt_impl.h"

#include "idas_impl.h"
#include "idas_sparse_impl.h"
#include "idas/idas_superlumt.h"

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* IDASSSUPERLUMT linit, lsetup, lsolve, and lfree routines */
 
static int IDASuperLUMTInit(IDAMem IDA_mem);

static int IDASuperLUMTSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
			     N_Vector rrp, N_Vector tmp1,
			     N_Vector tmp2, N_Vector tmp3);

static int IDASuperLUMTSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
			     N_Vector ycur, N_Vector ypcur, N_Vector rrcur);

static int IDASuperLUMTFree(IDAMem IDA_mem);

/* IDASUPERLUMT lfreeB function */

static int IDASuperLUMTFreeB(IDABMem IDAB_mem);


/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/*
 * -----------------------------------------------------------------
 * IDASuperLUMT
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the IDA / SuperLUMT linear solver module.  
 * IDASuperLUMT first calls the existing lfree routine if this is not NULL.
 * Then it sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDASuperLUMTInit, IDASuperLUMTSetup,
 * IDASuperLUMTSolve, NULL, and IDASuperLUMTFree, respectively.
 * It allocates memory for a structure of type IDAsluMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem) to TRUE, sets the d_jdata field
 * in the IDAsluMemRec structure to be the input parameter jdata,
 * and sets the d_jac field to be:
 *   (1) the input parameter djac, if djac != NULL, or                
 *   (2) throws an error, if djac == NULL.                             
 * Finally, it allocates memory for SuperLUMT.
 * The return value is IDASLS_SUCCESS = 0, IDASLS_LMEM_FAIL = -1,
 * or IDASLS_ILL_INPUT = -2.
 *
 * NOTE: The SuperLUMT linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, IDASuperLUMT will first 
 *       test for a compatible N_Vector internal representation
 *       by checking that the function N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */

int IDASuperLUMT(void *ida_mem, int num_threads, int n, int nnz)
{
  IDAMem IDA_mem;
  IDASlsMem idasls_mem;
  SLUMTData slumt_data;
  int *perm_c, *perm_r;
  int nrhs, panel_size, relax;
  double *bd;
  SuperMatrix *B;

  /* Return immediately if ida_mem is NULL. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDASuperLUMT", 
		    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the Direct solver */
  if (IDA_mem->ida_tempv1->ops->nvgetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", "IDASuperLUMT", 
		    MSGSP_BAD_NVECTOR);
    return(IDASLS_ILL_INPUT);
  }

  if (IDA_mem->ida_lfree != NULL) IDA_mem->ida_lfree(IDA_mem);

  /* Set five main function fields in IDA_mem. */
  IDA_mem->ida_linit  = IDASuperLUMTInit;
  IDA_mem->ida_lsetup = IDASuperLUMTSetup;
  IDA_mem->ida_lsolve = IDASuperLUMTSolve;
  IDA_mem->ida_lperf  = NULL;
  IDA_mem->ida_lfree  = IDASuperLUMTFree;

  /* Get memory for IDASlsMemRec. */
  idasls_mem = NULL;
  idasls_mem = (IDASlsMem) malloc(sizeof(struct IDASlsMemRec));
  if (idasls_mem == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDASuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }

  /* Get memory for SLUMT_data. */
  slumt_data = NULL;
  slumt_data = (SLUMTData)malloc(sizeof(struct SLUMTDataRec));
  if (slumt_data == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDASuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }

  IDA_mem->ida_setupNonNull = TRUE;

  /* Set default Jacobian routine and Jacobian data */
  idasls_mem->s_jaceval = NULL;
  idasls_mem->s_jacdata = IDA_mem->ida_user_data;

  /* Allocate memory for the sparse Jacobian */
  idasls_mem->s_JacMat = NULL;
  idasls_mem->s_JacMat = SparseNewMat(n, n, nnz, CSC_MAT);
  if (idasls_mem->s_JacMat == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDASuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }

  /* Set up memory for the permutations */
  perm_r = (int *)malloc(n*sizeof(int));
  if (perm_r == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDASuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }
  perm_c = (int *)malloc(n*sizeof(int));
  if (perm_c == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDASuperLUMT", 
		    MSGSP_MEM_FAIL);
    free(perm_r);
    return(IDASLS_MEM_FAIL);
  }
  slumt_data->perm_r = perm_r;
  slumt_data->perm_c = perm_c;

  idasls_mem->s_last_flag = IDASLS_SUCCESS;

  /* Set default parameters for SuperLU */
  slumt_data->num_threads = num_threads;
  slumt_data->diag_pivot_thresh = 1.0;

  /* Allocate structures for SuperLU */
  slumt_data->Gstat = (Gstat_t *)malloc(sizeof(Gstat_t));
  slumt_data->s_A = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  slumt_data->s_AC = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  slumt_data->s_L = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  slumt_data->s_U = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  slumt_data->s_A->Store  = NULL;
  slumt_data->s_AC->Store = NULL;
  slumt_data->s_L->Store  = NULL;
  slumt_data->s_U->Store  = NULL;
  slumt_data->superlumt_options = (superlumt_options_t *)malloc(sizeof(superlumt_options_t));

  dCreate_CompCol_Matrix(slumt_data->s_A, idasls_mem->s_JacMat->M, idasls_mem->s_JacMat->N, 
			 idasls_mem->s_JacMat->NNZ, idasls_mem->s_JacMat->data, 
			 idasls_mem->s_JacMat->indexvals, idasls_mem->s_JacMat->indexptrs, 
			 SLU_NC, SLU_D, SLU_GE);

  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  StatAlloc(idasls_mem->s_JacMat->N, num_threads, panel_size, relax, slumt_data->Gstat);
  
  /* Create RHS matrix */
  nrhs = 1;
  bd = NULL;
  B = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  B->Store = NULL;
  dCreate_Dense_Matrix(B, n, nrhs, bd, n, 
		       SLU_DN, SLU_D, SLU_GE);
  slumt_data->s_B = B;

  /* Set ordering to COLAMD as the idas default use.
     Users can set a different value with IDASuperLUMTSetOrdering,
     and the user-set value is loaded before any call to factorize the
     matrix in IDASuperLUMTSetup.  */
  slumt_data->s_ordering = 3;

  /* Attach linear solver memory to the integrator memory */
  idasls_mem->s_solver_data = (void *) slumt_data;
  IDA_mem->ida_lmem = idasls_mem;

  return(IDASLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDASuperLUMT interface functions
 * -----------------------------------------------------------------
 */

/*
  This routine does remaining initializations specific to the IDASuperLUMT
  linear solver module.  
  It returns 0 if successful.
*/

static int IDASuperLUMTInit(IDAMem IDA_mem)
{
  int num_threads, n;
  IDASlsMem idasls_mem;
  SLUMTData slumt_data;

  idasls_mem = (IDASlsMem)IDA_mem->ida_lmem;
  slumt_data = (SLUMTData) idasls_mem->s_solver_data;

  idasls_mem->s_nje = 0;
  idasls_mem->s_first_factorize = 1;

  /* ------------------------------------------------------------
     Allocate storage and initialize statistics variables. 
     ------------------------------------------------------------*/
  n = idasls_mem->s_JacMat->N;
  num_threads = slumt_data->num_threads;

  StatInit(n, num_threads, slumt_data->Gstat);

  idasls_mem->s_last_flag = 0;
  return(0);
}

/*
  This routine does the setup operations for the IDASuperLUMT linear 
  solver module.  It calls the Jacobian evaluation routine,
  updates counters, and calls the LU factorization routine.
  The return value is either
     IDASLS_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
*/

static int IDASuperLUMTSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
			     N_Vector rrp, N_Vector tmp1, N_Vector tmp2,
			     N_Vector tmp3)
{
  int retval, info;
  int nprocs, panel_size, relax, permc_spec, lwork;
  int *perm_r, *perm_c;
  realtype tn, cj;
  double diag_pivot_thresh, drop_tol;
  fact_t fact;
  trans_t trans;
  yes_no_t refact, usepr;
  IDASlsMem idasls_mem;
  IDASlsSparseJacFn jaceval;
  SuperMatrix *A, *AC, *L, *U;
  Gstat_t *Gstat;
  superlumt_options_t *superlumt_options;
  SLUMTData slumt_data;
  SlsMat JacMat;
  void *jacdata;
  void *work;
  
  idasls_mem = (IDASlsMem) (IDA_mem->ida_lmem);
  tn = IDA_mem->ida_tn; 
  cj = IDA_mem->ida_cj;

  slumt_data = (SLUMTData) idasls_mem->s_solver_data;
  jaceval = idasls_mem->s_jaceval;
  jacdata = idasls_mem->s_jacdata;
  JacMat = idasls_mem->s_JacMat;

  superlumt_options = slumt_data->superlumt_options;
  A = slumt_data->s_A;
  AC = slumt_data->s_AC;
  L = slumt_data->s_L;
  U = slumt_data->s_U;
  Gstat = slumt_data->Gstat;
  perm_r = slumt_data->perm_r;
  perm_c = slumt_data->perm_c;
  nprocs = slumt_data->num_threads;
  diag_pivot_thresh = slumt_data->diag_pivot_thresh;

  /* Set option values for SuperLU_MT */
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  fact = EQUILIBRATE;
  trans = NOTRANS;
  usepr = NO;
  drop_tol = 0.0;
  lwork = 0;
  work = NULL;

  /* Check that Jacobian eval routine is set */
  if (jaceval == NULL) {
    IDAProcessError(IDA_mem, IDASLS_JAC_NOSET, "IDASSLS", "IDASuperLUMTSetup", 
		    MSGSP_JAC_NOSET);
    free(idasls_mem); idasls_mem = NULL;
    return(IDASLS_JAC_NOSET);
  }

  /* Increment nje counter and call Jacobian eval routine. */
  idasls_mem->s_nje++;
  retval = jaceval(tn, cj, yyp, ypp, rrp, JacMat, jacdata, 
		   tmp1, tmp2, tmp3);

  if (retval < 0) {
    IDAProcessError(IDA_mem, IDASLS_JACFUNC_UNRECVR, "IDASSLS", 
		    "IDASuperLUMTSetup", MSGSP_JACFUNC_FAILED);
    idasls_mem->s_last_flag = IDASLS_JACFUNC_UNRECVR;
    return(IDASLS_JACFUNC_UNRECVR);
  }
  if (retval > 0) {
    idasls_mem->s_last_flag = IDASLS_JACFUNC_RECVR;
    return(+1);
  }

  if (idasls_mem->s_first_factorize) {
    /* ------------------------------------------------------------
       Get column permutation vector perm_c[], according to permc_spec:
       permc_spec = 3: approximate minimum degree for unsymmetric matrices
       ------------------------------------------------------------*/ 
    permc_spec = 3;
    get_perm_c(permc_spec, A, perm_c);
 
    refact= NO;
    idasls_mem->s_first_factorize = 0;
  }
  else {
    /* ------------------------------------------------------------
       Re-initialize statistics variables 
       ------------------------------------------------------------*/
    StatInit(JacMat->N, nprocs, Gstat);
    Destroy_CompCol_Permuted(AC);
    refact= YES;
  }

  /* ------------------------------------------------------------
     Initialize the option structure superlumt_options using the
     user-input parameters;  Subsequent calls will re-initialize
     options.
     Apply perm_c to the columns of original A to form AC.
     ------------------------------------------------------------*/
  pdgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
	       diag_pivot_thresh, usepr, drop_tol, perm_c, perm_r,
	       work, lwork, A, AC, superlumt_options, Gstat);
  /* ------------------------------------------------------------
     Compute the LU factorization of A.
     The following routine will create nprocs threads.
     ------------------------------------------------------------*/
  pdgstrf(superlumt_options, AC, perm_r, L, U, Gstat, &info);
    
  if (info != 0) {
    idasls_mem->s_last_flag = info;
    return(+1);
  }
  idasls_mem->s_last_flag = IDASLS_SUCCESS;

  return(0);
}

/*
  This routine handles the solve operation for the IDASuperLUMT linear
  solver module.  It calls the SuperLUMT solve routine, scales the
  solution vector according to cjratio, then returns IDASLU_SUCCESS = 0.
*/

static int IDASuperLUMTSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
			     N_Vector ycur, N_Vector ypcur, N_Vector rrcur)
{
  int info, trans;
  int *perm_r, *perm_c;
  double cjratio;
  IDASlsMem idasls_mem;
  SuperMatrix *L, *U, *B;
  Gstat_t *Gstat;
  DNformat *Bstore;
  SLUMTData slumt_data;
  realtype *bd;
  
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;
  cjratio = IDA_mem->ida_cjratio;

  slumt_data = (SLUMTData) idasls_mem->s_solver_data;
  L = slumt_data->s_L;
  U = slumt_data->s_U;
  perm_r = slumt_data->perm_r;
  perm_c = slumt_data->perm_c;
  Gstat = slumt_data->Gstat;
  B = slumt_data->s_B;
   
  bd = N_VGetArrayPointer(b);
  Bstore = (DNformat *) (B->Store);
  Bstore->nzval = bd;

  /* Call SuperLUMT to solve the linear system using L and U */
  trans = NOTRANS;
  dgstrs(trans, L, U, perm_r, perm_c, B, Gstat, &info);

  /* Scale the correction to account for change in cj. */
  if (cjratio != ONE) N_VScale(TWO/(ONE + cjratio), b, b);

  Bstore->nzval = NULL;

  idasls_mem->s_last_flag = IDASLS_SUCCESS;
  return(IDASLS_SUCCESS);
}

/*
  This routine frees memory specific to the IDASuperLUMT linear solver.
*/

static int IDASuperLUMTFree(IDAMem IDA_mem)
{
  IDASlsMem idasls_mem;
  SLUMTData slumt_data;
  
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;

  slumt_data = (SLUMTData) idasls_mem->s_solver_data;

  pxgstrf_finalize(slumt_data->superlumt_options, slumt_data->s_AC);

  free(slumt_data->perm_r);
  free(slumt_data->perm_c);
  free(slumt_data->superlumt_options);
  Destroy_SuperNode_SCP( (slumt_data->s_L) );
  Destroy_CompCol_NCP( (slumt_data->s_U) );
  StatFree( (slumt_data->Gstat) );
  free(slumt_data->Gstat);
  
  Destroy_SuperMatrix_Store(slumt_data->s_B);
  SUPERLU_FREE(slumt_data->s_A->Store);
  if (idasls_mem->s_JacMat) {
    SparseDestroyMat(idasls_mem->s_JacMat);
    idasls_mem->s_JacMat = NULL;
  }

  free(slumt_data->s_B);
  free(slumt_data->s_A);
  free(slumt_data->s_AC);
  free(slumt_data->s_L);
  free(slumt_data->s_U);

  free(slumt_data); 
  slumt_data = NULL;
 
  free(IDA_mem->ida_lmem); 
  IDA_mem->ida_lmem = NULL;

  return(IDASLS_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * IDASuperLUMTSetOrdering sets the ordering used by SuperLUMT for reducing fill.
 * Options are: 
 * 0 for natural ordering
 * 1 for minimal degree ordering on A'*A
 * 2 for minimal degree ordering on A'+A
 * 3 for approximate minimal degree ordering for unsymmetric matrices
 * The default used in SUNDIALS is 3 for COLAMD.
 * -----------------------------------------------------------------
 */

int IDASuperLUMTSetOrdering(void *ida_mem_v, int ordering_choice)
{
  IDAMem ida_mem;
  IDASlsMem idasls_mem;
  SLUMTData slumt_data;

 /* Return immediately if ida_mem is NULL */
  if (ida_mem_v == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASLS", "IDASuperLUMTSetOrdering",
		    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  ida_mem = (IDAMem) ida_mem_v;
  idasls_mem = (IDASlsMem) ida_mem->ida_lmem;

 /* Return if ordering choice argument is not valid */
  if ( (ordering_choice != 0) && (ordering_choice != 1) && 
       (ordering_choice != 2) && (ordering_choice != 3) ) {
    IDAProcessError(NULL, IDASLS_ILL_INPUT, "IDASLS", "IDASuperLUMTSetOrdering",
		    MSGSP_ILL_INPUT);
    return(IDASLS_ILL_INPUT);
  }

  slumt_data = (SLUMTData) idasls_mem->s_solver_data;

  slumt_data->s_ordering = ordering_choice;

  return(IDASLS_SUCCESS);
}



/* 
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */

/*
 * IDASuperLUMTB is a wrapper around IDASuperLUMT.
 */

int IDASuperLUMTB(void *ida_mem, int which, int num_threads, int n, int nnz)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  IDASlsMemB idaslsB_mem;
  void *ida_memB;
  int flag;
  
  /* Is ida_mem allright? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDASuperLUMTB", 
		    MSGSP_CAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASLS_NO_ADJ, "IDASSLS", "IDASuperLUMTB",  
		    MSGSP_NO_ADJ);
    return(IDASLS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", "IDASuperLUMTB", 
		    MSGSP_BAD_WHICH);
    return(IDASLS_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  /* Alloc memory for IDASlsMemRecB */
  idaslsB_mem = (IDASlsMemB) malloc(sizeof(struct IDASlsMemRecB));
  if (idaslsB_mem == NULL) {
    IDAProcessError(IDAB_mem->IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", 
		    "IDASuperLUMTB", MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  
  }

  /* set matrix type and initialize Jacob function. */
  idaslsB_mem->s_djacB = NULL;

  /* Attach lmemB data and lfreeB function. */
  IDAB_mem->ida_lmem  = idaslsB_mem;
  IDAB_mem->ida_lfree = IDASuperLUMTFreeB;

  /* Call IDASuperLUMT to the IDAS data of the backward problem. */
  ida_memB = (void *)IDAB_mem->IDA_mem;
  flag = IDASuperLUMT(ida_memB, num_threads, n, nnz);

  if (flag != IDASLS_SUCCESS) {
    free(idaslsB_mem);
    idaslsB_mem = NULL;
  }

  return(flag);
}


/*
 * IDASuperLUMTSetOrderingB is a wrapper around IDASuperLUMTSetOrdering.
 */
int IDASuperLUMTSetOrderingB(void *ida_mem, int which, int ordering_choiceB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  int flag;
  
  /* Is ida_mem allright? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDASuperLUMTB", 
		    MSGSP_CAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASLS_NO_ADJ, "IDASSLS", "IDASuperLUMTB",  
		    MSGSP_NO_ADJ);
    return(IDASLS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", "IDASuperLUMTB", 
		    MSGSP_BAD_WHICH);
    return(IDASLS_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  ida_memB = (void *)(IDAB_mem->IDA_mem);
  
  flag = IDASuperLUMTSetOrdering(ida_memB, ordering_choiceB);

  return(flag);
}




/*
 * IDASuperLUMTFreeB frees the linear solver's memory for that backward problem passed 
 * as argument. 
 */

static int IDASuperLUMTFreeB(IDABMem IDAB_mem)
{
  IDASlsMemB idaslsB_mem;

  idaslsB_mem = (IDASlsMemB) IDAB_mem->ida_lmem;

  free(idaslsB_mem);

  return(0);
}

