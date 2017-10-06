/*
 * -----------------------------------------------------------------
 * $Revision: 4924 $
 * $Date: 2016-09-19 14:36:05 -0700 (Mon, 19 Sep 2016) $
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
 * This is the implementation file for the KINSuperLUMT linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "sundials/sundials_math.h"
#include "sundials/sundials_superlumt_impl.h"

#include "kinsol/kinsol_superlumt.h"
#include "kinsol_impl.h"
#include "kinsol_sparse_impl.h"

/* Constants */

#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* KINSuperLUMT linit, lsetup, lsolve, and lfree routines */
 
static int kinSuperLUMTInit(KINMem kin_mem);
static int kinSuperLUMTSetup(KINMem kin_mem);
static int kinSuperLUMTSolve(KINMem kin_mem, N_Vector x, N_Vector b,
		       realtype *sJpnorm, realtype *sFdotJp);		       
static int kinSuperLUMTFree(KINMem kin_mem);

/*
 * -----------------------------------------------------------------
 * KINSuperLUMT
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the KINSOL / SuperLUMT linear solver module.  
 * KINSuperLUMT first calls the existing lfree routine if this is not NULL.
 * Then it sets the kin_linit, kin_lsetup, kin_lsolve, and
 * kin_lfree fields in (*kin_mem) to be kinSuperLUMTInit, kinSuperLUMTSetup,
 * kinSuperLUMTSolve, NULL, and kinSuperLUMTFree, respectively.
 * It allocates memory for a structure of type kinsluMemRec and sets
 * the kin_lmem field in (*kin_mem) to the address of this structure.
 * It sets setupNonNull in (*kin_mem) to TRUE.
 * Finally, it allocates memory for SuperLUMT.
 * The return value is KINSLS_SUCCESS = 0, KINSLS_LMEM_FAIL = -1,
 * or KINSLS_ILL_INPUT = -2.
 *
 * NOTE: The SuperLUMT linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINSuperLUMT will first 
 *       test for a compatible N_Vector internal representation
 *       by checking that the function N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */

int KINSuperLUMT(void *kin_mem_v, int num_threads, int n, int nnz)
{
  KINMem kin_mem;
  KINSlsMem kinsls_mem;
  SLUMTData slumt_data;
  int *perm_c, *perm_r;
  int nrhs, panel_size, relax;
  double *bd;
  SuperMatrix *B;

  /* Return immediately if kin_mem is NULL. */
  if (kin_mem_v == NULL) {
    KINProcessError(NULL, KINSLS_MEM_NULL, "KINSLS", "KINSuperLUMT", 
		    MSGSP_KINMEM_NULL);
    return(KINSLS_MEM_NULL);
  }
  kin_mem = (KINMem) kin_mem_v;

  /* Test if the NVECTOR package is compatible with the Direct solver */
  if (kin_mem->kin_vtemp1->ops->nvgetarraypointer == NULL) {
    KINProcessError(kin_mem, KINSLS_ILL_INPUT, "KINSLS", "KINSuperLUMT", 
		    MSGSP_BAD_NVECTOR);
    return(KINSLS_ILL_INPUT);
  }

  if (kin_mem->kin_lfree != NULL) kin_mem->kin_lfree(kin_mem);

  /* Set five main function fields in kin_mem. */
  kin_mem->kin_linit  = kinSuperLUMTInit;
  kin_mem->kin_lsetup = kinSuperLUMTSetup;
  kin_mem->kin_lsolve = kinSuperLUMTSolve;
  kin_mem->kin_lfree  = kinSuperLUMTFree;

  /* Get memory for kinSlsMemRec. */
  kinsls_mem = (KINSlsMem) malloc(sizeof(struct KINSlsMemRec));
  if (kinsls_mem == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "KINSuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }

  /* Get memory for SLUMTData. */
  slumt_data = (SLUMTData)malloc(sizeof(struct SLUMTDataRec));
  if (slumt_data == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "KINSuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }

  kin_mem->kin_setupNonNull = TRUE;

  /* Set default Jacobian routine and Jacobian data */
  kinsls_mem->s_jaceval = NULL;
  kinsls_mem->s_jacdata = kin_mem->kin_user_data;

  /* Allocate memory for the sparse Jacobian */
  kinsls_mem->s_JacMat = SparseNewMat(n, n, nnz, CSC_MAT);
  if (kinsls_mem->s_JacMat == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "KINSuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }

   /* Set up memory for the permutations */
  perm_r = (int *)malloc(n*sizeof(int));
  if (perm_r == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "kinSuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }
  perm_c = (int *)malloc(n*sizeof(int));
  if (perm_c == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "kinSuperLUMT", 
		    MSGSP_MEM_FAIL);
    free(perm_r);
    return(KINSLS_MEM_FAIL);
  }
  slumt_data->perm_r = perm_r;
  slumt_data->perm_c = perm_c;

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;

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

  dCreate_CompCol_Matrix(slumt_data->s_A, kinsls_mem->s_JacMat->M, kinsls_mem->s_JacMat->N, 
			 kinsls_mem->s_JacMat->NNZ, kinsls_mem->s_JacMat->data, 
			 kinsls_mem->s_JacMat->indexvals, kinsls_mem->s_JacMat->indexptrs, 
			 SLU_NC, SLU_D, SLU_GE);

  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  StatAlloc(kinsls_mem->s_JacMat->N, num_threads, panel_size, relax, slumt_data->Gstat);
  
  /* Create RHS matrix */
  nrhs = 1;
  bd = NULL;
  B = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  B->Store = NULL;
  dCreate_Dense_Matrix(B, n, nrhs, bd, n, 
		       SLU_DN, SLU_D, SLU_GE);
  slumt_data->s_B = B;

  /* Set ordering to COLAMD as the kinsol default use.
     Users can set a different value with KINSuperLUMTSetOrdering,
     and the user-set value is loaded before any call to factorize the
     matrix in kinSuperLUMTSetup.  */
  slumt_data->s_ordering = 3;

  /* This is a direct linear solver */
  kin_mem->kin_inexact_ls = FALSE;

  /* Attach linear solver memory to the nonlinear solver memory */
  kinsls_mem->s_solver_data = (void *) slumt_data;
  kin_mem->kin_lmem = kinsls_mem;

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;

  return(KINSLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINSuperLUMT interface functions
 * -----------------------------------------------------------------
 */

/*
  This routine does remaining initializations specific to the KINSuperLUMT
  linear solver module.  
  It returns 0 if successful.
*/

static int kinSuperLUMTInit(KINMem kin_mem)
{
  int num_threads, n;
  KINSlsMem kinsls_mem;
  SLUMTData slumt_data;

  kinsls_mem = (KINSlsMem)kin_mem->kin_lmem;
  slumt_data = (SLUMTData) kinsls_mem->s_solver_data;

  kinsls_mem->s_nje = 0;
  kinsls_mem->s_first_factorize = 1;

  /* ------------------------------------------------------------
     Allocate storage and initialize statistics variables. 
     ------------------------------------------------------------*/
  n = kinsls_mem->s_JacMat->N;
  num_threads = slumt_data->num_threads;

  StatInit(n, num_threads, slumt_data->Gstat);

  kinsls_mem->s_last_flag = 0;
  return(0);
}

/*
  This routine does the setup operations for the KINSuperLUMT linear 
  solver module.  It calls the Jacobian evaluation routine,
  updates counters, and calls the LU factorization routine.
  The return value is either
     KINSLS_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
*/

static int kinSuperLUMTSetup(KINMem kin_mem)
{
  int retval, info;
  int nprocs, panel_size, relax, permc_spec, lwork;
  int *perm_r, *perm_c;
  double diag_pivot_thresh, drop_tol;
  fact_t fact;
  trans_t trans;
  yes_no_t refact, usepr;
  KINSlsMem kinsls_mem;
  KINSlsSparseJacFn jaceval;
  SuperMatrix *A, *AC, *L, *U;
  Gstat_t *Gstat;
  superlumt_options_t *superlumt_options;
  SLUMTData slumt_data;
  SlsMat JacMat;
  void *jacdata;
  void *work;

  kinsls_mem = (KINSlsMem) (kin_mem->kin_lmem);

  slumt_data = (SLUMTData) kinsls_mem->s_solver_data;

  jaceval = kinsls_mem->s_jaceval;
  jacdata = kinsls_mem->s_jacdata;
  JacMat = kinsls_mem->s_JacMat;

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

  /* Set option values for SuperL_MT */
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
    KINProcessError(kin_mem, KINSLS_JAC_NOSET, "KINSLS", "kinSuperLUMTSetup", 
		    MSGSP_JAC_NOSET);
    free(kinsls_mem); kinsls_mem = NULL;
    return(KINSLS_JAC_NOSET);
  }

  /* Increment nje counter and call Jacobian eval routine. */
  kinsls_mem->s_nje++;
  retval = jaceval(kin_mem->kin_uu, kin_mem->kin_fval, JacMat, jacdata, 
		   kin_mem->kin_vtemp1, kin_mem->kin_vtemp2);

  if (retval < 0) {
    KINProcessError(kin_mem, KINSLS_JACFUNC_UNRECVR, "KINSLS", 
		    "kinSuperLUMTSetup", MSGSP_JACFUNC_FAILED);
    kinsls_mem->s_last_flag = KINSLS_JACFUNC_UNRECVR;
    return(KINSLS_JACFUNC_UNRECVR);
  }
  if (retval > 0) {
    kinsls_mem->s_last_flag = KINSLS_JACFUNC_RECVR;
    return(+1);
  }

  if (kinsls_mem->s_first_factorize) {
    /* Update the ordering option with any user-updated values from 
       calls to KINSuperLUMTSetOrdering */
    permc_spec = slumt_data->s_ordering;
    get_perm_c(permc_spec, A, perm_c);
 
    refact= NO;
    kinsls_mem->s_first_factorize = 0;
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
    kinsls_mem->s_last_flag = info;
    return(+1);
  }

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;

  return(0);
}

/*
  This routine handles the solve operation for the KINSuperLUMT linear
  solver module.  It calls the SuperLUMT solve routine, 
  then returns KINSLS_SUCCESS = 0.
*/

static int kinSuperLUMTSolve(KINMem kin_mem, N_Vector x, N_Vector b,
			     realtype *sJpnorm, realtype *sFdotJp)		       
{
  int info, trans;
  int *perm_r, *perm_c;
  KINSlsMem kinsls_mem;
  SuperMatrix *L, *U, *B;
  Gstat_t *Gstat;
  DNformat *Bstore;
  SLUMTData slumt_data;
  realtype *xd;
  
  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;
  slumt_data = (SLUMTData) kinsls_mem->s_solver_data;

  L = slumt_data->s_L;
  U = slumt_data->s_U;
  perm_r = slumt_data->perm_r;
  perm_c = slumt_data->perm_c;
  Gstat = slumt_data->Gstat;
  B = slumt_data->s_B;
   
  /* Copy the right-hand side into x */
  N_VScale(ONE, b, x);
  xd = N_VGetArrayPointer(x);
  Bstore = (DNformat *) (B->Store);
  Bstore->nzval = xd;

  /* Call SuperLUMT to solve the linear system using L and U */
  trans = NOTRANS;
  dgstrs(trans, L, U, perm_r, perm_c, B, Gstat, &info);

  Bstore->nzval = NULL;

  /* Compute the term sFdotJp for use in the linesearch routine.
     This term is subsequently corrected if the step is reduced by
     constraints or the linesearch.

     sFdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale.                            */

  N_VProd(b, kin_mem->kin_fscale, b);
  N_VProd(b, kin_mem->kin_fscale, b);
  *sFdotJp = N_VDotProd(kin_mem->kin_fval, b);

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;
  return(KINSLS_SUCCESS);
}

/*
  This routine frees memory specific to the KINSuperLUMT linear solver.
*/

static int kinSuperLUMTFree(KINMem kin_mem)
{
  KINSlsMem kinsls_mem;
  SLUMTData slumt_data;
  
  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;
  slumt_data = (SLUMTData) kinsls_mem->s_solver_data;

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
  if (kinsls_mem->s_JacMat) {
    SparseDestroyMat(kinsls_mem->s_JacMat);
    kinsls_mem->s_JacMat = NULL;
  }

  free(slumt_data->s_B);
  free(slumt_data->s_A);
  free(slumt_data->s_AC);
  free(slumt_data->s_L);
  free(slumt_data->s_U);

  free(slumt_data); 
  free(kin_mem->kin_lmem); 

  return(0);
}


/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * KINSuperLUMTSetOrdering sets the ordering used by SuperLUMT for reducing fill.
 * Options are: 
 * 0 for natural ordering
 * 1 for minimal degree ordering on A'*A
 * 2 for minimal degree ordering on A'+A
 * 3 for approximate minimal degree ordering for unsymmetric matrices
 * The default used in SUNDIALS is 3 for COLAMD.
 * -----------------------------------------------------------------
 */

int KINSuperLUMTSetOrdering(void *kin_mem_v, int ordering_choice)
{
  KINMem kin_mem;
  KINSlsMem kinsls_mem;
  SLUMTData slumt_data;

 /* Return immediately if kin_mem is NULL */
  if (kin_mem_v == NULL) {
    KINProcessError(NULL, KINSLS_MEM_NULL, "KINSLS", "KINSuperLUMTSetOrdering",
		    MSGSP_KINMEM_NULL);
    return(KINSLS_MEM_NULL);
  }
  kin_mem = (KINMem) kin_mem_v;
  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;

 /* Return if ordering choice argument is not valid */
  if ( (ordering_choice != 0) && (ordering_choice != 1) && 
       (ordering_choice != 2) && (ordering_choice != 3) ) {
    KINProcessError(NULL, KINSLS_ILL_INPUT, "KINSLS", "KINSuperLUMTSetOrdering",
		    MSGSP_ILL_INPUT);
    return(KINSLS_ILL_INPUT);
  }

  slumt_data = (SLUMTData) kinsls_mem->s_solver_data;

  slumt_data->s_ordering = ordering_choice;

  return(KINSLS_SUCCESS);
}

