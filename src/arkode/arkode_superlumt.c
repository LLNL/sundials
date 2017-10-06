/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Implementation file for the ARKSUPERLUMT linear solver module.
 *---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_superlumt.h"
#include "arkode/arkode_sparse.h"
#include "arkode_sparse_impl.h"
#include "arkode_impl.h"

#include <sundials/sundials_superlumt_impl.h>
#include <sundials/sundials_math.h>

/* Constants */
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)

/* ARKSUPERLUMT linit, lsetup, lsolve, and lfree routines */
static int arkSuperLUMTInit(ARKodeMem ark_mem);
static int arkSuperLUMTSetup(ARKodeMem ark_mem, int convfail, 
			     N_Vector ypred, N_Vector fpred, 
			     booleantype *jcurPtr, N_Vector tmp1, 
			     N_Vector tmp2, N_Vector tmp3);
static int arkSuperLUMTSolve(ARKodeMem ark_mem, N_Vector b, 
			     N_Vector weight, N_Vector ycur, 
			     N_Vector fcur);
static int arkSuperLUMTFree(ARKodeMem ark_mem);

/* ARKSUPERLUMT minit, msetup, msolve, mfree and mtimes routines */
static int arkMassSuperLUMTInit(ARKodeMem ark_mem);
static int arkMassSuperLUMTSetup(ARKodeMem ark_mem, N_Vector tmp1,
			     N_Vector tmp2, N_Vector tmp3);
static int arkMassSuperLUMTSolve(ARKodeMem ark_mem, N_Vector b, 
			     N_Vector weight);
static int arkMassSuperLUMTFree(ARKodeMem ark_mem);
static int arkMassSuperLUMTMultiply(N_Vector v, N_Vector Mv, 
				    realtype t, void *arkode_mem);


/*---------------------------------------------------------------
 ARKSuperLUMT

 This routine initializes the memory record and sets various 
 function fields specific to the ARKODE / SuperLUMT linear solver
 module.  ARKSuperLUMT first calls the existing lfree routine if
 this is not NULL.  Then it sets the ark_linit, ark_lsetup, 
 ark_lsolve and ark_lfree fields in (*ark_mem) to be 
 arkSuperLUMTInit, arkSuperLUMTSetup, arkSuperLUMTSolve, and 
 arkSuperLUMTFree, respectively.  It allocates memory for a 
 structure of type ARKSlsMemRec and sets the ark_lmem field in 
 (*arkode_mem) to the address of this structure.  It sets 
 setupNonNull in (*arkode_mem) to TRUE.  Finally, it allocates 
 memory for SuperLUMT.  The return value is ARKSLS_SUCCESS = 0, 
 ARKSLS_LMEM_FAIL = -1, or ARKSLS_ILL_INPUT = -2.

 NOTE: The SuperLUMT linear solver assumes a serial 
       implementation of the NVECTOR package. Therefore, 
       ARKSuperLUMT will first test for a compatible N_Vector 
       internal representation by checking that the function 
       N_VGetArrayPointer exists.
---------------------------------------------------------------*/
int ARKSuperLUMT(void *arkode_mem, int num_threads, int n, int nnz)
{
  ARKodeMem ark_mem;
  ARKSlsMem arksls_mem;
  SLUMTData slumt_data;
  int *perm_c, *perm_r;
  int nrhs, panel_size, relax;
  double *bd;
  SuperMatrix *B;

  /* Return immediately if ark_mem is NULL. */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKSuperLUMT", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKSLS_ILL_INPUT, "ARKSLS", 
		    "ARKSuperLUMT", MSGSP_BAD_NVECTOR);
    return(ARKSLS_ILL_INPUT);
  }

  if (ark_mem->ark_lfree != NULL) ark_mem->ark_lfree(ark_mem);

  /* Set four main function fields in ark_mem. */
  ark_mem->ark_linit  = arkSuperLUMTInit;
  ark_mem->ark_lsetup = arkSuperLUMTSetup;
  ark_mem->ark_lsolve = arkSuperLUMTSolve;
  ark_mem->ark_lfree  = arkSuperLUMTFree;
  ark_mem->ark_lsolve_type = 3;

  /* Get memory for ARKSlsMemRec. */
  arksls_mem = (ARKSlsMem) malloc(sizeof(struct ARKSlsMemRec));
  if (arksls_mem == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKSuperLUMT", MSGSP_MEM_FAIL);
    return(ARKSLS_MEM_FAIL);
  }

  /* Get memory for SLUMTData. */
  slumt_data = (SLUMTData) malloc(sizeof(struct SLUMTDataRec));
  if (slumt_data == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKSuperLUMT", MSGSP_MEM_FAIL);
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Initialize Jacobian-related data */
  arksls_mem->s_Jeval = NULL;
  arksls_mem->s_Jdata = ark_mem->ark_user_data;
  ark_mem->ark_setupNonNull = TRUE;
  arksls_mem->sparsetype = CSC_MAT;

  /* Initialize counters */
  arksls_mem->s_nje = 0;
  arksls_mem->s_first_factorize = 1;
  arksls_mem->s_nstlj = 0;

  /* Allocate memory for the sparse Jacobian */
  arksls_mem->s_A = NULL;
  arksls_mem->s_A = SparseNewMat(n, n, nnz, CSC_MAT);
  if (arksls_mem->s_A == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKSuperLUMT", MSGSP_MEM_FAIL);
    free(slumt_data); slumt_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Allocate memory for saved sparse Jacobian */
  arksls_mem->s_savedJ = NULL;
  arksls_mem->s_savedJ = SparseNewMat(n, n, nnz, CSC_MAT);
  if (arksls_mem->s_savedJ == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKSuperLUMT", MSGSP_MEM_FAIL);
    SparseDestroyMat(arksls_mem->s_A);
    free(slumt_data); slumt_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Set up memory for the permutations */
  perm_r = NULL;
  perm_r = (int *) malloc(n*sizeof(int));
  if (perm_r == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKSuperLUMT", MSGSP_MEM_FAIL);
    SparseDestroyMat(arksls_mem->s_A);
    SparseDestroyMat(arksls_mem->s_savedJ);
    free(slumt_data); slumt_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }
  perm_c = NULL;
  perm_c = (int *) malloc(n*sizeof(int));
  if (perm_c == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKSuperLUMT", MSGSP_MEM_FAIL);
    SparseDestroyMat(arksls_mem->s_A);
    SparseDestroyMat(arksls_mem->s_savedJ);
    free(slumt_data); slumt_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    free(perm_r); perm_r = NULL;
    return(ARKSLS_MEM_FAIL);
  }
  slumt_data->perm_r = perm_r;
  slumt_data->perm_c = perm_c;

  /* Set default parameters for SuperLU */
  slumt_data->num_threads = num_threads;
  slumt_data->diag_pivot_thresh = 1.0;

  /* Allocate structures for SuperLU */
  slumt_data->Gstat = (Gstat_t *) malloc(sizeof(Gstat_t));
  slumt_data->s_A   = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  slumt_data->s_AC  = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  slumt_data->s_L   = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  slumt_data->s_U   = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  slumt_data->s_A->Store  = NULL;
  slumt_data->s_AC->Store = NULL;
  slumt_data->s_L->Store  = NULL;
  slumt_data->s_U->Store  = NULL;
  slumt_data->superlumt_options = 
    (superlumt_options_t *) malloc(sizeof(superlumt_options_t));
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  StatAlloc(arksls_mem->s_A->N, num_threads, panel_size, 
	    relax, slumt_data->Gstat);
  
  /* Create RHS matrix */
  nrhs = 1;
  bd = NULL;
  B = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  B->Store = NULL;
  dCreate_Dense_Matrix(B, n, nrhs, bd, n, 
		       SLU_DN, SLU_D, SLU_GE);
  slumt_data->s_B = B;

  /* Set ordering to COLAMD as the arkode default use.
     Users can set a different value with ARKSuperLUMTSetOrdering,
     and the user-set value is loaded before any call to 
     factorize the matrix in arkSuperLUMTSetup.  */
  slumt_data->s_ordering = 3;

  /* Attach linear solver memory to the integrator memory */
  arksls_mem->s_solver_data = (void *) slumt_data;
  ark_mem->ark_lmem = arksls_mem;

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkSuperLUMTInit:

 This routine does remaining initializations specific to the 
 ARKSuperLUMT linear solver module.  It returns 0 if successful.
---------------------------------------------------------------*/
static int arkSuperLUMTInit(ARKodeMem ark_mem)
{
  ARKSlsMem arksls_mem;
  SLUMTData slumt_data;

  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

  arksls_mem->s_nje = 0;
  arksls_mem->s_first_factorize = 1;
  arksls_mem->s_nstlj = 0;

  /* Allocate storage and initialize statistics variables. */
  StatInit(arksls_mem->s_A->N, slumt_data->num_threads, 
	   slumt_data->Gstat);

  arksls_mem->s_Jdata = ark_mem->ark_user_data;
  arksls_mem->s_last_flag = 0;
  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkSuperLUMTSetup:

  This routine does the setup operations for the ARKSuperLUMT 
  linear solver module.  It calls the Jacobian evaluation 
  routine, updates counters, and calls the LU factorization 
  routine.  The return value is either
     ARKSLS_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
---------------------------------------------------------------*/
static int arkSuperLUMTSetup(ARKodeMem ark_mem, int convfail, 
			     N_Vector ypred, N_Vector fpred, 
			     booleantype *jcurPtr, N_Vector vtemp1, 
			     N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  int retval;
  int panel_size, relax, lwork;
  realtype dgamma;
  double drop_tol;
  fact_t fact;
  trans_t trans;
  yes_no_t refact, usepr;
  ARKSlsMem arksls_mem;
  ARKSlsMassMem arksls_mass_mem;
  SLUMTData slumt_data;
  void *work;
  
  arksls_mem = (ARKSlsMem) (ark_mem->ark_lmem);
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

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
  if (arksls_mem->s_Jeval == NULL) {
    arkProcessError(ark_mem, ARKSLS_JAC_NOSET, "ARKSLS", 
		    "arkSuperLUMTSetup", MSGSP_JAC_NOSET);
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_JAC_NOSET);
  }

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((ark_mem->ark_gamma/ark_mem->ark_gammap) - ONE);
  jbad = (ark_mem->ark_nst == 0) || 
    (ark_mem->ark_nst > arksls_mem->s_nstlj + ARKS_MSBJ) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKS_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
  
  /* If jok = TRUE, use saved copy of J */
  if (jok) {
    *jcurPtr = FALSE;
    SparseCopyMat(arksls_mem->s_savedJ, arksls_mem->s_A);

  /* If jok = FALSE, call jac routine for new J value */
  } else {
    arksls_mem->s_nje++;
    arksls_mem->s_nstlj = ark_mem->ark_nst;
    *jcurPtr = TRUE;
    SparseSetMatToZero(arksls_mem->s_A);

    retval = arksls_mem->s_Jeval(ark_mem->ark_tn, ypred, fpred, 
				 arksls_mem->s_A, arksls_mem->s_Jdata, 
				 vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      arkProcessError(ark_mem, ARKSLS_JACFUNC_UNRECVR, "ARKSLS", 
		      "arkSuperLUMTSetup", MSGSP_JACFUNC_FAILED);
      arksls_mem->s_last_flag = ARKSLS_JACFUNC_UNRECVR;
      return(ARKSLS_JACFUNC_UNRECVR);
    }
    if (retval > 0) {
      arksls_mem->s_last_flag = ARKSLS_JACFUNC_RECVR;
      return(ARKSLS_JACFUNC_RECVR);
    }

    SparseCopyMat(arksls_mem->s_A, arksls_mem->s_savedJ);
  }

  /* Scale J by -gamma */
  SparseScaleMat(-ark_mem->ark_gamma, arksls_mem->s_A);

  /* Add mass matrix to get A = M-gamma*J */
  if (ark_mem->ark_mass_matrix) {

    /* Compute mass matrix */
    arksls_mass_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
    SparseSetMatToZero(arksls_mass_mem->s_M);
    retval = arksls_mass_mem->s_Meval(ark_mem->ark_tn, 
				      arksls_mass_mem->s_M, 
				      arksls_mass_mem->s_Mdata, 
				      vtemp1, vtemp2, vtemp3);
    arksls_mass_mem->s_nme++;
    if (retval < 0) {
      arkProcessError(ark_mem, ARKSLS_MASSFUNC_UNRECVR, "ARKSLS", 
		      "arkSuperLUMTSetup",  MSGSP_MASSFUNC_FAILED);
      arksls_mem->s_last_flag = ARKSLS_MASSFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arksls_mem->s_last_flag = ARKSLS_MASSFUNC_RECVR;
      return(1);
    }
    
    /* add to A */
    retval = SparseAddMat(arksls_mem->s_A, arksls_mass_mem->s_M);
    if (retval < 0) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "arkSuperLUMTSetup",  
		      "Error in adding mass matrix to Jacobian");
      arksls_mem->s_last_flag = ARKSLS_PACKAGE_FAIL;
      return(retval);
    }
    if (retval > 0)  return(retval);
    
  } else {
    SparseAddIdentityMat(arksls_mem->s_A);
  }

  /* free and reallocate sparse matrix */
  if (slumt_data->s_A->Store) {
    SUPERLU_FREE(slumt_data->s_A->Store);
  }
  dCreate_CompCol_Matrix(slumt_data->s_A, arksls_mem->s_A->M, 
			 arksls_mem->s_A->N, arksls_mem->s_A->NNZ, 
			 arksls_mem->s_A->data, 
			 arksls_mem->s_A->indexvals, 
			 arksls_mem->s_A->indexptrs, 
			 SLU_NC, SLU_D, SLU_GE);

  /* On first decomposition, set up reusable pieces */ 
  if (arksls_mem->s_first_factorize) {

    /* Get column permutation vector perm_c[], according to s_ordering */
    get_perm_c(slumt_data->s_ordering, slumt_data->s_A, 
	       slumt_data->perm_c);
    refact= NO;
    arksls_mem->s_first_factorize = 0;

  } else {
    /* Re-initialize statistics variables */
    StatInit(arksls_mem->s_A->N, slumt_data->num_threads, 
	     slumt_data->Gstat);
    Destroy_CompCol_Permuted(slumt_data->s_AC);
    refact= YES;
  }

  /* Initialize the option structure superlumt_options using the
     user-input parameters. Subsequent calls will re-initialize
     options.  Apply perm_c to columns of original A to form AC */
  pdgstrf_init(slumt_data->num_threads, fact, trans, refact, 
	       panel_size, relax, slumt_data->diag_pivot_thresh, 
	       usepr, drop_tol, slumt_data->perm_c, 
	       slumt_data->perm_r, work, lwork, slumt_data->s_A, 
	       slumt_data->s_AC, slumt_data->superlumt_options, 
	       slumt_data->Gstat);

  /* Compute the LU factorization of A.
     The following routine will create num_threads threads. */
  pdgstrf(slumt_data->superlumt_options, slumt_data->s_AC, 
	  slumt_data->perm_r, slumt_data->s_L, slumt_data->s_U, 
	  slumt_data->Gstat, &retval);
  if (retval != 0) {
    arksls_mem->s_last_flag = retval;
    return(ARKSLS_PACKAGE_FAIL);
  }

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkSuperLUMTSolve:

 This routine handles the solve operation for the ARKSuperLUMT 
 linear solver module.  It calls the SuperLU_MT solve routine,
 then returns ARKSLS_SUCCESS = 0.
---------------------------------------------------------------*/
static int arkSuperLUMTSolve(ARKodeMem ark_mem, N_Vector b, 
			     N_Vector weight, N_Vector ycur, 
			     N_Vector fcur)
{
  int info;
  ARKSlsMem arksls_mem;
  DNformat *Bstore;
  SLUMTData slumt_data;
  realtype *bd;
  trans_t trans;
  
  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

  bd = N_VGetArrayPointer(b);
  Bstore = (DNformat *) (slumt_data->s_B->Store);
  Bstore->nzval = bd;

  /* Call SuperLUMT to solve the linear system using L and U */
  trans = NOTRANS;
  dgstrs(trans, slumt_data->s_L, slumt_data->s_U, 
	 slumt_data->perm_r, slumt_data->perm_c, 
	 slumt_data->s_B, slumt_data->Gstat, &info);

  /* Scale the correction to account for change in gamma. */
  if (ark_mem->ark_gamrat != ONE) 
    N_VScale(TWO/(ONE + ark_mem->ark_gamrat), b, b);
  Bstore->nzval = NULL;

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkSuperLUMTFree:

 This routine frees memory specific to the ARKSuperLUMT linear 
 solver.
---------------------------------------------------------------*/
static int arkSuperLUMTFree(ARKodeMem ark_mem)
{
  ARKSlsMem arksls_mem;
  SLUMTData slumt_data;
  
  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

  pxgstrf_finalize(slumt_data->superlumt_options, 
		   slumt_data->s_AC);
  free(slumt_data->perm_r);
  free(slumt_data->perm_c);
  free(slumt_data->superlumt_options);
  Destroy_SuperNode_SCP( (slumt_data->s_L) );
  Destroy_CompCol_NCP( (slumt_data->s_U) );
  StatFree(slumt_data->Gstat);
  free(slumt_data->Gstat);
  
  Destroy_SuperMatrix_Store(slumt_data->s_B);
  SUPERLU_FREE(slumt_data->s_A->Store);
  if (arksls_mem->s_A) {
    SparseDestroyMat(arksls_mem->s_A);
    arksls_mem->s_A = NULL;
  }

  if (arksls_mem->s_savedJ) {
    SparseDestroyMat(arksls_mem->s_savedJ);
    arksls_mem->s_savedJ = NULL;
  }

  free(slumt_data->s_B);
  free(slumt_data->s_A);
  free(slumt_data->s_AC);
  free(slumt_data->s_L);
  free(slumt_data->s_U);

  free(slumt_data); 
  slumt_data = NULL;
  free(arksls_mem); 
  ark_mem->ark_lmem = NULL;

  return(0);
}



/*---------------------------------------------------------------
 ARKMassSuperLUMT

 This routine initializes the memory record and sets various 
 function fields specific to the ARKODE / SuperLUMT mass matrix 
 linear solver module.  ARKMassSuperLUMT first calls the 
 existing mfree routine if this is not NULL.  Then it sets the 
 ark_minit, ark_msetup, ark_msolve and ark_mfree fields in 
 (*arkode_mem) to be arkMassSuperLUMTInit, arkMassSuperLUMTSetup, 
 arkMassSuperLUMTSolve, and arkMassSuperLUMTFree, respectively.  
 It allocates memory for a structure of type ARKSlsMassMemRec and
 sets the ark_mass_mem field in (*arkode_mem) to the address of 
 this structure.  It sets MassSetupNonNull in (*arkode_mem) to 
 TRUE.  Finally, it allocates memory for SuperLUMT.  The return
 value is ARKSLS_SUCCESS = 0, ARKSLS_LMEM_FAIL = -1, or 
 ARKSLS_ILL_INPUT = -2.

 NOTE: The SuperLUMT linear solver assumes a serial 
       implementation of the NVECTOR package. Therefore, 
       ARKMassSuperLUMT will first test for a compatible 
       N_Vector internal representation by checking that the 
       function N_VGetArrayPointer exists.
---------------------------------------------------------------*/
int ARKMassSuperLUMT(void *arkode_mem, int num_threads, 
		     int n, int nnz, ARKSlsSparseMassFn smass)
{
  ARKodeMem ark_mem;
  ARKSlsMassMem arksls_mem;
  SLUMTData slumt_data;
  int *perm_c, *perm_r;
  int nrhs, panel_size, relax;
  double *bd;
  SuperMatrix *B;

  /* Return immediately if arkode_mem is NULL. */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKMassSuperLUMT", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKSLS_ILL_INPUT, "ARKSLS", 
		    "ARKMassSuperLUMT", MSGSP_BAD_NVECTOR);
    return(ARKSLS_ILL_INPUT);
  }

  if (ark_mem->ark_mfree != NULL) ark_mem->ark_mfree(ark_mem);

  /* Set main function fields in ark_mem, enable mass matrix. */
  ark_mem->ark_mass_matrix = TRUE;
  ark_mem->ark_minit  = arkMassSuperLUMTInit;
  ark_mem->ark_msetup = arkMassSuperLUMTSetup;
  ark_mem->ark_msolve = arkMassSuperLUMTSolve;
  ark_mem->ark_mfree  = arkMassSuperLUMTFree;
  ark_mem->ark_mtimes = arkMassSuperLUMTMultiply;
  ark_mem->ark_mtimes_data = (void *) ark_mem;
  ark_mem->ark_msolve_type = 3;

  /* Get memory for ARKSlsMassMemRec. */
  arksls_mem = (ARKSlsMassMem) malloc(sizeof(struct ARKSlsMassMemRec));
  if (arksls_mem == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassSuperLUMT", MSGSP_MEM_FAIL);
    return(ARKSLS_MEM_FAIL);
  }
  
  /* Get memory for SLUMTData. */
  slumt_data = (SLUMTData) malloc(sizeof(struct SLUMTDataRec));
  if (slumt_data == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassSuperLUMT", MSGSP_MEM_FAIL);
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Initialize mass-matrix-related data */
  arksls_mem->s_nme = 0;
  arksls_mem->s_first_factorize = 1;
  arksls_mem->s_Meval = smass;
  arksls_mem->s_Mdata = ark_mem->ark_user_data;
  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  ark_mem->ark_MassSetupNonNull = TRUE;
  arksls_mem->sparsetype = CSC_MAT;

  /* Allocate memory for M and M_lu */
  arksls_mem->s_M = NULL;
  arksls_mem->s_M = SparseNewMat(n, n, nnz, CSC_MAT);
  if (arksls_mem->s_M == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassSuperLUMT", MSGSP_MEM_FAIL);
    free(slumt_data); slumt_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }
  arksls_mem->s_M_lu = NULL;
  arksls_mem->s_M_lu = SparseNewMat(n, n, nnz, CSC_MAT);
  if (arksls_mem->s_M_lu == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassSuperLUMT", MSGSP_MEM_FAIL);
    SparseDestroyMat(arksls_mem->s_M);
    free(slumt_data); slumt_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Set up memory for the permutations */
  perm_r = NULL;
  perm_r = (int *) malloc(n*sizeof(int));
  if (perm_r == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassSuperLUMT", MSGSP_MEM_FAIL);
    SparseDestroyMat(arksls_mem->s_M);
    SparseDestroyMat(arksls_mem->s_M_lu);
    free(slumt_data); slumt_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }
  perm_c = NULL;
  perm_c = (int *) malloc(n*sizeof(int));
  if (perm_c == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassSuperLUMT", MSGSP_MEM_FAIL);
    SparseDestroyMat(arksls_mem->s_M);
    SparseDestroyMat(arksls_mem->s_M_lu);
    free(slumt_data); slumt_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    free(perm_r); perm_r = NULL;
    return(ARKSLS_MEM_FAIL);
  }
  slumt_data->perm_r = perm_r;
  slumt_data->perm_c = perm_c;

  /* Set default parameters for SuperLU */
  slumt_data->num_threads = num_threads;
  slumt_data->diag_pivot_thresh = 1.0;

  /* Allocate structures for SuperLU */
  slumt_data->Gstat = (Gstat_t *) malloc(sizeof(Gstat_t));
  slumt_data->s_A   = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  slumt_data->s_AC  = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  slumt_data->s_L   = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  slumt_data->s_U   = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  slumt_data->s_A->Store  = NULL;
  slumt_data->s_AC->Store = NULL;
  slumt_data->s_L->Store  = NULL;
  slumt_data->s_U->Store  = NULL;
  slumt_data->superlumt_options = 
    (superlumt_options_t *) malloc(sizeof(superlumt_options_t));
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  StatAlloc(arksls_mem->s_M->N, num_threads, panel_size, 
	    relax, slumt_data->Gstat);
  
  /* Create RHS matrix */
  nrhs = 1;
  bd = NULL;
  B = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  B->Store = NULL;
  dCreate_Dense_Matrix(B, n, nrhs, bd, n, 
		       SLU_DN, SLU_D, SLU_GE);
  slumt_data->s_B = B;

  /* Set ordering to COLAMD as the arkode default use.
     Users can set a different value with ARKMassSuperLUMTSetOrdering,
     and the user-set value is loaded before any call to 
     factorize the matrix in arkMassSuperLUMTSetup.  */
  slumt_data->s_ordering = 3;

  /* Attach linear solver memory to the integrator memory */
  arksls_mem->s_solver_data = (void *) slumt_data;
  ark_mem->ark_mass_mem = arksls_mem;

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkMassSuperLUMTInit:

 This routine does remaining initializations specific to the 
 ARKSuperLUMT mass matrix linear solver module.  It returns 0 if
 successful.
---------------------------------------------------------------*/
static int arkMassSuperLUMTInit(ARKodeMem ark_mem)
{
  ARKSlsMassMem arksls_mem;
  SLUMTData slumt_data;

  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

  arksls_mem->s_nme = 0;
  arksls_mem->s_first_factorize = 1;

  /* Allocate storage and initialize statistics variables. */
  StatInit(arksls_mem->s_M->N, slumt_data->num_threads, 
	   slumt_data->Gstat);
  arksls_mem->s_Mdata = ark_mem->ark_user_data;
  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkMassSuperLUMTSetup:

  This routine does the setup operations for the ARKMassSuperLUMT 
  linear solver module.  It calls the mass matrix evaluation 
  routine, updates counters, and calls the LU factorization 
  routine.  The return value is either
     ARKSLS_SUCCESS = 0  if successful,
     +1  if the Meval routine failed recoverably or the
         LU factorization failed, or
     -1  if the Meval routine failed unrecoverably.
---------------------------------------------------------------*/
static int arkMassSuperLUMTSetup(ARKodeMem ark_mem, N_Vector vtemp1, 
				 N_Vector vtemp2, N_Vector vtemp3)
{
  int retval;
  int panel_size, relax, lwork;
  double drop_tol;
  fact_t fact;
  trans_t trans;
  yes_no_t refact, usepr;
  ARKSlsMassMem arksls_mem;
  SLUMTData slumt_data;
  void *work;
  
  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

  /* Set option values for SuperLU_MT */
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  fact = EQUILIBRATE;
  trans = NOTRANS;
  usepr = NO;
  drop_tol = 0.0;
  lwork = 0;
  work = NULL;

  /* Check that mass matrix eval routine is set */
  if (arksls_mem->s_Meval == NULL) {
    arkProcessError(ark_mem, ARKSLS_MASS_NOSET, "ARKSLS", 
		    "arkMassSuperLUMTSetup", MSGSP_MASS_NOSET);
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MASS_NOSET);
  }

  /* call Meval routine for new M matrix */
  SparseSetMatToZero(arksls_mem->s_M);
  retval = arksls_mem->s_Meval(ark_mem->ark_tn, arksls_mem->s_M, 
			       arksls_mem->s_Mdata, vtemp1, 
			       vtemp2, vtemp3);
  arksls_mem->s_nme++;
  if (retval < 0) {
    arkProcessError(ark_mem, ARKSLS_MASSFUNC_UNRECVR, "ARKSLS", 
		    "arkMassSuperLUMTSetup", MSGSP_MASSFUNC_FAILED);
    arksls_mem->s_last_flag = ARKSLS_MASSFUNC_UNRECVR;
    return(ARKSLS_MASSFUNC_UNRECVR);
  }
  if (retval > 0) {
    arksls_mem->s_last_flag = ARKSLS_MASSFUNC_RECVR;
    return(ARKSLS_MASSFUNC_RECVR);
  }

  /* Copy M into M_lu for LU decomposition */
  SparseCopyMat(arksls_mem->s_M, arksls_mem->s_M_lu);


  /* free and reallocate sparse matrix */
  if (slumt_data->s_A->Store) {
    SUPERLU_FREE(slumt_data->s_A->Store);
  }
  dCreate_CompCol_Matrix(slumt_data->s_A, arksls_mem->s_M->M, 
			 arksls_mem->s_M->N, arksls_mem->s_M->NNZ, 
			 arksls_mem->s_M->data, 
			 arksls_mem->s_M->indexvals, 
			 arksls_mem->s_M->indexptrs, 
			 SLU_NC, SLU_D, SLU_GE);

  /* On first decomposition, set up reusable pieces */ 
  if (arksls_mem->s_first_factorize) {

    /* Get column permutation vector perm_c[], according to s_ordering */
    get_perm_c(slumt_data->s_ordering, slumt_data->s_A, 
	       slumt_data->perm_c);
    refact= NO;
    arksls_mem->s_first_factorize = 0;

  } else {
    /* Re-initialize statistics variables */
    StatInit(arksls_mem->s_M->N, slumt_data->num_threads, 
	     slumt_data->Gstat);
    Destroy_CompCol_Permuted(slumt_data->s_AC);
    refact= YES;
  }

  /* Initialize the option structure superlumt_options using the
     user-input parameters. Subsequent calls will re-initialize
     options.  Apply perm_c to columns of original A to form AC */
  pdgstrf_init(slumt_data->num_threads, fact, trans, refact, 
	       panel_size, relax, slumt_data->diag_pivot_thresh, 
	       usepr, drop_tol, slumt_data->perm_c, 
	       slumt_data->perm_r, work, lwork, slumt_data->s_A, 
	       slumt_data->s_AC, slumt_data->superlumt_options, 
	       slumt_data->Gstat);

  /* Compute the LU factorization of A.
     The following routine will create num_threads threads. */
  pdgstrf(slumt_data->superlumt_options, slumt_data->s_AC, 
	  slumt_data->perm_r, slumt_data->s_L, slumt_data->s_U, 
	  slumt_data->Gstat, &retval);
  if (retval != 0) {
    arksls_mem->s_last_flag = retval;
    return(ARKSLS_PACKAGE_FAIL);
  }

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkMassSuperLUMTSolve:

 This routine handles the solve operation for the ARKSuperLUMT 
 mass matrix linear solver module.  It calls the SuperLU_MT 
 solve routine, then returns ARKSLS_SUCCESS = 0.
---------------------------------------------------------------*/
static int arkMassSuperLUMTSolve(ARKodeMem ark_mem, N_Vector b, 
				 N_Vector weight)
{
  int info;
  ARKSlsMassMem arksls_mem;
  DNformat *Bstore;
  SLUMTData slumt_data;
  realtype *bd;
  trans_t trans;
  
  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

  bd = N_VGetArrayPointer(b);
  Bstore = (DNformat *) (slumt_data->s_B->Store);
  Bstore->nzval = bd;

  /* Call SuperLUMT to solve the linear system using L and U */
  trans = NOTRANS;
  dgstrs(trans, slumt_data->s_L, slumt_data->s_U, 
	 slumt_data->perm_r, slumt_data->perm_c, 
	 slumt_data->s_B, slumt_data->Gstat, &info);

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkMassSuperLUMTFree:

 This routine frees memory specific to the ARKSuperLUMT mass 
 matrix linear solver.
---------------------------------------------------------------*/
static int arkMassSuperLUMTFree(ARKodeMem ark_mem)
{
  ARKSlsMassMem arksls_mem;
  SLUMTData slumt_data;
  
  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

  pxgstrf_finalize(slumt_data->superlumt_options, 
		   slumt_data->s_AC);
  free(slumt_data->perm_r);
  free(slumt_data->perm_c);
  free(slumt_data->superlumt_options);
  Destroy_SuperNode_SCP( (slumt_data->s_L) );
  Destroy_CompCol_NCP( (slumt_data->s_U) );
  StatFree(slumt_data->Gstat);
  free(slumt_data->Gstat);
  
  Destroy_SuperMatrix_Store(slumt_data->s_B);
  SUPERLU_FREE(slumt_data->s_A->Store);
  if (arksls_mem->s_M) {
    SparseDestroyMat(arksls_mem->s_M);
    arksls_mem->s_M = NULL;
  }

  if (arksls_mem->s_M_lu) {
    SparseDestroyMat(arksls_mem->s_M_lu);
    arksls_mem->s_M_lu = NULL;
  }

  free(slumt_data->s_B);
  free(slumt_data->s_A);
  free(slumt_data->s_AC);
  free(slumt_data->s_L);
  free(slumt_data->s_U);

  free(slumt_data); 
  slumt_data = NULL;
  free(arksls_mem); 
  ark_mem->ark_mass_mem = NULL;

  return(0);
}


/*===============================================================
 Utility Functions
===============================================================*/

/*---------------------------------------------------------------
 arkMassSuperLUMTMultiply:

 Multiplies the mass matrix by the vector v to fill in the 
 vector Mv.
---------------------------------------------------------------*/
static int arkMassSuperLUMTMultiply(N_Vector v, N_Vector Mv, 
				    realtype t, void *arkode_mem)
{

  /* extract the SlsMassMem structure from the arkode_mem pointer */
  ARKodeMem ark_mem;
  ARKSlsMassMem arksls_mem;
  realtype *vdata=NULL, *Mvdata=NULL;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "arkMassSuperLUMTMultiply", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;

  /* zero out the result */
  N_VConst(0.0, Mv);

  /* access the vector arrays (since they must be serial vectors) */
  vdata = N_VGetArrayPointer(v);
  Mvdata = N_VGetArrayPointer(Mv);
  if (vdata == NULL || Mvdata == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "arkMassSuperLUMTMultiply", 
		    "Vector data un-allocated.");
    return(ARKSLS_MEM_NULL);
  }

  /* perform matrix-vector product with arksls_mem->s_M and return */
  if (SparseMatvec(arksls_mem->s_M, vdata, Mvdata) != 0) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "arkMassSuperLUMTMultiply", 
		    "Mass matrix data un-allocated.");
    return(ARKSLS_MEM_NULL);
  }

  return(0);

}


/*===============================================================
 Optional Input Specification Functions
===============================================================*/


/*---------------------------------------------------------------
 ARKSuperLUMTSetOrdering:

 Sets the ordering used by SuperLUMT for reducing fill.
 Options are: 
     0 for natural ordering
     1 for minimal degree ordering on A'*A
     2 for minimal degree ordering on A'+A
     3 for AMD ordering for unsymmetric matrices
 The default used in SUNDIALS is 3 for COLAMD.
---------------------------------------------------------------*/
int ARKSuperLUMTSetOrdering(void *arkode_mem, int ordering_choice)
{
  ARKodeMem ark_mem;
  ARKSlsMem arksls_mem;
  SLUMTData slumt_data;

  /* Return immediately if ark_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKSuperLUMTSetOrdering", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Return if ordering choice argument is not valid */
  if ( (ordering_choice < 0) || (ordering_choice > 3) ) {
    arkProcessError(NULL, ARKSLS_ILL_INPUT, "ARKSLS", 
		    "ARKSuperLUMTSetOrdering", MSGSP_ILL_INPUT);
    return(ARKSLS_ILL_INPUT);
  }

  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

  slumt_data->s_ordering = ordering_choice;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKMassSuperLUMTSetOrdering:

 Sets the ordering used by SuperLUMT for reducing fill. Options 
 are: 
     0 for natural ordering
     1 for minimal degree ordering on A'*A
     2 for minimal degree ordering on A'+A
     3 for AMD ordering for unsymmetric matrices
 The default used in SUNDIALS is 3 for COLAMD.
---------------------------------------------------------------*/
int ARKMassSuperLUMTSetOrdering(void *arkode_mem, int ordering_choice)
{
  ARKodeMem     ark_mem;
  ARKSlsMassMem arksls_mem;
  SLUMTData     slumt_data;

  /* Return immediately if ark_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKMassSuperLUMTSetOrdering", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Return if ordering choice argument is not valid */
  if ( (ordering_choice < 0) || (ordering_choice > 3) ) {
    arkProcessError(NULL, ARKSLS_ILL_INPUT, "ARKSLS", 
		    "ARKMassSuperLUMTSetOrdering", MSGSP_ILL_INPUT);
    return(ARKSLS_ILL_INPUT);
  }

  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  slumt_data = (SLUMTData) arksls_mem->s_solver_data;

  slumt_data->s_ordering = ordering_choice;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
   EOF
---------------------------------------------------------------*/
