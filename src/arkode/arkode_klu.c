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
 * Implementation file for the ARKKLU linear solver module.
 *---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_klu.h"
#include "arkode/arkode_sparse.h"
#include "arkode_sparse_impl.h"
#include "arkode_impl.h"

#include <sundials/sundials_klu_impl.h>
#include <sundials/sundials_math.h>

/* Constants */
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)
#define TWOTHIRDS    RCONST(0.6666666666666667)

/* ARKKLU linit, lsetup, lsolve, and lfree routines */
static int arkKLUInit(ARKodeMem ark_mem);
static int arkKLUSetup(ARKodeMem ark_mem, int convfail, 
		       N_Vector ypred, N_Vector fpred, 
		       booleantype *jcurPtr, N_Vector tmp1, 
		       N_Vector tmp2, N_Vector tmp3);
static int arkKLUSolve(ARKodeMem ark_mem, N_Vector b, 
		       N_Vector weight, N_Vector ycur, 
		       N_Vector fcur);
static void arkKLUFree(ARKodeMem ark_mem);

/* ARKKLU minit, msetup, msolve, and mfree routines */
static int arkMassKLUInit(ARKodeMem ark_mem);
static int arkMassKLUSetup(ARKodeMem ark_mem, N_Vector tmp1, 
			   N_Vector tmp2, N_Vector tmp3);
static int arkMassKLUSolve(ARKodeMem ark_mem, N_Vector b, 
			   N_Vector weight);
static void arkMassKLUFree(ARKodeMem ark_mem);
static int arkMassKLUMultiply(N_Vector v, N_Vector Mv, 
			      realtype t, void *arkode_mem);


/*---------------------------------------------------------------
 ARKKLU

 This routine initializes the memory record and sets various 
 function fields specific to the ARKode / KLU linear solver 
 module.  ARKKLU first calls the existing lfree routine if this 
 is not NULL.  Then it sets the ark_linit, ark_lsetup, ark_lsolve
 and ark_lfree fields in (*arkode_mem) to be arkKLUInit, 
 arkKLUSetup, arkKLUSolve and arkKLUFree, respectively.   It 
 allocates memory for a structure of type ARKSlsMemRec and sets 
 the ark_lmem field in (*arkode_mem) to the address of this 
 structure.  It sets setupNonNull in (*arkode_mem) to TRUE.  
 Finally, it allocates memory for KLU.  The return value is 
 ARKSLS_SUCCESS = 0, ARKSLS_LMEM_FAIL = -1, or 
 ARKSLS_ILL_INPUT = -2.

 NOTE: The KLU linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKKLU will first 
       test for a compatible N_Vector internal representation
       by checking that the function N_VGetArrayPointer exists.
---------------------------------------------------------------*/
int ARKKLU(void *arkode_mem, int n, int nnz)
{
  ARKodeMem ark_mem;
  ARKSlsMem arksls_mem;
  KLUData klu_data;
  int flag;

  /* Return immediately if ark_mem is NULL. */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKKLU", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKSLS_ILL_INPUT, "ARKSLS", 
		    "ARKKLU", MSGSP_BAD_NVECTOR);
    return(ARKSLS_ILL_INPUT);
  }

  if (ark_mem->ark_lfree != NULL) ark_mem->ark_lfree(ark_mem);

  /* Set four main function fields in ark_mem. */
  ark_mem->ark_linit  = arkKLUInit;
  ark_mem->ark_lsetup = arkKLUSetup;
  ark_mem->ark_lsolve = arkKLUSolve;
  ark_mem->ark_lfree  = arkKLUFree;
  ark_mem->ark_lsolve_type = 3;

  /* Get memory for ARKSlsMemRec. */
  arksls_mem = (ARKSlsMem) malloc(sizeof(struct ARKSlsMemRec));
  if (arksls_mem == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKKLU", MSGSP_MEM_FAIL);
    return(ARKSLS_MEM_FAIL);
  }

  /* Get memory for KLUData. */
  klu_data = (KLUData) malloc(sizeof(struct KLUDataRec));
  if (klu_data == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKKLU", MSGSP_MEM_FAIL);
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Initialize Jacobian-related data */
  arksls_mem->s_Jeval = NULL;
  arksls_mem->s_Jdata = ark_mem->ark_user_data;
  ark_mem->ark_setupNonNull = TRUE;

  /* Initialize counters */
  arksls_mem->s_nje = 0;
  arksls_mem->s_first_factorize = 1;
  arksls_mem->s_nstlj = 0;

  /* Allocate memory for the sparse Jacobian */
  arksls_mem->s_A = NULL;
  arksls_mem->s_A = NewSparseMat(n, n, nnz);
  if (arksls_mem->s_A == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKKLU", MSGSP_MEM_FAIL);
    free(klu_data); klu_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Allocate memory for saved sparse Jacobian */
  arksls_mem->s_savedJ = NULL;
  arksls_mem->s_savedJ = NewSparseMat(n, n, nnz);
  if (arksls_mem->s_savedJ == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKKLU", MSGSP_MEM_FAIL);
    DestroySparseMat(arksls_mem->s_A);
    free(klu_data); klu_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Initialize KLU structures */
  klu_data->s_Symbolic = NULL;
  klu_data->s_Numeric = NULL;

  /* Set default parameters for KLU */
  flag = klu_defaults(&klu_data->s_Common);
  if (flag == 0) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKKLU", MSGSP_MEM_FAIL);
    klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
    free(klu_data->s_Numeric);  klu_data->s_Numeric = NULL;
    klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
    free(klu_data->s_Symbolic);  klu_data->s_Symbolic = NULL;
    DestroySparseMat(arksls_mem->s_A);
    DestroySparseMat(arksls_mem->s_savedJ);
    free(klu_data); klu_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Set ordering to COLAMD as the arkode default use.
     Users can set a different value with ARKKLUSetOrdering,
     and the user-set value is loaded before any call to klu_analyze in
     ARKKLUSetup.  */
  klu_data->s_ordering = 1;
  klu_data->s_Common.ordering = klu_data->s_ordering;

  /* Attach linear solver memory to the integrator memory */
  arksls_mem->s_solver_data = (void *) klu_data;
  ark_mem->ark_lmem = arksls_mem;

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
  ARKKLUReInit

  This routine reinitializes memory and flags for a new 
  factorization (symbolic and numeric) to be conducted at the 
  next solver setup call.  This routine is useful in the cases 
  where the number of nonzeroes has changed or if the structure 
  of the linear system has changed which would require a new 
  symbolic (and numeric) factorization.
 
  The reinit_type argument governs the level of reinitialization:
 
  reinit_type = 1: The Jacobian matrix will be destroyed and 
                   a new one will be allocated based on the nnz
                   value passed to this call. New symbolic and
                   numeric factorizations will be completed at the 
                   next solver setup.
 
  reinit_type = 2: Only symbolic and numeric factorizations will be 
                   completed.  It is assumed that the Jacobian size
                   has not exceeded the size of nnz given in the 
                   prior call to ARKMassKLU.
 
  This routine assumes no other changes to solver use are necessary.
 
  The return value is ARKSLS_SUCCESS = 0, ARKSLS_MEM_NULL = -1, 
  ARKSLS_MASSMEM_NULL = -8, ARKSLS_ILL_INPUT = -3, or 
  ARKSLS_MEM_FAIL = -4.
---------------------------------------------------------------*/
int ARKKLUReInit(void *arkode_mem, int n, int nnz, int reinit_type)
{
  ARKodeMem ark_mem;
  ARKSlsMem arksls_mem;
  KLUData klu_data;

  /* Return immediately if arkode_mem is NULL. */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", "ARKKLUReInit", 
		    MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Return immediately if ark_lmem is NULL. */
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(NULL, ARKSLS_LMEM_NULL, "ARKSLS", "ARKKLUReInit", 
		    MSGSP_LMEM_NULL);
    return(ARKSLS_LMEM_NULL);
  }

  arksls_mem = (ARKSlsMem) (ark_mem->ark_lmem);
  klu_data = (KLUData) arksls_mem->s_solver_data;

  /* Return if reinit_type is not valid */
  if ((reinit_type != 1) && (reinit_type != 2)) {
    arkProcessError(NULL, ARKSLS_ILL_INPUT, "ARKSLS", "ARKKLUReInit", 
		    MSGSP_ILL_INPUT);
    return(ARKSLS_ILL_INPUT);
  }

  if (reinit_type == 1) {

    /* Destroy previous Jacobian information */
    if (arksls_mem->s_A) {
      DestroySparseMat(arksls_mem->s_A);
    }

    /* Allocate memory for the sparse Jacobian */
    arksls_mem->s_A = NewSparseMat(n, n, nnz);
    if (arksls_mem->s_A == NULL) {
      arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", "ARKKLU", 
		    MSGSP_MEM_FAIL);
      return(ARKSLS_MEM_FAIL);
    }
  }

  /* Free the prior factorazation and reset for first factorization */
  if( klu_data->s_Symbolic != NULL)
    klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
  if( klu_data->s_Numeric != NULL)
    klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
  arksls_mem->s_first_factorize = 1;

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;

  return(0);
}

/*---------------------------------------------------------------
 arkKLUInit:

 This routine does remaining initializations specific to the 
 ARKKLU linear solver module.  It returns 0 if successful.
---------------------------------------------------------------*/
static int arkKLUInit(ARKodeMem ark_mem)
{
  ARKSlsMem arksls_mem;

  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;

  arksls_mem->s_nje = 0;
  arksls_mem->s_first_factorize = 1;
  arksls_mem->s_nstlj = 0;

  arksls_mem->s_last_flag = 0;
  return(0);
}


/*---------------------------------------------------------------
 arkKLUSetup:

  This routine does the setup operations for the ARKKLU linear 
  solver module.  It calls the Jacobian evaluation routine,
  updates counters, and calls the LU factorization routine.
  The return value is either
     ARKSLS_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
---------------------------------------------------------------*/
static int arkKLUSetup(ARKodeMem ark_mem, int convfail, 
		       N_Vector ypred, N_Vector fpred, 
		       booleantype *jcurPtr, N_Vector vtemp1, 
		       N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  ARKSlsMem arksls_mem;
  ARKSlsMassMem arksls_mass_mem;
  KLUData klu_data;
  int retval;

  realtype uround_twothirds;
  
  uround_twothirds = SUNRpowerR(ark_mem->ark_uround,TWOTHIRDS);

  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;
  klu_data = (KLUData) arksls_mem->s_solver_data;
  
  /* Check that Jacobian eval routine is set */
  if (arksls_mem->s_Jeval == NULL) {
    arkProcessError(ark_mem, ARKSLS_JAC_NOSET, "ARKSLS", 
		    "arkKLUSetup", MSGSP_JAC_NOSET);
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
    CopySparseMat(arksls_mem->s_savedJ, arksls_mem->s_A);

  /* If jok = FALSE, call jac routine for new J value */
  } else {
    arksls_mem->s_nje++;
    arksls_mem->s_nstlj = ark_mem->ark_nst;
    *jcurPtr = TRUE;
    SlsSetToZero(arksls_mem->s_A);

    retval = arksls_mem->s_Jeval(ark_mem->ark_tn, ypred, fpred, 
				 arksls_mem->s_A, arksls_mem->s_Jdata, 
				 vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      arkProcessError(ark_mem, ARKSLS_JACFUNC_UNRECVR, "ARKSLS", 
		      "arkKLUSetup", MSGSP_JACFUNC_FAILED);
      arksls_mem->s_last_flag = ARKSLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arksls_mem->s_last_flag = ARKSLS_JACFUNC_RECVR;
      return(1);
    }

    CopySparseMat(arksls_mem->s_A, arksls_mem->s_savedJ);
  }

  /* Scale J by -gamma */
  ScaleSparseMat(-ark_mem->ark_gamma, arksls_mem->s_A);

  /* Add mass matrix to get A = M-gamma*J */
  if (ark_mem->ark_mass_matrix) {

    /* Compute mass matrix */
    arksls_mass_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
    SlsSetToZero(arksls_mass_mem->s_M);
    retval = arksls_mass_mem->s_Meval(ark_mem->ark_tn, 
				      arksls_mass_mem->s_M, 
				      arksls_mass_mem->s_Mdata, 
				      vtemp1, vtemp2, vtemp3);
    arksls_mass_mem->s_nme++;
    if (retval < 0) {
      arkProcessError(ark_mem, ARKSLS_MASSFUNC_UNRECVR, "ARKSLS", 
		      "arkKLUSetup",  MSGSP_MASSFUNC_FAILED);
      arksls_mem->s_last_flag = ARKSLS_MASSFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arksls_mem->s_last_flag = ARKSLS_MASSFUNC_RECVR;
      return(1);
    }
    
    /* add to A */
    retval = SlsAddMat(arksls_mem->s_A, arksls_mass_mem->s_M);
    if (retval < 0) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "arkKLUSetup",  "Error in adding mass matrix to Jacobian");
      arksls_mem->s_last_flag = ARKSLS_PACKAGE_FAIL;
      return(retval);
    }
    if (retval > 0)  return(retval);
    
  } else {
    AddIdentitySparseMat(arksls_mem->s_A);
  }


  /* On first decomposition, get the symbolic factorization */ 
  if (arksls_mem->s_first_factorize) {

    /* Update the ordering option with user-updated values */
    klu_data->s_Common.ordering = klu_data->s_ordering;

    /* Perform symbolic analysis of sparsity structure */
    klu_data->s_Symbolic = klu_analyze(arksls_mem->s_A->N, 
				       arksls_mem->s_A->colptrs, 
				       arksls_mem->s_A->rowvals, 
				       &(klu_data->s_Common));
    if (klu_data->s_Symbolic == NULL) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "ARKKLUSetup", MSGSP_PACKAGE_FAIL);
      return(ARKSLS_PACKAGE_FAIL);
    }

    /* ------------------------------------------------------------
       Compute the LU factorization of  the Jacobian.
       ------------------------------------------------------------*/
    klu_data->s_Numeric = klu_factor(arksls_mem->s_A->colptrs, 
				     arksls_mem->s_A->rowvals, 
				     arksls_mem->s_A->data, 
				     klu_data->s_Symbolic, 
				     &(klu_data->s_Common));
    if (klu_data->s_Numeric == NULL) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "ARKKLUSetup", MSGSP_PACKAGE_FAIL);
      return(ARKSLS_PACKAGE_FAIL);
    }

    arksls_mem->s_first_factorize = 0;
  }
  else {

    retval = klu_refactor(arksls_mem->s_A->colptrs, 
			  arksls_mem->s_A->rowvals, 
			  arksls_mem->s_A->data, 
			  klu_data->s_Symbolic, klu_data->s_Numeric,
			  &(klu_data->s_Common));
    if (retval == 0) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "ARKKLUSetup", MSGSP_PACKAGE_FAIL);
      return(ARKSLS_PACKAGE_FAIL);
    }
    
    /*-----------------------------------------------------------
      Check if a cheap estimate of the reciprocal of the condition 
      number is getting too small.  If so, delete
      the prior numeric factorization and recompute it.
      -----------------------------------------------------------*/
    
    retval = klu_rcond(klu_data->s_Symbolic, klu_data->s_Numeric,
		       &(klu_data->s_Common));
    if (retval == 0) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "ARKKLUSetup", MSGSP_PACKAGE_FAIL);
      return(ARKSLS_PACKAGE_FAIL);
    }

    if ( (klu_data->s_Common.rcond)  < uround_twothirds ) {
      
      /* Condition number may be getting large.  
	 Compute more accurate estimate */
      retval = klu_condest(arksls_mem->s_A->colptrs, 
			   arksls_mem->s_A->data, 
			   klu_data->s_Symbolic, klu_data->s_Numeric,
			   &(klu_data->s_Common));
      if (retval == 0) {
	arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
			"ARKKLUSetup", MSGSP_PACKAGE_FAIL);
	return(ARKSLS_PACKAGE_FAIL);
      }
      
      if ( (klu_data->s_Common.condest) > 
	   (1.0/uround_twothirds) ) {

	/* More accurate estimate also says condition number is 
	   large, so recompute the numeric factorization */

	klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
	
	klu_data->s_Numeric = klu_factor(arksls_mem->s_A->colptrs, 
					 arksls_mem->s_A->rowvals, 
					 arksls_mem->s_A->data,
					 klu_data->s_Symbolic, 
					 &(klu_data->s_Common));

	if (klu_data->s_Numeric == NULL) {
	  arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
			  "ARKKLUSetup", MSGSP_PACKAGE_FAIL);
	  return(ARKSLS_PACKAGE_FAIL);
	}
      }
    }
  }

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkKLUSolve:

 This routine handles the solve operation for the ARKKLU linear
 solver module.  It calls the KLU solve routine, then returns 
 ARKSLS_SUCCESS = 0.
---------------------------------------------------------------*/
static int arkKLUSolve(ARKodeMem ark_mem, N_Vector b, 
		       N_Vector weight, N_Vector ycur, 
		       N_Vector fcur)
{
  int flag;
  ARKSlsMem arksls_mem;
  KLUData klu_data;
  realtype *bd;
  
  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;
  klu_data = (KLUData) arksls_mem->s_solver_data;

  bd = N_VGetArrayPointer(b);

  /* Call KLU to solve the linear system */
  flag = klu_solve(klu_data->s_Symbolic, klu_data->s_Numeric, 
		   arksls_mem->s_A->N, 1, bd, 
		   &(klu_data->s_Common));
  if (flag == 0) {
    arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		    "ARKKLUSolve", MSGSP_PACKAGE_FAIL);
    return(ARKSLS_PACKAGE_FAIL);
  }

  /* Scale the correction to account for change in gamma. */
  if (ark_mem->ark_gamrat != ONE)
    N_VScale(TWO/(ONE + ark_mem->ark_gamrat), b, b);

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkKLUFree:

 This routine frees memory specific to the ARKKLU linear solver.
---------------------------------------------------------------*/
static void arkKLUFree(ARKodeMem ark_mem)
{
  ARKSlsMem arksls_mem;
  KLUData klu_data;
  
  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;
  klu_data = (KLUData) arksls_mem->s_solver_data;

  klu_free_numeric(&(klu_data->s_Numeric), 
		   &(klu_data->s_Common));

  klu_free_symbolic(&(klu_data->s_Symbolic), 
		    &(klu_data->s_Common));

  if (arksls_mem->s_A) {
    DestroySparseMat(arksls_mem->s_A);
    arksls_mem->s_A = NULL;
  }

  if (arksls_mem->s_savedJ) {
    DestroySparseMat(arksls_mem->s_savedJ);
    arksls_mem->s_savedJ = NULL;
  }

  free(klu_data); 
  free(arksls_mem); 
  ark_mem->ark_lmem = NULL;
}



/*---------------------------------------------------------------
 ARKMassKLU

 This routine initializes the memory record and sets various 
 function fields specific to the ARKode / KLU mass matrix linear 
 solver module.  ARKMassKLU first calls the existing mfree 
 routine if this is not NULL.  Then it sets the ark_minit, 
 ark_msetup, ark_msolve and ark_mfree fields in (*arkode_mem) 
 to be arkMassKLUInit, arkMassKLUSetup, arkMassKLUSolve and 
 arkMassKLUFree, respectively.   It allocates memory for a 
 structure of type ARKSlsMassMemRec and sets the ark_mass_mem 
 field in (*arkode_mem) to the address of this structure.  It 
 sets MassSetupNonNull in (*arkode_mem) to TRUE.  Finally, it 
 allocates memory for KLU.  The return value is 
 ARKSLS_SUCCESS = 0, ARKSLS_LMEM_FAIL = -1, or 
 ARKSLS_ILL_INPUT = -2.

 NOTE: The KLU linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKMassKLU will first 
       test for a compatible N_Vector internal representation
       by checking that the function N_VGetArrayPointer exists.
---------------------------------------------------------------*/
int ARKMassKLU(void *arkode_mem, int n, int nnz, 
	       ARKSlsSparseMassFn smass)
{
  ARKodeMem ark_mem;
  ARKSlsMassMem arksls_mem;
  KLUData klu_data;
  int flag;

  /* Return immediately if arkode_mem is NULL. */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKMassKLU", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKSLS_ILL_INPUT, "ARKSLS", 
		    "ARKMassKLU", MSGSP_BAD_NVECTOR);
    return(ARKSLS_ILL_INPUT);
  }

  if (ark_mem->ark_mfree != NULL) ark_mem->ark_mfree(ark_mem);

  /* Set five main function fields in ark_mem, enable mass matrix. */
  ark_mem->ark_mass_matrix = TRUE;
  ark_mem->ark_minit  = arkMassKLUInit;
  ark_mem->ark_msetup = arkMassKLUSetup;
  ark_mem->ark_msolve = arkMassKLUSolve;
  ark_mem->ark_mfree  = arkMassKLUFree;
  ark_mem->ark_mtimes = arkMassKLUMultiply;
  ark_mem->ark_mtimes_data = (void *) ark_mem;
  ark_mem->ark_msolve_type = 3;

  /* Get memory for ARKSlsMassMemRec. */
  arksls_mem = (ARKSlsMassMem) malloc(sizeof(struct ARKSlsMassMemRec));
  if (arksls_mem == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassKLU", MSGSP_MEM_FAIL);
    return(ARKSLS_MEM_FAIL);
  }

  /* Get memory for KLUData. */
  klu_data = (KLUData) malloc(sizeof(struct KLUDataRec));
  if (klu_data == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassKLU", MSGSP_MEM_FAIL);
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

  /* Allocate memory for M and M_lu */
  arksls_mem->s_M = NULL;
  arksls_mem->s_M = NewSparseMat(n, n, nnz);
  if (arksls_mem->s_M == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassKLU", MSGSP_MEM_FAIL);
    free(klu_data); klu_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }
  arksls_mem->s_M_lu = NULL;
  arksls_mem->s_M_lu = NewSparseMat(n, n, nnz);
  if (arksls_mem->s_M_lu == NULL) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassKLU", MSGSP_MEM_FAIL);
    DestroySparseMat(arksls_mem->s_M);
    free(klu_data); klu_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Initialize KLU structures */
  klu_data->s_Symbolic = NULL;
  klu_data->s_Numeric = NULL;

  /* Set default parameters for KLU */
  flag = klu_defaults(&klu_data->s_Common);
  if (flag == 0) {
    arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", 
		    "ARKMassKLU", MSGSP_MEM_FAIL);
    klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
    free(klu_data->s_Numeric);  klu_data->s_Numeric = NULL;
    klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
    free(klu_data->s_Symbolic);  klu_data->s_Symbolic = NULL;
    DestroySparseMat(arksls_mem->s_M);
    DestroySparseMat(arksls_mem->s_M_lu);
    free(klu_data); klu_data = NULL;
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MEM_FAIL);
  }

  /* Set ordering to COLAMD as the arkode default use.
     Users can set a different value with ARKMassKLUSetOrdering,
     and the user-set value is loaded before any call to klu_analyze in
     ARKMassKLUSetup.  */
  klu_data->s_ordering = 1;
  klu_data->s_Common.ordering = klu_data->s_ordering;

  /* Attach linear solver memory to the integrator memory */
  arksls_mem->s_solver_data = (void *) klu_data;
  ark_mem->ark_mass_mem = arksls_mem;

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
  ARKMassKLUReInit

  This routine reinitializes memory and flags for a new 
  factorization (symbolic and numeric) to be conducted at the 
  next solver setup call.  This routine is useful in the cases 
  where the number of nonzeroes has changed or if the structure 
  of the linear system has changed which would require a new 
  symbolic (and numeric) factorization.
 
  The reinit_type argument governs the level of reinitialization:
 
  reinit_type = 1: The mass matrix will be destroyed and 
                   a new one will be allocated based on the nnz
                   value passed to this call. New symbolic and
                   numeric factorizations will be completed at the 
                   next solver setup.
 
  reinit_type = 2: Only symbolic and numeric factorizations will be 
                   completed.  It is assumed that the mass matrix
                   size has not exceeded the size of nnz given in 
                   the prior call to ARKMassKLU.
 
  This routine assumes no other changes to solver use are necessary.
 
  The return value is ARKSLS_SUCCESS = 0, ARKSLS_MEM_NULL = -1, 
  ARKSLS_MASSMEM_NULL = -8, ARKSLS_ILL_INPUT = -3, or 
  ARKSLS_MEM_FAIL = -4.
---------------------------------------------------------------*/
int ARKMassKLUReInit(void *arkode_mem, int n, int nnz, int reinit_type)
{
  ARKodeMem ark_mem;
  ARKSlsMassMem arksls_mem;
  KLUData klu_data;

  /* Return immediately if arkode_mem is NULL. */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", "ARKMassKLUReInit", 
		    MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Return immediately if ark_lmem is NULL. */
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MASSMEM_NULL, "ARKSLS", "ARKMassKLUReInit", 
		    MSGSP_MASSMEM_NULL);
    return(ARKSLS_MASSMEM_NULL);
  }

  arksls_mem = (ARKSlsMassMem) (ark_mem->ark_mass_mem);
  klu_data = (KLUData) arksls_mem->s_solver_data;

  /* Return if reinit_type is not valid */
  if ((reinit_type != 1) && (reinit_type != 2)) {
    arkProcessError(NULL, ARKSLS_ILL_INPUT, "ARKSLS", "ARKMassKLUReInit", 
		    MSGSP_ILL_INPUT);
    return(ARKSLS_ILL_INPUT);
  }

  if (reinit_type == 1) {

    /* Destroy previous Jacobian information */
    if (arksls_mem->s_M) {
      DestroySparseMat(arksls_mem->s_M);
    }

    /* Allocate memory for the sparse Jacobian */
    arksls_mem->s_M = NewSparseMat(n, n, nnz);
    if (arksls_mem->s_M == NULL) {
      arkProcessError(ark_mem, ARKSLS_MEM_FAIL, "ARKSLS", "ARKMassKLU", 
		    MSGSP_MEM_FAIL);
      return(ARKSLS_MEM_FAIL);
    }
  }

  /* Free the prior factorazation and reset for first factorization */
  if( klu_data->s_Symbolic != NULL)
    klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
  if( klu_data->s_Numeric != NULL)
    klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
  arksls_mem->s_first_factorize = 1;

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;

  return(0);
}

/*---------------------------------------------------------------
 arkMassKLUInit:

 This routine does remaining initializations specific to the 
 ARKKLU mass matrix linear solver module.  It returns 0 if 
 successful.
---------------------------------------------------------------*/
static int arkMassKLUInit(ARKodeMem ark_mem)
{
  ARKSlsMassMem arksls_mem;

  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  arksls_mem->s_nme = 0;
  arksls_mem->s_first_factorize = 1;

  /* Set mass matrix function data */
  arksls_mem->s_Mdata = ark_mem->ark_user_data;
  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkMassKLUSetup:

  This routine does the setup operations for the ARKMassKLU 
  linear solver module.  It calls the mass matrix evaluation 
  routine, updates counters, and calls the LU factorization 
  routine.  The return value is either
     ARkSLS_SUCCESS = 0  if successful,
     +1  if the Meval routine failed recoverably or the
         LU factorization failed, or
     -1  if the Meval routine failed unrecoverably.
---------------------------------------------------------------*/
static int arkMassKLUSetup(ARKodeMem ark_mem, N_Vector vtemp1, 
			   N_Vector vtemp2, N_Vector vtemp3)
{
  ARKSlsMassMem arksls_mem;
  KLUData klu_data;
  int retval;
  
  realtype uround_twothirds;
  
  uround_twothirds = SUNRpowerR(ark_mem->ark_uround,TWOTHIRDS);

  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  klu_data = (KLUData) arksls_mem->s_solver_data;
  
  /* Check that mass matrix eval routine is set */
  if (arksls_mem->s_Meval == NULL) {
    arkProcessError(ark_mem, ARKSLS_MASS_NOSET, "ARKSLS", 
		    "arkMassKLUSetup", MSGSP_MASS_NOSET);
    free(arksls_mem); arksls_mem = NULL;
    return(ARKSLS_MASS_NOSET);
  }

  /* call Meval routine for new M matrix */
  SlsSetToZero(arksls_mem->s_M);
  retval = arksls_mem->s_Meval(ark_mem->ark_tn, arksls_mem->s_M, 
			       arksls_mem->s_Mdata, vtemp1, 
			       vtemp2, vtemp3);
  arksls_mem->s_nme++;
  if (retval < 0) {
    arkProcessError(ark_mem, ARKSLS_MASSFUNC_UNRECVR, "ARKSLS", 
		    "arkMassKLUSetup", MSGSP_MASSFUNC_FAILED);
    arksls_mem->s_last_flag = ARKSLS_MASSFUNC_UNRECVR;
    return(-1);
  }
  if (retval > 0) {
    arksls_mem->s_last_flag = ARKSLS_MASSFUNC_RECVR;
    return(1);
  }

  /* Copy M into M_lu for LU decomposition */
  CopySparseMat(arksls_mem->s_M, arksls_mem->s_M_lu);

  /* On first decomposition, get the symbolic factorization */ 
  if (arksls_mem->s_first_factorize) {

    /* Update the ordering option with user-updated values */
    klu_data->s_Common.ordering = klu_data->s_ordering;

    /* Perform symbolic analysis of sparsity structure */
    klu_data->s_Symbolic = klu_analyze(arksls_mem->s_M_lu->N, 
				       arksls_mem->s_M_lu->colptrs, 
				       arksls_mem->s_M_lu->rowvals, 
				       &(klu_data->s_Common));
    if (klu_data->s_Symbolic == NULL) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "ARKMassKLUSetup", MSGSP_PACKAGE_FAIL);
      return(ARKSLS_PACKAGE_FAIL);
    }

    /* ------------------------------------------------------------
       Compute the LU factorization of  the Jacobian.
       ------------------------------------------------------------*/
    klu_data->s_Numeric = klu_factor(arksls_mem->s_M_lu->colptrs, 
				     arksls_mem->s_M_lu->rowvals, 
				     arksls_mem->s_M_lu->data, 
				     klu_data->s_Symbolic, 
				     &(klu_data->s_Common));
    if (klu_data->s_Numeric == NULL) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "ARKMassKLUSetup", MSGSP_PACKAGE_FAIL);
      return(ARKSLS_PACKAGE_FAIL);
    }
    
    arksls_mem->s_first_factorize = 0;
  }
  else {

    retval = klu_refactor(arksls_mem->s_M_lu->colptrs, 
			  arksls_mem->s_M_lu->rowvals, 
			  arksls_mem->s_M_lu->data,
			  klu_data->s_Symbolic, klu_data->s_Numeric,
			  &(klu_data->s_Common));
    if (retval == 0) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "ARKMassKLUSetup", MSGSP_PACKAGE_FAIL);
      return(ARKSLS_PACKAGE_FAIL);
    }
    
    /*-----------------------------------------------------------
      Check if a cheap estimate of the reciprocal of the condition 
      number is getting too small.  If so, delete
      the prior numeric factorization and recompute it.
      -----------------------------------------------------------*/
    
    retval = klu_rcond(klu_data->s_Symbolic, klu_data->s_Numeric,
		       &(klu_data->s_Common));
    if (retval == 0) {
      arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		      "ARKMassKLUSetup", MSGSP_PACKAGE_FAIL);
      return(ARKSLS_PACKAGE_FAIL);
    }

    if ( (klu_data->s_Common.rcond)  < uround_twothirds ) {
      
      /* Condition number may be getting large.  
	 Compute more accurate estimate */
      retval = klu_condest(arksls_mem->s_M_lu->colptrs, 
			   arksls_mem->s_M_lu->data,
			   klu_data->s_Symbolic, klu_data->s_Numeric,
			   &(klu_data->s_Common));
      if (retval == 0) {
	arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
			"ARKMassKLUSetup", MSGSP_PACKAGE_FAIL);
	return(ARKSLS_PACKAGE_FAIL);
      }
      
      if ( (klu_data->s_Common.condest) > 
	   (1.0/uround_twothirds) ) {

	/* More accurate estimate also says condition number is 
	   large, so recompute the numeric factorization */

	klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
	
	klu_data->s_Numeric = klu_factor(arksls_mem->s_M_lu->colptrs, 
					 arksls_mem->s_M_lu->rowvals, 
					 arksls_mem->s_M_lu->data, 
					 klu_data->s_Symbolic, 
					 &(klu_data->s_Common));

	if (klu_data->s_Numeric == NULL) {
	  arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
			  "ARKMassKLUSetup", MSGSP_PACKAGE_FAIL);
	  return(ARKSLS_PACKAGE_FAIL);
	}
      }
    }
  }

  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkMassKLUSolve:

  This routine handles the solve operation for the ARKKLU mass 
  matrix linear solver module.  It calls the KLU solve routine, 
  then returns ARKSLS_SUCCESS = 0.
---------------------------------------------------------------*/
static int arkMassKLUSolve(ARKodeMem ark_mem, N_Vector b, 
			   N_Vector weight)
{
  int flag;
  ARKSlsMassMem arksls_mem;
  KLUData klu_data;
  realtype *bd;
  
  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  klu_data = (KLUData) arksls_mem->s_solver_data;

  bd = N_VGetArrayPointer(b);

  /* Call KLU to solve the linear system */
  flag = klu_solve(klu_data->s_Symbolic, klu_data->s_Numeric, 
		   arksls_mem->s_M_lu->N, 1, bd, 
		   &(klu_data->s_Common));
  if (flag == 0) {
    arkProcessError(ark_mem, ARKSLS_PACKAGE_FAIL, "ARKSLS", 
		    "ARKMassKLUSolve", MSGSP_PACKAGE_FAIL);
    return(ARKSLS_PACKAGE_FAIL);
  }
  arksls_mem->s_last_flag = ARKSLS_SUCCESS;
  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkMassKLUFree:

 This routine frees memory specific to the ARKMassKLU mass 
 matrix linear solver.
---------------------------------------------------------------*/
static void arkMassKLUFree(ARKodeMem ark_mem)
{
  ARKSlsMassMem arksls_mem;
  KLUData klu_data;
  
  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  klu_data = (KLUData) arksls_mem->s_solver_data;

  klu_free_numeric(&(klu_data->s_Numeric), 
		   &(klu_data->s_Common));

  klu_free_symbolic(&(klu_data->s_Symbolic), 
		    &(klu_data->s_Common));

  if (arksls_mem->s_M) {
    DestroySparseMat(arksls_mem->s_M);
    arksls_mem->s_M = NULL;
  }

  if (arksls_mem->s_M_lu) {
    DestroySparseMat(arksls_mem->s_M_lu);
    arksls_mem->s_M_lu = NULL;
  }

  free(klu_data); 
  free(arksls_mem); 
  ark_mem->ark_mass_mem = NULL;
}


/*===============================================================
 Utility Functions
===============================================================*/

/*---------------------------------------------------------------
 arkMassKLUMultiply:

 Multiplies the mass matrix by the vector v to fill in the 
 vector Mv.
---------------------------------------------------------------*/
static int arkMassKLUMultiply(N_Vector v, N_Vector Mv, 
			      realtype t, void *arkode_mem)
{

  /* extract the SlsMassMem structure from the arkode_mem pointer */
  ARKodeMem ark_mem;
  ARKSlsMassMem arksls_mem;
  realtype *vdata=NULL, *Mvdata=NULL;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "arkMassKLUMultiply", MSGSP_ARKMEM_NULL);
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
		    "arkMassKLLUMultiply", 
		    "Vector data un-allocated.");
    return(ARKSLS_MEM_NULL);
  }

  /* perform matrix-vector product with arksls_mem->s_M and return */
  if (SlsMatvec(arksls_mem->s_M, vdata, Mvdata) != 0) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "arkMassKLUMultiply", 
		    "Mass matrix data un-allocated.");
    return(ARKSLS_MEM_NULL);
  }

  return(0);

}


/*===============================================================
 Optional Input Specification Functions
===============================================================*/


/*---------------------------------------------------------------
 ARKKLUSetOrdering:

 Sets the ordering used by KLU for reducing fill.  Options are: 
     0 for AMD, 
     1 for COLAMD, and 
     2 for the natural ordering.  
 The default used is 1 (COLAMD).
---------------------------------------------------------------*/
int ARKKLUSetOrdering(void *arkode_mem, int ordering_choice)
{
  ARKodeMem ark_mem;
  ARKSlsMem arksls_mem;
  KLUData klu_data;

  /* Return immediately if ark_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKKLUSetOrdering", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Return if ordering choice argument is not valid */
  if ( (ordering_choice != 0) && (ordering_choice != 1) && 
       (ordering_choice != 2) ) {
    arkProcessError(NULL, ARKSLS_ILL_INPUT, "ARKSLS", 
		    "ARKKLUSetOrdering", MSGSP_ILL_INPUT);
    return(ARKSLS_ILL_INPUT);
  }

  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;
  klu_data = (KLUData) arksls_mem->s_solver_data;

  klu_data->s_ordering = ordering_choice;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKMassKLUSetOrdering:

 Sets the ordering used by KLU for reducing fill.  Options are: 
     0 for AMD, 
     1 for COLAMD, and 
     2 for the natural ordering.  
 The default used is 1 (COLAMD).
---------------------------------------------------------------*/
int ARKMassKLUSetOrdering(void *arkode_mem, int ordering_choice)
{
  ARKodeMem ark_mem;
  ARKSlsMassMem arksls_mem;
  KLUData klu_data;

  /* Return immediately if ark_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKKLUSetOrdering", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Return if ordering choice argument is not valid */
  if ( (ordering_choice != 0) && (ordering_choice != 1) && 
       (ordering_choice != 2) ) {
    arkProcessError(NULL, ARKSLS_ILL_INPUT, "ARKSLS", 
		    "ARKKLUSetOrdering", MSGSP_ILL_INPUT);
    return(ARKSLS_ILL_INPUT);
  }

  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;
  klu_data = (KLUData) arksls_mem->s_solver_data;

  klu_data->s_ordering = ordering_choice;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
   EOF
---------------------------------------------------------------*/
