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
 * This is the implementation file for the KINKLU linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>

#include "kinsol/kinsol_klu.h"
#include "kinsol_impl.h"
#include "kinsol_sparse_impl.h"
#include "sundials/sundials_klu_impl.h"

/* Constants */

#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)
#define TWOTHIRDS    RCONST(0.6666666666666667)

/* KINKLU linit, lsetup, lsolve, and lfree routines */
 
static int kinKLUInit(KINMem kin_mem);
static int kinKLUSetup(KINMem kin_mem);
static int kinKLUSolve(KINMem kin_mem, N_Vector x, N_Vector b,
		       realtype *sJpnorm, realtype *sFdotJp);		       
static int kinKLUFree(KINMem kin_mem);

/*
 * -----------------------------------------------------------------
 * KINKLU
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the KINSOL / KLU linear solver module.  
 * KINKLU first calls the existing lfree routine if this is not NULL.
 * Then it sets the kin_linit, kin_lsetup, kin_lsolve, kin_lperf, and
 * kin_lfree fields in (*kin_mem) to be kinKLUInit, kinKLUSetup,
 * kinKLUSolve, NULL, and kinKLUFree, respectively.
 * It allocates memory for a structure of type kinkluMemRec and sets
 * the kin_lmem field in (*kin_mem) to the address of this structure.
 * It sets setupNonNull in (*kin_mem) to TRUE.
 * Finally, it allocates memory for KLU.
 * The return value is KINSLS_SUCCESS = 0, KINSLS_LMEM_FAIL = -1,
 * or KINSLS_ILL_INPUT = -2.
 *
 * NOTE: The KLU linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINKLU will first 
 *       test for a compatible N_Vector internal representation
 *       by checking that the function N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */

int KINKLU(void *kin_mem_v, int n, int nnz, int sparsetype)
{
  KINMem kin_mem;
  KINSlsMem kinsls_mem;
  KLUData klu_data;
  int flag;

  /* Return immediately if kin_mem is NULL. */
  if (kin_mem_v == NULL) {
    KINProcessError(NULL, KINSLS_MEM_NULL, "KINSLS", "KINKLU", 
                    MSGSP_KINMEM_NULL);
    return(KINSLS_MEM_NULL);
  }
  kin_mem = (KINMem) kin_mem_v;

  /* Test if the NVECTOR package is compatible with the Direct solver */
  if (kin_mem->kin_vtemp1->ops->nvgetarraypointer == NULL) {
    KINProcessError(kin_mem, KINSLS_ILL_INPUT, "KINSLS", "KINKLU", 
                    MSGSP_BAD_NVECTOR);
    return(KINSLS_ILL_INPUT);
  }

  if (kin_mem->kin_lfree != NULL) kin_mem->kin_lfree(kin_mem);

  /* Set five main function fields in kin_mem. */
  kin_mem->kin_linit  = kinKLUInit;
  kin_mem->kin_lsetup = kinKLUSetup;
  kin_mem->kin_lsolve = kinKLUSolve;
  kin_mem->kin_lfree  = kinKLUFree;

  /* Get memory for kinSlsMemRec. */
  kinsls_mem = (KINSlsMem) malloc(sizeof(struct KINSlsMemRec));
  if (kinsls_mem == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "KINKLU", 
                    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }

  /* Get memory for KLUData. */
  klu_data = (KLUData)malloc(sizeof(struct KLUDataRec));
  if (klu_data == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "KINKLU", 
                    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }

  kin_mem->kin_setupNonNull = TRUE;

  /* Set default Jacobian routine and Jacobian data */
  kinsls_mem->s_jaceval = NULL;
  kinsls_mem->s_jacdata = NULL;

  /* Allocate memory for the sparse Jacobian */
  kinsls_mem->s_JacMat = SparseNewMat(n, n, nnz, sparsetype);
  if (kinsls_mem->s_JacMat == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "KINKLU", 
                    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }

  /* Initialize KLU structures */
  switch (sparsetype) {
    case CSC_MAT:
      klu_data->sun_klu_solve = &klu_solve;
      break;
    case CSR_MAT:
      klu_data->sun_klu_solve = &klu_tsolve;
      break;
    default:
      SparseDestroyMat(kinsls_mem->s_JacMat);
      free(klu_data);
      free(kinsls_mem);
      return(KINSLS_ILL_INPUT);
  }
  klu_data->s_Symbolic = NULL;
  klu_data->s_Numeric = NULL;

  /* Set default parameters for KLU */
  flag = klu_defaults(&klu_data->s_Common);
  if (flag == 0) {
    KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "KINKLU", 
                    MSGSP_PACKAGE_FAIL);
    return(KINSLS_PACKAGE_FAIL);
  }
  /* Set ordering to COLAMD as the kinsol default use.
     Users can set a different value with KINKLUSetOrdering,
     and the user-set value is loaded before any call to klu_analyze in
     kinKLUSetup.  */
  klu_data->s_ordering = 1;
  klu_data->s_Common.ordering = klu_data->s_ordering;

  /* This is a direct linear solver */
  kin_mem->kin_inexact_ls = FALSE;

  /* Attach linear solver memory to the nonlinear solver memory */
  kinsls_mem->s_solver_data = (void *) klu_data;
  kin_mem->kin_lmem = kinsls_mem;

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;

  return(KINSLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINKLUReInit
 * -----------------------------------------------------------------
 * This routine reinitializes memory and flags for a new factorization 
 * (symbolic and numeric) to be conducted at the next solver setup
 * call.  This routine is useful in the cases where the number of nonzeroes 
 * has changed or if the structure of the linear system has changed
 * which would require a new symbolic (and numeric factorization).
 *
 * The reinit_type argumenmt governs the level of reinitialization:
 *
 * reinit_type = 1: The Jacobian matrix will be destroyed and 
 *                  a new one will be allocated based on the nnz
 *                  value passed to this call. New symbolic and
 *                  numeric factorizations will be completed at the next
 *                  solver setup.
 *
 * reinit_type = 2: Only symbolic and numeric factorizations will be 
 *                  completed.  It is assumed that the Jacobian size
 *                  has not exceeded the size of nnz given in the prior
 *                  call to KINKLU.
 *
 * This routine assumes no other changes to solver use are necessary.
 *
 * The return value is KINSLS_SUCCESS = 0, KINSLS_MEM_NULL = -1, 
 * KINSLS_LMEM_NULL = -2, KINSLS_ILL_INPUT = -3, or KINSLS_MEM_FAIL = -4.
 *
 * -----------------------------------------------------------------
 */

int KINKLUReInit(void *kin_mem_v, int n, int nnz, int reinit_type)
{
  KINMem kin_mem;
  KINSlsMem kinsls_mem;
  KLUData klu_data;

  /* Return immediately if kin_mem is NULL. */
  if (kin_mem_v == NULL) {
    KINProcessError(NULL, KINSLS_MEM_NULL, "KINSLS", "KINKLUReInit", 
                    MSGSP_KINMEM_NULL);
    return(KINSLS_MEM_NULL);
  }
  kin_mem = (KINMem) kin_mem_v;

  /* Return immediately if kin_lmem is NULL. */
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(NULL, KINSLS_LMEM_NULL, "KINSLS", "KINKLUReInit", 
                    MSGSP_LMEM_NULL);
    return(KINSLS_LMEM_NULL);
  }
  kinsls_mem = (KINSlsMem) (kin_mem->kin_lmem);
  klu_data = (KLUData) kinsls_mem->s_solver_data;

  /* Return if reinit_type is not valid */
  if ((reinit_type != 1) && (reinit_type != 2)) {
    KINProcessError(NULL, KINSLS_ILL_INPUT, "KINSLS", "KINKLUReInit", 
                    MSGSP_ILL_INPUT);
    return(KINSLS_ILL_INPUT);
  }

  if (reinit_type == 1) {

    /* Destroy previous Jacobian information */
    if (kinsls_mem->s_JacMat) {
      SparseDestroyMat(kinsls_mem->s_JacMat);
    }

    /* Allocate memory for the sparse Jacobian */
    kinsls_mem->s_JacMat = SparseNewMat(n, n, nnz, kinsls_mem->sparsetype);
    if (kinsls_mem->s_JacMat == NULL) {
      KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "KINKLU", 
                      MSGSP_MEM_FAIL);
      return(KINSLS_MEM_FAIL);
    }
  }

  /* Free the prior factorazation and reset for first factorization */
  if( klu_data->s_Symbolic != NULL)
    klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
  if( klu_data->s_Numeric != NULL)
    klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
  kinsls_mem->s_first_factorize = 1;

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * KINKLU interface functions
 * -----------------------------------------------------------------
 */

/*
  This routine does remaining initializations specific to the KINKLU
  linear solver module.  
  It returns 0 if successful.
*/

static int kinKLUInit(KINMem kin_mem)
{
  KINSlsMem kinsls_mem;

  kinsls_mem = (KINSlsMem)kin_mem->kin_lmem;

  kinsls_mem->s_jacdata = kin_mem->kin_user_data;

  kinsls_mem->s_nje = 0;
  /* This forces factorization for every call to KINSol */
  kinsls_mem->s_first_factorize = 1;

  kinsls_mem->s_last_flag = 0;
  return(0);
}

/*
  This routine does the setup operations for the KINKLU linear 
  solver module.  It calls the Jacobian evaluation routine,
  updates counters, and calls the LU factorization routine.
  The return value is either
     KINSLS_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
*/

static int kinKLUSetup(KINMem kin_mem)
{
  int retval;
  KINSlsMem kinsls_mem;
  KINSlsSparseJacFn jaceval;
  KLUData klu_data;
  SlsMat JacMat;
  void *jacdata;
  realtype uround_twothirds;

  uround_twothirds = SUNRpowerR(kin_mem->kin_uround,TWOTHIRDS);

  kinsls_mem = (KINSlsMem) (kin_mem->kin_lmem);

  klu_data = (KLUData) kinsls_mem->s_solver_data;

  jaceval = kinsls_mem->s_jaceval;
  jacdata = kinsls_mem->s_jacdata;
  JacMat = kinsls_mem->s_JacMat;

  /* Check that Jacobian eval routine is set */
  if (jaceval == NULL) {
    KINProcessError(kin_mem, KINSLS_JAC_NOSET, "KINSLS", "kinKLUSetup", 
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
		    "kinKLUSetup", MSGSP_JACFUNC_FAILED);
    kinsls_mem->s_last_flag = KINSLS_JACFUNC_UNRECVR;
    return(KINSLS_JACFUNC_UNRECVR);
  }
  if (retval == 1) {
    kinsls_mem->s_first_factorize = 1;
  }
  if (retval > 1) {
    kinsls_mem->s_last_flag = KINSLS_JACFUNC_RECVR;
    return(+1);
  }

  if (kinsls_mem->s_first_factorize) {
    /* ------------------------------------------------------------
       Get the symbolic factorization
       ------------------------------------------------------------*/ 
    /* Update the ordering option with any user-updated values from 
       calls to KINKLUSetOrdering */
    klu_data->s_Common.ordering = klu_data->s_ordering;

    if (klu_data->s_Symbolic != NULL) {
       klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
    }
    klu_data->s_Symbolic = klu_analyze(JacMat->NP, JacMat->indexptrs, 
				       JacMat->indexvals, &(klu_data->s_Common));
    if (klu_data->s_Symbolic == NULL) {
      KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "kinKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(KINSLS_PACKAGE_FAIL);
    }

    /* ------------------------------------------------------------
       Compute the LU factorization of the Jacobian.
       ------------------------------------------------------------*/
    /* If klu_factor previously called, free data */
    if( klu_data->s_Numeric != NULL) {
       klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
    }
    klu_data->s_Numeric = klu_factor(JacMat->indexptrs, JacMat->indexvals, 
				     JacMat->data, klu_data->s_Symbolic, 
				     &(klu_data->s_Common));

    if (klu_data->s_Numeric == NULL) {
      KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "kinKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(KINSLS_PACKAGE_FAIL);
    }

    kinsls_mem->s_first_factorize = 0;

  }
  else {
    
    retval = klu_refactor(JacMat->indexptrs, JacMat->indexvals, JacMat->data, 
			klu_data->s_Symbolic, klu_data->s_Numeric,
			&(klu_data->s_Common));
    if (retval == 0) {
      KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "kinKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(KINSLS_PACKAGE_FAIL);
    }

    /*-----------------------------------------------------------
      Check if a cheap estimate of the reciprocal of the condition 
      number is getting too small.  If so, delete
      the prior numeric factorization and recompute it.
      -----------------------------------------------------------*/
    
    retval = klu_rcond(klu_data->s_Symbolic, klu_data->s_Numeric,
		       &(klu_data->s_Common));
    if (retval == 0) {
      KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "kinKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(KINSLS_PACKAGE_FAIL);
    }

# if 1
    if ( (klu_data->s_Common.rcond)  < uround_twothirds ) {
      
      /* Condition number may be getting large.  
	 Compute more accurate estimate */
      retval = klu_condest(JacMat->indexptrs, JacMat->data, 
			   klu_data->s_Symbolic, klu_data->s_Numeric,
			   &(klu_data->s_Common));
      if (retval == 0) {
	KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "kinKLUSetup", 
			MSGSP_PACKAGE_FAIL);
	return(KINSLS_PACKAGE_FAIL);
      }
      
      if ( (klu_data->s_Common.condest) > 
	   (1.0/uround_twothirds) ) {

	/* More accurate estimate also says condition number is 
	   large, so recompute the numeric factorization */

	klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
	
	klu_data->s_Numeric = klu_factor(JacMat->indexptrs, JacMat->indexvals, 
					 JacMat->data, klu_data->s_Symbolic, 
					 &(klu_data->s_Common));

	if (klu_data->s_Numeric == NULL) {
	  KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", 
			  "kinKLUSetup", MSGSP_PACKAGE_FAIL);
	  return(KINSLS_PACKAGE_FAIL);
	}
      }
    }
#endif 

  }

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;

  return(0);
}

/*
  This routine handles the solve operation for the KINKLU linear
  solver module.  It calls the KLU solve routine, 
  then returns KINSLS_SUCCESS = 0.
*/

static int kinKLUSolve(KINMem kin_mem, N_Vector x, N_Vector b,
		       realtype *sJpnorm, realtype *sFdotJp)		       
{
  int flag;
  KINSlsMem kinsls_mem;
  KLUData klu_data;
  SlsMat JacMat;
  realtype *xd;
  
  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;
  JacMat = kinsls_mem->s_JacMat;

  klu_data = (KLUData) kinsls_mem->s_solver_data;

  /* Copy the right-hand side into x */
  N_VScale(ONE, b, x);
  xd = N_VGetArrayPointer(x);

  /* Call KLU to solve the linear system */
  flag = klu_data->sun_klu_solve(klu_data->s_Symbolic, klu_data->s_Numeric, JacMat->NP, 1, xd, 
                                 &(klu_data->s_Common));
  if (flag == 0) {
    KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "kinKLUSolve", 
		    MSGSP_PACKAGE_FAIL);
    return(KINSLS_PACKAGE_FAIL);
  }

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
  This routine frees memory specific to the KINKLU linear solver.
*/

static int kinKLUFree(KINMem kin_mem)
{
  KINSlsMem kinsls_mem;
  KLUData klu_data;
  
  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;
  klu_data = (KLUData) kinsls_mem->s_solver_data;

  if( klu_data->s_Numeric != NULL)
  {
     klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
  }
  if( klu_data->s_Symbolic != NULL)
  {
     klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
  }

  if (kinsls_mem->s_JacMat) {
    SparseDestroyMat(kinsls_mem->s_JacMat);
    kinsls_mem->s_JacMat = NULL;
  }

  free(klu_data); 
  free(kin_mem->kin_lmem); 

  return(0);
}


/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * KINKLUSetOrdering sets the ordering used by KLU for reducing fill.
 * Options are: 0 for AMD, 1 for COLAMD, and 2 for the natural ordering.
 * The default used in KINSOL is 1 for COLAMD.
 * -----------------------------------------------------------------
 */

int KINKLUSetOrdering(void *kin_mem_v, int ordering_choice)
{
  KINMem kin_mem;
  KINSlsMem kinsls_mem;
  KLUData klu_data;

 /* Return immediately if kin_mem is NULL */
  if (kin_mem_v == NULL) {
    KINProcessError(NULL, KINSLS_MEM_NULL, "KINSLS", "KINKLUSetOrdering",
		    MSGSP_KINMEM_NULL);
    return(KINSLS_MEM_NULL);
  }
  kin_mem = (KINMem) kin_mem_v;

 /* Return if ordering choice argument is not valid */
  if ( (ordering_choice != 0) && (ordering_choice != 1) && 
       (ordering_choice != 2) ) {
    KINProcessError(NULL, KINSLS_ILL_INPUT, "KINSLS", "KINKLUSetOrdering",
		    MSGSP_ILL_INPUT);
    return(KINSLS_ILL_INPUT);
  }

  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;
  klu_data = (KLUData) kinsls_mem->s_solver_data;

  klu_data->s_ordering = ordering_choice;

  return(KINSLS_SUCCESS);
}

