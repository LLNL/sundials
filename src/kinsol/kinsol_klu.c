/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * begincopyright(llns)
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * endcopyright(llns)
 * -----------------------------------------------------------------
 * This is the implementation file for the KINKLU linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>

#include "kinsol_impl.h"
#include "kinsol_sparse_impl.h"
#include "sundials/sundials_klu_impl.h"

/* Constants */

#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* kinKLU linit, lsetup, lsolve, and lfree routines */
 
static int kinKLUInit(KINMem kin_mem);
static int kinKLUSetup(KINMem kin_mem);
static int kinKLUSolve(KINMem kin_mem, N_Vector x, N_Vector b,
		       realtype *sJpnorm, realtype *sFdotJp);		       
static void kinKLUFree(KINMem kin_mem);

/*
 * -----------------------------------------------------------------
 * kinKLU
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the KINSOL / KLU linear solver module.  
 * kinKLU first calls the existing lfree routine if this is not NULL.
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
 *       of the NVECTOR package. Therefore, kinKLU will first 
 *       test for a compatible N_Vector internal
 *       representation by checking that the functions N_VGetArrayPointer
 *       and N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int kinKLU(void *kin_mem_v, int n, int nnz)
{
  KINMem kin_mem;
  KINSlsMem kinsls_mem;
  KLUData klu_data;
  int flag;

  /* Return immediately if kin_mem is NULL. */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINSLS_MEM_NULL, "KINSLS", "kinKLU", 
		    MSGSP_KINMEM_NULL);
    return(KINSLS_MEM_NULL);
  }
  kin_mem = (KINMem) kin_mem_v;

  /* Test if the NVECTOR package is compatible with the Direct solver */
  if(kin_mem->kin_vtemp1->ops->nvgetarraypointer == NULL ||
     kin_mem->kin_vtemp1->ops->nvsetarraypointer == NULL) {
    KINProcessError(kin_mem, KINSLS_ILL_INPUT, "KINSLS", "kinKLU", 
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
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "kinKLU", 
		    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }

  /* Get memory for KLUData. */
  klu_data = (KLUData)malloc(sizeof(struct KLUDataRec));
  if (klu_data == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "kinKLU", 
		    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }

  kin_mem->kin_setupNonNull = TRUE;

  /* Set default Jacobian routine and Jacobian data */
  kinsls_mem->s_jaceval = NULL;
  kinsls_mem->s_jacdata = kin_mem->kin_user_data;

  /* Allocate memory for the sparse Jacobian */
  kinsls_mem->s_JacMat = NewSparseMat(n, n, nnz);
  if (kinsls_mem->s_JacMat == NULL) {
    KINProcessError(kin_mem, KINSLS_MEM_FAIL, "KINSLS", "kinKLU", 
		    MSGSP_MEM_FAIL);
    return(KINSLS_MEM_FAIL);
  }

  /* Allocate structures for KLU */

  /* DO I ALLOCATE COMMON????*/

  klu_data->s_Symbolic = (klu_symbolic *)malloc(sizeof(klu_symbolic));
  klu_data->s_Numeric = (klu_numeric *)malloc(sizeof(klu_numeric));

  /* Set default parameters for KLU */
  flag = klu_defaults(&klu_data->s_Common);
  if (flag == 0) {
    KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "kinKLU", 
		    MSGSP_PACKAGE_FAIL);
    return(KINSLS_PACKAGE_FAIL);
  }
  /* Set ordering as COLAMD.  The default is AMD, ordering 0. */
  /* This should get changed to a user option */
  klu_data->s_Common.ordering = 1;

  /* Attach linear solver memory to the nonlinear solver memory */
  kinsls_mem->s_solver_data = (void *) klu_data;
  kin_mem->kin_lmem = kinsls_mem;

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;

  return(KINSLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * kinKLU interface functions
 * -----------------------------------------------------------------
 */

/*
  This routine does remaining initializations specific to the kinKLU
  linear solver module.  
  It returns 0 if successful.
*/

static int kinKLUInit(KINMem kin_mem)
{
  int retval, n;
  KINSlsMem kinsls_mem;
  KLUData klu_data;

  kinsls_mem = (KINSlsMem)kin_mem->kin_lmem;
  klu_data = (KLUData)kinsls_mem->s_solver_data;

  kinsls_mem->s_nje = 0;
  kinsls_mem->s_first_factorize = 1;

  /* ------------------------------------------------------------
     Allocate storage and initialize statistics variables. 
     ------------------------------------------------------------*/
  n = kinsls_mem->s_JacMat->N;

  kinsls_mem->s_last_flag = 0;
  return(0);
}

/*
  This routine does the setup operations for the kinKLU linear 
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
  if (retval > 0) {
    kinsls_mem->s_last_flag = KINSLS_JACFUNC_RECVR;
    return(+1);
  }

  if (kinsls_mem->s_first_factorize) {
    /* ------------------------------------------------------------
       Get the symbolic factorization
       ------------------------------------------------------------*/ 
    klu_data->s_Symbolic = klu_analyze(JacMat->N, JacMat->colptrs, 
				       JacMat->rowvals, &(klu_data->s_Common));
    if (klu_data->s_Symbolic == NULL) {
      KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "kinKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(KINSLS_PACKAGE_FAIL);
    }

    kinsls_mem->s_first_factorize = 0;
  }

  /* ------------------------------------------------------------
     Compute the LU factorization of  the Jacobian.
     ------------------------------------------------------------*/
  klu_data->s_Numeric = klu_factor(JacMat->colptrs, JacMat->rowvals, 
				   JacMat->data, klu_data->s_Symbolic, 
				   &(klu_data->s_Common));

  if (klu_data->s_Numeric == NULL) {
    KINProcessError(kin_mem, KINSLS_PACKAGE_FAIL, "KINSLS", "kinKLUSetup", 
		    MSGSP_PACKAGE_FAIL);
    return(KINSLS_PACKAGE_FAIL);
  }

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;

  return(0);
}

/*
  This routine handles the solve operation for the kinKLU linear
  solver module.  It calls the KLU solve routine, 
  then returns KINSLS_SUCCESS = 0.
*/

static int kinKLUSolve(KINMem kin_mem, N_Vector x, N_Vector b,
		       realtype *sJpnorm, realtype *sFdotJp)		       
{
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
  klu_solve(klu_data->s_Symbolic, klu_data->s_Numeric, JacMat->N, 1, xd, 
	    &(klu_data->s_Common));

  /* Compute the term sFdotJp for use in the linesearch routine.
     This term is subsequently corrected if the step is reduced by
     constraints or the linesearch.

     sFdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale.                            */

  /* CSW: The liens below mimic those in kinsol_dense.c but seem wrong.
     Need to check with Alan and likely will need a matvec for matrices in 
     CSC format.  */

  N_VProd(b, kin_mem->kin_fscale, b);
  N_VProd(b, kin_mem->kin_fscale, b);
  *sFdotJp = N_VDotProd(kin_mem->kin_fval, b);

  kinsls_mem->s_last_flag = KINSLS_SUCCESS;
  return(KINSLS_SUCCESS);
}

/*
  This routine frees memory specific to the kinKLU linear solver.
*/

static void kinKLUFree(KINMem kin_mem)
{
  KINSlsMem kinsls_mem;
  KLUData klu_data;
  
  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;
  klu_data = (KLUData) kinsls_mem->s_solver_data;

  klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
  klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));

  free(klu_data->s_Symbolic);

  if (kinsls_mem->s_JacMat) {
    DestroySparseMat(kinsls_mem->s_JacMat);
    kinsls_mem->s_JacMat = NULL;
  }

  free(klu_data); 
  free(kin_mem->kin_lmem); 
}

