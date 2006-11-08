/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:29 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a KINSOL dense linear solver
 * using BLAS and LAPACK functions.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "kinsol_lapack_impl.h"
#include "kinsol_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* KINLAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int kinLapackDenseInit(KINMem kin_mem);
static int kinLapackDenseSetup(KINMem kin_mem);
static int kinLapackDenseSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm);
static void kinLapackDenseFree(KINMem kin_mem);

/* KINLAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int kinLapackBandInit(KINMem kin_mem);
static int kinLapackBandSetup(KINMem kin_mem);
static int kinLapackBandSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm);
static void kinLapackBandFree(KINMem kin_mem);

/* KINLAPACK DENSE and BAND DQ integration Jacobian functions */
static int kinLapackDenseDQJac(int N,
                              N_Vector u, N_Vector fu, 
                              LapackMat Jac, void *jac_data,
                              N_Vector tmp1, N_Vector tmp2);
static int kinLapackBandDQJac(int N, int mupper, int mlower,
                             N_Vector u, N_Vector fu,
                             LapackMat Jac, void *jac_data,
                             N_Vector tmp1, N_Vector tmp2);

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define lrw1           (kin_mem->kin_lrw1)
#define liw1           (kin_mem->kin_liw1)
#define uround         (kin_mem->kin_uround)
#define func           (kin_mem->kin_func)
#define f_data         (kin_mem->kin_f_data)
#define printfl        (kin_mem->kin_printfl)
#define linit          (kin_mem->kin_linit)
#define lsetup         (kin_mem->kin_lsetup)
#define lsolve         (kin_mem->kin_lsolve)
#define lfree          (kin_mem->kin_lfree)
#define lmem           (kin_mem->kin_lmem)
#define inexact_ls     (kin_mem->kin_inexact_ls)
#define uu             (kin_mem->kin_uu)
#define fval           (kin_mem->kin_fval)
#define uscale         (kin_mem->kin_uscale)
#define fscale         (kin_mem->kin_fscale)
#define sqrt_relfunc   (kin_mem->kin_sqrt_relfunc)
#define sJpnorm        (kin_mem->kin_sJpnorm)
#define sfdotJp        (kin_mem->kin_sfdotJp)
#define errfp          (kin_mem->kin_errfp)
#define infofp         (kin_mem->kin_infofp)
#define setupNonNull   (kin_mem->kin_setupNonNull)
#define vtemp1         (kin_mem->kin_vtemp1)
#define vec_tmpl       (kin_mem->kin_vtemp1)
#define vtemp2         (kin_mem->kin_vtemp2)

#define mtype          (kinlapack_mem->l_mtype)
#define n              (kinlapack_mem->l_n)
#define ml             (kinlapack_mem->b_ml)
#define mu             (kinlapack_mem->b_mu)
#define smu            (kinlapack_mem->b_smu)
#define djac           (kinlapack_mem->d_jac)
#define bjac           (kinlapack_mem->b_jac)
#define J              (kinlapack_mem->l_J)
#define pivots         (kinlapack_mem->l_pivots)
#define nje            (kinlapack_mem->l_nje)
#define nfeDQ          (kinlapack_mem->l_nfeDQ)
#define J_data         (kinlapack_mem->l_J_data)
#define last_flag      (kinlapack_mem->l_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * KINLapackDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the linear solver module.  KINLapackDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the kin_linit, kin_lsetup, kin_lsolve, kin_lfree fields in (*kinmem)
 * to be kinLapackDenseInit, kinLapackDenseSetup, kinLapackDenseSolve, 
 * and kinLapackDenseFree, respectively.  It allocates memory for a 
 * structure of type KINLapackMemRec and sets the kin_lmem field in 
 * (*kinmem) to the address of this structure.  It sets lsetup_exists 
 * in (*kinmem) to TRUE, and the d_jac field to the default 
 * kinLapackDenseDQJac. Finally, it allocates memory for M, pivots, and 
 * (if needed) savedJ.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINLapackDense will first 
 *       test for a compatible N_Vector internal representation 
 *       by checking that N_VGetArrayPointer and N_VSetArrayPointer 
 *       exist.
 * -----------------------------------------------------------------
 */
int KINLapackDense(void *kinmem, int N)
{
  KINMem kin_mem;
  KINLapackMem kinlapack_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINLAPACK_MEM_NULL, "KINLAPACK", "KINLapackDense", MSGLS_KINMEM_NULL);
    return(KINLAPACK_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    KINProcessError(kin_mem, KINLAPACK_ILL_INPUT, "KINLAPACK", "KINLapackDense", MSGLS_BAD_NVECTOR);
    return(KINLAPACK_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */
  linit  = kinLapackDenseInit;
  lsetup = kinLapackDenseSetup;
  lsolve = kinLapackDenseSolve;
  lfree  = kinLapackDenseFree;

  /* Get memory for KINLapackMemRec */
  kinlapack_mem = NULL;
  kinlapack_mem = (KINLapackMem) malloc(sizeof(KINLapackMemRec));
  if (kinlapack_mem == NULL) {
    KINProcessError(kin_mem, KINLAPACK_MEM_FAIL, "KINLAPACK", "KINLapackDense", MSGLS_MEM_FAIL);
    return(KINLAPACK_MEM_FAIL);
  }

  /* Set default Jacobian routine and Jacobian data */
  djac = kinLapackDenseDQJac;
  J_data = kin_mem;
  last_flag = KINLAPACK_SUCCESS;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for J and pivot array */
  
  J = NULL;
  J = LapackAllocDenseMat(N, N);
  if (J == NULL) {
    KINProcessError(kin_mem, KINLAPACK_MEM_FAIL, "KINLAPACK", "KINLapackDense", MSGLS_MEM_FAIL);
    free(kinlapack_mem); kinlapack_mem = NULL;
    return(KINLAPACK_MEM_FAIL);
  }

  pivots = NULL;
  pivots = LapackAllocIntArray(N);
  if (pivots == NULL) {
    KINProcessError(kin_mem, KINLAPACK_MEM_FAIL, "KINLAPACK", "KINLapackDense", MSGLS_MEM_FAIL);
    LapackFreeMat(J);
    free(kinlapack_mem); kinlapack_mem = NULL;
    return(KINLAPACK_MEM_FAIL);
  }

  /* This is a direct linear solver */
  inexact_ls = FALSE;

  /* Attach linear solver memory to integrator memory */
  lmem = kinlapack_mem;

  return(KINLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINLapackBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module. It first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * kin_linit, kin_lsetup, kin_lsolve, and kin_lfree fields in (*kinmem)
 * to be kinLapackBandInit, kinLapackBandSetup, kinLapackBandSolve, 
 * and kinLapackBandFree, respectively.  It allocates memory for a 
 * structure of type KINLapackBandMemRec and sets the kin_lmem field in 
 * (*kinmem) to the address of this structure.  It sets lsetup_exists 
 * in (*kinmem) to be TRUE, mu to be mupper, ml to be mlower, and 
 * the jacE and jacI field to NULL.
 * Finally, it allocates memory for M, pivots, and (if needed) savedJ.  
 * The KINLapackBand return value is KINLAPACK_SUCCESS = 0, 
 * KINLAPACK_MEM_FAIL = -1, or KINLAPACK_ILL_INPUT = -2.
 *
 * NOTE: The KINLAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINLapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int KINLapackBand(void *kinmem, int N, int mupper, int mlower)
{
  KINMem kin_mem;
  KINLapackMem kinlapack_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINLAPACK_MEM_NULL, "KINLAPACK", "KINLapackBand", MSGLS_KINMEM_NULL);
    return(KINLAPACK_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL) {
    KINProcessError(kin_mem, KINLAPACK_ILL_INPUT, "KINLAPACK", "KINLapackBand", MSGLS_BAD_NVECTOR);
    return(KINLAPACK_ILL_INPUT);
  }

  if (lfree != NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */  
  linit  = kinLapackBandInit;
  lsetup = kinLapackBandSetup;
  lsolve = kinLapackBandSolve;
  lfree  = kinLapackBandFree;
  
  /* Get memory for KINLapackMemRec */
  kinlapack_mem = NULL;
  kinlapack_mem = (KINLapackMem) malloc(sizeof(KINLapackMemRec));
  if (kinlapack_mem == NULL) {
    KINProcessError(kin_mem, KINLAPACK_MEM_FAIL, "KINLAPACK", "KINLapackBand", MSGLS_MEM_FAIL);
    return(KINLAPACK_MEM_FAIL);
  }
  
  /* Set default Jacobian routine and Jacobian data */
  bjac = kinLapackBandDQJac;
  J_data = kinmem;
  last_flag = KINLAPACK_SUCCESS;

  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in kinlapack_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    KINProcessError(kin_mem, KINLAPACK_MEM_FAIL, "KINLAPACK", "KINLapackBand", MSGLS_MEM_FAIL);
    return(KINLAPACK_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(N-1, mu + ml);

  /* Allocate memory for J and pivot array */
  J = NULL;
  J = LapackAllocBandMat(N, mu, ml, smu);
  if (J == NULL) {
    KINProcessError(kin_mem, KINLAPACK_MEM_FAIL, "KINLAPACK", "KINLapackBand", MSGLS_MEM_FAIL);
    free(kinlapack_mem); kinlapack_mem = NULL;
    return(KINLAPACK_MEM_FAIL);
  }

  pivots = NULL;
  pivots = LapackAllocIntArray(N);
  if (pivots == NULL) {
    KINProcessError(kin_mem, KINLAPACK_MEM_FAIL, "KINLAPACK", "KINLapackBand", MSGLS_MEM_FAIL);
    LapackFreeMat(J);
    free(kinlapack_mem); kinlapack_mem = NULL;
    return(KINLAPACK_MEM_FAIL);
  }

  /* This is a direct linear solver */
  inexact_ls = FALSE;

  /* Attach linear solver memory to integrator memory */
  lmem = kinlapack_mem;

  return(KINLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINLapackSetJacFn
 * -----------------------------------------------------------------
 */

int KINLapackSetJacFn(void *kinmem, void *jac, void *jac_data)
{
  KINMem kin_mem;
  KINLapackMem kinlapack_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINLAPACK_MEM_NULL, "KINLAPACK", "KINLapackSetJacFn", MSGLS_KINMEM_NULL);
    return(KINLAPACK_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINLAPACK_LMEM_NULL, "KINLAPACK", "KINLapackSetJacFn", MSGLS_LMEM_NULL);
    return(KINLAPACK_LMEM_NULL);
  }
  kinlapack_mem = (KINLapackMem) lmem;

  if (mtype == LAPACK_DENSE)
    djac = (KINLapackDenseJacFn) jac;
  else if (mtype == LAPACK_BAND)
    bjac = (KINLapackBandJacFn) jac;

  J_data = jac_data;

  return(KINLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINBandGetWorkSpace
 * -----------------------------------------------------------------
 */

int KINBandGetWorkSpace(void *kinmem, long int *lenrwLS, long int *leniwLS)
{
  KINMem kin_mem;
  KINLapackMem kinlapack_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINLAPACK_MEM_NULL, "KINLAPACK", "KINBandGetWorkSpace", MSGLS_KINMEM_NULL);
    return(KINLAPACK_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINLAPACK_LMEM_NULL, "KINLAPACK", "KINBandGetWorkSpace", MSGLS_LMEM_NULL);
    return(KINLAPACK_LMEM_NULL);
  }
  kinlapack_mem = (KINLapackMem) lmem;

  if (mtype == LAPACK_DENSE) {
    *lenrwLS = n*n;
    *leniwLS = n;
  } else if (mtype == LAPACK_BAND) {
    *lenrwLS = n*(smu + mu + 2*ml + 2);
    *leniwLS = n;
  }

  return(KINLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINLapackGetNumJacEvals
 * -----------------------------------------------------------------
 */

int KINLapackGetNumJacEvals(void *kinmem, long int *njevals)
{
  KINMem kin_mem;
  KINLapackMem kinlapack_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINLAPACK_MEM_NULL, "KINLAPACK", "KINLapackGetNumJacEvals", MSGLS_KINMEM_NULL);
    return(KINLAPACK_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINLAPACK_LMEM_NULL, "KINLAPACK", "KINLapackGetNumJacEvals", MSGLS_LMEM_NULL);
    return(KINLAPACK_LMEM_NULL);
  }
  kinlapack_mem = (KINLapackMem) lmem;

  *njevals = nje;

  return(KINLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINLapackGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int KINLapackGetNumFuncEvals(void *kinmem, long int *nfevalsLS)
{
  KINMem kin_mem;
  KINLapackMem kinlapack_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINLAPACK_MEM_NULL, "KINLAPACK", "KINLapackGetNumFuncEvals", MSGLS_KINMEM_NULL);
    return(KINLAPACK_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINLAPACK_LMEM_NULL, "KINLAPACK", "KINLapackGetNumGuncEvals", MSGLS_LMEM_NULL);
    return(KINLAPACK_LMEM_NULL);
  }
  kinlapack_mem = (KINLapackMem) lmem;

  *nfevalsLS = nfeDQ;

  return(KINLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINLapackGetLastFlag
 * -----------------------------------------------------------------
 */

int KINLapackGetLastFlag(void *kinmem, int *flag)
{
  KINMem kin_mem;
  KINLapackMem kinlapack_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINLAPACK_MEM_NULL, "KINLAPACK", "KINLapackGetLastFlag", MSGLS_KINMEM_NULL);
    return(KINLAPACK_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINLAPACK_LMEM_NULL, "KINLAPACK", "KINLapackGetLastFlag", MSGLS_LMEM_NULL);
    return(KINLAPACK_LMEM_NULL);
  }
  kinlapack_mem = (KINLapackMem) lmem;

  *flag = last_flag;

  return(KINLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINLapackGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *KINLapackGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case KINLAPACK_SUCCESS:
    sprintf(name, "KINLAPACK_SUCCESS");
    break;
  case KINLAPACK_MEM_NULL:
    sprintf(name, "KINLAPACK_MEM_NULL");
    break;
  case KINLAPACK_LMEM_NULL:
    sprintf(name, "KINLAPACK_LMEM_NULL");
    break;
  case KINLAPACK_ILL_INPUT:
    sprintf(name, "KINLAPACK_ILL_INPUT");
    break;
  case KINLAPACK_MEM_FAIL:
    sprintf(name, "KINLAPACK_MEM_FAIL");
    break;
  default:
    sprintf(name, "NONE");
  }

  return(name);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR SOLUTION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * kinLapackDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int kinLapackDenseInit(KINMem kin_mem)
{
  KINLapackMem kinlapack_mem;

  kinlapack_mem = (KINLapackMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  
  if (djac == NULL) {
    djac = kinLapackDenseDQJac;
    J_data = kin_mem;
  }

  last_flag = KINLAPACK_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It calls the dense LU factorization routine.
 * -----------------------------------------------------------------
 */

static int kinLapackDenseSetup(KINMem kin_mem)
{
  KINLapackMem kinlapack_mem;
  int ier, retval;

  kinlapack_mem = (KINLapackMem) lmem;
 
  nje++;
  retval = djac(n, uu, fval, J, J_data, vtemp1, vtemp2);
  if (retval != 0) {
    last_flag = -1;
    return(-1);
  }

  /* Do LU factorization of J */
  dgetrf_f77(&n, &n, J->data, &(J->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return -1 */
  last_flag = ier;
  if (ier > 0) return(-1);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int kinLapackDenseSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINLapackMem kinlapack_mem;
  realtype *xd;
  int ier, one = 1;

  kinlapack_mem = (KINLapackMem) lmem;

  /* Copy the right-hand side into x */
  N_VScale(ONE, b, x);
  xd = N_VGetArrayPointer(x);

  /* Back-solve and get solution in x */
  dgetrs_f77("N", &n, &one, J->data, &(J->ldim), pivots, xd, &n, &ier, 1); 
  if (ier > 0) return(-1);

  /* Compute the terms Jpnorm and sfdotJp for use in the global strategy
   * routines and in KINForcingTerm. Both of these terms are subsequently
   * corrected if the step is reduced by constraints or the line search.
   *
   * sJpnorm is the norm of the scaled product (scaled by fscale) of
   * the current Jacobian matrix J and the step vector p.
   *
   * sfdotJp is the dot product of the scaled f vector and the scaled
   * vector J*p, where the scaling uses fscale. 
   */
  sJpnorm = N_VWL2Norm(b,fscale);
  N_VProd(b, fscale, b);
  N_VProd(b, fscale, b);
  sfdotJp = N_VDotProd(fval, b);

  last_flag = KINLAPACK_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void kinLapackDenseFree(KINMem kin_mem)
{
  KINLapackMem  kinlapack_mem;

  kinlapack_mem = (KINLapackMem) lmem;
  
  LapackFreeMat(J);
  LapackFreeArray(pivots);
  free(kinlapack_mem); kinlapack_mem = NULL;
}


/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR SOLUTION WITH BANDED JACOBIANS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * kinLapackBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int kinLapackBandInit(KINMem kin_mem)
{
  KINLapackMem kinlapack_mem;

  kinlapack_mem = (KINLapackMem) lmem;

  nje   = 0;
  nfeDQ = 0;

  if (bjac == NULL) {
    bjac = kinLapackBandDQJac;
    J_data = kin_mem;
  }

  last_flag = KINLAPACK_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the band LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int kinLapackBandSetup(KINMem kin_mem)
{
  KINLapackMem kinlapack_mem;
  int ier, retval;

  kinlapack_mem = (KINLapackMem) lmem;

  nje++;
  retval = jac(n, mu, ml, uu, fval, J, J_data, vtemp1, vtemp2);
  if (retval != 0) {
    last_flag = -1;
    return(-1);
  }
  
  /* Do LU factorization of J */
  dgbtrf_f77(&n, &n, &ml, &mu, J->data, &(J->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return -1 */
  last_flag = ier;
  if (ier > 0) return(-1);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackBandSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int kinLapackBandSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINLapackMem kinlapack_mem;
  realtype *xd;
  int ier, one = 1;

  kinlapack_mem = (KINLapackMem) lmem;

  /* Copy the right-hand side into x */
  N_VScale(ONE, b, x);
  xd = N_VGetArrayPointer(x);

  /* Back-solve and get solution in x */
  dgbtrs_f77("N", &n, &ml, &mu, &one, J->data, &(J->ldim), pivots, xd, &n, &ier, 1);
  if (ier > 0) return(-1);

  /* Compute the terms Jpnorm and sfdotJp for use in the global strategy
   * routines and in KINForcingTerm. Both of these terms are subsequently
   * corrected if the step is reduced by constraints or the line search.
   * 
   * sJpnorm is the norm of the scaled product (scaled by fscale) of
   * the current Jacobian matrix J and the step vector p.
   *
   * sfdotJp is the dot product of the scaled f vector and the scaled
   * vector J*p, where the scaling uses fscale. 
   */
  sJpnorm = N_VWL2Norm(b,fscale);
  N_VProd(b, fscale, b);
  N_VProd(b, fscale, b);
  sfdotJp = N_VDotProd(fval, b);

  last_flag = KINLAPACK_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static void kinLapackBandFree(KINMem kin_mem)
{
  KINLapackMem kinlapack_mem;

  kinlapack_mem = (KINLapackMem) lmem;

  LapackFreeMat(J);
  LapackFreeArray(pivots);
  free(kinlapack_mem); kinlapack_mem = NULL;
}




#undef n
#undef J




/*
 * -----------------------------------------------------------------
 * kinLapackDenseDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of F(u). It assumes that a dense matrix of type
 * LapackMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 *
 * The increment used in the finitie-difference approximation
 *   J_ij = ( F_i(u+sigma_j * e_j) - F_i(u)  ) / sigma_j
 * is
 *  sigma_j = max{|u_j|, |1/uscale_j|} * sqrt(uround)
 *
 * Note: uscale_j = 1/typ(u_j)
 *
 * NOTE: Any type of failure of the system function her leads to an
 *       unrecoverable failure of the Jacobian function and thus
 *       of the linear solver setup function, stopping KINSOL.
 * -----------------------------------------------------------------
 */

static int kinLapackDenseDQJac(int N,
                               N_Vector u, N_Vector fu,
                               LapackMat Jac, void *jac_data,
                               N_Vector tmp1, N_Vector tmp2)
{
  realtype inc, inc_inv, ujsaved, ujscale, sign;
  realtype *tmp2_data, *u_data, *uscale_data;
  N_Vector ftemp, jthCol;
  long int j;
  int retval;

  KINMem kin_mem;
  KINLapackMem  kinlapack_mem;

  /* jac_data points to kin_mem */
  kin_mem = (KINMem) jac_data;
  kinlapack_mem = (KINLapackMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for u and uscale */
  u_data   = N_VGetArrayPointer(u);
  uscale_data = N_VGetArrayPointer(uscale);

  /* This is the only for loop for 0..N-1 in KINSOL */

  for (j = 0; j < N; j++) {

    /* Generate the jth col of Jac(u) */

    N_VSetArrayPointer(LAPACK_DENSE_COL(Jac,j), jthCol);

    ujsaved = u_data[j];
    ujscale = ONE/uscale_data[j];
    sign = (ujsaved >= ZERO) ? ONE : -ONE;
    inc = sqrt_relfunc*MAX(ABS(ujsaved), ujscale)*sign;
    u_data[j] += inc;

    retval = func(u, ftemp, f_data);
    if (retval != 0) return(-1); 

    u_data[j] = ujsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fu, jthCol);

  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  /* Increment counter nfeDQ */
  nfeDQ += N;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackBandDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of F(u).  It assumes that a band matrix of type
 * BandMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 *
 * NOTE: Any type of failure of the system function her leads to an
 *       unrecoverable failure of the Jacobian function and thus
 *       of the linear solver setup function, stopping KINSOL.
 * -----------------------------------------------------------------
 */

static int kinLapackBandDQJac(int N, int mupper, int mlower,
                              N_Vector u, N_Vector fu,
                              LapackMat Jac, void *jac_data,
                              N_Vector tmp1, N_Vector tmp2)
{
  realtype inc, inc_inv;
  N_Vector futemp, utemp;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *fu_data, *futemp_data, *u_data, *utemp_data, *uscale_data;
  int retval;

  KINMem kin_mem;
  KINLapackMem kinlapack_mem;

  /* jac_dat points to kinmem */
  kin_mem = (KINMem) jac_data;
  kinlapack_mem = (KINLapackMem) lmem;

  /* Rename work vectors for use as temporary values of u and fu */
  futemp = tmp1;
  utemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, futemp, y, ytemp */
  fu_data    = N_VGetArrayPointer(fu);
  futemp_data = N_VGetArrayPointer(futemp);
  u_data     = N_VGetArrayPointer(u);
  uscale_data = N_VGetArrayPointer(uscale);
  utemp_data = N_VGetArrayPointer(utemp);

  /* Load utemp with u */
  N_VScale(ONE, u, utemp);

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);
  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all utemp components in group */
    for(j=group-1; j < N; j+=width) {
      inc = sqrt_relfunc*MAX(ABS(u_data[j]), ABS(uscale_data[j]));
      utemp_data[j] += inc;
    }

    /* Evaluate f with incremented u */
    retval = func(utemp, futemp, f_data);
    if (retval != 0) return(-1); 

    /* Restore utemp components, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      utemp_data[j] = u_data[j];
      col_j = LAPACK_BAND_COL(Jac,j);
      inc = sqrt_relfunc*MAX(ABS(u_data[j]), ABS(uscale_data[j]));
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        LAPACK_BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (futemp_data[i] - fu_data[i]);
    }
  }
  
  /* Increment counter nfeDQ */
  nfeDQ += ngroups;

  return(0);
}
