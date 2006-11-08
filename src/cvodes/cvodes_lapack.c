/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:19 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a CVODES dense linear solver
 * using BLAS and LAPACK functions.
 * -----------------------------------------------------------------
 */

/*
 * NOTE: the only operation that does not use Blas/Lapack functions
 *       is matrix plus identity (in calculating I-gamma*J in lsetup)
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_lapack_impl.h"
#include "cvodes_impl.h"
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

/* CVLAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int cvLapackDenseInit(CVodeMem cv_mem);
static int cvLapackDenseSetup(CVodeMem cv_mem, int convfail, 
                              N_Vector yP, N_Vector fctP, 
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cvLapackDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector fctC);
static void cvLapackDenseFree(CVodeMem cv_mem);
static void cvLapackDenseAddI(LapackMat A);

/* CVLAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int cvLapackBandInit(CVodeMem cv_mem);
static int cvLapackBandSetup(CVodeMem cv_mem, int convfail, 
                             N_Vector yP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cvLapackBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector fctC);
static void cvLapackBandFree(CVodeMem cv_mem);
static void cvLapackBandAddI(LapackMat A);


/* CVLAPACK DENSE and BAND DQ integration Jacobian functions */
static int cvLapackDenseDQJac(int N, realtype t,
                              N_Vector y, N_Vector fy, 
                              LapackMat Jac, void *jac_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int cvLapackBandDQJac(int N, int mupper, int mlower,
                             realtype t, N_Vector y, N_Vector fy, 
                             LapackMat Jac, void *jac_data,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define lmm            (cv_mem->cv_lmm)
#define f              (cv_mem->cv_f)
#define f_data         (cv_mem->cv_f_data)
#define uround         (cv_mem->cv_uround)
#define nst            (cv_mem->cv_nst)
#define tn             (cv_mem->cv_tn)
#define h              (cv_mem->cv_h)
#define gamma          (cv_mem->cv_gamma)
#define gammap         (cv_mem->cv_gammap)
#define gamrat         (cv_mem->cv_gamrat)
#define ewt            (cv_mem->cv_ewt)

#define linit          (cv_mem->cv_linit)
#define lsetup         (cv_mem->cv_lsetup)
#define lsolve         (cv_mem->cv_lsolve)
#define lfree          (cv_mem->cv_lfree)
#define lmem           (cv_mem->cv_lmem)
#define tempv          (cv_mem->cv_tempv)
#define setupNonNull   (cv_mem->cv_setupNonNull)

#define mtype          (cvlapack_mem->l_mtype)
#define n              (cvlapack_mem->l_n)
#define ml             (cvlapack_mem->b_ml)
#define mu             (cvlapack_mem->b_mu)
#define smu            (cvlapack_mem->b_smu)
#define djac           (cvlapack_mem->d_jac)
#define bjac           (cvlapack_mem->b_jac)
#define M              (cvlapack_mem->l_M)
#define savedJ         (cvlapack_mem->l_savedJ)
#define pivots         (cvlapack_mem->l_pivots)
#define nstlj          (cvlapack_mem->l_nstlj)
#define nje            (cvlapack_mem->l_nje)
#define nfeDQ          (cvlapack_mem->l_nfeDQ)
#define J_data         (cvlapack_mem->l_J_data)
#define last_flag      (cvlapack_mem->l_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * CVLapackDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the linear solver module.  CVLapackDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be cvLapackDenseInit, cvLapackDenseSetup, cvLapackDenseSolve, 
 * and cvLapackDenseFree, respectively.  It allocates memory for a 
 * structure of type CVLapackMemRec and sets the cv_lmem field in 
 * (*cvode_mem) to the address of this structure.  It sets setupNonNull 
 * in (*cvode_mem) to TRUE, and the d_jac field to the default 
 * cvLapackDenseDQJac. Finally, it allocates memory for M, pivots, and 
 * (if needed) savedJ.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVLapackDense will first 
 *       test for a compatible N_Vector internal representation 
 *       by checking that N_VGetArrayPointer and N_VSetArrayPointer 
 *       exist.
 * -----------------------------------------------------------------
 */
int CVLapackDense(void *cvode_mem, int N)
{
  CVodeMem cv_mem;
  CVLapackMem cvlapack_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVLAPACK_MEM_NULL, "CVLAPACK", "CVLapackDense", MSGLS_CVMEM_NULL);
    return(CVLAPACK_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    CVProcessError(cv_mem, CVLAPACK_ILL_INPUT, "CVLAPACK", "CVLapackDense", MSGLS_BAD_NVECTOR);
    return(CVLAPACK_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  linit  = cvLapackDenseInit;
  lsetup = cvLapackDenseSetup;
  lsolve = cvLapackDenseSolve;
  lfree  = cvLapackDenseFree;

  /* Get memory for CVLapackMemRec */
  cvlapack_mem = NULL;
  cvlapack_mem = (CVLapackMem) malloc(sizeof(CVLapackMemRec));
  if (cvlapack_mem == NULL) {
    CVProcessError(cv_mem, CVLAPACK_MEM_FAIL, "CVLAPACK", "CVLapackDense", MSGLS_MEM_FAIL);
    return(CVLAPACK_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = LAPACK_DENSE;

  /* Set default Jacobian routine and Jacobian data */
  djac   = NULL;
  J_data = NULL;

  last_flag = CVLAPACK_SUCCESS;
  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for M, pivot array, and (if needed) savedJ */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = LapackAllocDenseMat(N, N);
  if (M == NULL) {
    CVProcessError(cv_mem, CVLAPACK_MEM_FAIL, "CVLAPACK", "CVLapackDense", MSGLS_MEM_FAIL);
    free(cvlapack_mem);
    return(CVLAPACK_MEM_FAIL);
  }
  pivots = LapackAllocIntArray(N);
  if (pivots == NULL) {
    CVProcessError(cv_mem, CVLAPACK_MEM_FAIL, "CVLAPACK", "CVLapackDense", MSGLS_MEM_FAIL);
    LapackFreeMat(M);
    free(cvlapack_mem);
    return(CVLAPACK_MEM_FAIL);
  }
  savedJ = LapackAllocDenseMat(N, N);
  if (savedJ == NULL) {
    CVProcessError(cv_mem, CVLAPACK_MEM_FAIL, "CVLAPACK", "CVLapackDense", MSGLS_MEM_FAIL);
    LapackFreeMat(M);
    LapackFreeArray(pivots);
    free(cvlapack_mem);
    return(CVLAPACK_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cvlapack_mem;

  return(CVLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVLapackBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module. It first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 * to be cvLapackBandInit, cvLapackBandSetup, cvLapackBandSolve, 
 * and cvLapackBandFree, respectively.  It allocates memory for a 
 * structure of type CVLapackBandMemRec and sets the cv_lmem field in 
 * (*cvode_mem) to the address of this structure.  It sets setupNonNull 
 * in (*cvode_mem) to be TRUE, mu to be mupper, ml to be mlower, and 
 * the jacE and jacI field to NULL.
 * Finally, it allocates memory for M, pivots, and (if needed) savedJ.  
 * The CVLapackBand return value is CVLAPACK_SUCCESS = 0, 
 * CVLAPACK_MEM_FAIL = -1, or CVLAPACK_ILL_INPUT = -2.
 *
 * NOTE: The CVLAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVLapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int CVLapackBand(void *cvode_mem, int N, int mupper, int mlower)
{
  CVodeMem cv_mem;
  CVLapackMem cvlapack_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVLAPACK_MEM_NULL, "CVLAPACK", "CVLapackBand", MSGLS_CVMEM_NULL);
    return(CVLAPACK_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (tempv->ops->nvgetarraypointer == NULL) {
    CVProcessError(cv_mem, CVLAPACK_ILL_INPUT, "CVLAPACK", "CVLapackBand", MSGLS_BAD_NVECTOR);
    return(CVLAPACK_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */  
  linit  = cvLapackBandInit;
  lsetup = cvLapackBandSetup;
  lsolve = cvLapackBandSolve;
  lfree  = cvLapackBandFree;
  
  /* Get memory for CVLapackMemRec */
  cvlapack_mem = NULL;
  cvlapack_mem = (CVLapackMem) malloc(sizeof(CVLapackMemRec));
  if (cvlapack_mem == NULL) {
    CVProcessError(cv_mem, CVLAPACK_MEM_FAIL, "CVLAPACK", "CVLapackBand", MSGLS_MEM_FAIL);
    return(CVLAPACK_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = LAPACK_BAND;

  /* Set default Jacobian routine and Jacobian data */
  bjac  = NULL;
  J_data = NULL;

  last_flag = CVLAPACK_SUCCESS;
  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in cvlapack_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    CVProcessError(cv_mem, CVLAPACK_ILL_INPUT, "CVLAPACK", "CVLapackBand", MSGLS_BAD_SIZES);
    return(CVLAPACK_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(N-1, mu + ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = LapackAllocBandMat(N, mu, ml, smu);
  if (M == NULL) {
    CVProcessError(cv_mem, CVLAPACK_MEM_FAIL, "CVLAPACK", "CVLapackBand", MSGLS_MEM_FAIL);
    free(cvlapack_mem);
    return(CVLAPACK_MEM_FAIL);
  }  
  pivots = LapackAllocIntArray(N);
  if (pivots == NULL) {
    CVProcessError(cv_mem, CVLAPACK_MEM_FAIL, "CVLAPACK", "CVLapackBand", MSGLS_MEM_FAIL);
    LapackFreeMat(M);
    free(cvlapack_mem);
    return(CVLAPACK_MEM_FAIL);
  }
  savedJ = LapackAllocBandMat(N, mu, ml, smu);
  if (savedJ == NULL) {
    CVProcessError(cv_mem, CVLAPACK_MEM_FAIL, "CVLAPACK", "CVLapackBand", MSGLS_MEM_FAIL);
    LapackFreeMat(M);
    LapackFreeArray(pivots);
    free(cvlapack_mem);
    return(CVLAPACK_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cvlapack_mem;

  return(CVLAPACK_SUCCESS);
}




/*
 * -----------------------------------------------------------------
 * Optional I/O functions for 
 * -----------------------------------------------------------------
 */

/*
 * CVLapackSetJacFn specifies the (dense or band) Jacobian function.
 */
int CVLapackSetJacFn(void *cvode_mem, void *jac, void *jac_data)
{
  CVodeMem cv_mem;
  CVLapackMem cvlapack_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVLAPACK_MEM_NULL, "CVLAPACK", "CVLapackSetJacFn", MSGLS_CVMEM_NULL);
    return(CVLAPACK_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVLAPACK_LMEM_NULL, "CVLAPACK", "CVLapackSetJacFn", MSGLS_LMEM_NULL);
    return(CVLAPACK_LMEM_NULL);
  }
  cvlapack_mem = (CVLapackMem) lmem;

  if (mtype == LAPACK_DENSE)
    djac = (CVLapackDenseJacFn) jac;
  else if (mtype == LAPACK_BAND)
    bjac = (CVLapackBandJacFn) jac;

  J_data = jac_data;

  return(CVLAPACK_SUCCESS);
}

/*
 * CVLapackGetWorkSpace returns the length of workspace allocated for the
 * CVLAPACK linear solver.
 */
int CVLapackGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS)
{
  CVodeMem cv_mem;
  CVLapackMem cvlapack_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVLAPACK_MEM_NULL, "CVLAPACK", "CVLapackGetWorkSpace", MSGLS_CVMEM_NULL);
    return(CVLAPACK_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVLAPACK_LMEM_NULL, "CVLAPACK", "CVLapackGetWorkSpace", MSGLS_LMEM_NULL);
    return(CVLAPACK_LMEM_NULL);
  }
  cvlapack_mem = (CVLapackMem) lmem;

  if (mtype == LAPACK_DENSE) {
    *lenrwLS = 2*n*n;
    *leniwLS = n;
  } else if (mtype == LAPACK_BAND) {

  }

  return(CVLAPACK_SUCCESS);
}

/*
 * CVLapackGetNumJacEvals returns the number of Jacobian evaluations.
 */
int CVLapackGetNumJacEvals(void *cvode_mem, long int *njevals)
{
  CVodeMem cv_mem;
  CVLapackMem cvlapack_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVLAPACK_MEM_NULL, "CVLAPACK", "CVLapackGetNumJacEvals", MSGLS_CVMEM_NULL);
    return(CVLAPACK_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVLAPACK_LMEM_NULL, "CVLAPACK", "CVLapackGetNumJacEvals", MSGLS_LMEM_NULL);
    return(CVLAPACK_LMEM_NULL);
  }
  cvlapack_mem = (CVLapackMem) lmem;

  *njevals = nje;

  return(CVLAPACK_SUCCESS);
}

/*
 * CVLapackGetNumRhsEvals returns the number of calls to the ODE function
 * needed for the DQ Jacobian approximation.
 */
int CVLapackGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{
  CVodeMem cv_mem;
  CVLapackMem cvlapack_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVLAPACK_MEM_NULL, "CVLAPACK", "CVLapackGetNumRhsEvals", MSGLS_CVMEM_NULL);
    return(CVLAPACK_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVLAPACK_LMEM_NULL, "CVLAPACK", "CVLapackGetNumRhsEvals", MSGLS_LMEM_NULL);
    return(CVLAPACK_LMEM_NULL);
  }
  cvlapack_mem = (CVLapackMem) lmem;

  *nfevalsLS = nfeDQ;

  return(CVLAPACK_SUCCESS);
}

/*
 * CVLapackGetReturnFlagName returns the name associated with a CVLAPACK
 * return value.
 */
char *CVLapackGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CVLAPACK_SUCCESS:
    sprintf(name,"CVLAPACK_SUCCESS");
    break;   
  case CVLAPACK_MEM_NULL:
    sprintf(name,"CVLAPACK_MEM_NULL");
    break;
  case CVLAPACK_LMEM_NULL:
    sprintf(name,"CVLAPACK_LMEM_NULL");
    break;
  case CVLAPACK_ILL_INPUT:
    sprintf(name,"CVLAPACK_ILL_INPUT");
    break;
  case CVLAPACK_MEM_FAIL:
    sprintf(name,"CVLAPACK_MEM_FAIL");
    break;
  case CVLAPACK_JACFUNC_UNRECVR:
    sprintf(name,"CVLAPACK_JACFUNC_UNRECVR");
    break;
  case CVLAPACK_JACFUNC_RECVR:
    sprintf(name,"CVLAPACK_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * CVLapackGetLastFlag returns the last flag set in a CVLAPACK function.
 */
int CVLapackGetLastFlag(void *cvode_mem, int *flag)
{
  CVodeMem cv_mem;
  CVLapackMem cvlapack_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVLAPACK_MEM_NULL, "CVLAPACK", "CVLapackGetLastFlag", MSGLS_CVMEM_NULL);
    return(CVLAPACK_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVLAPACK_LMEM_NULL, "CVLAPACK", "CVLapackGetLastFlag", MSGLS_LMEM_NULL);
    return(CVLAPACK_LMEM_NULL);
  }
  cvlapack_mem = (CVLapackMem) lmem;

  *flag = last_flag;

  return(CVLAPACK_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * cvLapackDenseInit does remaining initializations specific to the dense
 * linear solver.
 */
static int cvLapackDenseInit(CVodeMem cv_mem)
{
  CVLapackMem cvlapack_mem;

  cvlapack_mem = (CVLapackMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;
  
  if (djac == NULL) {
    djac = cvLapackDenseDQJac;
    J_data = cv_mem;
  } 
  
  last_flag = CVLAPACK_SUCCESS;
  return(0);
}

/*
 * cvLapackDenseSetup does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy. In any case, it constructs the Newton matrix M = I - gamma*J
 * updates counters, and calls the dense LU factorization routine.
 */
static int cvLapackDenseSetup(CVodeMem cv_mem, int convfail,
                              N_Vector yP, N_Vector fctP,
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CVLapackMem cvlapack_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;

  cvlapack_mem = (CVLapackMem) lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + CVL_MSBJ) ||
    ((convfail == CV_FAIL_BAD_J) && (dgamma < CVL_DGMAX)) ||
    (convfail == CV_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {
    
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    dcopy_f77(&(savedJ->ldata), savedJ->data, &one, M->data, &one);
    
  } else {
    
    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    
    retval = djac(n, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      dcopy_f77(&(M->ldata), M->data, &one, savedJ->data, &one);
    } else if (retval < 0) {
      CVProcessError(cv_mem, CVLAPACK_JACFUNC_UNRECVR, "CVLAPACK", "cvLapackDenseSetup", MSGLS_JACFUNC_FAILED);
      last_flag = CVLAPACK_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CVLAPACK_JACFUNC_RECVR;
      return(1);
    }
    
  }

  /* Scale J by - gamma */
  fact = -gamma;
  dscal_f77(&(M->ldata), &fact, M->data, &one);
  
  /* Add identity to get M = I - gamma*J*/
  cvLapackDenseAddI(M);

  /* Do LU factorization of M */
  dgetrf_f77(&n, &n, M->data, &(M->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * cvLapackDenseAddI overwrites the dense matrix A with I+A
 */
static void cvLapackDenseAddI(LapackMat A)
{
  int j;
  realtype *col_j;
  for (j=0; j<A->N; j++) {
    col_j = A->cols[j];
    col_j[j] += ONE;
  }
}

/*
 * cvLapackDenseSolve handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.
 */
static int cvLapackDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector fctC)
{
  CVLapackMem cvlapack_mem;
  realtype *bd, fact;
  int ier, one = 1;

  cvlapack_mem = (CVLapackMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  dgetrs_f77("N", &n, &one, M->data, &(M->ldim), pivots, bd, &n, &ier, 1); 
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&n, &fact, bd, &one); 
  }
  
  last_flag = CVLAPACK_SUCCESS;
  return(0);
}

/*
 * cvLapackDenseFree frees memory specific to the dense linear solver.
 */
static void cvLapackDenseFree(CVodeMem cv_mem)
{
  CVLapackMem  cvlapack_mem;

  cvlapack_mem = (CVLapackMem) lmem;
  
  LapackFreeMat(M);
  LapackFreeArray(pivots);
  LapackFreeMat(savedJ);
  free(cvlapack_mem); 
  cvlapack_mem = NULL;
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH BAND JACOBIANS
 * =================================================================
 */

/*
 * cvLapackBandInit does remaining initializations specific to the band
 * linear solver.
 */
static int cvLapackBandInit(CVodeMem cv_mem)
{
  CVLapackMem cvlapack_mem;

  cvlapack_mem = (CVLapackMem) lmem;

  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  if (bjac == NULL) {
    bjac = cvLapackBandDQJac;
    J_data = cv_mem;
  } 
  
  last_flag = CVLAPACK_SUCCESS;
  return(0);
}

/*
 * cvLapackBandSetup does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy. In any case, it constructs the Newton matrix M = I - gamma*J, 
 * updates counters, and calls the band LU factorization routine.
 */
static int cvLapackBandSetup(CVodeMem cv_mem, int convfail, 
                             N_Vector yP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CVLapackMem cvlapack_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;

  cvlapack_mem = (CVLapackMem) lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + CVL_MSBJ) ||
    ((convfail == CV_FAIL_BAD_J) && (dgamma < CVL_DGMAX)) ||
    (convfail == CV_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {
    
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    dcopy_f77(&(savedJ->ldata), savedJ->data, &one, M->data, &one);
    
  } else {
    
    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    
    retval = bjac(n, mu, ml, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      dcopy_f77(&(M->ldata), M->data, &one, savedJ->data, &one);
    } else if (retval < 0) {
      CVProcessError(cv_mem, CVLAPACK_JACFUNC_UNRECVR, "CVLAPACK", "cvLapackBandSetup", MSGLS_JACFUNC_FAILED);
      last_flag = CVLAPACK_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CVLAPACK_JACFUNC_RECVR;
      return(1);
    }
    
  }
  
  /* Scale J by - gamma */
  fact = -gamma;
  dscal_f77(&(M->ldata), &fact, M->data, &one);
  
  /* Add identity to get M = I - gamma*J*/
  cvLapackBandAddI(M);
  
  /* Do LU factorization of M */
  dgbtrf_f77(&n, &n, &ml, &mu, M->data, &(M->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);

}

/*
 * cvLapackBandSolve handles the solve operation for the band linear solver
 * by calling the band backsolve routine.
 */
static int cvLapackBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector fctC)
{
  CVLapackMem cvlapack_mem;
  realtype *bd, fact;
  int ier, one = 1;

  cvlapack_mem = (CVLapackMem) lmem;

  bd = N_VGetArrayPointer(b);

  dgbtrs_f77("N", &n, &ml, &mu, &one, M->data, &(M->ldim), pivots, bd, &n, &ier, 1);
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&n, &fact, bd, &one); 
  }

  last_flag = CVLAPACK_SUCCESS;
  return(0);
}

/*
 * cvLapackBandFree frees memory specific to the band linear solver.
 */
static void cvLapackBandFree(CVodeMem cv_mem)
{
  CVLapackMem  cvlapack_mem;

  cvlapack_mem = (CVLapackMem) lmem;
  
  LapackFreeMat(M);
  LapackFreeArray(pivots);
  LapackFreeMat(savedJ);
  free(cvlapack_mem); 
  cvlapack_mem = NULL;
}

/*
 * cvLapackBandAddI overwrites the banded matrix A with I+A
 */
static void cvLapackBandAddI(LapackMat A)
{
  int j;
  realtype *col_j;
  for (j=0; j<A->N; j++) {
    col_j = (A->cols)[j] + A->storage_mu;
    col_j[0] += ONE;
  }
}

/* 
 * =================================================================
 *  DENSE DQ JACOBIAN APPROXIMATIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cvLapackDenseDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a dense matrix of type
 * LapackMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro LAPACK_DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */ 
static int cvLapackDenseDQJac(int N, realtype t,
                              N_Vector y, N_Vector fy, 
                              LapackMat Jac, void *jac_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype fnorm, minInc, inc, inc_inv, yjsaved, srur;
  realtype *tmp2_data, *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  int j;
  int retval = 0;

  CVodeMem cv_mem;
  CVLapackMem  cvlapack_mem;

  /* jac_data points to cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvlapack_mem = (CVLapackMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */

    N_VSetArrayPointer(LAPACK_DENSE_COL(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = MAX(srur*ABS(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = f(t, y, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;
    
    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

    LAPACK_DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}

/* 
 * =================================================================
 *  BAND DQ JACOBIAN APPROXIMATIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

static int cvLapackBandDQJac(int N, int mupper, int mlower,
                             realtype t, N_Vector y, N_Vector fy, 
                             LapackMat Jac, void *jac_data,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  N_Vector ftemp, ytemp;
  realtype fnorm, minInc, inc, inc_inv, srur;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int group, i, j, width, ngroups, i1, i2;
  int retval = 0;

  CVodeMem cv_mem;
  CVLapackMem  cvlapack_mem;

  /* jac_dat points to cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvlapack_mem = (CVLapackMem) lmem;

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */

    retval = f(tn, ytemp, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = LAPACK_BAND_COL(Jac,j);
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        LAPACK_BAND_COL_ELEM(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }
  
  return(retval);
}

