/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
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
 * This is the implementation file for a dense or banded CVODES 
 * linear solver using BLAS and LAPACK functions.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes_lapack.h>
#include "cvodes_direct_impl.h"
#include "cvodes_impl.h"

#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* CVSLAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int cvLapackDenseInit(CVodeMem cv_mem);
static int cvLapackDenseSetup(void* cur_state, N_Vector tmp1,
                              N_Vector tmp2, N_Vector tmp3);
static int cvLapackDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector fctC);
static int cvLapackDenseFree(CVodeMem cv_mem);

/* CVSLAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int cvLapackBandInit(CVodeMem cv_mem);
static int cvLapackBandSetup(void* cur_state, N_Vector tmp1,
                             N_Vector tmp2, N_Vector tmp3);
static int cvLapackBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector fctC);
static int cvLapackBandFree(CVodeMem cv_mem);

/* CVSLAPACK lfreeB functions */
static int cvLapackDenseFreeB(CVodeBMem cvB_mem);
static int cvLapackBandFreeB(CVodeBMem cvB_mem);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
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
 * structure of type CVDlsMemRec and sets the cv_lmem field in 
 * (*cvode_mem) to the address of this structure.  It sets setupNonNull 
 * in (*cvode_mem) to TRUE, and the d_jac field to the default 
 * cvDlsDenseDQJac. Finally, it allocates memory for M, pivots, and 
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
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVSLAPACK", "CVLapackDense", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the CVSLAPACK solver */
  if (cv_mem->cv_tempv->ops->nvgetarraypointer == NULL ||
      cv_mem->cv_tempv->ops->nvsetarraypointer == NULL) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackDense", MSGD_BAD_NVECTOR);
    return(CVDLS_ILL_INPUT);
  }

  if (cv_mem->cv_lfree !=NULL) cv_mem->cv_lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  cv_mem->cv_linit  = cvLapackDenseInit;
  cv_mem->cv_lsetup = cvLapackDenseSetup;
  cv_mem->cv_lsolve = cvLapackDenseSolve;
  cv_mem->cv_lfree  = cvLapackDenseFree;

  /* Get memory for CVDlsMemRec */
  cvdls_mem = NULL;
  cvdls_mem = (CVDlsMem) malloc(sizeof(struct CVDlsMemRec));
  if (cvdls_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDense", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* Set matrix type */
  cvdls_mem->d_type = SUNDIALS_DENSE;

  /* Initialize Jacobian-related data */
  cvdls_mem->d_jacDQ  = TRUE;
  cvdls_mem->d_djac   = NULL;
  cvdls_mem->d_J_data = NULL;

  cvdls_mem->d_last_flag = CVDLS_SUCCESS;
  cvDlsInitializeCounters(cvdls_mem);  
  cv_mem->cv_setupNonNull = TRUE;

  /* Set problem dimension */
  cvdls_mem->d_n = (sunindextype) N;

  /* Allocate memory for M, pivot array, and (if needed) savedJ */
  cvdls_mem->d_M = NULL;
  cvdls_mem->d_pivots = NULL;
  cvdls_mem->d_savedJ = NULL;

  cvdls_mem->d_M = NewDenseMat(cvdls_mem->d_n, cvdls_mem->d_n);
  if (cvdls_mem->d_M == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDense", MSGD_MEM_FAIL);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  cvdls_mem->d_pivots = NewIntArray(N);
  if (cvdls_mem->d_pivots == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDense", MSGD_MEM_FAIL);
    DestroyMat(cvdls_mem->d_M);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  cvdls_mem->d_savedJ = NewDenseMat(cvdls_mem->d_n, cvdls_mem->d_n);
  if (cvdls_mem->d_savedJ == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDense", MSGD_MEM_FAIL);
    DestroyMat(cvdls_mem->d_M);
    DestroyArray(cvdls_mem->d_pivots);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  cv_mem->cv_lmem = cvdls_mem;

  return(CVDLS_SUCCESS);
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
 * The CVLapackBand return value is CVDLS_SUCCESS = 0, 
 * CVDLS_MEM_FAIL = -1, or CVDLS_ILL_INPUT = -2.
 *
 * NOTE: The CVSLAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVLapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int CVLapackBand(void *cvode_mem, int N, int mupper, int mlower)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVSLAPACK", "CVLapackBand", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (cv_mem->cv_tempv->ops->nvgetarraypointer == NULL) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackBand", MSGD_BAD_NVECTOR);
    return(CVDLS_ILL_INPUT);
  }

  if (cv_mem->cv_lfree != NULL) cv_mem->cv_lfree(cv_mem);

  /* Set four main function fields in cv_mem */  
  cv_mem->cv_linit  = cvLapackBandInit;
  cv_mem->cv_lsetup = cvLapackBandSetup;
  cv_mem->cv_lsolve = cvLapackBandSolve;
  cv_mem->cv_lfree  = cvLapackBandFree;
  
  /* Get memory for CVDlsMemRec */
  cvdls_mem = NULL;
  cvdls_mem = (CVDlsMem) malloc(sizeof(struct CVDlsMemRec));
  if (cvdls_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBand", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* Set matrix type */
  cvdls_mem->d_type = SUNDIALS_BAND;

  /* Initialize Jacobian-related data */
  cvdls_mem->d_jacDQ  = TRUE;
  cvdls_mem->d_bjac   = NULL;
  cvdls_mem->d_J_data = NULL;

  cvdls_mem->d_last_flag = CVDLS_SUCCESS;
  cvDlsInitializeCounters(cvdls_mem);  
  cv_mem->cv_setupNonNull = TRUE;
  
  /* Load problem dimension */
  cvdls_mem->d_n = (sunindextype) N;

  /* Load half-bandwiths in cvdls_mem */
  cvdls_mem->d_ml = (sunindextype) mlower;
  cvdls_mem->d_mu = (sunindextype) mupper;

  /* Test ml and mu for legality */
  if ((cvdls_mem->d_ml < 0) || (cvdls_mem->d_mu < 0) ||
      (cvdls_mem->d_ml >= cvdls_mem->d_n) || (cvdls_mem->d_mu >= cvdls_mem->d_n)) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackBand", MSGD_BAD_SIZES);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  cvdls_mem->d_smu = cvdls_mem->d_mu + cvdls_mem->d_ml;

  /* Allocate memory for M, savedJ, and pivot arrays */
  cvdls_mem->d_M = NULL;
  cvdls_mem->d_pivots = NULL;
  cvdls_mem->d_savedJ = NULL;

  cvdls_mem->d_M = NewBandMat(cvdls_mem->d_n, cvdls_mem->d_mu,
                              cvdls_mem->d_ml, cvdls_mem->d_smu);
  if (cvdls_mem->d_M == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBand", MSGD_MEM_FAIL);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }  
  cvdls_mem->d_pivots = NewIntArray(N);
  if (cvdls_mem->d_pivots == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBand", MSGD_MEM_FAIL);
    DestroyMat(cvdls_mem->d_M);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  cvdls_mem->d_savedJ = NewBandMat(cvdls_mem->d_n, cvdls_mem->d_mu,
                                   cvdls_mem->d_ml, cvdls_mem->d_smu);
  if (cvdls_mem->d_savedJ == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBand", MSGD_MEM_FAIL);
    DestroyMat(cvdls_mem->d_M);
    DestroyArray(cvdls_mem->d_pivots);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  cv_mem->cv_lmem = cvdls_mem;

  return(CVDLS_SUCCESS);
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
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;
  
  cvDlsInitializeCounters(cvdls_mem);  

  /* Set Jacobian function and data, depending on jacDQ */
  if (cvdls_mem->d_jacDQ) {
    cvdls_mem->d_djac = cvDlsDenseDQJac;
    cvdls_mem->d_J_data = cv_mem;
  } else {
    cvdls_mem->d_J_data = cv_mem->cv_user_data;
  }

  cvdls_mem->d_last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * cvLapackDenseSetup does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy. In any case, it constructs the Newton matrix M = I - gamma*J
 * updates counters, and calls the dense LU factorization routine.
 */
static int cvLapackDenseSetup(void* current_state, N_Vector tmp1,
                              N_Vector tmp2, N_Vector tmp3)
{
  CVDlsMem cvdls_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;
  int intn, lenmat;
  CVodeMem cv_mem;
  cvLinPoint *cur_state;

  cur_state = (cvLinPoint*) current_state;
  cv_mem = cur_state->cv_mem;
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;
  intn = (int) cvdls_mem->d_n;
  lenmat = cvdls_mem->d_M->ldata;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((cv_mem->cv_gamma/cv_mem->cv_gammap) - ONE);
  jbad = (cv_mem->cv_nst == 0) ||
    (cv_mem->cv_nst > cvdls_mem->d_nstlj + CVD_MSBJ) ||
    ((cv_mem->cv_convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
    (cv_mem->cv_convfail == CV_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {
    
    /* If jok = TRUE, use saved copy of J */
    cv_mem->cv_jcur = FALSE;
    dcopy_f77(&lenmat, cvdls_mem->d_savedJ->data, &one,
              cvdls_mem->d_M->data, &one);
    
  } else {
    
    /* If jok = FALSE, call jac routine for new J value */
    cvdls_mem->d_nje++;
    cvdls_mem->d_nstlj = cv_mem->cv_nst;
    cv_mem->cv_jcur = TRUE;
    SetToZero(cvdls_mem->d_M);

    retval = cvdls_mem->d_djac(cvdls_mem->d_n, cur_state->t,
                               cur_state->y, cur_state->f, cvdls_mem->d_M, cvdls_mem->d_J_data,
                               tmp1, tmp2, tmp3);
    if (retval == 0) {
      dcopy_f77(&lenmat, cvdls_mem->d_M->data, &one, cvdls_mem->d_savedJ->data, &one);
    } else if (retval < 0) {
      cvProcessError(cv_mem, CVDLS_JACFUNC_UNRECVR, "CVSLAPACK", "cvLapackDenseSetup", MSGD_JACFUNC_FAILED);
      cvdls_mem->d_last_flag = CVDLS_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      cvdls_mem->d_last_flag = CVDLS_JACFUNC_RECVR;
      return(1);
    }
    
  }

  /* Scale J by - gamma */
  fact = -cv_mem->cv_gamma;
  dscal_f77(&lenmat, &fact, cvdls_mem->d_M->data, &one);
  
  /* Add identity to get M = I - gamma*J*/
  AddIdentity(cvdls_mem->d_M);

  /* Do LU factorization of M */
  dgetrf_f77(&intn, &intn, cvdls_mem->d_M->data,
             &intn, cvdls_mem->d_pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  cvdls_mem->d_last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * cvLapackDenseSolve handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.
 */
static int cvLapackDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector fctC)
{
  CVDlsMem cvdls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  intn = (int) cvdls_mem->d_n;

  bd = N_VGetArrayPointer(b);

  dgetrs_f77("N", &intn, &one, cvdls_mem->d_M->data, &intn,
             cvdls_mem->d_pivots, bd, &intn, &ier, 1); 
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((cv_mem->cv_lmm == CV_BDF) && (cv_mem->cv_gamrat != ONE)) {
    fact = TWO/(ONE + cv_mem->cv_gamrat);
    dscal_f77(&intn, &fact, bd, &one); 
  }
  
  cvdls_mem->d_last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * cvLapackDenseFree frees memory specific to the dense linear solver.
 */
static int cvLapackDenseFree(CVodeMem cv_mem)
{
  CVDlsMem  cvdls_mem;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;
  
  DestroyMat(cvdls_mem->d_M);
  DestroyArray(cvdls_mem->d_pivots);
  DestroyMat(cvdls_mem->d_savedJ);
  free(cvdls_mem); 
  cvdls_mem = NULL;

  return(0);
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
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  cvDlsInitializeCounters(cvdls_mem);  

  /* Set Jacobian function and data, depending on jacDQ */
  if (cvdls_mem->d_jacDQ) {
    cvdls_mem->d_bjac = cvDlsBandDQJac;
    cvdls_mem->d_J_data = cv_mem;
  } else {
    cvdls_mem->d_J_data = cv_mem->cv_user_data;
  }

  cvdls_mem->d_last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * cvLapackBandSetup does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy. In any case, it constructs the Newton matrix M = I - gamma*J, 
 * updates counters, and calls the band LU factorization routine.
 */
static int cvLapackBandSetup(void* current_state, N_Vector tmp1,
                             N_Vector tmp2, N_Vector tmp3)
{
  CVDlsMem cvdls_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;
  int intn, iml, imu, lenmat, ldmat;
  CVodeMem cv_mem;
  cvLinPoint *cur_state;

  cur_state = (cvLinPoint*) current_state;
  cv_mem = cur_state->cv_mem;
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  intn = (int) cvdls_mem->d_n;
  iml = (int) cvdls_mem->d_ml;
  imu = (int) cvdls_mem->d_mu;
  lenmat = cvdls_mem->d_M->ldata;
  ldmat = cvdls_mem->d_M->ldim;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((cv_mem->cv_gamma/cv_mem->cv_gammap) - ONE);
  jbad = (cv_mem->cv_nst == 0) ||
    (cv_mem->cv_nst > cvdls_mem->d_nstlj + CVD_MSBJ) ||
    ((cv_mem->cv_convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
    (cv_mem->cv_convfail == CV_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {
    
    /* If jok = TRUE, use saved copy of J */
    cv_mem->cv_jcur = FALSE;
    dcopy_f77(&lenmat, cvdls_mem->d_savedJ->data, &one,
              cvdls_mem->d_M->data, &one);
    
  } else {
    
    /* If jok = FALSE, call jac routine for new J value */
    cvdls_mem->d_nje++;
    cvdls_mem->d_nstlj = cv_mem->cv_nst;
    cv_mem->cv_jcur = TRUE;
    SetToZero(cvdls_mem->d_M); 

    retval = cvdls_mem->d_bjac(cvdls_mem->d_n, cvdls_mem->d_mu,
                               cvdls_mem->d_ml, cur_state->t,
                               cur_state->y, cur_state->f,
                               cvdls_mem->d_M, cvdls_mem->d_J_data,
                               tmp1, tmp2, tmp3);
    if (retval == 0) {
      dcopy_f77(&lenmat, cvdls_mem->d_M->data, &one,
                cvdls_mem->d_savedJ->data, &one);
    } else if (retval < 0) {
      cvProcessError(cv_mem, CVDLS_JACFUNC_UNRECVR, "CVSLAPACK", "cvLapackBandSetup", MSGD_JACFUNC_FAILED);
      cvdls_mem->d_last_flag = CVDLS_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      cvdls_mem->d_last_flag = CVDLS_JACFUNC_RECVR;
      return(1);
    }
    
  }
  
  /* Scale J by - gamma */
  fact = -cv_mem->cv_gamma;
  dscal_f77(&lenmat, &fact, cvdls_mem->d_M->data, &one);
  
  /* Add identity to get M = I - gamma*J*/
  AddIdentity(cvdls_mem->d_M);
  
  /* Do LU factorization of M */
  dgbtrf_f77(&intn, &intn, &iml, &imu, cvdls_mem->d_M->data,
             &ldmat, cvdls_mem->d_pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  cvdls_mem->d_last_flag = (long int) ier;
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
  CVDlsMem cvdls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn, iml, imu, ldmat;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  intn = (int) cvdls_mem->d_n;
  iml = (int) cvdls_mem->d_ml;
  imu = (int) cvdls_mem->d_mu;
  ldmat = cvdls_mem->d_M->ldim;

  bd = N_VGetArrayPointer(b);

  dgbtrs_f77("N", &intn, &iml, &imu, &one, cvdls_mem->d_M->data,
             &ldmat, cvdls_mem->d_pivots, bd, &intn, &ier, 1);
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((cv_mem->cv_lmm == CV_BDF) && (cv_mem->cv_gamrat != ONE)) {
    fact = TWO/(ONE + cv_mem->cv_gamrat);
    dscal_f77(&intn, &fact, bd, &one); 
  }

  cvdls_mem->d_last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * cvLapackBandFree frees memory specific to the band linear solver.
 */
static int cvLapackBandFree(CVodeMem cv_mem)
{
  CVDlsMem  cvdls_mem;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;
  
  DestroyMat(cvdls_mem->d_M);
  DestroyArray(cvdls_mem->d_pivots);
  DestroyMat(cvdls_mem->d_savedJ);
  free(cvdls_mem); 
  cvdls_mem = NULL;

  return(0);
}

/* 
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */

/*
 * CVLapackDenseB is a wrapper around CVLapackDense. It attaches the
 * dense CVSLAPACK linear solver to the backward problem memory block.
 */

int CVLapackDenseB(void *cvode_mem, int which, int nB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  CVDlsMemB cvdlsB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVSLAPACK", "CVLapackDenseB", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVDLS_NO_ADJ, "CVSLAPACK", "CVLapackDenseB", MSGD_NO_ADJ);
    return(CVDLS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackDenseB", MSGCV_BAD_WHICH);
    return(CVDLS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* Get memory for CVDlsMemRecB */
  cvdlsB_mem = (CVDlsMemB) malloc(sizeof(struct CVDlsMemRecB));
  if (cvdlsB_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDenseB", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* set matrix type */
  cvdlsB_mem->d_typeB = SUNDIALS_DENSE;

  /* initialize Jacobian function */
  cvdlsB_mem->d_djacB = NULL;

  /* attach lmemB and lfreeB */
  cvB_mem->cv_lmem = cvdlsB_mem;
  cvB_mem->cv_lfree = cvLapackDenseFreeB;

  flag = CVLapackDense(cvodeB_mem, nB);

  if (flag != CVDLS_SUCCESS) {
    free(cvdlsB_mem);
    cvdlsB_mem = NULL;
  }

  return(flag);
}

/*
 * cvLapackDenseFreeB frees the memory associated with the dense CVSLAPACK
 * linear solver for backward integration.
 */

static int cvLapackDenseFreeB(CVodeBMem cvB_mem)
{
  CVDlsMemB cvdlsB_mem;

  cvdlsB_mem = (CVDlsMemB) (cvB_mem->cv_lmem);

  free(cvdlsB_mem);

  return(0);
}

/*
 * CVLapackBandB is a wrapper around CVLapackBand. It attaches the band
 * CVSLAPACK linear solver to the backward problem memory block.
 */

int CVLapackBandB(void *cvode_mem, int which,
                  int nB, int mupperB, int mlowerB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  CVDlsMemB cvdlsB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVSLAPACK", "CVLapackBandB", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVDLS_NO_ADJ, "CVSLAPACK", "CVLapackBandB", MSGD_NO_ADJ);
    return(CVDLS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackBandB", MSGCV_BAD_WHICH);
    return(CVDLS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* Get memory for CVDlsMemRecB */
  cvdlsB_mem = (CVDlsMemB) malloc(sizeof(struct CVDlsMemRecB));
  if (cvdlsB_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBandB", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* set matrix type */
  cvdlsB_mem->d_typeB = SUNDIALS_BAND;

  /* initialize Jacobian function */
  cvdlsB_mem->d_bjacB = NULL;

  /* attach lmemB and lfreeB */
  cvB_mem->cv_lmem = cvdlsB_mem;
  cvB_mem->cv_lfree = cvLapackBandFreeB;

  flag = CVLapackBand(cvodeB_mem, nB, mupperB, mlowerB);

  if (flag != CVDLS_SUCCESS) {
    free(cvdlsB_mem);
    cvdlsB_mem = NULL;
  }

  return(flag);
}

/*
 * cvLapackBandFreeB frees the memory associated with the band CVSLAPACK
 * linear solver for backward integration.
 */

static int cvLapackBandFreeB(CVodeBMem cvB_mem)
{
  CVDlsMemB cvdlsB_mem;

  cvdlsB_mem = (CVDlsMemB) (cvB_mem->cv_lmem);

  free(cvdlsB_mem);

  return(0);
}
