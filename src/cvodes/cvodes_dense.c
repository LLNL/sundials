/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-11-22 00:12:49 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CVSDENSE linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes_dense.h>
#include "cvodes_direct_impl.h"
#include "cvodes_impl.h"

#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* CVSDENSE linit, lsetup, lsolve, and lfree routines */
static int cvDenseInit(CVodeMem cv_mem);
static int cvDenseSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, 
                        N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
static int cvDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                        N_Vector ycur, N_Vector fcur);
static void cvDenseFree(CVodeMem cv_mem);

/* CVSDENSE lfreeB function */
static void cvDenseFreeB(CVadjMem ca_mem);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */


/* Readability Replacements */

#define lmm       (cv_mem->cv_lmm)
#define f         (cv_mem->cv_f)
#define f_data    (cv_mem->cv_f_data)
#define uround    (cv_mem->cv_uround)
#define nst       (cv_mem->cv_nst)
#define tn        (cv_mem->cv_tn)
#define h         (cv_mem->cv_h)
#define gamma     (cv_mem->cv_gamma)
#define gammap    (cv_mem->cv_gammap)
#define gamrat    (cv_mem->cv_gamrat)
#define ewt       (cv_mem->cv_ewt)
#define linit     (cv_mem->cv_linit)
#define lsetup    (cv_mem->cv_lsetup)
#define lsolve    (cv_mem->cv_lsolve)
#define lfree     (cv_mem->cv_lfree)
#define lmem      (cv_mem->cv_lmem)
#define vec_tmpl     (cv_mem->cv_tempv)
#define setupNonNull (cv_mem->cv_setupNonNull)

#define mtype     (cvdls_mem->d_type)
#define n         (cvdls_mem->d_n)
#define jac       (cvdls_mem->d_djac)
#define M         (cvdls_mem->d_M)
#define pivots    (cvdls_mem->d_pivots)
#define savedJ    (cvdls_mem->d_savedJ)
#define nstlj     (cvdls_mem->d_nstlj)
#define nje       (cvdls_mem->d_nje)
#define nfeDQ     (cvdls_mem->d_nfeDQ)
#define J_data    (cvdls_mem->d_J_data)
#define last_flag (cvdls_mem->d_last_flag)

/*
 * -----------------------------------------------------------------
 * CVDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module.  CVDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be cvDenseInit, cvDenseSetup, cvDenseSolve, and cvDenseFree,
 * respectively.  It allocates memory for a structure of type
 * CVDlsMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure.  It sets setupNonNull in (*cvode_mem) to
 * TRUE, and the d_jac field to the default CVDenseDQJac.
 * Finally, it allocates memory for M, savedJ, and pivots.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int CVDense(void *cvode_mem, int N)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_MEM_NULL, "CVSDENSE", "CVDense", MSGD_CVMEM_NULL);
    return(CVDIRECT_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    CVProcessError(cv_mem, CVDIRECT_ILL_INPUT, "CVSDENSE", "CVDense", MSGD_BAD_NVECTOR);
    return(CVDIRECT_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  linit  = cvDenseInit;
  lsetup = cvDenseSetup;
  lsolve = cvDenseSolve;
  lfree  = cvDenseFree;

  /* Get memory for CVDlsMemRec */
  cvdls_mem = NULL;
  cvdls_mem = (CVDlsMem) malloc(sizeof(CVDlsMemRec));
  if (cvdls_mem == NULL) {
    CVProcessError(cv_mem, CVDIRECT_MEM_FAIL, "CVSDENSE", "CVDense", MSGD_MEM_FAIL);
    return(CVDIRECT_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Set default Jacobian routine and Jacobian data */
  jac = cvDlsDenseDQJac;
  J_data = cvode_mem;
  last_flag = CVDIRECT_SUCCESS;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for M, savedJ, and pivot array */

  M = NULL;
  M = NewDenseMat(N, N);
  if (M == NULL) {
    CVProcessError(cv_mem, CVDIRECT_MEM_FAIL, "CVSDENSE", "CVDense", MSGD_MEM_FAIL);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDIRECT_MEM_FAIL);
  }
  savedJ = NULL;
  savedJ = NewDenseMat(N, N);
  if (savedJ == NULL) {
    CVProcessError(cv_mem, CVDIRECT_MEM_FAIL, "CVSDENSE", "CVDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDIRECT_MEM_FAIL);
  }
  pivots = NULL;
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    CVProcessError(cv_mem, CVDIRECT_MEM_FAIL, "CVSDENSE", "CVDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    DestroyMat(savedJ);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDIRECT_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cvdls_mem;

  return(CVDIRECT_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * cvDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cvDenseInit(CVodeMem cv_mem)
{
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;
  
  if (jac == NULL) {
    jac = cvDlsDenseDQJac;
    J_data = cv_mem;
  }

  last_flag = CVDIRECT_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the dense LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int cvDenseSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, 
                        N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  CVDlsMem cvdls_mem;
  booleantype jbad, jok;
  realtype dgamma;
  int ier, retval;

  cvdls_mem = (CVDlsMem) lmem;
 
  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
 
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + CVD_MSBJ) ||
         ((convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
         (convfail == CV_FAIL_OTHER);
  jok = !jbad;
 
  if (jok) {

    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    DenseCopy(savedJ, M);

  } else {

    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    DenseZero(M);

    retval = jac(n, tn, ypred, fpred, M, J_data, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      CVProcessError(cv_mem, CVDIRECT_JACFUNC_UNRECVR, "CVSDENSE", "cvDenseSetup", MSGD_JACFUNC_FAILED);
      last_flag = CVDIRECT_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      last_flag = CVDIRECT_JACFUNC_RECVR;
      return(1);
    }

    DenseCopy(M, savedJ);

  }
  
  /* Scale and add I to get M = I - gamma*J */
  DenseScale(-gamma, M);
  DenseAddI(M);

  /* Do LU factorization of M */
  ier = DenseGETRF(M, pivots); 

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int cvDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                        N_Vector ycur, N_Vector fcur)
{
  CVDlsMem cvdls_mem;
  realtype *bd;

  cvdls_mem = (CVDlsMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseGETRS(M, pivots, bd);

  /* If CV_BDF, scale the correction to account for change in gamma */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }
  
  last_flag = CVDIRECT_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void cvDenseFree(CVodeMem cv_mem)
{
  CVDlsMem  cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;
  
  DestroyMat(M);
  DestroyMat(savedJ);
  DestroyArray(pivots);
  free(cvdls_mem); cvdls_mem = NULL;
}

/* 
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */

/*
 * CVDenseB is a wraper around CVDense. It attaches the CVSDENSE linear solver
 * to the adjoint memory block.
 */

int CVDenseB(void *cvadj_mem, int nB)
{
  CVadjMem ca_mem;
  CVDlsMemB cvdlsB_mem;
  CVodeMem cvB_mem;
  int flag;

  if (cvadj_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_ADJMEM_NULL, "CVSDENSE", "CVDenseB", MSGD_CAMEM_NULL);
    return(CVDIRECT_ADJMEM_NULL);
  }
  ca_mem = (CVadjMem) cvadj_mem;

  cvB_mem = ca_mem->cvb_mem;

  /* Get memory for CVDlsMemRecB */
  cvdlsB_mem = (CVDlsMemB) malloc(sizeof(CVDlsMemRecB));
  if (cvdlsB_mem == NULL) {
    CVProcessError(cvB_mem, CVDIRECT_MEM_FAIL, "CVSDENSE", "CVDenseB", MSGD_MEM_FAIL);
    return(CVDIRECT_MEM_FAIL);
  }

  /* set matrix type */
  cvdlsB_mem->d_typeB = SUNDIALS_DENSE;

  /* initialize Jacobian function */
  cvdlsB_mem->d_djacB = NULL;
  cvdlsB_mem->d_jac_dataB = NULL;

  /* attach lmemB and lfreeB */
  ca_mem->ca_lmemB = cvdlsB_mem;
  ca_mem->ca_lfreeB = cvDenseFreeB;

  flag = CVDense(cvB_mem, nB);

  if (flag != CVDIRECT_SUCCESS) {
    free(cvdlsB_mem);
    cvdlsB_mem = NULL;
  }

  return(flag);
}

/*
 * cvDenseFreeB frees the memory associated with the CVSDENSE linear
 * solver for backward integration.
 */

static void cvDenseFreeB(CVadjMem ca_mem)
{
  CVDlsMemB cvdlsB_mem;

  cvdlsB_mem = (CVDlsMemB) ca_mem->ca_lmemB;

  free(cvdlsB_mem);
}

