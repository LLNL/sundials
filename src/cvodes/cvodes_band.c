/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2007-03-21 18:56:33 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CVSBAND linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes_band.h>
#include "cvodes_direct_impl.h"
#include "cvodes_impl.h"

#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* CVSBAND linit, lsetup, lsolve, and lfree routines */
static int cvBandInit(CVodeMem cv_mem);
static int cvBandSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3);
static int cvBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur);
static void cvBandFree(CVodeMem cv_mem);

/* CVSBAND lfreeB function */
static void cvBandFreeB(CVodeBMem cvB_mem);

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
#define nfe       (cv_mem->cv_nfe)
#define linit     (cv_mem->cv_linit)
#define lsetup    (cv_mem->cv_lsetup)
#define lsolve    (cv_mem->cv_lsolve)
#define lfree     (cv_mem->cv_lfree)
#define lmem      (cv_mem->cv_lmem)
#define vec_tmpl      (cv_mem->cv_tempv)
#define setupNonNull  (cv_mem->cv_setupNonNull)

#define mtype      (cvdls_mem->d_type)
#define n          (cvdls_mem->d_n)
#define jac        (cvdls_mem->d_bjac)
#define M          (cvdls_mem->d_M)
#define mu         (cvdls_mem->d_mu)
#define ml         (cvdls_mem->d_ml)
#define smu        (cvdls_mem->d_smu)
#define pivots     (cvdls_mem->d_pivots)
#define savedJ     (cvdls_mem->d_savedJ)
#define nstlj      (cvdls_mem->d_nstlj)
#define nje        (cvdls_mem->d_nje)
#define nfeDQ       (cvdls_mem->d_nfeDQ)
#define J_data     (cvdls_mem->d_J_data)
#define last_flag  (cvdls_mem->d_last_flag)

/*
 * -----------------------------------------------------------------
 * CVBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module.  CVBand first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 * to be cvBandInit, cvBandSetup, cvBandSolve, and cvBandFree,
 * respectively.  It allocates memory for a structure of type
 * CVDlsMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure.  It sets setupNonNull in (*cvode_mem) to be
 * TRUE, b_mu to be mupper, b_ml to be mlower, and the b_jac field to be 
 * CVBandDQJac.
 * Finally, it allocates memory for M, savedJ, and pivot.  The CVBand
 * return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.
 *
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */
                  
int CVBand(void *cvode_mem, int N, int mupper, int mlower)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_MEM_NULL, "CVSBAND", "CVBand", MSGD_CVMEM_NULL);
    return(CVDIRECT_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL) {
    CVProcessError(cv_mem, CVDIRECT_ILL_INPUT, "CVSBAND", "CVBand", MSGD_BAD_NVECTOR);
    return(CVDIRECT_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */  
  linit  = cvBandInit;
  lsetup = cvBandSetup;
  lsolve = cvBandSolve;
  lfree  = cvBandFree;
  
  /* Get memory for CVDlsMemRec */
  cvdls_mem = NULL;
  cvdls_mem = (CVDlsMem) malloc(sizeof(CVDlsMemRec));
  if (cvdls_mem == NULL) {
    CVProcessError(cv_mem, CVDIRECT_MEM_FAIL, "CVSBAND", "CVBand", MSGD_MEM_FAIL);
    return(CVDIRECT_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_BAND;
  
  /* Set default Jacobian routine and Jacobian data */
  jac = cvDlsBandDQJac;
  J_data = cvode_mem;
  last_flag = CVDIRECT_SUCCESS;

  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in cvdls_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    CVProcessError(cv_mem, CVDIRECT_ILL_INPUT, "CVSBAND", "CVBand", MSGD_BAD_SIZES);
    return(CVDIRECT_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(N-1, mu + ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  M = NULL;
  M = NewBandMat(N, mu, ml, smu);
  if (M == NULL) {
    CVProcessError(cv_mem, CVDIRECT_MEM_FAIL, "CVSBAND", "CVBand", MSGD_MEM_FAIL);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDIRECT_MEM_FAIL);
  }
  savedJ = NULL;
  savedJ = NewBandMat(N, mu, ml, mu);
  if (savedJ == NULL) {
    CVProcessError(cv_mem, CVDIRECT_MEM_FAIL, "CVSBAND", "CVBand", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDIRECT_MEM_FAIL);
  }
  pivots = NULL;
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    CVProcessError(cv_mem, CVDIRECT_MEM_FAIL, "CVSBAND", "CVBand", MSGD_MEM_FAIL);
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
 * cvBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cvBandInit(CVodeMem cv_mem)
{
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;

  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  if (jac == NULL) {
    jac = cvDlsBandDQJac;
    J_data = cv_mem;
  }

  last_flag = CVDIRECT_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the band LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int cvBandSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3)
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
    BandCopy(savedJ, M, mu, ml);

  } else {

    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    BandZero(M); 

    retval = jac(n, mu, ml, tn, ypred, fpred, M, J_data, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      CVProcessError(cv_mem, CVDIRECT_JACFUNC_UNRECVR, "CVSBAND", "cvBandSetup", MSGD_JACFUNC_FAILED);
      last_flag = CVDIRECT_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      last_flag = CVDIRECT_JACFUNC_RECVR;
      return(1);
    }

    BandCopy(M, savedJ, mu, ml);

  }
  
  /* Scale and add I to get M = I - gamma*J */
  BandScale(-gamma, M);
  BandAddI(M);

  /* Do LU factorization of M */
  ier = BandGBTRF(M, pivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) {
    last_flag = ier;
    return(1);
  }
  last_flag = CVDIRECT_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvBandSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int cvBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur)
{
  CVDlsMem cvdls_mem;
  realtype *bd;

  cvdls_mem = (CVDlsMem) lmem;

  bd = N_VGetArrayPointer(b);

  BandGBTRS(M, pivots, bd);

  /* If CV_BDF, scale the correction to account for change in gamma */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }

  last_flag = CVDIRECT_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static void cvBandFree(CVodeMem cv_mem)
{
  CVDlsMem cvdls_mem;

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
 * CVBandB is a wraper around CVBand. It attaches the CVSBAND linear solver
 * to the backward problem memory block.
 */

int CVBandB(void *cvb_mem, int nB, int mupperB, int mlowerB)
{
  CVodeBMem cvB_mem;
  void *cvode_mem;
  CVDlsMemB cvdlsB_mem;
  int flag;

  if (cvb_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_ADJMEM_NULL, "CVSBAND", "CVBandB", MSGD_CAMEM_NULL);
    return(CVDIRECT_ADJMEM_NULL);
  }
  cvB_mem = (CVodeBMem) cvb_mem;

  cvode_mem = (void *) (cvB_mem->cv_mem);

  /* Get memory for CVDlsMemRecB */
  cvdlsB_mem = (CVDlsMemB) malloc(sizeof(CVDlsMemRecB));
  if (cvdlsB_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_MEM_FAIL, "CVSBAND", "CVBandB", MSGD_MEM_FAIL);
    return(CVDIRECT_MEM_FAIL);
  }

  /* set matrix type */
  cvdlsB_mem->d_typeB = SUNDIALS_BAND;

  /* initialize Jacobian function */
  cvdlsB_mem->d_bjacB = NULL;
  cvdlsB_mem->d_jac_dataB = NULL;

  /* attach lmemB and lfreeB */
  cvB_mem->cv_lmem = cvdlsB_mem;
  cvB_mem->cv_lfree = cvBandFreeB;

  flag = CVBand(cvode_mem, nB, mupperB, mlowerB);

  if (flag != CVDIRECT_SUCCESS) {
    free(cvdlsB_mem);
    cvdlsB_mem = NULL;
  }

  return(flag);
}

/*
 * cvBandFreeB frees the memory associated with the CVSBAND linear
 * solver for backward integration.
 */

static void cvBandFreeB(CVodeBMem cvB_mem)
{
  CVDlsMemB cvdlsB_mem;

  cvdlsB_mem = (CVDlsMemB) (cvB_mem->cv_lmem);

  free(cvdlsB_mem);
}

