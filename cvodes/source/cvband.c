/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2004-10-08 16:35:27 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the CVBAND linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvband_impl.h"
#include "cvodes_impl.h"

#include "sundialsmath.h"

/* Error Messages */

#define CVBAND                "CVBand-- "
  
#define MSG_MEM_FAIL          CVBAND "A memory request failed.\n\n"

#define MSG_BAD_SIZES_1       CVBAND "Illegal bandwidth parameter(s) "
#define MSG_BAD_SIZES_2       "ml = %ld, mu = %ld.\n"
#define MSG_BAD_SIZES_3       "Must have 0 <=  ml, mu <= N-1=%ld.\n\n"
#define MSG_BAD_SIZES         MSG_BAD_SIZES_1 MSG_BAD_SIZES_2 MSG_BAD_SIZES_3

#define MSG_BAD_NVECTOR       CVBAND "A required vector operation is not implemented.\n\n"

#define MSG_CVMEM_NULL        CVBAND "Integrator memory is NULL.\n\n"

#define MSG_SETGET_CVMEM_NULL "CVBandSet*/CVBandGet*-- Integrator memory is NULL. \n\n"

#define MSG_SETGET_LMEM_NULL  "CVBandSet*/CVBandGet*-- cvband memory is NULL. \n\n"

/* Other Constants */

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* CVBAND linit, lsetup, lsolve, and lfree routines */

static int CVBandInit(CVodeMem cv_mem);

static int CVBandSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3);

static int CVBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur);

static void CVBandFree(CVodeMem cv_mem);

/* CVBAND DQJac routine */

static void CVBandDQJac(long int n, long int mupper, long int mlower,
                        BandMat J, realtype t,
                        N_Vector y, N_Vector fy, void *jac_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

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
#define errfp     (cv_mem->cv_errfp)
#define linit     (cv_mem->cv_linit)
#define lsetup    (cv_mem->cv_lsetup)
#define lsolve    (cv_mem->cv_lsolve)
#define lfree     (cv_mem->cv_lfree)
#define lmem      (cv_mem->cv_lmem)
#define vec_tmpl      (cv_mem->cv_tempv)
#define setupNonNull  (cv_mem->cv_setupNonNull)

#define n          (cvband_mem->b_n)
#define jac        (cvband_mem->b_jac)
#define M          (cvband_mem->b_M)
#define mu         (cvband_mem->b_mu)
#define ml         (cvband_mem->b_ml)
#define storage_mu (cvband_mem->b_storage_mu)
#define pivots     (cvband_mem->b_pivots)
#define savedJ     (cvband_mem->b_savedJ)
#define nstlj      (cvband_mem->b_nstlj)
#define nje        (cvband_mem->b_nje)
#define nfeB       (cvband_mem->b_nfeB)
#define J_data     (cvband_mem->b_J_data)
#define last_flag  (cvband_mem->b_last_flag)

/*
 * -----------------------------------------------------------------
 * CVBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module.  CVBand first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 * to be CVBandInit, CVBandSetup, CVBandSolve, and CVBandFree,
 * respectively.  It allocates memory for a structure of type
 * CVBandMemRec and sets the cv_lmem field in (*cvode_mem) to the
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
                  
int CVBand(void *cvode_mem, long int N,
           long int mupper, long int mlower)
{
  CVodeMem cv_mem;
  CVBandMem cvband_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_CVMEM_NULL);
    return(CVBAND_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_NVECTOR);
    return(CVBAND_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */  
  linit  = CVBandInit;
  lsetup = CVBandSetup;
  lsolve = CVBandSolve;
  lfree  = CVBandFree;
  
  /* Get memory for CVBandMemRec */
  cvband_mem = (CVBandMem) malloc(sizeof(CVBandMemRec));
  if (cvband_mem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    return(CVBAND_MEM_FAIL);
  }
  
  /* Set default Jacobian routine and Jacobian data */
  jac = CVBandDQJac;
  J_data = cvode_mem;
  last_flag = CVBAND_SUCCESS;

  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in cvband_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_SIZES, ml, mu, N-1);
    return(CVBAND_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  storage_mu = MIN(N-1, mu + ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  M = BandAllocMat(N, mu, ml, storage_mu);
  if (M == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    return(CVBAND_MEM_FAIL);
  }
  savedJ = BandAllocMat(N, mu, ml, mu);
  if (savedJ == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    BandFreeMat(M);
    return(CVBAND_MEM_FAIL);
  }
  pivots = BandAllocPiv(N);
  if (pivots == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    BandFreeMat(M);
    BandFreeMat(savedJ);
    return(CVBAND_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cvband_mem;

  return(CVBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVBandSetJacFn
 * -----------------------------------------------------------------
 */

int CVBandSetJacFn(void *cvode_mem, CVBandJacFn bjac)
{
  CVodeMem cv_mem;
  CVBandMem cvband_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_SETGET_CVMEM_NULL);
    return(CVBAND_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(CVBAND_LMEM_NULL);
  }
  cvband_mem = (CVBandMem) lmem;

  jac = bjac;

  return(CVBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVBandSetJacData
 * -----------------------------------------------------------------
 */

int CVBandSetJacData(void *cvode_mem, void *jac_data)
{
  CVodeMem cv_mem;
  CVBandMem cvband_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_SETGET_CVMEM_NULL);
    return(CVBAND_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(CVBAND_LMEM_NULL);
  }
  cvband_mem = (CVBandMem) lmem;

  J_data = jac_data;

  return(CVBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVBandGetWorkSpace
 * -----------------------------------------------------------------
 */

int CVBandGetWorkSpace(void *cvode_mem, long int *lenrwB, long int *leniwB)
{
  CVodeMem cv_mem;
  CVBandMem cvband_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_SETGET_CVMEM_NULL);
    return(CVBAND_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(CVBAND_LMEM_NULL);
  }
  cvband_mem = (CVBandMem) lmem;

  *lenrwB = n*(storage_mu + mu + 2*ml + 2);
  *leniwB = n;

  return(CVBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVBandGetNumJacEvals
 * -----------------------------------------------------------------
 */

int CVBandGetNumJacEvals(void *cvode_mem, long int *njevalsB)
{
  CVodeMem cv_mem;
  CVBandMem cvband_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_SETGET_CVMEM_NULL);
    return(CVBAND_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(CVBAND_LMEM_NULL);
  }
  cvband_mem = (CVBandMem) lmem;

  *njevalsB = nje;

  return(CVBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVBandGetNumRhsEvals
 * -----------------------------------------------------------------
 */

int CVBandGetNumRhsEvals(void *cvode_mem, long int *nfevalsB)
{
  CVodeMem cv_mem;
  CVBandMem cvband_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_SETGET_CVMEM_NULL);
    return(CVBAND_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(CVBAND_LMEM_NULL);
  }
  cvband_mem = (CVBandMem) lmem;

  *nfevalsB = nfeB;

  return(CVBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVBandGetLastFlag
 * -----------------------------------------------------------------
 */

int CVBandGetLastFlag(void *cvode_mem, int *flag)
{
  CVodeMem cv_mem;
  CVBandMem cvband_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_SETGET_CVMEM_NULL);
    return(CVBAND_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(CVBAND_LMEM_NULL);
  }
  cvband_mem = (CVBandMem) lmem;

  *flag = last_flag;

  return(CVBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int CVBandInit(CVodeMem cv_mem)
{
  CVBandMem cvband_mem;

  cvband_mem = (CVBandMem) lmem;

  nje   = 0;
  nfeB  = 0;
  nstlj = 0;

  if (jac == NULL) {
    jac = CVBandDQJac;
    J_data = cv_mem;
  }

  last_flag = CVBAND_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * CVBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the band LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int CVBandSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  long int ier;
  CVBandMem cvband_mem;
  
  cvband_mem = (CVBandMem) lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */

  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + CVB_MSBJ) ||
         ((convfail == CV_FAIL_BAD_J) && (dgamma < CVB_DGMAX)) ||
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
    jac(n, mu, ml, M, tn, ypred, fpred, J_data, vtemp1, vtemp2, vtemp3);
    BandCopy(M, savedJ, mu, ml);
  }
  
  /* Scale and add I to get M = I - gamma*J */
  BandScale(-gamma, M);
  BandAddI(M);

  /* Do LU factorization of M */
  ier = BandFactor(M, pivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) {
    last_flag = ier;
    return(1);
  }
  last_flag = CVBAND_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * CVBandSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int CVBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur)
{
  CVBandMem cvband_mem;
  realtype *bd;

  cvband_mem = (CVBandMem) lmem;

  bd = N_VGetArrayPointer(b);

  BandBacksolve(M, pivots, bd);

  /* If CV_BDF, scale the correction to account for change in gamma */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }

  last_flag = CVBAND_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * CVBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static void CVBandFree(CVodeMem cv_mem)
{
  CVBandMem cvband_mem;

  cvband_mem = (CVBandMem) lmem;

  BandFreeMat(M);
  BandFreeMat(savedJ);
  BandFreePiv(pivots);
  free(cvband_mem);
}

/*
 * -----------------------------------------------------------------
 * CVBandDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y).  It assumes that a band matrix of type
 * BandMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

static void CVBandDQJac(long int N, long int mupper, long int mlower,
                        BandMat J, realtype t,
                        N_Vector y, N_Vector fy, void *jac_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype fnorm, minInc, inc, inc_inv, srur;
  N_Vector ftemp, ytemp;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;

  CVodeMem cv_mem;
  CVBandMem cvband_mem;

  /* jac_dat points to cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvband_mem = (CVBandMem) lmem;

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
  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */

    f(tn, ytemp, ftemp, f_data);

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(J,j);
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }
  
  /* Increment counter nfeB */
  nfeB += ngroups;
}
