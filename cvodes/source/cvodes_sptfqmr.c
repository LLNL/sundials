/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-11 21:13:51 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the CVSPTFQMR linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_sptfqmr_impl.h"
#include "cvodes_impl.h"
#include "sundials_math.h"

/* Other Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* CVSPTFQMR linit, lsetup, lsolve, and lfree routines */

static int CVSptfqmrInit(CVodeMem cv_mem);

static int CVSptfqmrSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
			  N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
			  N_Vector vtemp2, N_Vector vtemp3);

static int CVSptfqmrSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
			  N_Vector ynow, N_Vector fnow);

static void CVSptfqmrFree(CVodeMem cv_mem);

/* CVSPTFQMR Atimes and PSolve routines called by generic SPTFQMR solver */

static int CVSptfqmrAtimes(void *cv_mem, N_Vector v, N_Vector z);

static int CVSptfqmrPSolve(void *cv_mem, N_Vector r, N_Vector z, int lr);

/* CVSPTFQMR difference quotient routine for J*v */

static int CVSptfqmrDQJtimes(N_Vector v, N_Vector Jv, realtype t,
			     N_Vector y, N_Vector fy, void *jac_data,
			     N_Vector work);

/* Readability Replacements */

#define lrw1         (cv_mem->cv_lrw1)
#define liw1         (cv_mem->cv_liw1)
#define tq           (cv_mem->cv_tq)
#define nst          (cv_mem->cv_nst)
#define tn           (cv_mem->cv_tn)
#define gamma        (cv_mem->cv_gamma)
#define gammap       (cv_mem->cv_gammap)
#define f            (cv_mem->cv_f)
#define f_data       (cv_mem->cv_f_data)
#define ewt          (cv_mem->cv_ewt)
#define errfp        (cv_mem->cv_errfp)
#define mnewt        (cv_mem->cv_mnewt)
#define linit        (cv_mem->cv_linit)
#define lsetup       (cv_mem->cv_lsetup)
#define lsolve       (cv_mem->cv_lsolve)
#define lfree        (cv_mem->cv_lfree)
#define lmem         (cv_mem->cv_lmem)
#define vec_tmpl     (cv_mem->cv_tempv)
#define setupNonNull (cv_mem->cv_setupNonNull)

#define sqrtN       (cvsptfqmr_mem->q_sqrtN)   
#define ytemp       (cvsptfqmr_mem->q_ytemp)
#define x           (cvsptfqmr_mem->q_x)
#define ycur        (cvsptfqmr_mem->q_ycur)
#define fcur        (cvsptfqmr_mem->q_fcur)
#define delta       (cvsptfqmr_mem->q_delta)
#define deltar      (cvsptfqmr_mem->q_deltar)
#define npe         (cvsptfqmr_mem->q_npe)
#define nli         (cvsptfqmr_mem->q_nli)
#define nps         (cvsptfqmr_mem->q_nps)
#define ncfl        (cvsptfqmr_mem->q_ncfl)
#define nstlpre     (cvsptfqmr_mem->q_nstlpre)
#define njtimes     (cvsptfqmr_mem->q_njtimes)
#define nfeSQ       (cvsptfqmr_mem->q_nfeSQ)
#define sptfqmr_mem (cvsptfqmr_mem->q_sptfqmr_mem)
#define last_flag   (cvsptfqmr_mem->q_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmr
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the Sptfqmr linear solver module. CVSptfqmr first
 * calls the existing lfree routine if this is not NULL. It then sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be CVSptfqmrInit, CVSptfqmrSetup, CVSptfqmrSolve, and CVSptfqmrFree,
 * respectively. It allocates memory for a structure of type
 * CVSptfqmrMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure. It sets setupNonNull in (*cvode_mem),
 * and sets the following fields in the CVSptfqmrMemRec structure:
 *
 *   q_pretype   = pretype
 *   q_maxl      = CVSPTFQMR_MAXL  if maxl <= 0
 *               = maxl            if maxl >  0
 *   q_delt      = CVSPTFQMR_DELT
 *   q_P_data    = NULL
 *   q_pset      = NULL
 *   q_psolve    = NULL
 *   q_jtimes    = CVSptfqmrDQJtimes
 *   q_j_data    = cvode_mem
 *   q_last_flag = CVSPTFQMR_SUCCESS
 *
 * Finally, CVSptfqmr allocates memory for ytemp and x, and calls
 * SptfqmrMalloc to allocate memory for the Sptfqmr solver.
 * -----------------------------------------------------------------
 */

int CVSptfqmr(void *cvode_mem, int pretype, int maxl)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Check if N_VDotProd is present */
  if (vec_tmpl->ops->nvdotprod == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_BAD_NVECTOR);
    return(CVSPTFQMR_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  linit  = CVSptfqmrInit;
  lsetup = CVSptfqmrSetup;
  lsolve = CVSptfqmrSolve;
  lfree  = CVSptfqmrFree;

  /* Get memory for CVSptfqmrMemRec */
  cvsptfqmr_mem = (CVSptfqmrMem) malloc(sizeof(CVSptfqmrMemRec));
  if (cvsptfqmr_mem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_MEM_FAIL);
    return(CVSPTFQMR_MEM_FAIL);
  }

  /* Set Sptfqmr parameters that have been passed in call sequence */
  cvsptfqmr_mem->q_pretype = pretype;
  mxl = cvsptfqmr_mem->q_maxl = (maxl <= 0) ? CVSPTFQMR_MAXL : maxl;

  /* Set default values for the rest of the Sptfqmr parameters */
  cvsptfqmr_mem->q_delt      = CVSPTFQMR_DELT;
  cvsptfqmr_mem->q_P_data    = NULL;
  cvsptfqmr_mem->q_pset      = NULL;
  cvsptfqmr_mem->q_psolve    = NULL;
  cvsptfqmr_mem->q_jtimes    = CVSptfqmrDQJtimes;
  cvsptfqmr_mem->q_j_data    = cvode_mem;
  cvsptfqmr_mem->q_last_flag = CVSPTFQMR_SUCCESS;

  setupNonNull = FALSE;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    if (errfp != NULL) 
      fprintf(errfp, MSGTFQMR_BAD_PRETYPE);
    return(CVSPTFQMR_ILL_INPUT);
  }

  /* Allocate memory for ytemp and x */
  ytemp = N_VClone(vec_tmpl);
  if (ytemp == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_MEM_FAIL);
    return(CVSPTFQMR_MEM_FAIL);
  }
  x = N_VClone(vec_tmpl);
  if (x == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_MEM_FAIL);
    N_VDestroy(ytemp);
    return(CVSPTFQMR_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  sqrtN = RSqrt(N_VDotProd(ytemp, ytemp));

  /* Call SptfqmrMalloc to allocate workspace for Sptfqmr */
  sptfqmr_mem = SptfqmrMalloc(mxl, vec_tmpl);
  if (sptfqmr_mem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(x);
    return(CVSPTFQMR_MEM_FAIL);
  }
  
  /* Attach linear solver memory to integrator memory */
  lmem = cvsptfqmr_mem;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrSetPrecType
 * -----------------------------------------------------------------
 */

int CVSptfqmrSetPrecType(void *cvode_mem, int pretype)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    if (errfp != NULL) 
      fprintf(errfp, MSGTFQMR_SET_BAD_PRETYPE);
    return(CVSPTFQMR_ILL_INPUT);
  }

  cvsptfqmr_mem->q_pretype = pretype;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrSetMaxl
 * -----------------------------------------------------------------
 */

int CVSptfqmrSetMaxl(void *cvode_mem, int maxl)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  mxl = (maxl <= 0) ? CVSPTFQMR_MAXL : maxl;
  cvsptfqmr_mem->q_maxl = mxl;
  sptfqmr_mem->l_max  = mxl;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrSetDelt
 * -----------------------------------------------------------------
 */

int CVSptfqmrSetDelt(void *cvode_mem, realtype delt)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  /* Check for legal delt */
  if (delt < ZERO) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SET_BAD_DELT);
    return(CVSPTFQMR_ILL_INPUT);
  }

  cvsptfqmr_mem->q_delt = (delt == ZERO) ? CVSPTFQMR_DELT : delt;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrSetPreconditioner
 * -----------------------------------------------------------------
 */

int CVSptfqmrSetPreconditioner(void *cvode_mem,
			       CVSpilsPrecSetupFn pset,
			       CVSpilsPrecSolveFn psolve,
			       void *P_data)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  cvsptfqmr_mem->q_pset = pset;
  cvsptfqmr_mem->q_psolve = psolve;
  if (psolve != NULL) cvsptfqmr_mem->q_P_data = P_data;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrSetJacTimesVecFn
 * -----------------------------------------------------------------
 */

int CVSptfqmrSetJacTimesVecFn(void *cvode_mem, 
			      CVSpilsJacTimesVecFn jtimes,
			      void *jac_data)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  cvsptfqmr_mem->q_jtimes = jtimes;
  if (jtimes != NULL) cvsptfqmr_mem->q_j_data = jac_data;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrGetWorkSpace
 * -----------------------------------------------------------------
 */

int CVSptfqmrGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS)
{
  CVodeMem cv_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }

#ifdef DEBUG
  *lenrwLS = lrw1*12;
  *leniwLS = liw1*12;
#else
  *lenrwLS = lrw1*11;
  *leniwLS = liw1*11;
#endif

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrGetNumPrecEvals
 * -----------------------------------------------------------------
 */

int CVSptfqmrGetNumPrecEvals(void *cvode_mem, long int *npevals)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  *npevals = npe;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrGetNumPrecSolves
 * -----------------------------------------------------------------
 */

int CVSptfqmrGetNumPrecSolves(void *cvode_mem, long int *npsolves)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  *npsolves = nps;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrGetNumLinIters
 * -----------------------------------------------------------------
 */

int CVSptfqmrGetNumLinIters(void *cvode_mem, long int *nliters)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  *nliters = nli;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrGetNumConvFails
 * -----------------------------------------------------------------
 */

int CVSptfqmrGetNumConvFails(void *cvode_mem, long int *nlcfails)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  *nlcfails = ncfl;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrGetNumJtimesEvals
 * -----------------------------------------------------------------
 */

int CVSptfqmrGetNumJtimesEvals(void *cvode_mem, long int *njvevals)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  *njvevals = njtimes;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrGetNumRhsEvals
 * -----------------------------------------------------------------
 */

int CVSptfqmrGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  *nfevalsLS = nfeSQ;

  return(CVSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrGetLastFlag
 * -----------------------------------------------------------------
 */

int CVSptfqmrGetLastFlag(void *cvode_mem, int *flag)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_CVMEM_NULL);
    return(CVSPTFQMR_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(CVSPTFQMR_LMEM_NULL);
  }
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  *flag = last_flag;

  return(CVSPTFQMR_SUCCESS);
}


/* Additional readability replacements */

#define pretype (cvsptfqmr_mem->q_pretype)
#define delt    (cvsptfqmr_mem->q_delt)
#define maxl    (cvsptfqmr_mem->q_maxl)
#define psolve  (cvsptfqmr_mem->q_psolve)
#define pset    (cvsptfqmr_mem->q_pset)
#define P_data  (cvsptfqmr_mem->q_P_data)
#define jtimes  (cvsptfqmr_mem->q_jtimes)
#define j_data  (cvsptfqmr_mem->q_j_data)

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the Sptfqmr
 * linear solver.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrInit(CVodeMem cv_mem)
{
  CVSptfqmrMem cvsptfqmr_mem;
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  /* Initialize counters */
  npe = nli = nps = ncfl = nstlpre = 0;
  njtimes = nfeSQ = 0;

  /* Check for legal combination pretype - psolve */
  if ((pretype != PREC_NONE) && (psolve == NULL)) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_PSOLVE_REQ);
    last_flag = -1;
    return(-1);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning
     (pretype != PREC_NONE)  and there is a preconditioning
     setup phase (pset != NULL) */
  setupNonNull = (pretype != PREC_NONE) && (pset != NULL);

  /* If jtimes is NULL at this time, set it to DQ */
  if (jtimes == NULL) {
    jtimes = CVSptfqmrDQJtimes;
    j_data = cv_mem;
  }

  last_flag = CVSPTFQMR_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the Sptfqmr linear solver.
 * It makes a decision as to whether or not to signal for reevaluation
 * of Jacobian data in the pset routine, based on various state
 * variables, then it calls pset. If we signal for reevaluation,
 * then we reset jcur = *jcurPtr to TRUE, regardless of the pset output.
 * In any case, if jcur == TRUE, we increment npe and save nst in nstlpre.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
			  N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
			  N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  ier;
  CVSptfqmrMem cvsptfqmr_mem;

  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlpre + CVSPTFQMR_MSBPRE) ||
      ((convfail == CV_FAIL_BAD_J) && (dgamma < CVSPTFQMR_DGMAX)) ||
      (convfail == CV_FAIL_OTHER);
  *jcurPtr = jbad;
  jok = !jbad;

  /* Call pset routine and possibly reset jcur */
  ier = pset(tn, ypred, fpred, jok, jcurPtr, gamma, P_data, 
             vtemp1, vtemp2, vtemp3);
  if (jbad) *jcurPtr = TRUE;

  /* If jcur = TRUE, increment npe and save nst value */
  if (*jcurPtr) {
    npe++;
    nstlpre = nst;
  }

  /* Return the same value ier that pset returned */
  last_flag = ier;
  return(ier);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic solver SptfqmrSolve
 * for the solution of the linear system Ax = b with the SPTFQMR method.
 * The solution x is returned in the vector b.
 *
 * If the WRMS norm of b is small, we return x = b (if this is the first
 * Newton iteration) or x = 0 (if a later Newton iteration).
 *
 * Otherwise, we set the tolerance parameter and initial guess (x = 0),
 * call SptfqmrSolve, and copy the solution x into b. The x-scaling and
 * b-scaling arrays are both equal to weight.
 *
 * The counters nli, nps, and ncfl are incremented, and the return value
 * is set according to the success of SptfqmrSolve. The success flag is
 * returned if SptfqmrSolve converged, or if this is the first Newton
 * iteration and the residual norm was reduced below its initial value.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
			  N_Vector ynow, N_Vector fnow)
{
  realtype bnorm, res_norm;
  CVSptfqmrMem cvsptfqmr_mem;
  int nli_inc, nps_inc, ier;
  
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  /* Test norm(b); if small, return x = 0 or x = b */
  deltar = delt * tq[4]; 

  bnorm = N_VWrmsNorm(b, weight);
  if (bnorm <= deltar) {
    if (mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
  ycur = ynow;
  fcur = fnow;

  /* Set inputs delta and initial guess x = 0 to SptfqmrSolve */  
  delta = deltar * sqrtN;
  N_VConst(ZERO, x);
  
  /* Call SptfqmrSolve and copy x to b */
  ier = SptfqmrSolve(sptfqmr_mem, cv_mem, x, b, pretype, delta,
                   cv_mem, weight, weight, CVSptfqmrAtimes, CVSptfqmrPSolve,
                   &res_norm, &nli_inc, &nps_inc);

  N_VScale(ONE, x, b);
  
  /* Increment counters nli, nps, and ncfl */
  nli += nli_inc;
  nps += nps_inc;
  if (ier != 0) ncfl++;

  /* Set return value to -1, 0, or 1 */
  last_flag = ier;

  if (ier < 0) return(-1);  

  if ((ier == SPTFQMR_SUCCESS) || 
      ((ier == SPTFQMR_RES_REDUCED) && (mnewt == 0)))
    return(0);

  return(1);  
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the Sptfqmr linear solver.
 * -----------------------------------------------------------------
 */

static void CVSptfqmrFree(CVodeMem cv_mem)
{
  CVSptfqmrMem cvsptfqmr_mem;

  cvsptfqmr_mem = (CVSptfqmrMem) lmem;
  
  N_VDestroy(ytemp);
  N_VDestroy(x);
  SptfqmrFree(sptfqmr_mem);
  free(cvsptfqmr_mem);

  return;
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrAtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = Mv, where
 * M = I - gamma*J. The product J*v is obtained by calling the jtimes
 * routine. It is then scaled by -gamma and added to v to obtain M*v.
 * The return value is the same as the value returned by jtimes --
 * 0 if successful, nonzero otherwise.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrAtimes(void *cvode_mem, N_Vector v, N_Vector z)
{
  CVodeMem   cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;
  int jtflag;

  cv_mem = (CVodeMem) cvode_mem;
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  jtflag = jtimes(v, z, tn, ycur, fcur, j_data, ytemp);
  njtimes++;
  if (jtflag != 0) return(jtflag);

  N_VLinearSum(ONE, v, -gamma, z, z);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic SptfqmrSolve routine and
 * the user's psolve routine. It passes to psolve all required state
 * information from cvode_mem. Its return value is the same as that
 * returned by psolve. Note that the generic SPTFQMR solver guarantees
 * that CVSptfqmrPSolve will not be called in the case in which
 * preconditioning is not done. This is the only case in which the
 * user's psolve routine is allowed to be NULL.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrPSolve(void *cvode_mem, N_Vector r, N_Vector z, int lr)
{
  CVodeMem   cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;
  int ier;

  cv_mem = (CVodeMem) cvode_mem;
  cvsptfqmr_mem = (CVSptfqmrMem)lmem;

  ier = psolve(tn, ycur, fcur, r, z, gamma, delta, lr, P_data, ytemp);
  /* This call is counted in nps within the CVSptfqmrSolve routine */

  return(ier);     
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrDQJtimes
 * -----------------------------------------------------------------
 * This routine generates a difference quotient approximation to
 * the Jacobian times vector f_y(t,y) * v. The approximation is
 * Jv = vnrm[f(y + v/vnrm) - f(y)], where vnrm = (WRMS norm of v) is
 * input, i.e. the WRMS norm of v/vnrm is 1.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrDQJtimes(N_Vector v, N_Vector Jv, realtype t, 
			     N_Vector y, N_Vector fy,
			     void *jac_data, N_Vector work)
{
  CVodeMem cv_mem;
  CVSptfqmrMem cvsptfqmr_mem;
  realtype vnrm;

  /* jac_data is cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvsptfqmr_mem = (CVSptfqmrMem) lmem;

  /* Evaluate norm of v */
  vnrm = N_VWrmsNorm(v, ewt);

  /* Set work = y + (1/vnrm) v */
  N_VLinearSum(ONE/vnrm, v, ONE, y, work);

  /* Set Jv = f(tn, work) */
  f(t, work, Jv, f_data); 
  nfeSQ++;

  /* Replace Jv by vnrm*(Jv - fy) */
  N_VLinearSum(ONE, Jv, -ONE, fy, Jv);
  N_VScale(vnrm, Jv, Jv);

  return(0);
}
