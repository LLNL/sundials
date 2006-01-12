/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-01-12 20:24:07 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the CVSPBCG linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_spbcgs_impl.h"
#include "cvodes_spils_impl.h"
#include "cvodes_impl.h"
#include "cvodea_impl.h"

#include "sundials_math.h"

/* Other Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* CVSPBCG linit, lsetup, lsolve, and lfree routines */

static int CVSpbcgInit(CVodeMem cv_mem);

static int CVSpbcgSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                        N_Vector vtemp2, N_Vector vtemp3);

static int CVSpbcgSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                        N_Vector ynow, N_Vector fnow);

static void CVSpbcgFree(CVodeMem cv_mem);

/* CVSPBCG Atimes and PSolve routines called by generic SPBCG solver */

static int CVSpbcgAtimes(void *cv_mem, N_Vector v, N_Vector z);

static int CVSpbcgPSolve(void *cv_mem, N_Vector r, N_Vector z, int lr);

/* CVSPBCG difference quotient routine for J*v */

static int CVSpbcgDQJtimes(N_Vector v, N_Vector Jv, realtype t,
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

#define sqrtN     (cvspbcg_mem->b_sqrtN)   
#define ytemp     (cvspbcg_mem->b_ytemp)
#define x         (cvspbcg_mem->b_x)
#define ycur      (cvspbcg_mem->b_ycur)
#define fcur      (cvspbcg_mem->b_fcur)
#define delta     (cvspbcg_mem->b_delta)
#define deltar    (cvspbcg_mem->b_deltar)
#define npe       (cvspbcg_mem->b_npe)
#define nli       (cvspbcg_mem->b_nli)
#define nps       (cvspbcg_mem->b_nps)
#define ncfl      (cvspbcg_mem->b_ncfl)
#define nstlpre   (cvspbcg_mem->b_nstlpre)
#define njtimes   (cvspbcg_mem->b_njtimes)
#define nfeSB     (cvspbcg_mem->b_nfeSB)
#define spbcg_mem (cvspbcg_mem->b_spbcg_mem)
#define last_flag (cvspbcg_mem->b_last_flag)

/* 
 * =================================================================
 * PART I - forward problems
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcg
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the Spbcg linear solver module. CVSpbcg first
 * calls the existing lfree routine if this is not NULL. It then sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be CVSpbcgInit, CVSpbcgSetup, CVSpbcgSolve, and CVSpbcgFree,
 * respectively. It allocates memory for a structure of type
 * CVSpbcgMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure. It sets setupNonNull in (*cvode_mem),
 * and sets the following fields in the CVSpbcgMemRec structure:
 *
 *   b_pretype   = pretype
 *   b_maxl      = CVSPBCG_MAXL  if maxl <= 0
 *               = maxl          if maxl >  0
 *   b_delt      = CVSPBCG_DELT
 *   b_P_data    = NULL
 *   b_pset      = NULL
 *   b_psolve    = NULL
 *   b_jtimes    = CVSpbcgDQJtimes
 *   b_j_data    = cvode_mem
 *   b_last_flag = CVSPBCG_SUCCESS
 *
 * Finally, CVSpbcg allocates memory for ytemp and x, and calls
 * SpbcgMalloc to allocate memory for the Spbcg solver.
 * -----------------------------------------------------------------
 */

int CVSpbcg(void *cvode_mem, int pretype, int maxl)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Check if N_VDotProd is present */
  if (vec_tmpl->ops->nvdotprod == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_BAD_NVECTOR);
    return(CVSPBCG_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  linit  = CVSpbcgInit;
  lsetup = CVSpbcgSetup;
  lsolve = CVSpbcgSolve;
  lfree  = CVSpbcgFree;

  /* Get memory for CVSpbcgMemRec */
  cvspbcg_mem = (CVSpbcgMem) malloc(sizeof(CVSpbcgMemRec));
  if (cvspbcg_mem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_MEM_FAIL);
    return(CVSPBCG_MEM_FAIL);
  }

  /* Set Spbcg parameters that have been passed in call sequence */
  cvspbcg_mem->b_pretype = pretype;
  mxl = cvspbcg_mem->b_maxl = (maxl <= 0) ? CVSPBCG_MAXL : maxl;

  /* Set default values for the rest of the Spbcg parameters */
  cvspbcg_mem->b_delt      = CVSPBCG_DELT;
  cvspbcg_mem->b_P_data    = NULL;
  cvspbcg_mem->b_pset      = NULL;
  cvspbcg_mem->b_psolve    = NULL;
  cvspbcg_mem->b_jtimes    = CVSpbcgDQJtimes;
  cvspbcg_mem->b_j_data    = cvode_mem;
  cvspbcg_mem->b_last_flag = CVSPBCG_SUCCESS;

  setupNonNull = FALSE;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    if (errfp != NULL) 
      fprintf(errfp, MSGBCG_BAD_PRETYPE);
    return(CVSPBCG_ILL_INPUT);
  }

  /* Allocate memory for ytemp and x */
  ytemp = N_VClone(vec_tmpl);
  if (ytemp == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_MEM_FAIL);
    return(CVSPBCG_MEM_FAIL);
  }
  x = N_VClone(vec_tmpl);
  if (x == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_MEM_FAIL);
    N_VDestroy(ytemp);
    return(CVSPBCG_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  sqrtN = RSqrt(N_VDotProd(ytemp, ytemp));

  /* Call SpbcgMalloc to allocate workspace for Spbcg */
  spbcg_mem = SpbcgMalloc(mxl, vec_tmpl);
  if (spbcg_mem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(x);
    return(CVSPBCG_MEM_FAIL);
  }
  
  /* Attach linear solver memory to integrator memory */
  lmem = cvspbcg_mem;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgSetPrecType
 * -----------------------------------------------------------------
 */

int CVSpbcgSetPrecType(void *cvode_mem, int pretype)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    if (errfp != NULL) 
      fprintf(errfp, MSGBCG_SET_BAD_PRETYPE);
    return(CVSPBCG_ILL_INPUT);
  }

  cvspbcg_mem->b_pretype = pretype;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgSetMaxl
 * -----------------------------------------------------------------
 */

int CVSpbcgSetMaxl(void *cvode_mem, int maxl)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  mxl = (maxl <= 0) ? CVSPBCG_MAXL : maxl;
  cvspbcg_mem->b_maxl = mxl;
  spbcg_mem->l_max  = mxl;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgSetDelt
 * -----------------------------------------------------------------
 */

int CVSpbcgSetDelt(void *cvode_mem, realtype delt)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  /* Check for legal delt */
  if (delt < ZERO) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SET_BAD_DELT);
    return(CVSPBCG_ILL_INPUT);
  }

  cvspbcg_mem->b_delt = (delt == ZERO) ? CVSPBCG_DELT : delt;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgSetPreconditioner
 * -----------------------------------------------------------------
 */

int CVSpbcgSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn pset,
                             CVSpilsPrecSolveFn psolve, void *P_data)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  cvspbcg_mem->b_pset = pset;
  cvspbcg_mem->b_psolve = psolve;
  if (psolve != NULL) cvspbcg_mem->b_P_data = P_data;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgSetJacTimesVecFn
 * -----------------------------------------------------------------
 */

int CVSpbcgSetJacTimesVecFn(void *cvode_mem, 
                            CVSpilsJacTimesVecFn jtimes, void *jac_data)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  cvspbcg_mem->b_jtimes = jtimes;
  if (jtimes != NULL) cvspbcg_mem->b_j_data = jac_data;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgGetWorkSpace
 * -----------------------------------------------------------------
 */

int CVSpbcgGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS)
{
  CVodeMem cv_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }

  *lenrwLS = lrw1 * 9;
  *leniwLS = liw1 * 9;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgGetNumPrecEvals
 * -----------------------------------------------------------------
 */

int CVSpbcgGetNumPrecEvals(void *cvode_mem, long int *npevals)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  *npevals = npe;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgGetNumPrecSolves
 * -----------------------------------------------------------------
 */

int CVSpbcgGetNumPrecSolves(void *cvode_mem, long int *npsolves)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  *npsolves = nps;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgGetNumLinIters
 * -----------------------------------------------------------------
 */

int CVSpbcgGetNumLinIters(void *cvode_mem, long int *nliters)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  *nliters = nli;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgGetNumConvFails
 * -----------------------------------------------------------------
 */

int CVSpbcgGetNumConvFails(void *cvode_mem, long int *nlcfails)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  *nlcfails = ncfl;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgGetNumJtimesEvals
 * -----------------------------------------------------------------
 */

int CVSpbcgGetNumJtimesEvals(void *cvode_mem, long int *njvevals)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  *njvevals = njtimes;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgGetNumRhsEvals
 * -----------------------------------------------------------------
 */

int CVSpbcgGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  *nfevalsLS = nfeSB;

  return(CVSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgGetLastFlag
 * -----------------------------------------------------------------
 */

int CVSpbcgGetLastFlag(void *cvode_mem, int *flag)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_CVMEM_NULL);
    return(CVSPBCG_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(CVSPBCG_LMEM_NULL);
  }
  cvspbcg_mem = (CVSpbcgMem) lmem;

  *flag = last_flag;

  return(CVSPBCG_SUCCESS);
}


/* 
 * =================================================================
 * PART II - backward problems
 * =================================================================
 */

/* Additional readability replacements */

#define pset_B      (ca_mem->ca_psetB)
#define psolve_B    (ca_mem->ca_psolveB)
#define jtimes_B    (ca_mem->ca_jtimesB)
#define P_data_B    (ca_mem->ca_P_dataB)
#define jac_data_B  (ca_mem->ca_jac_dataB)

/*
 * CVSpbcgB and CVSpbcgSet*B
 *
 * Wrappers for the backward phase around the corresponding 
 * CVODES functions
 */

int CVSpbcgB(void *cvadj_mem, int pretypeB, int maxlB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CVSPBCG_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;
  
  flag = CVSpbcg(cvode_mem, pretypeB, maxlB);

  return(flag);
}

int CVSpbcgSetPrecTypeB(void *cvadj_mem, int pretypeB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CVSPBCG_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpbcgSetPrecType(cvode_mem, pretypeB);

  return(flag);
}

int CVSpbcgSetDeltB(void *cvadj_mem, realtype deltB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CVSPBCG_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpbcgSetDelt(cvode_mem,deltB);

  return(flag);
}

int CVSpbcgSetPreconditionerB(void *cvadj_mem, CVSpilsPrecSetupFnB psetB,
                              CVSpilsPrecSolveFnB psolveB, void *P_dataB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CVSPBCG_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  pset_B   = psetB;
  psolve_B = psolveB;
  P_data_B = P_dataB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpbcgSetPreconditioner(cvode_mem, 
                                  CVAspilsPrecSetup, CVAspilsPrecSolve, cvadj_mem);

  return(flag);
}

int CVSpbcgSetJacTimesVecFnB(void *cvadj_mem, 
                             CVSpilsJacTimesVecFnB jtimesB, void *jac_dataB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CVSPBCG_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  jtimes_B   = jtimesB;
  jac_data_B = jac_dataB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpbcgSetJacTimesVecFn(cvode_mem, CVAspilsJacTimesVec, cvadj_mem);

  return(flag);
}


/* 
 * =================================================================
 * PART III - private functions
 * =================================================================
 */



/* Additional readability replacements */

#define pretype (cvspbcg_mem->b_pretype)
#define delt    (cvspbcg_mem->b_delt)
#define maxl    (cvspbcg_mem->b_maxl)
#define psolve  (cvspbcg_mem->b_psolve)
#define pset    (cvspbcg_mem->b_pset)
#define P_data  (cvspbcg_mem->b_P_data)
#define jtimes  (cvspbcg_mem->b_jtimes)
#define j_data  (cvspbcg_mem->b_j_data)

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the Spbcg
 * linear solver.
 * -----------------------------------------------------------------
 */

static int CVSpbcgInit(CVodeMem cv_mem)
{
  CVSpbcgMem cvspbcg_mem;
  cvspbcg_mem = (CVSpbcgMem) lmem;

  /* Initialize counters */
  npe = nli = nps = ncfl = nstlpre = 0;
  njtimes = nfeSB = 0;

  /* Check for legal combination pretype - psolve */
  if ((pretype != PREC_NONE) && (psolve == NULL)) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_PSOLVE_REQ);
    last_flag = -1;
    return(-1);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning
     (pretype != PREC_NONE)  and there is a preconditioning
     setup phase (pset != NULL) */
  setupNonNull = (pretype != PREC_NONE) && (pset != NULL);

  /* If jtimes is NULL at this time, set it to DQ */
  if (jtimes == NULL) {
    jtimes = CVSpbcgDQJtimes;
    j_data = cv_mem;
  }

  last_flag = CVSPBCG_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the Spbcg linear solver.
 * It makes a decision as to whether or not to signal for reevaluation
 * of Jacobian data in the pset routine, based on various state
 * variables, then it calls pset. If we signal for reevaluation,
 * then we reset jcur = *jcurPtr to TRUE, regardless of the pset output.
 * In any case, if jcur == TRUE, we increment npe and save nst in nstlpre.
 * -----------------------------------------------------------------
 */

static int CVSpbcgSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                        N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  ier;
  CVSpbcgMem cvspbcg_mem;

  cvspbcg_mem = (CVSpbcgMem) lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlpre + CVSPBCG_MSBPRE) ||
    ((convfail == CV_FAIL_BAD_J) && (dgamma < CVSPBCG_DGMAX)) ||
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
 * Function : CVSpbcgSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic solver SpbcgSolve
 * for the solution of the linear system Ax = b with the SPBCG method.
 * The solution x is returned in the vector b.
 *
 * If the WRMS norm of b is small, we return x = b (if this is the first
 * Newton iteration) or x = 0 (if a later Newton iteration).
 *
 * Otherwise, we set the tolerance parameter and initial guess (x = 0),
 * call SpbcgSolve, and copy the solution x into b. The x-scaling and
 * b-scaling arrays are both equal to weight.
 *
 * The counters nli, nps, and ncfl are incremented, and the return value
 * is set according to the success of SpbcgSolve. The success flag is
 * returned if SpbcgSolve converged, or if this is the first Newton
 * iteration and the residual norm was reduced below its initial value.
 * -----------------------------------------------------------------
 */

static int CVSpbcgSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                        N_Vector ynow, N_Vector fnow)
{
  realtype bnorm, res_norm;
  CVSpbcgMem cvspbcg_mem;
  int nli_inc, nps_inc, ier;
  
  cvspbcg_mem = (CVSpbcgMem) lmem;

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

  /* Set inputs delta and initial guess x = 0 to SpbcgSolve */  
  delta = deltar * sqrtN;
  N_VConst(ZERO, x);
  
  /* Call SpbcgSolve and copy x to b */
  ier = SpbcgSolve(spbcg_mem, cv_mem, x, b, pretype, delta,
                   cv_mem, weight, weight, CVSpbcgAtimes, CVSpbcgPSolve,
                   &res_norm, &nli_inc, &nps_inc);

  N_VScale(ONE, x, b);
  
  /* Increment counters nli, nps, and ncfl */
  nli += nli_inc;
  nps += nps_inc;
  if (ier != 0) ncfl++;

  /* Set return value to -1, 0, or 1 */
  last_flag = ier;

  if (ier < 0) return(-1);  

  if ((ier == SPBCG_SUCCESS) || 
      ((ier == SPBCG_RES_REDUCED) && (mnewt == 0)))
    return(0);

  return(1);  
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the Spbcg linear solver.
 * -----------------------------------------------------------------
 */

static void CVSpbcgFree(CVodeMem cv_mem)
{
  CVSpbcgMem cvspbcg_mem;

  cvspbcg_mem = (CVSpbcgMem) lmem;
  
  N_VDestroy(ytemp);
  N_VDestroy(x);
  SpbcgFree(spbcg_mem);
  free(cvspbcg_mem);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgAtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = Mv, where
 * M = I - gamma*J. The product J*v is obtained by calling the jtimes
 * routine. It is then scaled by -gamma and added to v to obtain M*v.
 * The return value is the same as the value returned by jtimes --
 * 0 if successful, nonzero otherwise.
 * -----------------------------------------------------------------
 */

static int CVSpbcgAtimes(void *cvode_mem, N_Vector v, N_Vector z)
{
  CVodeMem   cv_mem;
  CVSpbcgMem cvspbcg_mem;
  int jtflag;

  cv_mem = (CVodeMem) cvode_mem;
  cvspbcg_mem = (CVSpbcgMem) lmem;

  jtflag = jtimes(v, z, tn, ycur, fcur, j_data, ytemp);
  njtimes++;
  if (jtflag != 0) return(jtflag);

  N_VLinearSum(ONE, v, -gamma, z, z);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic SpbcgSolve routine and
 * the user's psolve routine. It passes to psolve all required state
 * information from cvode_mem. Its return value is the same as that
 * returned by psolve. Note that the generic SPBCG solver guarantees
 * that CVSpbcgPSolve will not be called in the case in which
 * preconditioning is not done. This is the only case in which the
 * user's psolve routine is allowed to be NULL.
 * -----------------------------------------------------------------
 */

static int CVSpbcgPSolve(void *cvode_mem, N_Vector r, N_Vector z, int lr)
{
  CVodeMem   cv_mem;
  CVSpbcgMem cvspbcg_mem;
  int ier;

  cv_mem = (CVodeMem) cvode_mem;
  cvspbcg_mem = (CVSpbcgMem)lmem;

  ier = psolve(tn, ycur, fcur, r, z, gamma, delta, lr, P_data, ytemp);
  /* This call is counted in nps within the CVSpbcgSolve routine */

  return(ier);     
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgDQJtimes
 * -----------------------------------------------------------------
 * This routine generates a difference quotient approximation to
 * the Jacobian times vector f_y(t,y) * v. The approximation is
 * Jv = vnrm[f(y + v/vnrm) - f(y)], where vnrm = (WRMS norm of v) is
 * input, i.e. the WRMS norm of v/vnrm is 1.
 * -----------------------------------------------------------------
 */

static int CVSpbcgDQJtimes(N_Vector v, N_Vector Jv, realtype t, 
                           N_Vector y, N_Vector fy,
                           void *jac_data, N_Vector work)
{
  CVodeMem cv_mem;
  CVSpbcgMem cvspbcg_mem;
  realtype vnrm;

  /* jac_data is cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvspbcg_mem = (CVSpbcgMem) lmem;

  /* Evaluate norm of v */
  vnrm = N_VWrmsNorm(v, ewt);

  /* Set work = y + (1/vnrm) v */
  N_VLinearSum(ONE/vnrm, v, ONE, y, work);

  /* Set Jv = f(tn, work) */
  f(t, work, Jv, f_data); 
  nfeSB++;

  /* Replace Jv by vnrm*(Jv - fy) */
  N_VLinearSum(ONE, Jv, -ONE, fy, Jv);
  N_VScale(vnrm, Jv, Jv);

  return(0);
}
