/*
 * -----------------------------------------------------------------
 * $Revision: 1.1.2.1 $
 * $Date: 2005-01-19 16:44:00 $
 * -----------------------------------------------------------------
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the optional input and output
 * functions for the CVODES solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_impl.h"
#include "sundialstypes.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* 
 * =================================================================
 * CVODE optional input functions
 * =================================================================
 */

/* 
 * CVodeSetErrFile
 *
 * Specifies the FILE pointer for output (NULL means no messages)
 */

int CVodeSetErrFile(void *cvode_mem, FILE *errfp)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_errfp = errfp;

  return(CV_SUCCESS);
}

#define errfp (cv_mem->cv_errfp)

/* 
 * CVodeSetIterType
 *
 * Specifies the iteration type (CV_FUNCTIONAL or CV_NEWTON)
 */

int CVodeSetIterType(void *cvode_mem, int iter)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if ((iter != CV_FUNCTIONAL) && (iter != CV_NEWTON)) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_BAD_ITER);
    return (CV_ILL_INPUT);
  }

  cv_mem->cv_iter = iter;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetFdata
 *
 * Specifies the user data pointer for f
 */

int CVodeSetFdata(void *cvode_mem, void *f_data)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_f_data = f_data;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetGdata
 *
 * Specifies the user data pointer for g
 */

int CVodeSetGdata(void *cvode_mem, void *g_data)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_g_data = g_data;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMaxOrd
 *
 * Specifies the maximum method order
 */

int CVodeSetMaxOrd(void *cvode_mem, int maxord)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (maxord <= 0) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NEG_MAXORD);
    return(CV_ILL_INPUT);
  }
  
  if (maxord > cv_mem->cv_qmax) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_BAD_MAXORD);
    return(CV_ILL_INPUT);
  }

  cv_mem->cv_qmax = maxord;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMaxNumSteps
 *
 * Specifies the maximum number of integration steps
 */

int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (mxsteps<=0) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NEG_MXSTEPS);
    return(CV_ILL_INPUT);
  }

  cv_mem->cv_mxstep = mxsteps;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMaxHnilWarns
 *
 * Specifies the maximum number of warnings for small h
 */

int CVodeSetMaxHnilWarns(void *cvode_mem, int mxhnil)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_mxhnil = mxhnil;

  return(CV_SUCCESS);
}

/* 
 *CVodeSetStabLimDet
 *
 * Turns on/off the stability limit detection algorithm
 */

int CVodeSetStabLimDet(void *cvode_mem, booleantype sldet)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if(cv_mem->cv_lmm != CV_BDF) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_SLDET);
    return(CV_ILL_INPUT);
  }

  cv_mem->cv_sldeton = sldet;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetInitStep
 *
 * Specifies the initial step size
 */

int CVodeSetInitStep(void *cvode_mem, realtype hin)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_hin = hin;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMinStep
 *
 * Specifies the minimum step size
 */

int CVodeSetMinStep(void *cvode_mem, realtype hmin)
{
  realtype hmax;
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (hmin<=0) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NEG_HMIN);
    return(CV_ILL_INPUT);
  }

  if (hmin * cv_mem->cv_hmax_inv > ONE) {
    hmax = ONE/cv_mem->cv_hmax_inv;
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_BAD_HMIN_HMAX);
    return(CV_ILL_INPUT);
  }

  cv_mem->cv_hmin = hmin;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMaxStep
 *
 * Specifies the maximum step size
 */

int CVodeSetMaxStep(void *cvode_mem, realtype hmax)
{
  realtype hmax_inv;
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (hmax <= 0) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NEG_HMAX);
    return(CV_ILL_INPUT);
  }

  hmax_inv = ONE/hmax;
  if (hmax_inv * cv_mem->cv_hmin > ONE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_BAD_HMIN_HMAX);
    return(CV_ILL_INPUT);
  }

  cv_mem->cv_hmax_inv = hmax_inv;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetStopTime
 *
 * Specifies the time beyond which the integration is not to
 * proceed
 */

int CVodeSetStopTime(void *cvode_mem, realtype tstop)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_tstop = tstop;
  cv_mem->cv_tstopset = TRUE;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMaxErrTestFails
 *
 * Specifies the maximum number of error test failures during one
 * step try.
 */

int CVodeSetMaxErrTestFails(void *cvode_mem, int maxnef)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_maxnef = maxnef;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMaxConvFails
 *
 * Specifies the maximum number of nonlinear convergence failures 
 * during one step try.
 */

int CVodeSetMaxConvFails(void *cvode_mem, int maxncf)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_maxncf = maxncf;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMaxNonlinIters
 *
 * Specifies the maximum number of nonlinear iterations during
 * one solve.
 */

int CVodeSetMaxNonlinIters(void *cvode_mem, int maxcor)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_maxcor = maxcor;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetNonlinConvCoef
 *
 * Specifies the coeficient in the nonlinear solver convergence
 * test
 */

int CVodeSetNonlinConvCoef(void *cvode_mem, realtype nlscoef)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_nlscoef = nlscoef;

  return(CV_SUCCESS);
}

/*
 * CVodeSetTolerances
 *
 * Changes te integration tolerances between calls to CVode()
 */

int CVodeSetTolerances(void *cvode_mem, 
                       int itol, realtype *reltol, void *abstol)
{
  CVodeMem cv_mem;
  booleantype neg_abstol;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if ((itol != CV_SS) && (itol != CV_SV)) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_BAD_ITOL);
    return(CV_ILL_INPUT);
  }

  if (*reltol < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_BAD_RELTOL);
    return(CV_ILL_INPUT);
  }

  if (abstol == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_ABSTOL_NULL);
    return(CV_ILL_INPUT);
  }

  if (itol == CV_SS) {
    neg_abstol = (*((realtype *)abstol) < ZERO);
  } else {
    neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);
  }
  if (neg_abstol) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_BAD_ABSTOL);
    return(CV_ILL_INPUT);
  }

  cv_mem->cv_itol   = itol;
  cv_mem->cv_reltol = reltol;      
  cv_mem->cv_abstol = abstol;

  return(CV_SUCCESS);
}

/* 
 * =================================================================
 * Quadrature optional input functions
 * =================================================================
 */

int CVodeSetQuadFdata(void *cvode_mem, void *fQ_data)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_fQ_data = fQ_data;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetQuadErrCon(void *cvode_mem, booleantype errconQ)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_errconQ = errconQ;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetQuadTolerances(void *cvode_mem, int itolQ, 
                           realtype *reltolQ, void *abstolQ)
{
  CVodeMem cv_mem;
  booleantype neg_abstol;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }
  
  cv_mem = (CVodeMem) cvode_mem;

  if ((itolQ != CV_SS) && (itolQ != CV_SV)) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_ITOLQ);
    return(CV_ILL_INPUT);
  }

  if (reltolQ == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_RELTOLQ_NULL);
    return(CV_ILL_INPUT);
  }
  
  if (*reltolQ < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_RELTOLQ);
    return(CV_ILL_INPUT);
  }

  if (abstolQ == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_ABSTOLQ_NULL);
    return(CV_ILL_INPUT);
  }

  if (itolQ == CV_SS) {
    neg_abstol = (*((realtype *)abstolQ) < ZERO);
  } else {
    neg_abstol = (N_VMin((N_Vector)abstolQ) < ZERO);
  }
  if (neg_abstol) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_ABSTOLQ);
    return(CV_ILL_INPUT);
  }

  cv_mem->cv_itolQ    = itolQ;
  cv_mem->cv_reltolQ  = reltolQ;
  cv_mem->cv_abstolQ  = abstolQ;
  
  return(CV_SUCCESS);
}

/* 
 * =================================================================
 * FSA optional input functions
 * =================================================================
 */


int CVodeSetSensRhsFn(void *cvode_mem, CVSensRhsFn fS)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_ifS  = CV_ALLSENS;

  if (fS != NULL) {
    cv_mem->cv_fS      = fS;
    cv_mem->cv_fSDQ    = FALSE;
  } else {
    cv_mem->cv_fS      = CVSensRhsDQ;
    cv_mem->cv_fS_data = cvode_mem;
    cv_mem->cv_fSDQ    = TRUE;
  }

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensRhs1Fn(void *cvode_mem, CVSensRhs1Fn fS1)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;
  
  cv_mem->cv_ifS  = CV_ONESENS;

  if(fS1 != NULL) {
    cv_mem->cv_fS1     = fS1;
    cv_mem->cv_fSDQ    = FALSE;
  } else {
    cv_mem->cv_fS1     = CVSensRhs1DQ;
    cv_mem->cv_fS_data = cvode_mem;
    cv_mem->cv_fSDQ    = TRUE;
  }

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensFdata(void *cvode_mem, void *fS_data)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_fS_data = fS_data;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensRho(void *cvode_mem, realtype rho)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_rhomax = rho;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensPbar(void *cvode_mem, realtype *pbar)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_pbar = pbar;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensErrCon(void *cvode_mem, booleantype errconS)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_errconS = errconS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensMaxNonlinIters(void *cvode_mem, int maxcorS)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_maxcorS = maxcorS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensTolerances(void *cvode_mem, int itolS,
                           realtype *reltolS, void *abstolS)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if ((itolS != CV_SS) && (itolS != CV_SV) && (itolS != CV_EE)) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_ITOLS);
    return(CV_ILL_INPUT);
  }

  if (itolS == CV_EE) {

    /* CVODES will set tolerances */
    cv_mem->cv_setSensTol = TRUE;
    cv_mem->cv_testSensTol = FALSE;

  } else {

    /* Test user-supplied tolerances */
    if (reltolS == NULL) {
      if(errfp!=NULL) fprintf(errfp, MSGCVS_RELTOLS_NULL);
      return(CV_ILL_INPUT);
    }
    
    if (abstolS == NULL) {
      if(errfp!=NULL) fprintf(errfp, MSGCVS_ABSTOLS_NULL);
      return(CV_ILL_INPUT);
    }
    
    cv_mem->cv_itolS    = itolS;
    cv_mem->cv_reltolS  = reltolS;
    cv_mem->cv_abstolS  = abstolS;

    cv_mem->cv_setSensTol = FALSE;
    cv_mem->cv_testSensTol = TRUE;

  }
    
  return(CV_SUCCESS);
}

/* 
 * =================================================================
 * CVODE optional output functions
 * =================================================================
 */

/* 
 * Readability constants
 */

#define nst            (cv_mem->cv_nst)
#define nfe            (cv_mem->cv_nfe)
#define ncfn           (cv_mem->cv_ncfn)
#define netf           (cv_mem->cv_netf)
#define nni            (cv_mem->cv_nni)
#define nsetups        (cv_mem->cv_nsetups)
#define q              (cv_mem->cv_q)
#define next_q         (cv_mem->cv_next_q)
#define ewt            (cv_mem->cv_ewt)  
#define h              (cv_mem->cv_h)
#define next_h         (cv_mem->cv_next_h)
#define h0u            (cv_mem->cv_h0u)
#define tolsf          (cv_mem->cv_tolsf)  
#define acor           (cv_mem->cv_acor)
#define lrw            (cv_mem->cv_lrw)
#define liw            (cv_mem->cv_liw)
#define nge            (cv_mem->cv_nge)
#define iroots         (cv_mem->cv_iroots)
#define nor            (cv_mem->cv_nor)
#define sldeton        (cv_mem->cv_sldeton)
#define tn             (cv_mem->cv_tn)

/*
 * CVodeGetNumSteps
 *
 * Returns the current number of integration steps
 */

int CVodeGetNumSteps(void *cvode_mem, long int *nsteps)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nsteps = nst;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetNumRhsEvals
 *
 * Returns the current number of calls to f
 */

int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevals)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nfevals = nfe;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetNumLinSolvSetups
 *
 * Returns the current number of calls to the linear solver setup routine
 */

int CVodeGetNumLinSolvSetups(void *cvode_mem, long int *nlinsetups)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nlinsetups = nsetups;

  return(CV_SUCCESS);
}

/*
 * CVodeGetNumErrTestFails
 *
 * Returns the current number of error test failures
 */

int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *netfails = netf;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetLastOrder
 *
 * Returns the order on the last succesful step
 */

int CVodeGetLastOrder(void *cvode_mem, int *qlast)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *qlast = q;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetCurrentOrder
 *
 * Returns the order to be attempted on the next step
 */

int CVodeGetCurrentOrder(void *cvode_mem, int *qcur)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *qcur = next_q;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetNumStabLimOrderReds
 *
 * Returns the number of order reductions triggered by the stability
 * limit detection algorithm
 */

int CVodeGetNumStabLimOrderReds(void *cvode_mem, long int *nslred)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sldeton==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SLDET);
    return(CV_NO_SLDET);
  }

  *nslred = nor;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetActualInitStep
 *
 * Returns the step size used on the first step
 */

int CVodeGetActualInitStep(void *cvode_mem, realtype *hinused)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *hinused = h0u;

  return(CV_SUCCESS);
}

/*
 * CVodeGetLastStep
 *
 * Returns the step size used on the last successful step
 */

int CVodeGetLastStep(void *cvode_mem, realtype *hlast)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *hlast = h;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetCurrentStep
 *
 * Returns the step size to be attempted on the next step
 */

int CVodeGetCurrentStep(void *cvode_mem, realtype *hcur)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;
  
  *hcur = next_h;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetCurrentTime
 *
 * Returns the current value of the independent variable
 */

int CVodeGetCurrentTime(void *cvode_mem, realtype *tcur)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *tcur = tn;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetTolScaleFactor
 *
 * Returns a suggested factor for scaling tolerances
 */

int CVodeGetTolScaleFactor(void *cvode_mem, realtype *tolsfact)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *tolsfact = tolsf;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetErrWeights
 *
 * This routine returns the current weight vector for y in weight.
 * Note that weight need not be allocated by the user.
 */

int CVodeGetErrWeights(void *cvode_mem, N_Vector *eweight)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *eweight = ewt;

  return(CV_SUCCESS);
}

/*
 * CVodeGetEstLocalErrors
 *
 * Returns an estimate of the local error
 */

int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector *ele)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *ele = acor;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetWorkSpace
 *
 * Returns integrator work space requirements
 */

int CVodeGetWorkSpace(void *cvode_mem, long int *lenrw, long int *leniw)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *leniw = liw;
  *lenrw = lrw;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetIntegratorStats
 *
 * Returns integrator statistics
 */

int CVodeGetIntegratorStats(void *cvode_mem, long int *nsteps, long int *nfevals, 
                            long int *nlinsetups, long int *netfails, int *qlast, 
                            int *qcur, realtype *hinused, realtype *hlast, 
                            realtype *hcur, realtype *tcur)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nsteps = nst;
  *nfevals = nfe;
  *nlinsetups = nsetups;
  *netfails = netf;
  *qlast = q;
  *qcur = next_q;
  *hinused = h0u;
  *hlast = h;
  *hcur = next_h;
  *tcur = tn;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetNumGEvals
 *
 * Returns the current number of calls to g (for rootfinding)
 */

int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *ngevals = nge;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetRootInfo
 *
 * Returns pointer to array rootsfound showing roots found
 */

int CVodeGetRootInfo(void *cvode_mem, int **rootsfound)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *rootsfound = iroots;

  return(CV_SUCCESS);
}


/* 
 * CVodeGetNumNonlinSolvIters
 *
 * Returns the current number of iterations in the nonlinear solver
 */

int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nniters = nni;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetNumNonlinSolvConvFails
 *
 * Returns the current number of convergence failures in the
 * nonlinear solver
 */

int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nncfails = ncfn;

  return(CV_SUCCESS);
}

/* 
 * CVodeGetNonlinSolvStats
 *
 * Returns nonlinear solver statistics
 */

int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters, 
                            long int *nncfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nniters = nni;
  *nncfails = ncfn;

  return(CV_SUCCESS);
}

/* 
 * =================================================================
 * Quadrature optional output functions
 * =================================================================
 */

/* 
 * Readability constants
 */

#define quadr          (cv_mem->cv_quadr)
#define nfQe           (cv_mem->cv_nfQe)
#define netfQ          (cv_mem->cv_netfQ)
#define ewtQ           (cv_mem->cv_ewtQ)
#define errconQ        (cv_mem->cv_errconQ)

/*-----------------------------------------------------------------*/

int CVodeGetQuadNumRhsEvals(void *cvode_mem, long int *nfQevals)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (quadr==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_QUAD);
    return(CV_NO_QUAD);
  }

  *nfQevals = nfQe;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetQuadNumErrTestFails(void *cvode_mem, long int *nQetfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (quadr==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_QUAD);
    return(CV_NO_QUAD);
  }

  *nQetfails = netfQ;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetQuadErrWeights(void *cvode_mem, N_Vector *eQweight)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (quadr==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_QUAD);
    return(CV_NO_QUAD);
  }

  if(errconQ) *eQweight = ewtQ;
  else        *eQweight = NULL;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetQuadStats(void *cvode_mem, long int *nfQevals, long int *nQetfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (quadr==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_QUAD);
    return(CV_NO_QUAD);
  }

  *nfQevals = nfQe;
  *nQetfails = netfQ;

  return(CV_SUCCESS);
}

/* 
 * =================================================================
 * FSA optional output functions
 * =================================================================
 */

/* 
 * Readability constants
 */

#define sensi          (cv_mem->cv_sensi)
#define ism            (cv_mem->cv_ism)
#define ewtS           (cv_mem->cv_ewtS)
#define nfSe           (cv_mem->cv_nfSe)
#define nfeS           (cv_mem->cv_nfeS)
#define nniS           (cv_mem->cv_nniS)
#define ncfnS          (cv_mem->cv_ncfnS)
#define netfS          (cv_mem->cv_netfS)
#define nsetupsS       (cv_mem->cv_nsetupsS)
#define nniS1          (cv_mem->cv_nniS1)
#define ncfnS1         (cv_mem->cv_ncfnS1)
#define ncfS1          (cv_mem->cv_ncfS1)

/*-----------------------------------------------------------------*/

int CVodeGetNumSensRhsEvals(void *cvode_mem, long int *nfSevals)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  *nfSevals = nfSe;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumRhsEvalsSens(void *cvode_mem, long int *nfevalsS)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  *nfevalsS = nfeS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumSensErrTestFails(void *cvode_mem, long int *nSetfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  *nSetfails = netfS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumSensLinSolvSetups(void *cvode_mem, long int *nlinsetupsS)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  *nlinsetupsS = nsetupsS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensErrWeights(void *cvode_mem, N_Vector_S *eSweight)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  *eSweight = ewtS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensStats(void *cvode_mem, long int *nfSevals, long int *nfevalsS, 
                      long int *nSetfails, long int *nlinsetupsS)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  *nfSevals = nfSe;
  *nfevalsS = nfeS;
  *nSetfails = netfS;
  *nlinsetupsS = nsetupsS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumSensNonlinSolvIters(void *cvode_mem, long int *nSniters)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  *nSniters = nniS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumSensNonlinSolvConvFails(void *cvode_mem, long int *nSncfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  *nSncfails = ncfnS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumStgrSensNonlinSolvIters(void *cvode_mem, long int *nSTGR1niters)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  if(ism==CV_STAGGERED1) nSTGR1niters = nniS1;
  else                nSTGR1niters = NULL;
    

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumStgrSensNonlinSolvConvFails(void *cvode_mem, long int *nSTGR1ncfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  if(ism==CV_STAGGERED1) nSTGR1ncfails = ncfnS1;
  else                nSTGR1ncfails = NULL;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensNonlinSolvStats(void *cvode_mem, long int *nSniters, 
                                long int *nSncfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  *nSniters = nniS;
  *nSncfails = ncfnS;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetStgrSensNonlinSolvStats(void *cvode_mem, long int *nSTGR1niters, 
                                    long int *nSTGR1ncfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  if(ism==CV_STAGGERED1) {
    nSTGR1niters  = nniS1;
    nSTGR1ncfails = ncfnS1;
  } else {
    nSTGR1niters  = NULL;
    nSTGR1ncfails = NULL;
  }
  return(CV_SUCCESS);
}

