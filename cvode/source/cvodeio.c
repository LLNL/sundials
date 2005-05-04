/*
 * -----------------------------------------------------------------
 * $Revision: 1.1.2.5 $
 * $Date: 2005-05-04 22:39:01 $
 * -----------------------------------------------------------------
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the optional input and output
 * functions for the CVODE solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "sundialstypes.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

#define lrw (cv_mem->cv_lrw)
#define liw (cv_mem->cv_liw)
#define lrw1 (cv_mem->cv_lrw1)
#define liw1 (cv_mem->cv_liw1)

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
    fprintf(stderr, MSGCV_SET_NO_MEM);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if ((iter != CV_FUNCTIONAL) && (iter != CV_NEWTON)) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_BAD_ITER);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_f_data = f_data;

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
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (maxord <= 0) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_NEG_MAXORD);
    return(CV_ILL_INPUT);
  }
  
  if (maxord > cv_mem->cv_qmax) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_BAD_MAXORD);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (mxsteps < 0) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_NEG_MXSTEPS);
    return(CV_ILL_INPUT);
  }

  /* Passing 0 sets the default */
  if (mxsteps == 0)
    cv_mem->cv_mxstep = MXSTEP_DEFAULT;
  else
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if(cv_mem->cv_lmm != CV_BDF) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_SLDET);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (hmin<0) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_NEG_HMIN);
    return(CV_ILL_INPUT);
  }

  /* Passing 0 sets hmin = zero */
  if (hmin == ZERO) {
    cv_mem->cv_hmin = HMIN_DEFAULT;
    return(CV_SUCCESS);
  }

  if (hmin * cv_mem->cv_hmax_inv > ONE) {
    hmax = ONE/cv_mem->cv_hmax_inv;
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_BAD_HMIN_HMAX);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (hmax < 0) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_NEG_HMAX);
    return(CV_ILL_INPUT);
  }

  /* Passing 0 sets hmax = infinity */
  if (hmax == ZERO) {
    cv_mem->cv_hmax_inv = HMAX_INV_DEFAULT;
    return(CV_SUCCESS);
  }

  hmax_inv = ONE/hmax;
  if (hmax_inv * cv_mem->cv_hmin > ONE) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_BAD_HMIN_HMAX);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
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
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_nlscoef = nlscoef;

  return(CV_SUCCESS);
}

/*
 * CVodeSetTolerances
 *
 * Changes the integration tolerances between calls to CVode()
 */

int CVodeSetTolerances(void *cvode_mem, 
                       int itol, realtype reltol, void *abstol)
{
  CVodeMem cv_mem;
  booleantype neg_abstol;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  /* Check if cvode_mem was allocated */

  if (cv_mem->cv_MallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_NO_MALLOC);
    return(CV_NO_MALLOC);
  }

  /* Check inputs */

  if ( (itol != CV_SS) && (itol != CV_SV) ) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_BAD_ITOL);
    return(CV_ILL_INPUT);
  }

  if (abstol == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_ABSTOL_NULL);
    return(CV_ILL_INPUT);
  }

  if (reltol < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_BAD_RELTOL);
    return(CV_ILL_INPUT);
  }

  if (itol == CV_SS)
    neg_abstol = (*((realtype *)abstol) < ZERO);
  else
    neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);
    
  if (neg_abstol) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_SET_BAD_ABSTOL);
    return(CV_ILL_INPUT);
  }

  /* Copy tolerances into memory */

  if ( (itol != CV_SV) && (cv_mem->cv_VabstolMallocDone) ) {
    N_VDestroy(cv_mem->cv_Vabstol);
    lrw -= lrw1;
    liw -= liw1;
    cv_mem->cv_VabstolMallocDone = FALSE;
  }

  if ( (itol == CV_SV) && !(cv_mem->cv_VabstolMallocDone) ) {
    cv_mem->cv_Vabstol = N_VClone(cv_mem->cv_ewt);
    lrw += lrw1;
    liw += liw1;
    cv_mem->cv_VabstolMallocDone = TRUE;
  }

  cv_mem->cv_itol   = itol;
  cv_mem->cv_reltol = reltol;
  if (itol == CV_SS)
    cv_mem->cv_Sabstol = *((realtype *)abstol);
  else
    N_VScale(ONE, (N_Vector)abstol, cv_mem->cv_Vabstol);

  cv_mem->cv_efun = CVEwtSet;
  cv_mem->cv_e_data = cvode_mem;


  return(CV_SUCCESS);
}

/* 
 * CVodeSetEwtFn
 *
 * Specifies the user-provide EwtSet function and data pointer for e
 */

int CVodeSetEwtFn(void *cvode_mem, CVEwtFn efun, void *e_data)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCV_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if ( cv_mem->cv_VabstolMallocDone ) {
    N_VDestroy(cv_mem->cv_Vabstol);
    lrw -= lrw1;
    liw -= liw1;
    cv_mem->cv_VabstolMallocDone = FALSE;
  }

  cv_mem->cv_itol = CV_WF;
  cv_mem->cv_efun = efun;
  cv_mem->cv_e_data = e_data;

  return(CV_SUCCESS);
}

/* 
 * =================================================================
 * Readability constants
 * =================================================================
 */

#define nst            (cv_mem->cv_nst)
#define nfe            (cv_mem->cv_nfe)
#define ncfn           (cv_mem->cv_ncfn)
#define netf           (cv_mem->cv_netf)
#define nni            (cv_mem->cv_nni)
#define nsetups        (cv_mem->cv_nsetups)
#define qu             (cv_mem->cv_qu)
#define next_q         (cv_mem->cv_next_q)
#define ewt            (cv_mem->cv_ewt)  
#define hu             (cv_mem->cv_hu)
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
 * =================================================================
 * CVODE optional output functions
 * =================================================================
 */

/*
 * CVodeGetNumSteps
 *
 * Returns the current number of integration steps
 */

int CVodeGetNumSteps(void *cvode_mem, long int *nsteps)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *qlast = qu;

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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (sldeton==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_GET_NO_SLDET);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *hlast = hu;

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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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

int CVodeGetErrWeights(void *cvode_mem, N_Vector eweight)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCV_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  N_VScale(ONE, ewt, eweight);

  return(CV_SUCCESS);
}

/*
 * CVodeGetEstLocalErrors
 *
 * Returns an estimate of the local error
 */

int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCV_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  N_VScale(ONE, acor, ele);

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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nsteps = nst;
  *nfevals = nfe;
  *nlinsetups = nsetups;
  *netfails = netf;
  *qlast = qu;
  *qcur = next_q;
  *hinused = h0u;
  *hlast = hu;
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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

int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
{
  CVodeMem cv_mem;
  int i, nrt;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCV_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  nrt = cv_mem->cv_nrtfn;

  for (i=0; i<nrt; i++) rootsfound[i] = iroots[i];

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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
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
    fprintf(stderr, MSGCV_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nniters = nni;
  *nncfails = ncfn;

  return(CV_SUCCESS);
}

