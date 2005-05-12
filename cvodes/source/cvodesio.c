/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2005-05-12 21:03:17 $
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
#include "sundialsmath.h"
#include "sundialstypes.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

#define lrw   (cv_mem->cv_lrw)
#define liw   (cv_mem->cv_liw)
#define lrw1  (cv_mem->cv_lrw1)
#define liw1  (cv_mem->cv_liw1)
#define lrw1Q (cv_mem->cv_lrw1Q)
#define liw1Q (cv_mem->cv_liw1Q)

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

  if (mxsteps < 0) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NEG_MXSTEPS);
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

  if(sldet && (cv_mem->cv_lmm != CV_BDF) ) {
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
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (hmin < 0) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NEG_HMIN);
    return(CV_ILL_INPUT);
  }

  /* Passing 0 sets hmin = zero */
  if (hmin == ZERO) {
    cv_mem->cv_hmin = HMIN_DEFAULT;
    return(CV_SUCCESS);
  }

  if (hmin * cv_mem->cv_hmax_inv > ONE) {
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

  if (hmax < 0) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NEG_HMAX);
    return(CV_ILL_INPUT);
  }

  /* Passing 0 sets hmax = infinity */
  if (hmax == ZERO) {
    cv_mem->cv_hmax_inv = HMAX_INV_DEFAULT;
    return(CV_SUCCESS);
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
 * Changes the integration tolerances between calls to CVode()
 */

int CVodeSetTolerances(void *cvode_mem, 
                       int itol, realtype reltol, void *abstol)
{
  CVodeMem cv_mem;
  booleantype neg_abstol;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  /* Check if cvode_mem was allocated */

  if (cv_mem->cv_MallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NO_MALLOC);
    return(CV_NO_MALLOC);
  }

  /* Check inputs */

  if ( (itol != CV_SS) && (itol != CV_SV) ) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_BAD_ITOL);
    return(CV_ILL_INPUT);
  }

  if (abstol == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_ABSTOL_NULL);
    return(CV_ILL_INPUT);
  }

  if (reltol < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_RELTOL);
    return (CV_ILL_INPUT);
  }

  if (itol == CV_SS)
    neg_abstol = (*((realtype *)abstol) < ZERO);
  else
    neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);

  if (neg_abstol) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_BAD_ABSTOL);
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
    cv_mem->cv_Vabstol = N_VClone(cv_mem->cv_tempv);
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
    fprintf(stderr, MSGCVS_SET_NO_MEM);
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

int CVodeSetQuadErrCon(void *cvode_mem, booleantype errconQ, 
                       int itolQ, realtype reltolQ, void *abstolQ)
{
  CVodeMem cv_mem;
  booleantype neg_abstol;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }
  
  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_errconQ = errconQ;

  /* Ckeck if quadrature was initialized? */

  if (cv_mem->cv_quadMallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NO_QUAD);
    return(CV_NO_QUAD);
  }

  /* Check inputs */

  if(errconQ == FALSE) {
    if (cv_mem->cv_VabstolQMallocDone) {
      N_VDestroy(cv_mem->cv_VabstolQ);
      lrw -= lrw1Q;
      liw -= liw1Q;
      cv_mem->cv_VabstolQMallocDone = FALSE;
    }
    return(CV_SUCCESS);
  }
  
  if ((itolQ != CV_SS) && (itolQ != CV_SV)) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_ITOLQ);
    return(CV_ILL_INPUT);
  }

  if (abstolQ == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_ABSTOLQ_NULL);
    return(CV_ILL_INPUT);
  }

  if (reltolQ < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_RELTOLQ);
    return(CV_ILL_INPUT);
  }

  if (itolQ == CV_SS)
    neg_abstol = (*((realtype *)abstolQ) < ZERO);
  else
    neg_abstol = (N_VMin((N_Vector)abstolQ) < ZERO);

  if (neg_abstol) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_ABSTOLQ);
    return(CV_ILL_INPUT);
  }

  /* See if we need to free or allocate memory */

  if ( (itolQ != CV_SV) && (cv_mem->cv_VabstolQMallocDone) ) {
    N_VDestroy(cv_mem->cv_VabstolQ);
    lrw -= lrw1Q;
    liw -= liw1Q;
    cv_mem->cv_VabstolQMallocDone = FALSE;
  }

  if ( (itolQ == CV_SV) && !(cv_mem->cv_VabstolQMallocDone) ) {
    cv_mem->cv_VabstolQ = N_VClone(cv_mem->cv_tempvQ);
    lrw += lrw1Q;
    liw += liw1Q;
    cv_mem->cv_VabstolMallocDone = TRUE;
  }

  /* Copy tolerances into memory */

  cv_mem->cv_itolQ    = itolQ;
  cv_mem->cv_reltolQ  = reltolQ;

  if (itolQ == CV_SS)
    cv_mem->cv_SabstolQ = *((realtype *)abstolQ);
  else
    N_VScale(ONE, (N_Vector)abstolQ, cv_mem->cv_VabstolQ);
  
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

int CVodeSetSensParams(void *cvode_mem, realtype *p, realtype *pbar, int *plist)
{
  CVodeMem cv_mem;
  int is, Ns;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  Ns = cv_mem->cv_Ns;

  /* Parameters */

  cv_mem->cv_p = p;

  /* pbar */

  if (pbar != NULL)
    for (is=0; is<Ns; is++) {
      if (pbar[is] == ZERO) {
        if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_PBAR);
        return(CV_ILL_INPUT);
      }
      cv_mem->cv_pbar[is] = ABS(pbar[is]);
    }
  else
    for (is=0; is<Ns; is++)
      cv_mem->cv_pbar[is] = ONE;

  /* plist */

  if (plist != NULL)
    for (is=0; is<Ns; is++) {
      if ( plist[is] <1 ) {
        if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_PLIST);
        return(CV_ILL_INPUT);
      }
      cv_mem->cv_plist[is] = plist[is];
    }
  else
    for (is=0; is<Ns; is++)
      cv_mem->cv_plist[is] = is+1;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensTolerances(void *cvode_mem, int itolS,
                           realtype reltolS, void *abstolS)
{
  CVodeMem cv_mem;
  booleantype neg_abstol;
  realtype *atolSS;
  N_Vector *atolSV;
  int is, Ns;

  atolSS = NULL;
  atolSV = NULL;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_SET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  /* Was sensitivity initialized? */

  if (cv_mem->cv_sensMallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_SET_NO_SENSI);
    return(CV_NO_SENS);
  } 

  /* Check inputs */

  Ns = cv_mem->cv_Ns;

  if ((itolS != CV_SS) && (itolS != CV_SV) && (itolS != CV_EE)) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_ITOLS);
    return(CV_ILL_INPUT);
  }

  if (itolS != CV_EE) {

    /* Test user-supplied tolerances */
    
    if (reltolS < ZERO) {
      if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_RELTOLS);
      return(CV_ILL_INPUT);
    }

    if (abstolS == NULL) {
      if(errfp!=NULL) fprintf(errfp, MSGCVS_ABSTOLS_NULL);
      return(CV_ILL_INPUT);
    }

    neg_abstol = FALSE;

    if (itolS == CV_SS) {
      atolSS = (realtype *) abstolS;
      for (is=0; is<Ns; is++)
        if (atolSS[is] < ZERO) {neg_abstol = TRUE; break;}
    } else {
      atolSV = (N_Vector *) abstolS;
      for (is=0; is<Ns; is++) 
        if (N_VMin(atolSV[is]) < ZERO) {neg_abstol = TRUE; break;}
    }

    if (neg_abstol) {
      if(errfp!=NULL) fprintf(errfp, MSGCVS_BAD_ABSTOLS);
      return(CV_ILL_INPUT);
    }
    
  }

  /* See if we should release some memory */

  if ( (itolS != CV_SV) && (cv_mem->cv_VabstolSMallocDone) ) {
    N_VDestroyVectorArray(cv_mem->cv_VabstolS, Ns);
    lrw -= Ns*lrw1;
    liw -= Ns*liw1;
    cv_mem->cv_VabstolSMallocDone = FALSE;
  }

  if ( (itolS != CV_SS) && (cv_mem->cv_SabstolSMallocDone) ) {
    free(cv_mem->cv_SabstolS);
    lrw -= Ns;
    cv_mem->cv_SabstolSMallocDone = FALSE;
  }

  /* If tolerances will be estimated, return now */

  if (itolS == CV_EE) return(CV_SUCCESS);

  /* See if we need to allocate some memory */

  if ( (itolS == CV_SV) && !(cv_mem->cv_VabstolSMallocDone) ) {
    cv_mem->cv_VabstolS = N_VCloneVectorArray(Ns, cv_mem->cv_tempv);
    lrw += Ns*lrw1;
    liw += Ns*liw1;
    cv_mem->cv_VabstolSMallocDone = TRUE;
  }

  if ( (itolS == CV_SS) && !(cv_mem->cv_SabstolSMallocDone) ) {
    cv_mem->cv_SabstolS = (realtype *)malloc(Ns*sizeof(realtype));
    lrw += Ns;
    cv_mem->cv_SabstolSMallocDone = TRUE;
  }

  /* Copy tolerances into memory */

  cv_mem->cv_itolS   = itolS;
  cv_mem->cv_reltolS = reltolS;

  if (itolS == CV_SS)
    for (is=0; is<Ns; is++)
      cv_mem->cv_SabstolS[is] = atolSS[is];
  else
    for (is=0; is<Ns; is++)    
      N_VScale(ONE, atolSV[is], cv_mem->cv_VabstolS[is]);

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
#define efun           (cv_mem->cv_efun)
#define e_data         (cv_mem->cv_e_data) 

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
 * This routine returns the error weight vector for y in eweight.
 */

int CVodeGetErrWeights(void *cvode_mem, N_Vector yy, N_Vector eweight)
{
  CVodeMem cv_mem;
  int ewtsetOK;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  ewtsetOK = efun(yy, eweight, e_data);

  if (ewtsetOK != 0) {
    fprintf(stderr, MSGCVS_GET_EWT_BAD);
    return(CV_ILL_INPUT);
  }

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
    fprintf(stderr, MSGCVS_GET_NO_MEM);
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

int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
{
  CVodeMem cv_mem;
  int i, nrt;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
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

int CVodeGetQuadErrWeights(void *cvode_mem, N_Vector eQweight)
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

  if(errconQ) N_VScale(ONE, ewtQ, eQweight);

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
  int is, Ns;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  Ns = cv_mem->cv_Ns;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  if(ism==CV_STAGGERED1) 
    for(is=0; is<Ns; is++) nSTGR1niters[is] = nniS1[is];

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumStgrSensNonlinSolvConvFails(void *cvode_mem, long int *nSTGR1ncfails)
{
  CVodeMem cv_mem;
  int is, Ns;

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCVS_GET_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  Ns = cv_mem->cv_Ns;

  if (sensi==FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCVS_GET_NO_SENSI);
    return(CV_NO_SENS);
  }

  if(ism==CV_STAGGERED1) 
    for(is=0; is<Ns; is++) nSTGR1ncfails[is] = ncfnS1[is];

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


