/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2005-06-21 23:58:06 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/ida/LICENSE
 * -----------------------------------------------------------------
 * This is the implementation file for the optional inputs and     
 * outputs for the IDA solver.                                    
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "ida_impl.h"

#define ZERO    RCONST(0.0)
#define HALF    RCONST(0.5)
#define ONE     RCONST(1.0)
#define TWOPT5  RCONST(2.5)

#define lrw  (IDA_mem->ida_lrw)
#define liw  (IDA_mem->ida_liw)
#define lrw1 (IDA_mem->ida_lrw1)
#define liw1 (IDA_mem->ida_liw1)

/* 
 * =================================================================
 * IDA optional input functions
 * =================================================================
 */

int IDASetErrFile(void *ida_mem, FILE *errfp)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_errfp = errfp;

  return(IDA_SUCCESS);
}

#define errfp (IDA_mem->ida_errfp)

/*-----------------------------------------------------------------*/

int IDASetRdata(void *ida_mem, void *res_data)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rdata = res_data;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxOrd(void *ida_mem, int maxord)
{
  IDAMem IDA_mem;
  int maxord_alloc;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxord <= 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_MAXORD);
    return(IDA_ILL_INPUT);
  }

  /* Cannot increase maximum order beyond the value that
     was used when allocating memory */
  maxord_alloc = IDA_mem->ida_maxord_alloc;

  if (maxord > maxord_alloc) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_MAXORD);
    return(IDA_ILL_INPUT);
  }  

  IDA_mem->ida_maxord = maxord;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumSteps(void *ida_mem, long int mxsteps)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (mxsteps < 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_MXSTEPS);
    return(IDA_ILL_INPUT);
  }

  /* Passing 0 sets the default */
  if (mxsteps == 0)
    IDA_mem->ida_mxstep = MXSTEP_DEFAULT;
  else
    IDA_mem->ida_mxstep = mxsteps;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetInitStep(void *ida_mem, realtype hin)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_hin = hin;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxStep(void *ida_mem, realtype hmax)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (hmax < 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_HMAX);
    return(IDA_ILL_INPUT);
  }

  /* Passing 0 sets hmax = infinity */
  if (hmax == ZERO) {
    IDA_mem->ida_hmax_inv = HMAX_INV_DEFAULT;
    return(IDA_SUCCESS);
  }

  IDA_mem->ida_hmax_inv = ONE/hmax;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetStopTime(void *ida_mem, realtype tstop)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_tstop = tstop;
  IDA_mem->ida_tstopset = TRUE;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetNonlinConvCoef(void *ida_mem, realtype epcon)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (epcon < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_EPCON);
    return(IDA_ILL_INPUT);
  }

  IDA_mem->ida_epcon = epcon;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxErrTestFails(void *ida_mem, int maxnef)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_maxnef = maxnef;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxConvFails(void *ida_mem, int maxncf)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_maxncf = maxncf;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNonlinIters(void *ida_mem, int maxcor)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_maxcor = maxcor;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_suppressalg = suppressalg;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetId(void *ida_mem, N_Vector id)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (id == NULL) {
    if (IDA_mem->ida_idMallocDone) {
      N_VDestroy(IDA_mem->ida_id);
      lrw -= lrw1;
      liw -= liw1;
    }
    IDA_mem->ida_idMallocDone = FALSE;    
    return(IDA_SUCCESS);
  }

  if ( !(IDA_mem->ida_idMallocDone) ) {
    IDA_mem->ida_id = N_VClone(id);
    lrw += lrw1;
    liw += liw1;
    IDA_mem->ida_idMallocDone = TRUE;
  }

  /* Load the id vector */

  N_VScale(ONE, id, IDA_mem->ida_id);

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetConstraints(void *ida_mem, N_Vector constraints)
{
  IDAMem IDA_mem;
  realtype temptest;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (constraints == NULL) {
    if (IDA_mem->ida_constraintsMallocDone) {
      N_VDestroy(IDA_mem->ida_constraints);
      lrw -= lrw1;
      liw -= liw1;
    }
    IDA_mem->ida_constraintsMallocDone = FALSE;
    IDA_mem->ida_constraintsSet = FALSE;
    return(IDA_SUCCESS);
  }

  /* Test if required vector ops. are defined */

  if (constraints->ops->nvdiv         == NULL ||
      constraints->ops->nvmaxnorm     == NULL ||
      constraints->ops->nvcompare     == NULL ||
      constraints->ops->nvconstrmask  == NULL ||
      constraints->ops->nvminquotient == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_NVECTOR);
    return(IDA_ILL_INPUT);
  }

  /*  Check the constraints vector */

  temptest = N_VMaxNorm(constraints);
  if((temptest > TWOPT5) || (temptest < HALF)){ 
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_CONSTRAINTS); 
    return(IDA_ILL_INPUT); 
  }

  if ( !(IDA_mem->ida_constraintsMallocDone) ) {
    IDA_mem->ida_constraints = N_VClone(constraints);
    lrw += lrw1;
    liw += liw1;
    IDA_mem->ida_constraintsMallocDone = TRUE;
  }

  /* Load the constraints vector */

  N_VScale(ONE, constraints, IDA_mem->ida_constraints);

  IDA_mem->ida_constraintsSet = TRUE;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetTolerances(void *ida_mem, 
                     int itol, realtype rtol, void *atol)
{
  IDAMem IDA_mem;
  booleantype neg_atol;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  /* Check if ida_mem was allocated */

  if (IDA_mem->ida_MallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NO_MALLOC);
    return(IDA_NO_MALLOC);
  }

  /* Check inputs */

  if ((itol != IDA_SS) && (itol != IDA_SV)) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_ITOL);
    return(IDA_ILL_INPUT);
  }

  if (atol == NULL) { 
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_ATOL_NULL); 
    return(IDA_ILL_INPUT); 
  }

  if (rtol < ZERO) { 
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_RTOL); 
    return(IDA_ILL_INPUT); 
  }

    
  if (itol == IDA_SS) { 
    neg_atol = (*((realtype *)atol) < ZERO); 
  } else { 
    neg_atol = (N_VMin((N_Vector)atol) < ZERO); 
  }

  if (neg_atol) { 
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_ATOL); 
    return(IDA_ILL_INPUT); 
  }

  /* Copy tolerances into memory */

  if ( (itol != IDA_SV) && (IDA_mem->ida_VatolMallocDone) ) {
    N_VDestroy(IDA_mem->ida_Vatol);
    lrw -= lrw1;
    liw -= liw1;
    IDA_mem->ida_VatolMallocDone = FALSE;
  }

  if ( (itol == IDA_SV) && !(IDA_mem->ida_VatolMallocDone) ) {
    IDA_mem->ida_Vatol = N_VClone(IDA_mem->ida_ewt);
    lrw += lrw1;
    liw += liw1;
    IDA_mem->ida_VatolMallocDone = TRUE;
  }

  IDA_mem->ida_itol = itol;
  IDA_mem->ida_rtol = rtol;      
  if (itol == IDA_SS)
    IDA_mem->ida_Satol = *((realtype *)atol);
  else 
    N_VScale(ONE, (N_Vector)atol, IDA_mem->ida_Vatol);

  IDA_mem->ida_efun = IDAEwtSet;
  IDA_mem->ida_edata = ida_mem;

  return(IDA_SUCCESS);
}

/* 
 * IDASetEwtFn
 *
 * Specifies the user-provide EwtSet function and data pointer for e
 */

int IDASetEwtFn(void *ida_mem, IDAEwtFn efun, void *edata)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if ( IDA_mem->ida_VatolMallocDone ) {
    N_VDestroy(IDA_mem->ida_Vatol);
    lrw -= lrw1;
    liw -= liw1;
    IDA_mem->ida_VatolMallocDone = FALSE;
  }

  IDA_mem->ida_itol = IDA_WF;
  IDA_mem->ida_efun = efun;
  IDA_mem->ida_edata = edata;

  return(IDA_SUCCESS);
}


/* 
 * =================================================================
 * IDA IC optional input functions
 * =================================================================
 */

int IDASetNonlinConvCoefIC(void *ida_mem, realtype epiccon)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (epiccon < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_EPICCON);
    return(IDA_ILL_INPUT);
  }

  IDA_mem->ida_epiccon = epiccon;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumStepsIC(void *ida_mem, int maxnh)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxnh < 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_MAXNH);
    return(IDA_ILL_INPUT);
  }

  IDA_mem->ida_maxnh = maxnh;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumJacsIC(void *ida_mem, int maxnj)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

   if (maxnj < 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_MAXNJ);
    return(IDA_ILL_INPUT);
  } 

  IDA_mem->ida_maxnj = maxnj;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumItersIC(void *ida_mem, int maxnit)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxnit < 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_MAXNIT);
    return(IDA_ILL_INPUT);
  }

  IDA_mem->ida_maxnit = maxnit;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_lsoff = lsoff;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetStepToleranceIC(void *ida_mem, realtype steptol)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (steptol < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_STEPTOL);
    return(IDA_ILL_INPUT);
  }

  IDA_mem->ida_steptol = steptol;

  return(IDA_SUCCESS);
}

/* 
 * =================================================================
 * Readability constants
 * =================================================================
 */

#define ewt         (IDA_mem->ida_ewt)
#define kk          (IDA_mem->ida_kk)
#define hh          (IDA_mem->ida_hh)
#define h0u         (IDA_mem->ida_h0u)
#define tn          (IDA_mem->ida_tn)
#define nbacktr     (IDA_mem->ida_nbacktr)
#define nst         (IDA_mem->ida_nst)
#define nre         (IDA_mem->ida_nre)
#define ncfn        (IDA_mem->ida_ncfn)
#define netf        (IDA_mem->ida_netf)
#define nni         (IDA_mem->ida_nni)
#define nsetups     (IDA_mem->ida_nsetups)
#define lrw         (IDA_mem->ida_lrw)
#define liw         (IDA_mem->ida_liw)
#define kused       (IDA_mem->ida_kused)          
#define hused       (IDA_mem->ida_hused)         
#define tolsf       (IDA_mem->ida_tolsf) 
#define efun        (IDA_mem->ida_efun)
#define edata       (IDA_mem->ida_edata)
#define nge         (IDA_mem->ida_nge)
#define iroots      (IDA_mem->ida_iroots)

/* 
 * =================================================================
 * IDA optional input functions
 * =================================================================
 */

int IDAGetNumSteps(void *ida_mem, long int *nsteps)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nsteps = nst;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumResEvals(void *ida_mem, long int *nrevals)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nrevals = nre;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumLinSolvSetups(void *ida_mem, long int *nlinsetups)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nlinsetups = nsetups;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumErrTestFails(void *ida_mem, long int *netfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *netfails = netf;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumBacktrackOps(void *ida_mem, long int *nbacktracks)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nbacktracks = nbacktr;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetLastOrder(void *ida_mem, int *klast)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *klast = kused;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentOrder(void *ida_mem, int *kcur)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *kcur = kk;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetActualInitStep(void *ida_mem, realtype *hinused)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *hinused = h0u;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetLastStep(void *ida_mem, realtype *hlast)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *hlast = hused;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentStep(void *ida_mem, realtype *hcur)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *hcur = hh;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentTime(void *ida_mem, realtype *tcur)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *tcur = tn;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetTolScaleFactor(void *ida_mem, realtype *tolsfact)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *tolsfact = tolsf;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetErrWeights(void *ida_mem, N_Vector eweight)
{
  IDAMem IDA_mem;
  
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem; 

  N_VScale(ONE, ewt, eweight);

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetErrWeightsAtY(void *ida_mem, N_Vector y, N_Vector eweight)
{
  IDAMem IDA_mem;
  int ewtsetOK;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem; 

  ewtsetOK = efun(y, eweight, edata);

  if (ewtsetOK != 0) {
    fprintf(stderr, MSG_IDAG_EWT_BAD);
    return(IDA_ILL_INPUT);
  }

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetWorkSpace(void *ida_mem, long int *lenrw, long int *leniw)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *leniw = liw;
  *lenrw = lrw;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetIntegratorStats(void *ida_mem, long int *nsteps, long int *nrevals, 
                          long int *nlinsetups, long int *netfails,
                          int *klast, int *kcur, realtype *hlast, 
                          realtype *hcur, realtype *tcur)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nsteps     = nst;
  *nrevals    = nre;
  *nlinsetups = nsetups;
  *netfails   = netf;
  *klast      = kused;
  *kcur       = kk;
  *hlast      = hused;
  *hcur       = hh;  
  *tcur       = tn;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumGEvals(void *ida_mem, long int *ngevals)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *ngevals = nge;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetRootInfo(void *ida_mem, int *rootsfound)
{
  IDAMem IDA_mem;
  int i, nrt;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  nrt = IDA_mem->ida_nrtfn;

  for (i=0; i<nrt; i++) rootsfound[i] = iroots[i];

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumNonlinSolvIters(void *ida_mem, long int *nniters)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nniters = nni;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumNonlinSolvConvFails(void *ida_mem, long int *nncfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nncfails = ncfn;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNonlinSolvStats(void *ida_mem, long int *nniters, long int *nncfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nniters  = nni;
  *nncfails = ncfn;

  return(IDA_SUCCESS);
}

