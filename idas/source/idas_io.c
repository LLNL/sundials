/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-11 21:13:57 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * This is the implementation file for the optional inputs and     
 * outputs for the IDAS solver.                                    
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "idas_impl.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

extern int IDASensResDQ(int Ns, realtype t, 
                        N_Vector yy, N_Vector yp, N_Vector resval,
                        N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                        void *rdataS,
                        N_Vector ytemp, N_Vector yptemp, N_Vector restemp);

extern int IDASensRes1DQ(int Ns, realtype t, 
                         N_Vector yy, N_Vector yp, N_Vector resval,
                         int iS,
                         N_Vector yyS, N_Vector ypS, N_Vector resvalS,
                         void *rdataS,
                         N_Vector ytemp, N_Vector yptemp, N_Vector restemp);

/*=================================================================*/
/*BEGIN        INTEGRATOR OPTIONAL INPUT FUNCTIONS                 */
/*=================================================================*/

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

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxord <= 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_MAXORD);
    return(IDA_ILL_INPUT);
  }

  if (maxord > IDA_mem->ida_maxord) {
    if(errfp!=NULL) 
      fprintf(errfp, MSG_IDAS_BAD_MAXORD);
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

  if (mxsteps <= 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_MXSTEPS);
    return(IDA_ILL_INPUT);
  }

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

  if (hmax <= 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_HMAX);
    return(IDA_ILL_INPUT);
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

  IDA_mem->ida_id = id;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetConstraints(void *ida_mem, N_Vector constraints)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_constraints = constraints;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetTolerances(void *ida_mem, 
                     int itol, realtype *rtol, void *atol)
{
  IDAMem IDA_mem;
  booleantype neg_atol;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if ((itol != IDA_SS) && (itol != IDA_SV)) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_ITOL);
    return(IDA_ILL_INPUT);
  }

  if (rtol == NULL) { 
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_RTOL_NULL); 
    return(IDA_ILL_INPUT); 
  }

  if (*rtol < ZERO) { 
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_RTOL); 
    return(IDA_ILL_INPUT); 
  }
   
  if (atol == NULL) { 
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_ATOL_NULL); 
    return(IDA_ILL_INPUT); 
  }

  /* Test absolute tolerances */
  if (itol == IDA_SS) { 
    neg_atol = (*((realtype *)atol) < ZERO); 
  } else { 
    neg_atol = (N_VMin((N_Vector)atol) < ZERO); 
  }
  if (neg_atol) { 
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_ATOL); 
    return(IDA_ILL_INPUT); 
  }

  IDA_mem->ida_itol   = itol;
  IDA_mem->ida_reltol = rtol;      
  IDA_mem->ida_abstol = atol;  

  return(IDA_SUCCESS);
}

/*=================================================================*/
/*       INITIAL CONDITION CALCULATION OPTIONAL INPUT FUNCTIONS    */
/*=================================================================*/

int IDASetNonlinConvFactorIC(void *ida_mem, realtype epiccon)
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

/*=================================================================*/
/*             QUADRATURE OPTIONAL INPUT FUNCTIONS                 */
/*=================================================================*/

int IDASetQuadErrCon(void *ida_mem, booleantype errconQ)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_errconQ = errconQ;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetQuadTolerances(void *ida_mem, int itolQ, 
                         realtype *reltolQ, void *abstolQ)
{
  IDAMem IDA_mem;
  booleantype neg_abstol;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if ((itolQ != IDA_SS) && (itolQ != IDA_SV)) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_ITOLQ);
    return(IDA_ILL_INPUT);
  }

  if (reltolQ == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_RTOLQ_NULL);
    return(IDA_ILL_INPUT);
  }
  
  if (*reltolQ < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_RTOLQ);
    return(IDA_ILL_INPUT);
  }

  if (abstolQ == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_ATOLQ_NULL);
    return(IDA_ILL_INPUT);
  }

  if (itolQ == IDA_SS) {
    neg_abstol = (*((realtype *)abstolQ) < ZERO);
  } else {
    neg_abstol = (N_VMin((N_Vector)abstolQ) < ZERO);
  }
  if (neg_abstol) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_BAD_ATOLQ);
    return(IDA_ILL_INPUT);
  }

  IDA_mem->ida_itolQ   = itolQ;
  IDA_mem->ida_reltolQ = reltolQ;
  IDA_mem->ida_abstolQ = abstolQ;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetQuadRdata(void *ida_mem, void *rdataQ)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rdataQ = rdataQ;

  return(IDA_SUCCESS);
}

/*=================================================================*/
/*             SENSITIVITY OPTIONAL INPUT FUNCTIONS                */
/*=================================================================*/

int IDASetSensResFn(void *ida_mem, IDASensResFn resS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_iresS = IDA_ALLSENS;

  if (resS != NULL) {
    IDA_mem->ida_resS   = resS;
    IDA_mem->ida_resSDQ = FALSE;
  } else {
    IDA_mem->ida_resS   = IDASensResDQ;
    IDA_mem->ida_rdataS = ida_mem;
    IDA_mem->ida_resSDQ = TRUE;
  }

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensRes1Fn(void *ida_mem, IDASensRes1Fn resS1)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_iresS = IDA_ONESENS;

  if (resS1 != NULL) {
    IDA_mem->ida_resS1  = resS1;
    IDA_mem->ida_resSDQ = FALSE;
  } else {
    IDA_mem->ida_resS1  = IDASensRes1DQ;
    IDA_mem->ida_rdataS = ida_mem;
    IDA_mem->ida_resSDQ = TRUE;
  }

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensErrCon(void *ida_mem, booleantype errconS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_errconS = errconS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensRho(void *ida_mem, realtype rho)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rhomax = rho;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensPbar(void *ida_mem, realtype *pbar)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_pbar = pbar;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensTolerances(void *ida_mem, int itolS, 
                         realtype *reltolS, void *abstolS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if ((itolS != IDA_SS) && (itolS != IDA_SV)) {
    if(errfp!=NULL) 
      fprintf(errfp, MSG_IDAS_BAD_ITOLS);
    return(IDA_ILL_INPUT);
  }

  if (itolS == IDA_EE) {

    /* IDAS will set tolerances */
    IDA_mem->ida_setSensTol = TRUE;
    IDA_mem->ida_testSensTol = FALSE;
    
  } else {

    /* Test user-supplied tolerances */
    if (reltolS == NULL) {
      if(errfp!=NULL) fprintf(errfp, MSG_IDAS_RTOLS_NULL);
      return(IDA_ILL_INPUT);
    }
    
    if (abstolS == NULL) {
      if(errfp!=NULL) fprintf(errfp, MSG_IDAS_ATOLS_NULL);
      return(IDA_ILL_INPUT);
    }
    
    IDA_mem->ida_itolS   = itolS;
    IDA_mem->ida_reltolS = reltolS;
    IDA_mem->ida_abstolS = abstolS;

    IDA_mem->ida_setSensTol = FALSE;
    IDA_mem->ida_testSensTol = TRUE;
  }

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensRdata(void *ida_mem, void *rdataS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rdataS = rdataS;

  return(IDA_SUCCESS);
}


/*-----------------------------------------------------------------*/

int IDASetSensMaxNonlinIters(void *ida_mem, int maxcorS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stderr, MSG_IDAS_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_maxcorS = maxcorS;

  return(IDA_SUCCESS);
}

#define lrw      (IDA_mem->ida_lrw)
#define liw      (IDA_mem->ida_liw)
#define nst      (IDA_mem->ida_nst)
#define nre      (IDA_mem->ida_nre)
#define ncfn     (IDA_mem->ida_ncfn)
#define netf     (IDA_mem->ida_netf)
#define nni      (IDA_mem->ida_nni)
#define nsetups  (IDA_mem->ida_nsetups)
#define kk       (IDA_mem->ida_kk)
#define hh       (IDA_mem->ida_hh)
#define h0u      (IDA_mem->ida_h0u)
#define kused    (IDA_mem->ida_kused)          
#define hused    (IDA_mem->ida_hused)         
#define tn       (IDA_mem->ida_tn)
#define ewt      (IDA_mem->ida_ewt)  
#define tolsf    (IDA_mem->ida_tolsf)  
#define nbacktr  (IDA_mem->ida_nbacktr)

#define quad     (IDA_mem->ida_quad)
#define nrQe     (IDA_mem->ida_nrQe)
#define netfQ    (IDA_mem->ida_netfQ)
#define ewtQ     (IDA_mem->ida_ewtQ)

#define sensi    (IDA_mem->ida_sensi)
#define ism      (IDA_mem->ida_ism)
#define ewtS     (IDA_mem->ida_ewtS)
#define nrSe     (IDA_mem->ida_nrSe)
#define nreS     (IDA_mem->ida_nreS)
#define nniS     (IDA_mem->ida_nniS)
#define ncfnS    (IDA_mem->ida_ncfnS)
#define netfS    (IDA_mem->ida_netfS)
#define nsetupsS (IDA_mem->ida_nsetupsS)
#define netfS1   (IDA_mem->ida_netfS1)
#define ncfnS1   (IDA_mem->ida_ncfnS1)
#define nniS1    (IDA_mem->ida_nniS1)

/*=================================================================*/
/*BEGIN        INTEGRATOR OPTIONAL INPUT FUNCTIONS                 */
/*=================================================================*/

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

int IDAGetErrWeights(void *ida_mem, N_Vector *eweight)
{
  IDAMem IDA_mem;
  
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem; 

  *eweight = ewt;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetWorkSpace(void *ida_mem,long int *lenrw, long int *leniw)
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

int IDAGetNonlinSolvStats(void *ida_mem, long int *nniters, 
                          long int *nncfails)
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

/*=================================================================*/
/*END          INTEGRATOR OPTIONAL INPUT FUNCTIONS                 */
/*=================================================================*/

/*=================================================================*/
/*BEGIN        QUADRATURE OPTIONAL INPUT FUNCTIONS                 */
/*=================================================================*/

int IDAGetQuadNumRhsEvals(void *ida_mem, long int *nrhsQevals)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if(quad != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_QUAD);
    return (IDA_NO_QUAD);
  }

  *nrhsQevals = nrQe;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetQuadNumErrTestFails(void *ida_mem, long int *nQetfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if(quad != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_QUAD);
    return (IDA_NO_QUAD);
  }

  *nQetfails = netfQ;

  return(IDA_SUCCESS);

}

/*-----------------------------------------------------------------*/

int IDAGetQuadErrWeights(void *ida_mem, N_Vector *eQweight)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if(quad != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_QUAD);
    return (IDA_NO_QUAD);
  }

  if(IDA_mem->ida_errconQ) *eQweight = ewtQ;
  else                     *eQweight = NULL;

  return(IDA_SUCCESS);

}

/*-----------------------------------------------------------------*/

int IDAGetQuadStats(void *ida_mem, long int *nrhsQevals, 
                    long int *nQetfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if(quad != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_QUAD);
    return (IDA_NO_QUAD);
  }

  *nrhsQevals = nrQe;
  *nQetfails  = netfQ;

  return(IDA_SUCCESS);

}

/*=================================================================*/
/*END          QUADRATURE OPTIONAL OUTPUT FUNCTIONS                */
/*=================================================================*/

/*=================================================================*/
/*BEGIN        SENSITIVITY OPTIONAL OUTPUT FUNCTIONS               */
/*=================================================================*/

int IDAGetSensNumRhsEvals(void *ida_mem, long int *nresSevals)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  *nresSevals = nrSe;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumRhsEvalsSens(void *ida_mem, long int *nresevalsS)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  *nresevalsS = nreS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensNumErrTestFails(void *ida_mem, long int *nSetfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  *nSetfails = netfS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensNumLinSolvSetups(void *ida_mem, long int *nlinsetupsS)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  *nlinsetupsS = nsetupsS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensErrWeights(void *ida_mem, N_Vector_S *eSweight)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  *eSweight = ewtS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensStats(void *ida_mem, long int *nresSevals, 
                    long int *nresevalsS, 
                    long int *nSetfails, long int *nlinsetupsS)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  *nresSevals  = nrSe;
  *nresevalsS  = nreS;
  *nSetfails   = netfS;
  *nlinsetupsS = nsetupsS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensNumNonlinSolvIters(void *ida_mem, long int *nSniters)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  *nSniters = nniS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensNumNonlinSolvConvFails(void *ida_mem, long int *nSncfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  *nSncfails = ncfnS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensNonlinSolvStats(void *ida_mem, long int *nSniters, 
                              long int *nSncfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  *nSniters  = nniS;
  *nSncfails = ncfnS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensNumStgrErrTestFails(void *ida_mem, long int *nSTGR1etfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  if (ism != IDA_STAGGERED1) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_STGR1);
    return (IDA_NO_STGR1);
  }

  nSTGR1etfails = netfS1;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensNumStgrNonlinSolvIters(void *ida_mem, long int *nSTGR1niters)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  if (ism != IDA_STAGGERED1) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_STGR1);
    return (IDA_NO_STGR1);
  }

  nSTGR1niters = nniS1;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensNumStgrNonlinSolvConvFails(void *ida_mem, long int *nSTGR1ncfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (sensi != TRUE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_SENS);
    return (IDA_NO_SENS);
  }

  if (ism != IDA_STAGGERED1) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_NO_STGR1);
    return (IDA_NO_STGR1);
  }

  nSTGR1ncfails = ncfnS1;  

  return(IDA_SUCCESS);
}

/*=================================================================*/
/*END          SENSITIVITY OPTIONAL OUTPUT FUNCTIONS               */
/*=================================================================*/

/*=================================================================*/
/*END          EXPORTED FUNCTIONS IMPLEMENTATION                   */
/*=================================================================*/
