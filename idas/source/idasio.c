/*******************************************************************
 * File          : idasio.c                                        *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 24 September 2003                               *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/idas/LICENSE                          *
 *-----------------------------------------------------------------*
 * This is the implementation file for the optional inputs and     *
 * outputs for the IDAS solver.                                    *
 *                                                                 *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "idas.h"
#include "sundialstypes.h"
#include "nvector.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/*=================================================================*/
/*BEGIN        IDAS Error Messages                                 */
/*=================================================================*/

/* IDASet* error messages */

#define MSG_IDAS_NO_MEM      "ida_mem=NULL in an IDASet routine illegal. \n\n"

#define MSG_IDAS_NEG_MAXORD  "IDASetMaxOrd-- maxord<=0 illegal. \n\n"

#define MSG_IDAS_BAD_MAXORD1 "IDASetMaxOrd-- Illegal attempt to increase "
#define MSG_IDAS_BAD_MAXORD2 "maximum method order from %d to %d.\n\n"
#define MSG_IDAS_BAD_MAXORD  MSG_IDAS_BAD_MAXORD1 MSG_IDAS_BAD_MAXORD2 

#define MSG_IDAS_NEG_MXSTEPS "IDASetMaxNumSteps-- mxsteps<=0 illegal. \n\n"

#define MSG_IDAS_NEG_HMAX    "IDASetMaxStep-- hmax<=0 illegal. \n\n"

#define MSG_IDAS_NEG_EPCON   "IDASetNlinConvCoef-- epcon < 0.0 illegal. \n\n"

#define MSG_IDAS_BAD_EPICCON "IDASetNlinConvCoefIC-- epiccon < 0.0 illegal.\n\n"

#define MSG_IDAS_BAD_MAXNH   "IDASetMaxNumStepsIC-- maxnh < 0 illegal.\n\n"

#define MSG_IDAS_BAD_MAXNJ   "IDASetMaxNumJacsIC-- maxnj < 0 illegal.\n\n"

#define MSG_IDAS_BAD_MAXNIT  "IDASetMaxNumItersIC-- maxnit < 0 illegal.\n\n"

#define MSG_IDAS_BAD_STEPTOL "IDASetLineSearchOffIC-- steptol < 0.0 illegal.\n\n"

#define MSG_BAD_ITOLQ1      "IDASetQuadTolerances-- itolQ=%d illegal.\n"
#define MSG_BAD_ITOLQ2      "The legal values are SS=%d and SV=%d.\n\n"
#define MSG_BAD_ITOLQ       MSG_BAD_ITOLQ1 MSG_BAD_ITOLQ2

#define MSG_BAD_ITOLS1      "IDASetSensTolerances-- itolS=%d illegal.\n"
#define MSG_BAD_ITOLS2      "The legal values are SS=%d and SV=%d.\n\n"
#define MSG_BAD_ITOLS       MSG_BAD_ITOLS1 MSG_BAD_ITOLS2

/* IDAGet* Error Messages */

#define MSG_IDAG_NO_MEM    "ida_mem=NULL in an IDAGet routine illegal. \n\n"

#define MSG_IDAG_NO_QUAD1  "IDAGetQuad*-- Illegal attempt to call before "
#define MSG_IDAG_NO_QUAD2  "calling IDAQuadMalloc.\n\n"
#define MSG_IDAG_NO_QUAD   MSG_IDAG_NO_QUAD1 MSG_IDAG_NO_QUAD2

#define MSG_IDAG_NO_SENSI1 "IDAGetSens*-- Illegal attempt to call before "
#define MSG_IDAG_NO_SENSI2 "calling IDASensMalloc.\n\n"
#define MSG_IDAG_NO_SENSI  MSG_IDAG_NO_SENSI1 MSG_IDAG_NO_SENSI2

/*=================================================================*/
/*END          IDAS Error Messages                                 */
/*=================================================================*/

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
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_errfp = errfp;

  return(SUCCESS);
}

#define errfp (IDA_mem->ida_errfp)

/*-----------------------------------------------------------------*/

int IDASetRdata(void *ida_mem, void *rdata)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rdata = rdata;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxOrd(void *ida_mem, int maxord)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxord <= 0) {
    fprintf(errfp, MSG_IDAS_NEG_MAXORD);
    return(IDAS_ILL_INPUT);
  }

  if (maxord > IDA_mem->ida_maxord) {
    fprintf(errfp, MSG_IDAS_BAD_MAXORD, IDA_mem->ida_maxord, maxord);
    return(IDAS_ILL_INPUT);
  }  

  IDA_mem->ida_maxord = maxord;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumSteps(void *ida_mem, int mxsteps)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (mxsteps <= 0) {
    fprintf(errfp, MSG_IDAS_NEG_MXSTEPS);
    return(IDAS_ILL_INPUT);
  }

  IDA_mem->ida_mxstep = mxsteps;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetInitStep(void *ida_mem, realtype hin)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_hin = hin;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxStep(void *ida_mem, realtype hmax)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (hmax <= 0) {
    fprintf(errfp, MSG_IDAS_NEG_HMAX);
    return(IDAS_ILL_INPUT);
  }

  IDA_mem->ida_hmax_inv = ONE/hmax;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetStopTime(void *ida_mem, realtype tstop)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_tstop = tstop;
  IDA_mem->ida_tstopset = TRUE;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetNlinConvCoef(void *ida_mem, realtype epcon)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (epcon < ZERO) {
    fprintf(errfp, MSG_IDAS_NEG_EPCON);
    return(IDAS_ILL_INPUT);
  }

  IDA_mem->ida_epcon = epcon;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxErrTestFails(void *ida_mem, int maxnef)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return (IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_maxnef = maxnef;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxConvFails(void *ida_mem, int maxncf)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return (IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_maxncf = maxncf;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNonlinIters(void *ida_mem, int maxcor)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return (IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_maxcor = maxcor;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_suppressalg = suppressalg;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetID(void *ida_mem, N_Vector id)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_id = id;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetConstraints(void *ida_mem, N_Vector constraints)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_constraints = constraints;

  return(SUCCESS);
}

/*=================================================================*/
/*END        INTEGRATOR OPTIONAL INPUT FUNCTIONS                   */
/*=================================================================*/

/*=================================================================*/
/*BEGIN  INITIAL CONDITION CALCULATION OPTIONAL INPUT FUNCTIONS    */
/*=================================================================*/

int IDASetNlinConvFactorIC(void *ida_mem, realtype epiccon)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (epiccon < ZERO) {
    fprintf(errfp, MSG_IDAS_BAD_EPICCON);
    return(IDAS_ILL_INPUT);
  }

  IDA_mem->ida_epiccon = epiccon;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumStepsIC(void *ida_mem, int maxnh)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxnh < 0) {
    fprintf(errfp, MSG_IDAS_BAD_MAXNH);
    return(IDAS_ILL_INPUT);
  }

  IDA_mem->ida_maxnh = maxnh;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumJacsIC(void *ida_mem, int maxnj)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

   if (maxnj < 0) {
    fprintf(errfp, MSG_IDAS_BAD_MAXNJ);
    return(IDAS_ILL_INPUT);
  } 

  IDA_mem->ida_maxnj = maxnj;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumItersIC(void *ida_mem, int maxnit)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxnit < 0) {
    fprintf(errfp, MSG_IDAS_BAD_MAXNIT);
    return(IDAS_ILL_INPUT);
  }

  IDA_mem->ida_maxnit = maxnit;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_lsoff = lsoff;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetStepToleranceIC(void *ida_mem, realtype steptol)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (steptol < ZERO) {
    fprintf(errfp, MSG_IDAS_BAD_STEPTOL);
    return(IDAS_ILL_INPUT);
  }

  IDA_mem->ida_steptol = steptol;

  return(SUCCESS);
}

/*=================================================================*/
/*END  INITIAL CONDITION CALCULATION OPTIONAL INPUT FUNCTIONS      */
/*=================================================================*/

/*=================================================================*/
/*BEGIN        QUADRATURE OPTIONAL INPUT FUNCTIONS                 */
/*=================================================================*/

int IDASetQuadErrCon(void *ida_mem, booleantype errconQ)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_errconQ = errconQ;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetQuadTolerances(void *ida_mem, int itolQ, 
                         realtype *reltolQ, void *abstolQ)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  if ((itolQ != SS) && (itolQ != SV)) {
    fprintf(errfp, MSG_BAD_ITOLQ, itolQ, SS, SV);
    return(IDAS_ILL_INPUT);
  }

  IDA_mem->ida_itolQ   = itolQ;
  IDA_mem->ida_reltolQ = reltolQ;
  IDA_mem->ida_abstolQ = abstolQ;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetQuadRdata(void *ida_mem, void *rdataQ)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rdataQ = rdataQ;

  return(SUCCESS);
}

/*=================================================================*/
/*END        QUADRATURE OPTIONAL INPUT FUNCTIONS                   */
/*=================================================================*/

/*=================================================================*/
/*BEGIN        SENSITIVITY OPTIONAL INPUT FUNCTIONS                */
/*=================================================================*/

int IDASetSensResFn(void *ida_mem, SensResFn resS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_iresS = ALLSENS;

  if (resS != NULL) {
    IDA_mem->ida_resS   = resS;
    IDA_mem->ida_resSDQ = FALSE;
  } else {
    IDA_mem->ida_resS   = IDASensResDQ;
    IDA_mem->ida_rdataS = ida_mem;
    IDA_mem->ida_resSDQ = TRUE;
  }

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensRes1Fn(void *ida_mem, SensRes1Fn resS1)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_iresS = ONESENS;

  if (resS1 != NULL) {
    IDA_mem->ida_resS1  = resS1;
    IDA_mem->ida_resSDQ = FALSE;
  } else {
    IDA_mem->ida_resS1  = IDASensRes1DQ;
    IDA_mem->ida_rdataS = ida_mem;
    IDA_mem->ida_resSDQ = TRUE;
  }

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensErrCon(void *ida_mem, booleantype errconS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_errconS = errconS;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensRho(void *ida_mem, realtype rho)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rhomax = rho;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensPbar(void *ida_mem, realtype *pbar)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_pbar = pbar;

  return(SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensTolerances(void *ida_mem, int itolS, 
                         realtype *reltolS, void *abstolS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  if ((itolS != SS) && (itolS != SV)) {
    fprintf(errfp, MSG_BAD_ITOLS, itolS, SS, SV);
    return(IDAS_ILL_INPUT);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_itolS   = itolS;
  IDA_mem->ida_reltolS = reltolS;
  IDA_mem->ida_abstolS = abstolS;

  return(SUCCESS);
}


/*-----------------------------------------------------------------*/

int IDASetSensRdata(void *ida_mem, void *rdataS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return(IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rdataS = rdataS;

  return(SUCCESS);
}


/*-----------------------------------------------------------------*/

int IDASetSensMaxNonlinIters(void *ida_mem, int maxcorS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAS_NO_MEM);
    return (IDAS_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_maxcorS = maxcorS;

  return(SUCCESS);
}

/*=================================================================*/
/*END        SENSITIVITY OPTIONAL INPUT FUNCTIONS                  */
/*=================================================================*/

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

/*-----------------------------------------------------------------*/

int IDAGetIntWorkSpace(void *ida_mem, long int *leniw)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *leniw = liw;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetRealWorkSpace(void *ida_mem, long int *lenrw)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *lenrw = lrw;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNumSteps(void *ida_mem, int *nsteps)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nsteps = nst;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNumResEvals(void *ida_mem, int *nrevals)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nrevals = nre;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNumLinSolvSetups(void *ida_mem, int *nlinsetups)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nlinsetups = nsetups;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNumErrTestFails(void *ida_mem, int *netfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *netfails = netf;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNumBacktrackOps(void *ida_mem, int *nbacktracks)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nbacktracks = nbacktr;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetLastOrder(void *ida_mem, int *klast)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *klast = kused;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentOrder(void *ida_mem, int *kcur)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *kcur = kk;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetActualInitStep(void *ida_mem, realtype *hinused)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *hinused = h0u;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetLastStep(void *ida_mem, realtype *hlast)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *hlast = hused;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentStep(void *ida_mem, realtype *hcur)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *hcur = hh;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentTime(void *ida_mem, realtype *tcur)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *tcur = tn;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetTolScaleFactor(void *ida_mem, realtype *tolsfact)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *tolsfact = tolsf;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetErrWeights(void *ida_mem, N_Vector *eweight)
{
  IDAMem IDA_mem;
  
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return (IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem; 

  *eweight = ewt;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetWorkSpace(void *ida_mem,long int *leniw, long int *lenrw)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *leniw = liw;
  *lenrw = lrw;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetIntegratorStats(void *ida_mem, int *nsteps, int *nrevals, 
                          int *nlinsetups, int *netfails,
                          int *klast, int *kcur, realtype *hlast, 
                          realtype *hcur, realtype *tcur)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
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

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNumNonlinSolvIters(void *ida_mem, int *nniters)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nniters = nni;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNumNonlinSolvConvFails(void *ida_mem, int *nncfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nncfails = ncfn;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNonlinSolvStats(void *ida_mem, int *nniters, int *nncfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return(IDAG_NO_MEM);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nniters  = nni;
  *nncfails = ncfn;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNumQuadRhsEvals(void *ida_mem, int *nrhsQevals)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return (IDAG_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if(quad != TRUE) {
    fprintf(errfp, MSG_IDAG_NO_QUAD);
    return (IDAG_NO_QUAD);
  }

  *nrhsQevals = nrQe;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int IDAGetNumQuadErrTestFails(void *ida_mem, int *nQetfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return (IDAG_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if(quad != TRUE) {
    fprintf(errfp, MSG_IDAG_NO_QUAD);
    return (IDAG_NO_QUAD);
  }

  *nQetfails = netfQ;

  return(OKAY);

}

/*-----------------------------------------------------------------*/

int IDAGetQuadErrWeights(void *ida_mem, N_Vector *eQweight)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return (IDAG_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if(quad != TRUE) {
    fprintf(errfp, MSG_IDAG_NO_QUAD);
    return (IDAG_NO_QUAD);
  }

  if(IDA_mem->ida_errconQ) *eQweight = ewtQ;
  else                     *eQweight = NULL;

  return(OKAY);

}

/*-----------------------------------------------------------------*/

int IDAGetQuadStats(void *ida_mem, int *nrhsQevals, int *nQetfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    fprintf(stdout, MSG_IDAG_NO_MEM);
    return (IDAG_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if(quad != TRUE) {
    fprintf(errfp, MSG_IDAG_NO_QUAD);
    return (IDAG_NO_QUAD);
  }

  *nrhsQevals = nrQe;
  *nQetfails = netfQ;

  return(OKAY);

}

/*=================================================================*/
/*END          EXPORTED FUNCTIONS IMPLEMENTATION                   */
/*=================================================================*/
