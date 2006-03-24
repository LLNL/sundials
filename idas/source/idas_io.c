/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-03-24 15:57:25 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
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
#include "sundials_math.h"
#include "sundials_types.h"

#define ZERO    RCONST(0.0)
#define HALF    RCONST(0.5)
#define ONE     RCONST(1.0)
#define TWOPT5  RCONST(2.5)

/* 
 * =================================================================
 * IDA optional input functions
 * =================================================================
 */

/* 
 * Readability constants
 */

#define lrw  (IDA_mem->ida_lrw)
#define liw  (IDA_mem->ida_liw)
#define lrw1 (IDA_mem->ida_lrw1)
#define liw1 (IDA_mem->ida_liw1)

int IDASetErrHandlerFn(void *ida_mem, IDAErrHandlerFn ehfun, void *eh_data)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetErrHandlerFn", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_ehfun = ehfun;
  IDA_mem->ida_eh_data = eh_data;

  return(IDA_SUCCESS);
}


int IDASetErrFile(void *ida_mem, FILE *errfp)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetErrFile", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_errfp = errfp;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetRdata(void *ida_mem, void *res_data)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetRdata", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxOrd", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxord <= 0) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxOrd", MSG_NEG_MAXORD);
    return(IDA_ILL_INPUT);
  }

  /* Cannot increase maximum order beyond the value that
     was used when allocating memory */
  maxord_alloc = IDA_mem->ida_maxord_alloc;

  if (maxord > maxord_alloc) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxOrd", MSG_BAD_MAXORD);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNumSteps", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (mxsteps < 0) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxNumSteps", MSG_NEG_MXSTEPS);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetInitStep", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxStep", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (hmax < 0) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxStep", MSG_NEG_HMAX);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetStopTime", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetNonlinConvCoef", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (epcon < ZERO) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetNonlinConvCoef", MSG_NEG_EPCON);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxErrTestFails", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxConvFails", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNonlinIters", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetSuppressAlg", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetId", MSG_NO_MEM);
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
    IDA_mem->ida_id = NULL;
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetConstraints", MSG_NO_MEM);
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
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetConstraints", MSG_BAD_NVECTOR);
    return(IDA_ILL_INPUT);
  }

  /*  Check the constraints vector */

  temptest = N_VMaxNorm(constraints);
  if((temptest > TWOPT5) || (temptest < HALF)){ 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetConstraints", MSG_BAD_CONSTR);
    return(IDA_ILL_INPUT); 
  }

  if ( !(IDA_mem->ida_constraintsMallocDone) ) {
    IDA_mem->ida_constraints = NULL;
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetTolerances", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  /* Check if ida_mem was allocated */

  if (IDA_mem->ida_MallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_MALLOC, "IDA", "IDASetTolerances", MSG_NO_MALLOC);
    return(IDA_NO_MALLOC);
  }

  /* Check inputs */

  if ((itol != IDA_SS) && (itol != IDA_SV)) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetTolerances", MSG_BAD_ITOL);
    return(IDA_ILL_INPUT);
  }

  if (atol == NULL) { 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetTolerances", MSG_ATOL_NULL);
    return(IDA_ILL_INPUT); 
  }

  if (rtol < ZERO) { 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetTolerances", MSG_BAD_RTOL);
    return(IDA_ILL_INPUT); 
  }

    
  if (itol == IDA_SS) { 
    neg_atol = (*((realtype *)atol) < ZERO); 
  } else { 
    neg_atol = (N_VMin((N_Vector)atol) < ZERO); 
  }

  if (neg_atol) { 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetTolerances", MSG_BAD_ATOL);
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
    IDA_mem->ida_Vatol = NULL;
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetEwtFn", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetNonlinConvCoefIC", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (epiccon < ZERO) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetNonlinConvCoefIC", MSG_BAD_EPICCON);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNumStepsIC", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxnh < 0) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxNumStepsIC", MSG_BAD_MAXNH);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNumJacsIC", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

   if (maxnj < 0) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxNumJacsIC", MSG_BAD_MAXNJ);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNumItersIC", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (maxnit < 0) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxNumItersIC", MSG_BAD_MAXNIT);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetLineSearchOffIC", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetStepToleranceIC", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (steptol < ZERO) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetStepToleranceIC", MSG_BAD_STEPTOL);
    return(IDA_ILL_INPUT);
  }

  IDA_mem->ida_steptol = steptol;

  return(IDA_SUCCESS);
}

/* 
 * =================================================================
 * Quadrature optional input functions
 * =================================================================
 */

/* 
 * Readability constants
 */

#define lrw1Q (IDA_mem->ida_lrw1Q)
#define liw1Q (IDA_mem->ida_liw1Q)

/*-----------------------------------------------------------------*/

int IDASetQuadRdata(void *ida_mem, void *rhsQ_data)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDASetQuadRdata", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rdataQ = rhsQ_data;

  return(IDA_SUCCESS);
}

int IDASetQuadErrCon(void *ida_mem, booleantype errconQ, 
                       int itolQ, realtype rtolQ, void *atolQ)
{
  IDAMem IDA_mem;
  booleantype neg_atol;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDASetQuadErrCon", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }
  
  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_errconQ = errconQ;

  /* Ckeck if quadrature was initialized? */

  if (IDA_mem->ida_quadMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_QUAD, "IDAS", "IDASetQuadErrCon", MSG_NO_QUAD); 
    return(IDA_NO_QUAD);
  }

  /* Check inputs */

  if(errconQ == FALSE) {
    if (IDA_mem->ida_VatolQMallocDone) {
      N_VDestroy(IDA_mem->ida_VatolQ);
      lrw -= lrw1Q;
      liw -= liw1Q;
      IDA_mem->ida_VatolQMallocDone = FALSE;
    }
    return(IDA_SUCCESS);
  }
  
  if ((itolQ != IDA_SS) && (itolQ != IDA_SV)) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetQuadErrCon", MSG_BAD_ITOLQ);
    return(IDA_ILL_INPUT);
  }

  if (atolQ == NULL) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetQuadErrCon", MSG_NULL_ATOLQ);
    return(IDA_ILL_INPUT);
  }

  if (rtolQ < ZERO) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetQuadErrCon", MSG_BAD_RTOLQ);
    return(IDA_ILL_INPUT);
  }

  if (itolQ == IDA_SS)
    neg_atol = (*((realtype *)atolQ) < ZERO);
  else
    neg_atol = (N_VMin((N_Vector)atolQ) < ZERO);

  if (neg_atol) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetQuadErrCon", MSG_BAD_ATOLQ);
    return(IDA_ILL_INPUT);
  }

  /* See if we need to free or allocate memory */

  if ( (itolQ != IDA_SV) && (IDA_mem->ida_VatolQMallocDone) ) {
    N_VDestroy(IDA_mem->ida_VatolQ);
    lrw -= lrw1Q;
    liw -= liw1Q;
    IDA_mem->ida_VatolQMallocDone = FALSE;
  }

  if ( (itolQ == IDA_SV) && !(IDA_mem->ida_VatolQMallocDone) ) {
    IDA_mem->ida_VatolQ = NULL;
    IDA_mem->ida_VatolQ = N_VClone(IDA_mem->ida_ypQ);
    lrw += lrw1Q;
    liw += liw1Q;
    IDA_mem->ida_VatolQMallocDone = TRUE;
  }

  /* Copy tolerances into memory */

  IDA_mem->ida_itolQ = itolQ;
  IDA_mem->ida_rtolQ = rtolQ;

  if (itolQ == IDA_SS)
    IDA_mem->ida_SatolQ = *((realtype *)atolQ);
  else
    N_VScale(ONE, (N_Vector)atolQ, IDA_mem->ida_VatolQ);
  
  return(IDA_SUCCESS);
}


/* 
 * =================================================================
 * FSA optional input functions
 * =================================================================
 */

int IDASetSensResFn(void *ida_mem, IDASensResFn resS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDASetSensResFn", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (resS != NULL) {
    IDA_mem->ida_resS    = resS;
    IDA_mem->ida_resSDQ  = FALSE;
  } else {
    IDA_mem->ida_resS    = IDASensResDQ;
    IDA_mem->ida_rdataS  = ida_mem;
    IDA_mem->ida_resSDQ  = TRUE;
  }

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensRdata(void *ida_mem, void *resS_data)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDASetSensRdata", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rdataS = resS_data;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensRho(void *ida_mem, realtype rho)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDASetSensRho", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_rhomax = rho;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensErrCon(void *ida_mem, booleantype errconS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDASetSensErrCon", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_errconS = errconS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensMaxNonlinIters(void *ida_mem, int maxcorS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDASetSensMaxNonlinIters", MSG_NO_MEM);    
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  IDA_mem->ida_maxcorS = maxcorS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensParams(void *ida_mem, realtype *p, realtype *pbar, int *plist)
{
  IDAMem IDA_mem;
  int is, Ns;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDASetSensParams", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  /* Was sensitivity initialized? */

  if (IDA_mem->ida_sensMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDASetSensParams", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  Ns = IDA_mem->ida_Ns;

  /* Parameters */

  IDA_mem->ida_p = p;

  /* pbar */

  if (pbar != NULL)
    for (is=0; is<Ns; is++) {
      if (pbar[is] == ZERO) {
        IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetSensParams", MSG_BAD_PBAR);
        return(IDA_ILL_INPUT);
      }
      IDA_mem->ida_pbar[is] = ABS(pbar[is]);
    }
  else
    for (is=0; is<Ns; is++)
      IDA_mem->ida_pbar[is] = ONE;

  /* plist */

  if (plist != NULL)
    for (is=0; is<Ns; is++) {
      if ( plist[is] < 0 ) {
        IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetSensParams", MSG_BAD_PLIST);
        return(IDA_ILL_INPUT);
      }
      IDA_mem->ida_plist[is] = plist[is];
    }
  else
    for (is=0; is<Ns; is++)
      IDA_mem->ida_plist[is] = is;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetSensTolerances(void *ida_mem, int itolS,
                         realtype rtolS, void *atolS)
{
  IDAMem IDA_mem;
  booleantype neg_atol;
  realtype *atolSS;
  N_Vector *atolSV;
  int is, Ns;

  atolSS = NULL;
  atolSV = NULL;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDASetSensTolerances", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  /* Was sensitivity initialized? */

  if (IDA_mem->ida_sensMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDASetSensTolerances", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  } 

  /* Check inputs */

  Ns = IDA_mem->ida_Ns;

  if ((itolS != IDA_SS) && (itolS != IDA_SV) && (itolS != IDA_EE)) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetSensTolerances", MSG_BAD_ITOLS);
    return(IDA_ILL_INPUT);
  }

  if (itolS != IDA_EE) {

    /* Test user-supplied tolerances */
    
    if (rtolS < ZERO) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetSensTolerances", MSG_BAD_RTOLS);
      return(IDA_ILL_INPUT);
    }

    if (atolS == NULL) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetSensTolerances", MSG_NULL_ATOLS);
      return(IDA_ILL_INPUT);
    }

    neg_atol = FALSE;

    if (itolS == IDA_SS) {
      atolSS = (realtype *) atolS;
      for (is=0; is<Ns; is++)
        if (atolSS[is] < ZERO) {neg_atol = TRUE; break;}
    } else {
      atolSV = (N_Vector *) atolS;
      for (is=0; is<Ns; is++) 
        if (N_VMin(atolSV[is]) < ZERO) {neg_atol = TRUE; break;}
    }

    if (neg_atol) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAS", "IDASetSensTolerances", MSG_BAD_ATOLS);
      return(IDA_ILL_INPUT);
    }
    
  }

  /* See if we should release some memory */

  if ( (itolS != IDA_SV) && (IDA_mem->ida_VatolSMallocDone) ) {
    N_VDestroyVectorArray(IDA_mem->ida_VatolS, Ns);
    lrw -= Ns*lrw1;
    liw -= Ns*liw1;
    IDA_mem->ida_VatolSMallocDone = FALSE;
  }

  if ( (itolS != IDA_SS) && (IDA_mem->ida_SatolSMallocDone) ) {
    free(IDA_mem->ida_SatolS); IDA_mem->ida_SatolS = NULL;
    lrw -= Ns;
    IDA_mem->ida_SatolSMallocDone = FALSE;
  }

  /* If tolerances will be estimated, return now */

  if (itolS == IDA_EE) return(IDA_SUCCESS);

  /* See if we need to allocate some memory */

  if ( (itolS == IDA_SV) && !(IDA_mem->ida_VatolSMallocDone) ) {
    IDA_mem->ida_VatolS = NULL;
    IDA_mem->ida_VatolS = N_VCloneVectorArray(Ns, IDA_mem->ida_tempv1);
    lrw += Ns*lrw1;
    liw += Ns*liw1;
    IDA_mem->ida_VatolSMallocDone = TRUE;
  }

  if ( (itolS == IDA_SS) && !(IDA_mem->ida_SatolSMallocDone) ) {
    IDA_mem->ida_SatolS = NULL;
    IDA_mem->ida_SatolS = (realtype *)malloc(Ns*sizeof(realtype));
    lrw += Ns;
    IDA_mem->ida_SatolSMallocDone = TRUE;
  }

  /* Copy tolerances into memory */

  IDA_mem->ida_itolS = itolS;
  IDA_mem->ida_rtolS = rtolS;

  if (itolS == IDA_SS)
    for (is=0; is<Ns; is++)
      IDA_mem->ida_SatolS[is] = atolSS[is];
  else
    for (is=0; is<Ns; is++)    
      N_VScale(ONE, atolSV[is], IDA_mem->ida_VatolS[is]);

  return(IDA_SUCCESS);
}


/* 
 * =================================================================
 * IDA optional output functions
 * =================================================================
 */

/* 
 * Readability constants
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
#define ee          (IDA_mem->ida_ee)

int IDAGetNumSteps(void *ida_mem, long int *nsteps)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumSteps", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumResEvals", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumLinSolvSetups", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumErrTestFails", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumBacktrackOps", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetLastOrder", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetCurrentOrder", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetActualInitStep", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetLastStep", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetCurrentStep", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetCurrentTime", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetTolScaleFactor", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetErrWeights", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem; 

  N_VScale(ONE, ewt, eweight);

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetEstLocalErrors(void *ida_mem, N_Vector ele)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetEstLocalErrors", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  N_VScale(ONE, ee, ele);

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetWorkSpace(void *ida_mem, long int *lenrw, long int *leniw)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetWorkSpace", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetIntegratorStats", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumGEvals", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetRootInfo", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumNonlinSolvIters", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumNonlinSolvConvFails", MSG_NO_MEM);
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
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNonlinSolvStats", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  *nniters  = nni;
  *nncfails = ncfn;

  return(IDA_SUCCESS);
}

/* 
 * =================================================================
 * Quadrature optional output functions
 * =================================================================
 */

/* 
 * Readability constants
 */

#define quadr          (IDA_mem->ida_quadr)
#define nrQe           (IDA_mem->ida_nrQe)
#define netfQ          (IDA_mem->ida_netfQ)
#define ewtQ           (IDA_mem->ida_ewtQ)
#define errconQ        (IDA_mem->ida_errconQ)

/*-----------------------------------------------------------------*/

int IDAGetQuadNumRhsEvals(void *ida_mem, long int *nrQevals)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetQuadNumRhsEvals", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (quadr==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_QUAD, "IDAS", "IDAGetQuadNumRhsEvals", MSG_NO_QUAD); 
    return(IDA_NO_QUAD);
  }

  *nrQevals = nrQe;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetQuadNumErrTestFails(void *ida_mem, long int *nQetfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetQuadNumErrTestFails", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (quadr==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_QUAD, "IDAS", "IDAGetQuadNumErrTestFails", MSG_NO_QUAD); 
    return(IDA_NO_QUAD);
  }

  *nQetfails = netfQ;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetQuadErrWeights(void *ida_mem, N_Vector eQweight)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetQuadErrWeights", MSG_NO_MEM); 
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (quadr==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_QUAD, "IDAS", "IDAGetQuadErrWeights", MSG_NO_QUAD); 
    return(IDA_NO_QUAD);
  }

  if(errconQ) N_VScale(ONE, ewtQ, eQweight);

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetQuadStats(void *ida_mem, long int *nrQevals, long int *nQetfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetQuadStats", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (quadr==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_QUAD, "IDAS", "IDAGetQuadStats", MSG_NO_QUAD); 
    return(IDA_NO_QUAD);
  }

  *nrQevals = nrQe;
  *nQetfails = netfQ;

  return(IDA_SUCCESS);
}


/* 
 * =================================================================
 * FSA optional output functions
 * =================================================================
 */

/* 
 * Readability constants
 */

#define sensi          (IDA_mem->ida_sensi)
#define ism            (IDA_mem->ida_ism)
#define ewtS           (IDA_mem->ida_ewtS)
#define nrSe           (IDA_mem->ida_nrSe)
#define nreS           (IDA_mem->ida_nreS)
#define nniS           (IDA_mem->ida_nniS)
#define ncfnS          (IDA_mem->ida_ncfnS)
#define netfS          (IDA_mem->ida_netfS)
#define nsetupsS       (IDA_mem->ida_nsetupsS)

/*-----------------------------------------------------------------*/

int IDAGetNumSensResEvals(void *ida_mem, long int *nrSevals)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetNumSensResEvals", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (sensi==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDAGetNumSensResEvals", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  *nrSevals = nrSe;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumResEvalsSens(void *ida_mem, long int *nrevalsS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetNumResEvalsSens", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (sensi==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDAGetNumResEvalsSens", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  *nrevalsS = nreS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumSensErrTestFails(void *ida_mem, long int *nSetfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetNumSensErrTestFails", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (sensi==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDAGetNumSensErrTestFails", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  *nSetfails = netfS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumSensLinSolvSetups(void *ida_mem, long int *nlinsetupsS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetNumSensLinSolvSetups", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (sensi==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDAGetNumSensLinSolvSetups", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  *nlinsetupsS = nsetupsS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensErrWeights(void *ida_mem, N_Vector_S eSweight)
{
  IDAMem IDA_mem;
  int is, Ns;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetSensErrWeights", MSG_NO_MEM);    
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (sensi==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDAGetSensErrWeights", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  Ns = IDA_mem->ida_Ns;

  for (is=0; is<Ns; is++)
    N_VScale(ONE, ewtS[is], eSweight[is]);

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensStats(void *ida_mem, long int *nrSevals, long int *nrevalsS, 
                      long int *nSetfails, long int *nlinsetupsS)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetSensStats", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (sensi==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDAGetSensStats", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  *nrSevals = nrSe;
  *nrevalsS = nreS;
  *nSetfails = netfS;
  *nlinsetupsS = nsetupsS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumSensNonlinSolvIters(void *ida_mem, long int *nSniters)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetNumSensNonlinSolvIters", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (sensi==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDAGetNumSensNonlinSolvIters", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  *nSniters = nniS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumSensNonlinSolvConvFails(void *ida_mem, long int *nSncfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetNumSensNonlinSolvConvFails", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (sensi==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDAGetNumSensNonlinSolvConvFails", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  *nSncfails = ncfnS;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetSensNonlinSolvStats(void *ida_mem, long int *nSniters, long int *nSncfails)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAS", "IDAGetSensNonlinSolvstats", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  if (sensi==FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_SENS, "IDAS", "IDAGetSensNonlinSolvStats", MSG_NO_SENSI);
    return(IDA_NO_SENS);
  }

  *nSniters = nniS;
  *nSncfails = ncfnS;

  return(IDA_SUCCESS);
}

/* 
 * =================================================================
 * IDAGetReturnFlagName
 * =================================================================
 */


char *IDAGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(24*sizeof(char));

  switch(flag) {
  case IDA_SUCCESS:
    sprintf(name,"IDA_SUCCESS");
    break;
  case IDA_TSTOP_RETURN:
    sprintf(name,"IDA_TSTOP_RETURN");
    break;
  case IDA_ROOT_RETURN:
    sprintf(name,"IDA_ROOT_RETURN");
    break;
  case IDA_TOO_MUCH_WORK:
    sprintf(name,"IDA_TOO_MUCH_WORK");
    break;
  case IDA_TOO_MUCH_ACC:
    sprintf(name,"IDA_TOO_MUCH_ACC");
    break;
  case IDA_ERR_FAIL:
    sprintf(name,"IDA_ERR_FAIL");
    break;
  case IDA_CONV_FAIL:
    sprintf(name,"IDA_CONV_FAIL");
    break;
  case IDA_LINIT_FAIL:
    sprintf(name,"IDA_LINIT_FAIL");
    break;
  case IDA_LSETUP_FAIL:
    sprintf(name,"IDA_LSETUP_FAIL");
    break;
  case IDA_LSOLVE_FAIL:
    sprintf(name,"IDA_LSOLVE_FAIL");
    break;
  case IDA_CONSTR_FAIL:
    sprintf(name,"IDA_CONSTR_FAIL");
    break;
  case IDA_RES_FAIL:
    sprintf(name,"IDA_RES_FAIL");
    break;
  case IDA_FIRST_RES_FAIL:
    sprintf(name,"IDA_FIRST_RES_FAIL");
    break;
  case IDA_REP_RES_ERR:
    sprintf(name,"IDA_REP_RES_ERR");
    break;
  case IDA_RTFUNC_FAIL:
    sprintf(name,"IDA_RTFUNC_FAIL");
    break;
  case IDA_MEM_FAIL:
    sprintf(name,"IDA_MEM_FAIL");
    break;
  case IDA_MEM_NULL:
    sprintf(name,"IDA_MEM_NULL");
    break;
  case IDA_ILL_INPUT:
    sprintf(name,"IDA_ILL_INPUT");
    break;
  case IDA_NO_MALLOC:
    sprintf(name,"IDA_NO_MALLOC");
    break;
  case IDA_BAD_T:
    sprintf(name,"IDA_BAD_T");
    break;
  case IDA_BAD_EWT:
    sprintf(name,"IDA_BAD_EWT");
    break;
  case IDA_NO_RECOVERY:
    sprintf(name,"IDA_NO_RECOVERY");
    break;
  case IDA_LINESEARCH_FAIL:
    sprintf(name,"IDA_LINESEARCH_FAIL");
    break;

  default:
    sprintf(name,"NONE");
  }

  return(name);
}

