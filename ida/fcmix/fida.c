/*
 * -----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2005-08-15 18:06:46 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the IDA package. See fida.h for usage.
 * NOTE: Some routines are necessarily stored elsewhere to avoid
 * linking problems.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "idaband.h"        /* prototypes for IDABAND interface routines      */
#include "idadense.h"       /* prototypes for IDADENSE interface routines     */
#include "ida.h"            /* IDA constants and prototypes                   */
#include "idasptfqmr.h"     /* prototypes for IDASPTFQMR interface routines   */
#include "idaspbcg.h"       /* prototypes for IDASPBCG interface routines     */
#include "idaspgmr.h"       /* prototypes for IDASPGMR interface routines     */
#include "fida.h"           /* actual function names, prototypes and global
			       variables                                      */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */

/*************************************************/

/* Definitions for global variables shared amongst various routines */

N_Vector F2C_IDA_ypvec, F2C_IDA_ewtvec;

void *IDA_idamem;
long int *IDA_iout;
realtype *IDA_rout;
int IDA_ls;
int IDA_nrtfn;

/*************************************************/

/* private constant(s) */
#define ZERO RCONST(0.0)

/*************************************************/

/* Prototype of user-supplied Fortran routine (IDAResFn) */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FIDA_RESFUN(realtype*, realtype*, realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_MALLOC(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 long int *iout, realtype *rout, 
                 int *ier)
{
  int itol;
  N_Vector Vatol;
  void *atolptr;

  /* Check for required vector operations */
  if ((F2C_IDA_vec->ops->nvgetarraypointer == NULL) ||
      (F2C_IDA_vec->ops->nvsetarraypointer == NULL)) {
    *ier = -1;
    printf("A required vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize all pointers to NULL */
  IDA_idamem = NULL;
  Vatol = NULL;
  atolptr = NULL;
  F2C_IDA_ypvec = F2C_IDA_ewtvec = NULL;

  /* Create IDA object */
  IDA_idamem = IDACreate();
  if (IDA_idamem == NULL) {
    *ier = -1;
    return;
  }

  /* Attach user's yy0 to F2C_IDA_vec */
  N_VSetArrayPointer(yy0, F2C_IDA_vec);

  /* Create F2C_IDA_ypvec and attach user's yp0 to it */
  F2C_IDA_ypvec = N_VCloneEmpty(F2C_IDA_vec);
  N_VSetArrayPointer(yp0, F2C_IDA_ypvec);

  /* Treat absolute tolerances */
  itol = -1;
  switch (*iatol) {
  case 1:
    itol = IDA_SS;
    atolptr = (void *) atol;
    break;
  case 2:
    itol = IDA_SV;
    Vatol= N_VCloneEmpty(F2C_IDA_vec);
    N_VSetArrayPointer(atol, Vatol);
    atolptr = (void *) Vatol;
    break;
  case 3:
    itol = IDA_WF;
    break;
  }

  /* Call IDAMalloc */
  *ier = 0;
  *ier = IDAMalloc(IDA_idamem, (IDAResFn) FIDAresfn, *t0, F2C_IDA_vec, F2C_IDA_ypvec,
		   itol, *rtol, atolptr);

  /* destroy Vatol if allocated */
  if (itol == IDA_SV) N_VDestroy(Vatol);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  /* On failure, clean-up and exit */
  if (*ier != IDA_SUCCESS) {
    N_VDestroy(F2C_IDA_ypvec);
    *ier = -1;
    return;
  }

  /* Grab optional output arrays and store them in global variables */
  IDA_iout = iout;
  IDA_rout = rout;

  /* Store the unit roundoff in rout for user access */
  IDA_rout[5] = UNIT_ROUNDOFF;

  return;
}

/*************************************************/

void FIDA_REINIT(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 int *ier)
{
  int itol;
  N_Vector Vatol;
  void *atolptr;

  /* Initialize all pointers to NULL */
  atolptr = NULL;
  Vatol = NULL;

  /* Attach user's yy0 to F2C_IDA_vec */
  N_VSetArrayPointer(yy0, F2C_IDA_vec);

  /* Attach user's yp0 to F2C_IDA_ypvec */
  N_VSetArrayPointer(yp0, F2C_IDA_ypvec);

  /* Treat absolute tolerances */
  itol = -1;
  switch (*iatol) {
  case 1:
    itol = IDA_SS;
    atolptr = atol;
    break;
  case 2:
    Vatol = N_VCloneEmpty(F2C_IDA_vec);
    N_VSetArrayPointer(atol, Vatol);
    itol = IDA_SV;
    atolptr = (void *) Vatol;
    break;
  case 3:
    itol = IDA_WF;
    break;
  }

  /* Call IDAReInit */
  *ier = 0;
  *ier = IDAReInit(IDA_idamem, (IDAResFn) FIDAresfn, *t0, F2C_IDA_vec, F2C_IDA_ypvec,
		   itol, *rtol, atolptr);

  /* destroy Vatol if allocated */
  if (itol == IDA_SV) N_VDestroy(Vatol);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  /* On failure, clean-up and exit */
  if (*ier != IDA_SUCCESS) {
    N_VDestroy(F2C_IDA_ypvec);
    *ier = -1;
    return;
  }

  return;
}

/*************************************************/

void FIDA_SETIIN(char key_name[], long int *ival, int *ier, int key_len)
{
  if (!strncmp(key_name,"MAX_ORD", (size_t)key_len)) 
    *ier = IDASetMaxOrd(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NSTEPS", (size_t)key_len)) 
    *ier = IDASetMaxNumSteps(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_ERRFAIL", (size_t)key_len)) 
    *ier = IDASetMaxErrTestFails(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NITERS", (size_t)key_len)) 
    *ier = IDASetMaxNonlinIters(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_CONVFAIL", (size_t)key_len)) 
    *ier = IDASetMaxConvFails(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"SUPPRESS_ALG", (size_t)key_len)) 
    *ier = IDASetSuppressAlg(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NSTEPS_IC", (size_t)key_len)) 
    *ier = IDASetMaxNumStepsIC(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NITERS_IC", (size_t)key_len)) 
    *ier = IDASetMaxNumItersIC(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NJE_IC", (size_t)key_len)) 
    *ier = IDASetMaxNumJacsIC(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"LS_OFF_IC", (size_t)key_len)) 
    *ier = IDASetLineSearchOffIC(IDA_idamem, (int) *ival);
  else {
    *ier = -99;
    printf("FIDASETIIN: Unrecognized key.\n\n");
  }

}

/***************************************************************************/

void FIDA_SETRIN(char key_name[], realtype *rval, int *ier, int key_len)
{

  if (!strncmp(key_name,"INIT_STEP", (size_t)key_len)) 
    *ier = IDASetInitStep(IDA_idamem, *rval);
  else if (!strncmp(key_name,"MAX_STEP", (size_t)key_len)) 
    *ier = IDASetMaxStep(IDA_idamem, *rval);
  else if (!strncmp(key_name,"STOP_TIME", (size_t)key_len)) 
    *ier = IDASetStopTime(IDA_idamem, *rval);
  else if (!strncmp(key_name,"NLCONV_COEF", (size_t)key_len)) 
    *ier = IDASetNonlinConvCoef(IDA_idamem, *rval);
  else if (!strncmp(key_name,"NLCONV_COEF_IC", (size_t)key_len)) 
    *ier = IDASetNonlinConvCoefIC(IDA_idamem, *rval);
  else if (!strncmp(key_name,"STEP_TOL_IC", (size_t)key_len)) 
    *ier = IDASetStepToleranceIC(IDA_idamem, *rval);
  else {
    *ier = -99;
    printf("FIDASETRIN: Unrecognized key.\n\n");
  }

}

/*************************************************/

void FIDA_SETVIN(char key_name[], realtype *vval, int *ier, int key_len)
{
  N_Vector Vec;

  if (!strncmp(key_name,"ID_VEC", (size_t)key_len)) {
    Vec = N_VCloneEmpty(F2C_IDA_vec);
    N_VSetArrayPointer(vval, Vec);
    IDASetId(IDA_idamem, Vec);
    N_VDestroy(Vec);
  } else if (!strncmp(key_name,"CONSTR_VEC", (size_t)key_len)) {
    Vec = N_VCloneEmpty(F2C_IDA_vec);
    N_VSetArrayPointer(vval, Vec);
    IDASetConstraints(IDA_idamem, Vec);
    N_VDestroy(Vec);
  } else {
    *ier = -99;
    printf("FIDASETVIN: Unrecognized key.\n\n");
  }

}

/*************************************************/

void FIDA_TOLREINIT(int *iatol, realtype *rtol, realtype *atol, int *ier)
{
  int itol;
  N_Vector Vatol;
  void *atolptr;

  atolptr = NULL;

  itol = -1;
  if (*iatol == 1) {
    itol = IDA_SS;
    atolptr = atol;
  }
  else {
    itol = IDA_SV;
    Vatol = N_VCloneEmpty(F2C_IDA_vec);
    N_VSetArrayPointer(atol, Vatol);
    atolptr = (void *) Vatol; 
  }

  *ier = 0;
  *ier = IDASetTolerances(IDA_idamem, itol, *rtol, atolptr);

  if (itol == IDA_SV) N_VDestroy(Vatol);

  return;
}

/*************************************************/

void FIDA_CALCIC(realtype *t0, realtype *yy0, realtype *yp0,
		 int *icopt, realtype *tout1, int *ier)
{
  /* Attach user data to vectors */
  N_VSetArrayPointer(yy0, F2C_IDA_vec);
  N_VSetArrayPointer(yp0, F2C_IDA_ypvec);

  *ier = 0;
  *ier = IDACalcIC(IDA_idamem, *t0, F2C_IDA_vec, F2C_IDA_ypvec, *icopt, *tout1);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  return;
}

/*************************************************/

void FIDA_SPTFQMR(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier)
{

  *ier = 0;

  *ier = IDASptfqmr(IDA_idamem, *maxl);
  if (*ier != IDASPTFQMR_SUCCESS) return;

  if (*eplifac != ZERO) {
    *ier = IDASptfqmrSetEpsLin(IDA_idamem, *eplifac);
    if (*ier != IDASPTFQMR_SUCCESS) return;
  }

  if (*dqincfac != ZERO) {
    *ier = IDASptfqmrSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPTFQMR_SUCCESS) return;
  }

  IDA_ls = IDA_LS_SPTFQMR;

  return;
}

/*************************************************/

void FIDA_SPBCG(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier)
{

  *ier = 0;

  *ier = IDASpbcg(IDA_idamem, *maxl);
  if (*ier != IDASPBCG_SUCCESS) return;

  if (*eplifac != ZERO) {
    *ier = IDASpbcgSetEpsLin(IDA_idamem, *eplifac);
    if (*ier != IDASPBCG_SUCCESS) return;
  }

  if (*dqincfac != ZERO) {
    *ier = IDASpbcgSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPBCG_SUCCESS) return;
  }

  IDA_ls = IDA_LS_SPBCG;

  return;
}

/*************************************************/

void FIDA_SPGMR(int *maxl, int *gstype, int *maxrs,
		realtype *eplifac, realtype *dqincfac, int *ier)
{

  *ier = 0;

  *ier = IDASpgmr(IDA_idamem, *maxl);
  if (*ier != IDASPGMR_SUCCESS) return;

  if (*gstype != 0) {
    *ier = IDASpgmrSetGSType(IDA_idamem, *gstype);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  if (*maxrs != 0) {
    *ier = IDASpgmrSetMaxRestarts(IDA_idamem, *maxrs);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  if (*eplifac != ZERO) {
    *ier = IDASpgmrSetEpsLin(IDA_idamem, *eplifac);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  if (*dqincfac != ZERO) {
    *ier = IDASpgmrSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  IDA_ls = IDA_LS_SPGMR;

  return;
}

/*************************************************/

void FIDA_DENSE(long int *neq, int *ier)
{

  *ier = 0;

  *ier = IDADense(IDA_idamem, *neq);

  IDA_ls = IDA_LS_DENSE;

  return;
}

/*************************************************/

void FIDA_BAND(long int *neq, long int *mupper, long int *mlower, int *ier)
{

  *ier = 0;

  *ier = IDABand(IDA_idamem, *neq, *mupper, *mlower);

  IDA_ls = IDA_LS_BAND;

  return;
}

/*************************************************/

void FIDA_SPTFQMRREINIT(realtype *eplifac, realtype *dqincfac, int *ier)
{

  *ier = 0;

  if (*eplifac != ZERO) {
    *ier = IDASptfqmrSetEpsLin(IDA_idamem, *eplifac);
    if (*ier != IDASPTFQMR_SUCCESS) return;
  }

  if (*dqincfac != ZERO) {
    *ier = IDASptfqmrSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPTFQMR_SUCCESS) return;
  }

  IDA_ls = IDA_LS_SPTFQMR;

  return;
}

/*************************************************/

void FIDA_SPBCGREINIT(realtype *eplifac, realtype *dqincfac, int *ier)
{

  *ier = 0;

  if (*eplifac != ZERO) {
    *ier = IDASpbcgSetEpsLin(IDA_idamem, *eplifac);
    if (*ier != IDASPBCG_SUCCESS) return;
  }

  if (*dqincfac != ZERO) {
    *ier = IDASpbcgSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPBCG_SUCCESS) return;
  }

  IDA_ls = IDA_LS_SPBCG;

  return;
}

/*************************************************/

void FIDA_SPGMRREINIT(int *gstype, int *maxrs, realtype *eplifac,
		      realtype *dqincfac, int *ier)
{

  *ier = 0;

  if (*gstype != 0) {
    *ier = IDASpgmrSetGSType(IDA_idamem, *gstype);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  if (*maxrs != 0) {
    *ier = IDASpgmrSetMaxRestarts(IDA_idamem, *maxrs);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  if (*eplifac != ZERO) {
    *ier = IDASpgmrSetEpsLin(IDA_idamem, *eplifac);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  if (*dqincfac != ZERO) {
    *ier = IDASpgmrSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  IDA_ls = IDA_LS_SPGMR;

  return;
}

/*************************************************/

void FIDA_SOLVE(realtype *tout, realtype *tret, realtype *yret,
		realtype *ypret, int *itask, int *ier)
{

  *ier = 0;

  /* Attach user data to vectors */
  N_VSetArrayPointer(yret, F2C_IDA_vec);
  N_VSetArrayPointer(ypret, F2C_IDA_ypvec);

  *ier = IDASolve(IDA_idamem, *tout, tret, F2C_IDA_vec, F2C_IDA_ypvec, *itask);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  /* Set optional outputs */

  IDAGetWorkSpace(IDA_idamem,
                  &IDA_iout[0],                 /* LENRW */
                  &IDA_iout[1]);                /* LENIW */

  IDAGetIntegratorStats(IDA_idamem,
                        &IDA_iout[2],           /* NST */
                        &IDA_iout[3],           /* NRE */
                        &IDA_iout[7],           /* NSETUPS */
                        &IDA_iout[4],           /* NETF */
                        (int *) &IDA_iout[8],   /* KLAST */
                        (int *) &IDA_iout[9],   /* KCUR */
                        &IDA_rout[1],           /* HLAST */
                        &IDA_rout[2],           /* HCUR */
                        &IDA_rout[3]);          /* TCUR */
  IDAGetNonlinSolvStats(IDA_idamem,
                        &IDA_iout[6],           /* NNI */
                        &IDA_iout[5]);          /* NCFN */
  IDAGetNumBacktrackOps(IDA_idamem, 
                        &IDA_iout[10]);         /* NBCKTRK */
  IDAGetActualInitStep(IDA_idamem, 
                       &IDA_rout[0]);           /* HINUSED */
  IDAGetTolScaleFactor(IDA_idamem,    
                       &IDA_rout[4]);           /* TOLSFAC */
  
  /* Root finding is on */
  if (IDA_nrtfn != 0)
    IDAGetNumGEvals(IDA_idamem, &IDA_iout[11]); /* NGE */
  
  switch(IDA_ls) {
  case IDA_LS_DENSE:
    IDADenseGetWorkSpace(IDA_idamem, &IDA_iout[12], &IDA_iout[13]);   /* LRW and LIW */
    IDADenseGetLastFlag(IDA_idamem, (int *) &IDA_iout[14]);           /* LSTF */
    IDADenseGetNumResEvals(IDA_idamem, &IDA_iout[15]);                /* NRE */
    IDADenseGetNumJacEvals(IDA_idamem, &IDA_iout[16]);                /* NJE */
    break;
  case IDA_LS_BAND:
    IDABandGetWorkSpace(IDA_idamem, &IDA_iout[12], &IDA_iout[13]);    /* LRW and LIW */
    IDABandGetLastFlag(IDA_idamem, (int *) &IDA_iout[14]);            /* LSTF */
    IDABandGetNumResEvals(IDA_idamem, &IDA_iout[15]);                 /* NRE */
    IDABandGetNumJacEvals(IDA_idamem, &IDA_iout[16]);                 /* NJE */
    break;
  case IDA_LS_SPGMR:
    IDASpgmrGetWorkSpace(IDA_idamem, &IDA_iout[12], &IDA_iout[13]);   /* LRW and LIW */
    IDASpgmrGetLastFlag(IDA_idamem, (int *) &IDA_iout[14]);           /* LSTF */
    IDASpgmrGetNumResEvals(IDA_idamem, &IDA_iout[15]);                /* NRE */
    IDASpgmrGetNumJtimesEvals(IDA_idamem, &IDA_iout[16]);             /* NJE */
    IDASpgmrGetNumPrecEvals(IDA_idamem, &IDA_iout[17]);               /* NPE */
    IDASpgmrGetNumPrecSolves(IDA_idamem, &IDA_iout[18]);              /* NPS */
    IDASpgmrGetNumLinIters(IDA_idamem, &IDA_iout[19]);                /* NLI */
    IDASpgmrGetNumConvFails(IDA_idamem, &IDA_iout[20]);               /* NCFL */
    break;
  case IDA_LS_SPBCG:
    IDASpbcgGetWorkSpace(IDA_idamem, &IDA_iout[12], &IDA_iout[13]);   /* LRW and LIW */
    IDASpbcgGetLastFlag(IDA_idamem, (int *) &IDA_iout[14]);           /* LSTF */
    IDASpbcgGetNumResEvals(IDA_idamem, &IDA_iout[15]);                /* NRE */
    IDASpbcgGetNumJtimesEvals(IDA_idamem, &IDA_iout[16]);             /* NJE */
    IDASpbcgGetNumPrecEvals(IDA_idamem, &IDA_iout[17]);               /* NPE */
    IDASpbcgGetNumPrecSolves(IDA_idamem, &IDA_iout[18]);              /* NPS */
    IDASpbcgGetNumLinIters(IDA_idamem, &IDA_iout[19]);                /* NLI */
    IDASpbcgGetNumConvFails(IDA_idamem, &IDA_iout[20]);               /* NCFL */
    break;
  case IDA_LS_SPTFQMR:
    IDASptfqmrGetWorkSpace(IDA_idamem, &IDA_iout[12], &IDA_iout[13]); /* LRW and LIW */
    IDASptfqmrGetLastFlag(IDA_idamem, (int *) &IDA_iout[14]);         /* LSTF */
    IDASptfqmrGetNumResEvals(IDA_idamem, &IDA_iout[15]);              /* NRE */
    IDASptfqmrGetNumJtimesEvals(IDA_idamem, &IDA_iout[16]);           /* NJE */
    IDASptfqmrGetNumPrecEvals(IDA_idamem, &IDA_iout[17]);             /* NPE */
    IDASptfqmrGetNumPrecSolves(IDA_idamem, &IDA_iout[18]);            /* NPS */
    IDASptfqmrGetNumLinIters(IDA_idamem, &IDA_iout[19]);              /* NLI */
    IDASptfqmrGetNumConvFails(IDA_idamem, &IDA_iout[20]);             /* NCFL */
    break;
  }


  return;
}

/*************************************************/

void FIDA_GETSOL(realtype *t, realtype *yret, realtype *ypret, int *ier)
{
  /* Attach user data to vectors */
  N_VSetArrayPointer(yret, F2C_IDA_vec);
  N_VSetArrayPointer(ypret, F2C_IDA_ypvec);

  *ier = 0;
  *ier = IDAGetSolution(IDA_idamem, *t, F2C_IDA_vec, F2C_IDA_ypvec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  return;
}

/*************************************************/

void FIDA_GETERRWEIGHTS(realtype *y, realtype *eweight, int *ier)
{
  /* Attach user data to vectors */
  N_VSetArrayPointer(y, F2C_IDA_vec);

  *ier = 0;
  *ier = IDAGetErrWeights(IDA_idamem, F2C_IDA_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);

  return;
}

/*************************************************/

void FIDA_GETESTLOCALERR(realtype *ele, int *ier)
{
  /* Attach user data to vector */
  N_VSetArrayPointer(ele, F2C_IDA_vec);

  *ier = 0;
  *ier = IDAGetEstLocalErrors(IDA_idamem, F2C_IDA_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);

  return;
}

/*************************************************/

void FIDA_FREE(void)
{
  if (IDA_idamem != NULL) IDAFree(IDA_idamem);

  /* Detach user data from F2C_IDA_vec */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);

  /* Free F2C_IDA_ypvec */
  if (F2C_IDA_ypvec != NULL) {
    N_VSetArrayPointer(NULL, F2C_IDA_ypvec);
    N_VDestroy(F2C_IDA_ypvec);
  }

  /* Free F2C_IDA_ewtvec */
  if (F2C_IDA_ewtvec != NULL) N_VDestroy(F2C_IDA_ewtvec);

  return;
}

/*************************************************/

int FIDAresfn(realtype t, N_Vector yy, N_Vector yp,
	      N_Vector rr, void *res_data)
{
  realtype *yy_data, *yp_data, *rr_data;
  int ier;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);

  /* Call user-supplied routine */
  FIDA_RESFUN(&t, yy_data, yp_data, rr_data, &ier);

  return(ier);
}
