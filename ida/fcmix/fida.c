/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-05-11 23:10:54 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the IDA package. See fida.h for usage.
 * NOTE: Some routines are necessarily stored elsewhere to avoid
 * linking problems.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "idaband.h"        /* prototypes for IDABAND interface routines      */
#include "idadense.h"       /* prototypes for IDADENSE interface routines     */
#include "ida.h"            /* IDA constants and prototypes                   */
#include "idaspbcg.h"       /* prototypes for IDASPBCG interface routines     */
#include "idaspgmr.h"       /* prototypes for IDASPGMR interface routines     */
#include "fida.h"           /* actual function names, prototypes and global
			       variables                                      */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */

/*************************************************/

/* Definitions for global variables shared amongst various routines */

N_Vector F2C_IDA_ypvec, F2C_IDA_atolvec, F2C_IDA_ewtvec;

void *IDA_idamem;
booleantype IDA_optin;
long int *IDA_iopt;
realtype *IDA_ropt;
int IDA_ls;

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
		 realtype *id, realtype *constr,
                 int *optin, long int *iopt, realtype *ropt,
                 int *ier)
{
  int itol;
  void *atolptr;

  /* Initialize itol to avoid compiler warning message */
  itol = -1;

  /* Check for required vector operations */
  if ((F2C_IDA_vec->ops->nvgetarraypointer == NULL) ||
      (F2C_IDA_vec->ops->nvsetarraypointer == NULL)) {
    *ier = -1;
    printf("A required vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize all pointers to NULL */
  IDA_idamem = NULL;
  atolptr = NULL;
  F2C_IDA_ypvec = F2C_IDA_atolvec = F2C_IDA_ewtvec = NULL;

  /* Attach user data to F2C_IDA_vec */
  N_VSetArrayPointer(yy0, F2C_IDA_vec);

  switch (*iatol) {
  case 1:
    itol = IDA_SS;
    atolptr = (void *) atol;
    break;
  case 2:
    itol = IDA_SV;
    F2C_IDA_atolvec = N_VCloneEmpty(F2C_IDA_vec);
    N_VSetArrayPointer(atol, F2C_IDA_atolvec);
    atolptr = (void *) F2C_IDA_atolvec;
    break;
  case 3:
    itol = IDA_WF;
    break;
  }

  *ier = 0;

  IDA_idamem = IDACreate();

  /* If IDACreate() fails, then deallocate memory and exit */
  if (IDA_idamem == NULL) {
    N_VSetArrayPointer(NULL, F2C_IDA_vec);
    if (F2C_IDA_atolvec != NULL) {
      N_VSetArrayPointer(NULL, F2C_IDA_atolvec);
      N_VDestroy(F2C_IDA_atolvec);
      F2C_IDA_atolvec = NULL;
    }
    *ier = -1;
    return;
  }

  /* Create F2C_IDA_ypvec */
  F2C_IDA_ypvec = N_VCloneEmpty(F2C_IDA_vec);
  /* Attach user data to F2C_IDA_ypvec */
  N_VSetArrayPointer(yp0, F2C_IDA_ypvec);

  *ier = IDAMalloc(IDA_idamem, (IDAResFn) FIDAresfn, *t0, F2C_IDA_vec, F2C_IDA_ypvec,
		   itol, *rtol, atolptr);

  /* If IDAMalloc() fails, then deallocate memory and exit */
  if (*ier != IDA_SUCCESS) {
    N_VSetArrayPointer(NULL, F2C_IDA_vec);
    N_VSetArrayPointer(NULL, F2C_IDA_ypvec);
    N_VDestroy(F2C_IDA_ypvec);
    if (F2C_IDA_atolvec != NULL) {
      N_VSetArrayPointer(NULL, F2C_IDA_atolvec);
      N_VDestroy(F2C_IDA_atolvec);
      F2C_IDA_atolvec = NULL;
    }
    *ier = -1;
    return;
  }

  /* Set optional inputs */
  if (*optin == 1) {
    IDA_optin = TRUE;
    if (iopt[0] > 0) IDASetMaxOrd(IDA_idamem, (int) iopt[0]);
    if (iopt[1] > 0) IDASetMaxNumSteps(IDA_idamem, iopt[1]);
    if (iopt[2] > 0) IDASetMaxErrTestFails(IDA_idamem, (int) iopt[2]);
    if (iopt[3] > 0) IDASetMaxNonlinIters(IDA_idamem, (int) iopt[3]);
    if (iopt[4] > 0) IDASetMaxConvFails(IDA_idamem, (int) iopt[4]);
    if (iopt[5] > 0) IDASetSuppressAlg(IDA_idamem, TRUE);
    if (iopt[6] > 0) {
      /* NOTE: Reusing F2C_IDA_atolvec to save memory */
      if (F2C_IDA_atolvec == NULL) F2C_IDA_atolvec = N_VCloneEmpty(F2C_IDA_vec);
      N_VSetArrayPointer(id, F2C_IDA_atolvec);
      IDASetId(IDA_idamem, F2C_IDA_atolvec);
    }
    if (iopt[7] > 0) {
      /* NOTE: Reusing F2C_ypvec to save memory */
      if (F2C_IDA_atolvec == NULL) F2C_IDA_atolvec = N_VCloneEmpty(F2C_IDA_vec);
      N_VSetArrayPointer(constr, F2C_IDA_atolvec);
      IDASetConstraints(IDA_idamem, F2C_IDA_atolvec);
    }
    if (iopt[8]  > 0) IDASetMaxNumStepsIC(IDA_idamem, (int) iopt[8]);
    if (iopt[9]  > 0) IDASetMaxNumJacsIC(IDA_idamem, (int) iopt[9]);
    if (iopt[10] > 0) IDASetMaxNumItersIC(IDA_idamem, (int) iopt[10]);
    if (iopt[11] > 0) IDASetLineSearchOffIC(IDA_idamem, TRUE);

    if (ropt[0] != ZERO) IDASetInitStep(IDA_idamem, ropt[0]);
    if (ropt[1] > ZERO)  IDASetMaxStep(IDA_idamem, ropt[1]);
    if (ropt[2] != ZERO) IDASetStopTime(IDA_idamem, ropt[2]);
    if (ropt[3] > ZERO)  IDASetNonlinConvCoef(IDA_idamem, ropt[3]);
    if (ropt[4] > ZERO)  IDASetNonlinConvCoefIC(IDA_idamem, ropt[4]);
    if (ropt[5] > ZERO)  IDASetStepToleranceIC(IDA_idamem, ropt[5]);
  }
  else IDA_optin = FALSE;

  /* NOTE: Private vectors are NOT deallocated until the call to
     the FIDFREE routine */

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);
  N_VSetArrayPointer(NULL, F2C_IDA_atolvec);

  atolptr = NULL;

  /* Store the unit roundoff in ropt[] for user access */
  ropt[14] = UNIT_ROUNDOFF;

  IDA_iopt = iopt;
  IDA_ropt = ropt;

  return;
}

/*************************************************/

void FIDA_REINIT(realtype *t0, realtype *yy0, realtype *yp0,
		 int *iatol, realtype *rtol, realtype *atol,
		 realtype *id, realtype *constr,
		 int *optin, long int *iopt, realtype *ropt,
		 int *ier)
{
  int itol;
  void *atolptr;

  /* Initialize itol to avoid compiler warning message */
  itol = -1;

  /* Initialize all pointers to NULL */
  atolptr = NULL;

  /* Attach user data to F2C_IDA_vec */
  N_VSetArrayPointer(yy0, F2C_IDA_vec);

  switch (*iatol) {
  case 1:
    itol = IDA_SS;
    atolptr = atol;
    break;
  case 2:
    F2C_IDA_atolvec = N_VCloneEmpty(F2C_IDA_vec);
    N_VSetArrayPointer(atol, F2C_IDA_atolvec);
    itol = IDA_SV;
    atolptr = (void *) F2C_IDA_atolvec;
    break;
  case 3:
    itol = IDA_WF;
    break;
  }

  /* Attach user data to F2C_IDA_ypvec */
  N_VSetArrayPointer(yp0, F2C_IDA_ypvec);

  *ier = 0;

  *ier = IDAReInit(IDA_idamem, (IDAResFn) FIDAresfn, *t0, F2C_IDA_vec, F2C_IDA_ypvec,
		   itol, *rtol, atolptr);

  /* If IDAReInit() fails, then deallocate memory and exit */
  if (*ier != IDA_SUCCESS) {
    N_VSetArrayPointer(NULL, F2C_IDA_vec);
    N_VSetArrayPointer(NULL, F2C_IDA_ypvec);
    N_VDestroy(F2C_IDA_ypvec);
    if (F2C_IDA_atolvec != NULL) {
      N_VSetArrayPointer(NULL, F2C_IDA_atolvec);
      N_VDestroy(F2C_IDA_atolvec);
      F2C_IDA_atolvec = NULL;
    }
    *ier = -1;
    return;
  }

  /* Set optional inputs */
  if (*optin == 1) {
    IDA_optin = TRUE;
    if (iopt[0] > 0) IDASetMaxOrd(IDA_idamem, (int) iopt[0]);
    if (iopt[1] > 0) IDASetMaxNumSteps(IDA_idamem, iopt[1]);
    if (iopt[2] > 0) IDASetMaxErrTestFails(IDA_idamem, (int) iopt[2]);
    if (iopt[3] > 0) IDASetMaxNonlinIters(IDA_idamem, (int) iopt[3]);
    if (iopt[4] > 0) IDASetMaxConvFails(IDA_idamem, (int) iopt[4]);
    if (iopt[5] > 0) IDASetSuppressAlg(IDA_idamem, TRUE);
    if (iopt[6] > 0) {
      /* NOTE: Reusing F2C_IDA_atolvec to save memory */
      if (F2C_IDA_atolvec == NULL) N_VCloneEmpty(F2C_IDA_vec);
      N_VSetArrayPointer(id, F2C_IDA_atolvec);
      IDASetId(IDA_idamem, F2C_IDA_atolvec);
    }
    if (iopt[7] > 0) {
      /* NOTE: Reusing F2C_IDA_atolvec to save memory */
      if (F2C_IDA_atolvec == NULL) N_VCloneEmpty(F2C_IDA_vec);
      N_VSetArrayPointer(constr, F2C_IDA_atolvec);
      IDASetConstraints(IDA_idamem, F2C_IDA_atolvec);
    }
    if (iopt[8]  > 0) IDASetMaxNumStepsIC(IDA_idamem, (int) iopt[8]);
    if (iopt[9]  > 0) IDASetMaxNumJacsIC(IDA_idamem, (int) iopt[9]);
    if (iopt[10]  > 0) IDASetMaxNumItersIC(IDA_idamem, (int) iopt[10]);
    if (iopt[11] > 0) IDASetLineSearchOffIC(IDA_idamem, TRUE);

    if (ropt[0] != ZERO) IDASetInitStep(IDA_idamem, ropt[0]);
    if (ropt[1] > ZERO)  IDASetMaxStep(IDA_idamem, ropt[1]);
    if (ropt[2] != ZERO) IDASetStopTime(IDA_idamem, ropt[2]);
    if (ropt[3] > ZERO)  IDASetNonlinConvCoef(IDA_idamem, ropt[3]);
    if (ropt[4] > ZERO)  IDASetNonlinConvCoefIC(IDA_idamem, ropt[4]);
    if (ropt[5] > ZERO)  IDASetStepToleranceIC(IDA_idamem, ropt[5]);
  }
  else IDA_optin = FALSE;

  /* NOTE: Private vectors are NOT deallocated until the call to
     the FIDFREE routine */

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);
  N_VSetArrayPointer(NULL, F2C_IDA_atolvec);

  atolptr = NULL;

  IDA_iopt = iopt;
  IDA_ropt = ropt;

  return;
}

/*************************************************/

void FIDA_TOLREINIT(int *iatol, realtype *rtol, realtype *atol, int *ier)
{
  int itol;
  void *atolptr;

  atolptr = NULL;

  if (*iatol == 1) {
    itol = IDA_SS;
    atolptr = atol;
  }
  else {
    N_VSetArrayPointer(atol, F2C_IDA_atolvec);
    itol = IDA_SV;
    atolptr = (void *) F2C_IDA_atolvec;
  }

  *ier = 0;

  *ier = IDASetTolerances(IDA_idamem, itol, *rtol, atolptr);

  /* Reset pointers */
  if (itol == IDA_SV) N_VSetArrayPointer(NULL, F2C_IDA_atolvec);
  else if (itol == IDA_SS) atolptr = NULL;

  if (*ier != IDA_SUCCESS) *ier = -1;

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

  if (*ier != IDA_SUCCESS) {
    N_VSetArrayPointer(NULL, F2C_IDA_vec);
    N_VSetArrayPointer(NULL, F2C_IDA_ypvec);
    *ier = -1;
  }

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  return;
}

/*************************************************/

void FIDA_SPBCG(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier)
{
  *ier = IDASpbcg(IDA_idamem, *maxl);
  if (*ier != IDASPBCG_SUCCESS) return;

  *ier = IDASpbcgSetEpsLin(IDA_idamem, *eplifac);
  if (*ier != IDASPBCG_SUCCESS) return;

  if (*dqincfac != ZERO) {
    *ier = IDASpbcgSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPBCG_SUCCESS) return;
  }

  IDA_ls = IDA_SPBCG;

  return;
}

/*************************************************/

void FIDA_SPGMR(int *maxl, int *gstype, int *maxrs,
		realtype *eplifac, realtype *dqincfac, int *ier)
{
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

  *ier = IDASpgmrSetEpsLin(IDA_idamem, *eplifac);
  if (*ier != IDASPGMR_SUCCESS) return;

  if (*dqincfac != ZERO) {
    *ier = IDASpgmrSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  IDA_ls = IDA_SPGMR;

  return;
}

/*************************************************/

void FIDA_DENSE(long int *neq, int *ier)
{
  *ier = 0;

  *ier = IDADense(IDA_idamem, *neq);

  IDA_ls = IDA_DENSE;

  return;
}

/*************************************************/

void FIDA_BAND(long int *neq, long int *mupper, long int *mlower, int *ier)
{
  *ier = IDABand(IDA_idamem, *neq, *mupper, *mlower);

  IDA_ls = IDA_BAND;

  return;
}

/*************************************************/

void FIDA_SPBCGREINIT(realtype *eplifac, realtype *dqincfac, int *ier)
{
  *ier = IDASpbcgSetEpsLin(IDA_idamem, *eplifac);
  if (*ier != IDASPBCG_SUCCESS) return;

  if (*dqincfac != ZERO) {
    *ier = IDASpbcgSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPBCG_SUCCESS) return;
  }

  IDA_ls = IDA_SPBCG;

  return;
}

/*************************************************/

void FIDA_SPGMRREINIT(int *gstype, realtype *eplifac,
		      realtype *dqincfac, int *ier)
{
  if (*gstype != 0) {
    *ier = IDASpgmrSetGSType(IDA_idamem, *gstype);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  *ier = IDASpgmrSetEpsLin(IDA_idamem, *eplifac);
  if (*ier != IDASPGMR_SUCCESS) return;

  if (*dqincfac != ZERO) {
    *ier = IDASpgmrSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  IDA_ls = IDA_SPGMR;

  return;
}

/*************************************************/

void FIDA_SOLVE(realtype *tout, realtype *tret, realtype *yret,
		realtype *ypret, int *itask, int *ier)
{
  int qlast, qcur, last_flag;

  *ier = 0;

  /* Attach user data to vectors */
  N_VSetArrayPointer(yret, F2C_IDA_vec);
  N_VSetArrayPointer(ypret, F2C_IDA_ypvec);

  *ier = IDASolve(IDA_idamem, *tout, tret, F2C_IDA_vec, F2C_IDA_ypvec, *itask);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  /* Set optional outputs */
  if ((IDA_iopt != NULL) && (IDA_ropt != NULL)) {

    /* NOTE: ropt[14] = UNIT_ROUNDOFF */
    IDAGetIntegratorStats(IDA_idamem,
			  &IDA_iopt[14],
			  &IDA_iopt[15],
			  &IDA_iopt[16],
			  &IDA_iopt[17],
			  &qlast,
			  &qcur,
			  &IDA_ropt[15],
			  &IDA_ropt[16],
			  &IDA_ropt[17]);

    IDA_iopt[18] = (long int) qlast;
    IDA_iopt[19] = (long int) qcur;

    IDAGetNonlinSolvStats(IDA_idamem,
			  &IDA_iopt[20],
			  &IDA_iopt[21]);

    IDAGetNumBacktrackOps(IDA_idamem, &IDA_iopt[22]);

    IDAGetActualInitStep(IDA_idamem, &IDA_ropt[18]);
    IDAGetTolScaleFactor(IDA_idamem, &IDA_ropt[19]);

    IDAGetWorkSpace(IDA_idamem,
		    &IDA_iopt[23],
		    &IDA_iopt[24]);

    switch(IDA_ls) {
    case IDA_SPGMR:
      IDASpgmrGetWorkSpace(IDA_idamem, &IDA_iopt[25], &IDA_iopt[26]);
      IDASpgmrGetNumPrecEvals(IDA_idamem, &IDA_iopt[27]);
      IDASpgmrGetNumPrecSolves(IDA_idamem, &IDA_iopt[28]);
      IDASpgmrGetNumLinIters(IDA_idamem, &IDA_iopt[29]);
      IDASpgmrGetNumConvFails(IDA_idamem, &IDA_iopt[30]);
      IDASpgmrGetNumJtimesEvals(IDA_idamem, &IDA_iopt[31]);
      IDASpgmrGetNumResEvals(IDA_idamem, &IDA_iopt[32]);
      IDASpgmrGetLastFlag(IDA_idamem, &last_flag);
      IDA_iopt[33] = (long int) last_flag;
      break;
    case IDA_SPBCG:
      IDASpbcgGetWorkSpace(IDA_idamem, &IDA_iopt[25], &IDA_iopt[26]);
      IDASpbcgGetNumPrecEvals(IDA_idamem, &IDA_iopt[27]);
      IDASpbcgGetNumPrecSolves(IDA_idamem, &IDA_iopt[28]);
      IDASpbcgGetNumLinIters(IDA_idamem, &IDA_iopt[29]);
      IDASpbcgGetNumConvFails(IDA_idamem, &IDA_iopt[30]);
      IDASpbcgGetNumJtimesEvals(IDA_idamem, &IDA_iopt[31]);
      IDASpbcgGetNumResEvals(IDA_idamem, &IDA_iopt[32]);
      IDASpbcgGetLastFlag(IDA_idamem, &last_flag);
      IDA_iopt[33] = (long int) last_flag;
      break;
    case IDA_DENSE:
      IDADenseGetWorkSpace(IDA_idamem, &IDA_iopt[25], &IDA_iopt[26]);
      IDADenseGetNumJacEvals(IDA_idamem, &IDA_iopt[27]);
      IDADenseGetNumResEvals(IDA_idamem, &IDA_iopt[28]);
      IDADenseGetLastFlag(IDA_idamem, &last_flag);
      IDA_iopt[29] = (long int) last_flag;
      break;
    case IDA_BAND:
      IDABandGetWorkSpace(IDA_idamem, &IDA_iopt[25], &IDA_iopt[26]);
      IDABandGetNumJacEvals(IDA_idamem, &IDA_iopt[27]);
      IDABandGetNumResEvals(IDA_idamem, &IDA_iopt[28]);
      IDABandGetLastFlag(IDA_idamem, &last_flag);
      IDA_iopt[29] = (long int) last_flag;
      break;
    }
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
  N_VSetArrayPointer(eweight, F2C_IDA_ypvec);

  *ier = 0;

  *ier = IDAGetErrWeights(IDA_idamem, F2C_IDA_vec, F2C_IDA_ypvec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

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
