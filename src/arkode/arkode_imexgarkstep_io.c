/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * Based on arkode_arkstep_io.c written by Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for the optional input and
 * output functions for the ARKode IMEXGARKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here,
 * with slightly different names.  The code transition will be
 * minimal, but the documentation changes will be significant.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_imexgarkstep_impl.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_types.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  IMEXGARKStep Optional input functions (wrappers for generic ARKode
  utility routines)
  ===============================================================*/

/*---------------------------------------------------------------
  IMEXGARKStepSetDenseOrder: Specifies the polynomial order for dense
  output.  Positive values are sent to the interpolation module;
  negative values imply to use the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetDenseOrder(void *arkode_mem, int dord)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetDenseOrder", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetDenseOrder(ark_mem, dord));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetErrHandlerFn: Specifies the error handler function
  ---------------------------------------------------------------*/
int IMEXGARKStepSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                                void *eh_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetErrHandlerFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetErrHandlerFn(ark_mem, ehfun, eh_data));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetErrFile: Specifies the FILE pointer for output (NULL
  means no messages)
  ---------------------------------------------------------------*/
int IMEXGARKStepSetErrFile(void *arkode_mem, FILE *errfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetErrFile", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetErrFile(ark_mem, errfp));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetUserData: Specifies the user data pointer for f
  ---------------------------------------------------------------*/
int IMEXGARKStepSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetUserData", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetUserData(ark_mem, user_data));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetDiagnostics: Specifies to enable solver diagnostics,
  and specifies the FILE pointer for output (diagfp==NULL
  disables output)
  ---------------------------------------------------------------*/
int IMEXGARKStepSetDiagnostics(void *arkode_mem, FILE *diagfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetDiagnostics", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetDiagnostics(ark_mem, diagfp));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetMaxNumSteps: Specifies the maximum number of
  integration steps
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxNumSteps(void *arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxNumSteps(ark_mem, mxsteps));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetMaxHnilWarns: Specifies the maximum number of warnings
  for small h
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxHnilWarns(void *arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxHnilWarns", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxHnilWarns(ark_mem, mxhnil));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetInitStep: Specifies the initial step size
  ---------------------------------------------------------------*/
int IMEXGARKStepSetInitStep(void *arkode_mem, realtype hin)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetInitStep(ark_mem, hin));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetMinStep: Specifies the minimum step size
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMinStep(void *arkode_mem, realtype hmin)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMinStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMinStep(ark_mem, hmin));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetMaxStep: Specifies the maximum step size
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxStep(void *arkode_mem, realtype hmax)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxStep(ark_mem, hmax));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetStopTime: Specifies the time beyond which the
  integration is not to proceed.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetStopTime(void *arkode_mem, realtype tstop)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetStopTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetStopTime(ark_mem, tstop));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetFixedStep: Specifies to use a fixed time step size
  instead of performing any form of temporal adaptivity.  ARKStep
  will use this step size for all steps (unless tstop is set, in
  which case it may need to modify that last step approaching
  tstop.  If any solver failure occurs in the timestepping
  module, ARKStep will typically immediately return with an error
  message indicating that the selected step size cannot be used.

  Any nonzero argument will result in the use of that fixed step
  size; an argument of 0 will re-enable temporal adaptivity.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetFixedStep(void *arkode_mem, realtype hfixed)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "ARKStepSetFixedStep",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* allocate or free adaptivity memory as needed */
  if (hfixed != ZERO) {
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
    }
  } else if (step_mem->hadapt_mem == NULL) {
    step_mem->hadapt_mem = arkAdaptInit();
    if (step_mem->hadapt_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                      "ARKStepSetFixedStep",
                      "Allocation of Step Adaptivity Structure Failed");
      return(ARK_MEM_FAIL);
    }
  }

  return(arkSetFixedStep(ark_mem, hfixed));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetRootDirection: Specifies the direction of zero-crossings
  to be monitored.  The default is to monitor both crossings.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetRootDirection(void *arkode_mem, int *rootdir)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetRootDirection", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetRootDirection(ark_mem, rootdir));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetNoInactiveRootWarn:  Disables issuing a warning if
  some root function appears to be identically zero at the
  beginning of the integration
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNoInactiveRootWarn(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetNoInactiveRootWarn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetNoInactiveRootWarn(ark_mem));
}

/*---------------------------------------------------------------
  IMEXGARKStepSetPostprocessStepFn:  Specifies a user-provided step
  postprocessing function having type ARKPostProcessStepFn.  A
  NULL input function disables step postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetPostprocessStepFn(void *arkode_mem,
                                     ARKPostProcessStepFn ProcessStep)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetPostprocessStepFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetPostprocessStepFn(ark_mem, ProcessStep));
}

/*---------------------------------------------------------------
  These wrappers for ARKLs module 'set' routines all are
  documented in arkode_arkstep.h.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetLinearSolver(void *arkode_mem, SUNLinearSolver LS,
                                SUNMatrix A) {
  return(arkLSSetLinearSolver(arkode_mem, LS, A)); }
int IMEXGARKStepSetMassLinearSolver(void *arkode_mem, SUNLinearSolver LS,
                                    SUNMatrix M, booleantype time_dep) {
  return(arkLSSetMassLinearSolver(arkode_mem, LS, M, time_dep)); }
int IMEXGARKStepSetJacFn(void *arkode_mem, ARKLsJacFn jac) {
  return(arkLSSetJacFn(arkode_mem, jac)); }
int IMEXGARKStepSetMassFn(void *arkode_mem, ARKLsMassFn mass) {
  return(arkLSSetMassFn(arkode_mem, mass)); }
int IMEXGARKStepSetMaxStepsBetweenJac(void *arkode_mem, long int msbj) {
  return(arkLSSetMaxStepsBetweenJac(arkode_mem, msbj)); }
int IMEXGARKStepSetEpsLin(void *arkode_mem, realtype eplifac) {
  return(arkLSSetEpsLin(arkode_mem, eplifac)); }
int IMEXGARKStepSetMassEpsLin(void *arkode_mem, realtype eplifac) {
  return(arkLSSetMassEpsLin(arkode_mem, eplifac)); }
int IMEXGARKStepSetPreconditioner(void *arkode_mem, ARKLsPrecSetupFn psetup,
                                  ARKLsPrecSolveFn psolve) {
  return(arkLSSetPreconditioner(arkode_mem, psetup, psolve)); }
int IMEXGARKStepSetMassPreconditioner(void *arkode_mem, ARKLsMassPrecSetupFn psetup,
                                      ARKLsMassPrecSolveFn psolve) {
  return(arkLSSetMassPreconditioner(arkode_mem, psetup, psolve)); }
int IMEXGARKStepSetJacTimes(void *arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                            ARKLsJacTimesVecFn jtimes) {
  return(arkLSSetJacTimes(arkode_mem, jtsetup, jtimes)); }
int IMEXGARKStepSetMassTimes(void *arkode_mem, ARKLsMassTimesSetupFn msetup,
                             ARKLsMassTimesVecFn mtimes, void *mtimes_data) {
  return(arkLSSetMassTimes(arkode_mem, msetup, mtimes, mtimes_data)); }

int IMEXGARKStepSetLinSysFn(void *arkode_mem, ARKLsLinSysFn linsys)
{
  return(arkLSSetLinSysFn(arkode_mem, linsys));
}


/*===============================================================
  IMEXGARKStep Optional output functions (wrappers for generic ARKode
  utility routines)
  ===============================================================*/

/*---------------------------------------------------------------
  IMEXGARKStepGetNumSteps:  Returns the current number of integration
  steps
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetNumSteps(ark_mem, nsteps));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetActualInitStep: Returns the step size used on the
  first step
  ---------------------------------------------------------------*/
int IMEXGARKStepGetActualInitStep(void *arkode_mem, realtype *hinused)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetActualInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetActualInitStep(ark_mem, hinused));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetLastStep: Returns the step size used on the last
  successful step
  ---------------------------------------------------------------*/
int IMEXGARKStepGetLastStep(void *arkode_mem, realtype *hlast)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetLastStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetLastStep(ark_mem, hlast));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetCurrentStep: Returns the step size to be attempted on
  the next step
  ---------------------------------------------------------------*/
int IMEXGARKStepGetCurrentStep(void *arkode_mem, realtype *hcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetCurrentStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetCurrentStep(ark_mem, hcur));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetCurrentTime: Returns the current value of the
  independent variable
  ---------------------------------------------------------------*/
int IMEXGARKStepGetCurrentTime(void *arkode_mem, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetCurrentTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetCurrentTime(ark_mem, tcur));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetCurrentY: Returns the current value of the
  dependent variable
  ---------------------------------------------------------------*/
int IMEXGARKStepGetCurrentState(void *arkode_mem, N_Vector *ycur)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepGetCurrentY", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  *ycur = ark_mem->ycur;
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  IMEXGARKStepGetCurrentGamma: Returns the current value for gamma
  ---------------------------------------------------------------*/
int IMEXGARKStepGetCurrentGamma(void *arkode_mem, realtype *gamma)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  retval = imexgarkStep_AccessStepMem(arkode_mem, NULL, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);
  *gamma = step_mem->gamma;
  return(retval);
}

/*---------------------------------------------------------------
  IMEXGARKStepGetTolScaleFactor: Returns a suggested factor for scaling
  tolerances
  ---------------------------------------------------------------*/
int IMEXGARKStepGetTolScaleFactor(void *arkode_mem, realtype *tolsfact)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetTolScaleFactor", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetTolScaleFactor(ark_mem, tolsfact));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetErrWeights: This routine returns the current error
  weight vector.
  ---------------------------------------------------------------*/
int IMEXGARKStepGetErrWeights(void *arkode_mem, N_Vector eweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetErrWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetErrWeights(ark_mem, eweight));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetResWeights: This routine returns the current residual
  weight vector.
  ---------------------------------------------------------------*/
int IMEXGARKStepGetResWeights(void *arkode_mem, N_Vector rweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetResWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetResWeights(ark_mem, rweight));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetWorkSpace: Returns integrator work space requirements
  ---------------------------------------------------------------*/
int IMEXGARKStepGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetWorkSpace", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetWorkSpace(ark_mem, lenrw, leniw));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetNumGEvals: Returns the current number of calls to g
  (for rootfinding)
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumGEvals(void *arkode_mem, long int *ngevals)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetNumGEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetNumGEvals(ark_mem, ngevals));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetRootInfo: Returns pointer to array rootsfound showing
  roots found
  ---------------------------------------------------------------*/
int IMEXGARKStepGetRootInfo(void *arkode_mem, int *rootsfound)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetRootInfo", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetRootInfo(ark_mem, rootsfound));
}

/*---------------------------------------------------------------
  ARKStepGetStepStats: Returns step statistics
  ---------------------------------------------------------------*/
int IMEXGARKStepGetStepStats(void *arkode_mem, long int *nsteps,
                             realtype *hinused, realtype *hlast,
                             realtype *hcur, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetStepStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetStepStats(ark_mem, nsteps, hinused, hlast, hcur, tcur));
}

/*---------------------------------------------------------------
  IMEXGARKStepGetReturnFlagName: translates from return flags IDs to
  names
  ---------------------------------------------------------------*/
char *IMEXGARKStepGetReturnFlagName(long int flag)
{ return(arkGetReturnFlagName(flag)); }

/*---------------------------------------------------------------
  These wrappers for ARKLs module 'get' routines all are
  documented in arkode_arkstep.h.
  ---------------------------------------------------------------*/
int IMEXGARKStepGetLinWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS) {
  return(arkLSGetWorkSpace(arkode_mem, lenrwLS, leniwLS)); }
int IMEXGARKStepGetNumJacEvals(void *arkode_mem, long int *njevals) {
  return(arkLSGetNumJacEvals(arkode_mem, njevals)); }
int IMEXGARKStepGetNumPrecEvals(void *arkode_mem, long int *npevals) {
  return(arkLSGetNumPrecEvals(arkode_mem, npevals)); }
int IMEXGARKStepGetNumPrecSolves(void *arkode_mem, long int *npsolves) {
  return(arkLSGetNumPrecSolves(arkode_mem, npsolves)); }
int IMEXGARKStepGetNumLinIters(void *arkode_mem, long int *nliters) {
  return(arkLSGetNumLinIters(arkode_mem, nliters)); }
int IMEXGARKStepGetNumLinConvFails(void *arkode_mem, long int *nlcfails) {
  return(arkLSGetNumConvFails(arkode_mem, nlcfails)); }
int IMEXGARKStepGetNumJTSetupEvals(void *arkode_mem, long int *njtsetups) {
  return(arkLSGetNumJTSetupEvals(arkode_mem, njtsetups)); }
int IMEXGARKStepGetNumJtimesEvals(void *arkode_mem, long int *njvevals) {
  return(arkLSGetNumJtimesEvals(arkode_mem, njvevals)); }
int IMEXGARKStepGetNumLinRhsEvals(void *arkode_mem, long int *nfevalsLS) {
  return(arkLSGetNumRhsEvals(arkode_mem, nfevalsLS)); } 
int IMEXGARKStepGetLastLinFlag(void *arkode_mem, long int *flag) {
  return(arkLSGetLastFlag(arkode_mem, flag)); }

int IMEXGARKStepGetMassWorkSpace(void *arkode_mem, long int *lenrwMLS, long int *leniwMLS) {
  return(arkLSGetMassWorkSpace(arkode_mem, lenrwMLS, leniwMLS)); }
int IMEXGARKStepGetNumMassSetups(void *arkode_mem, long int *nmsetups) {
  return(arkLSGetNumMassSetups(arkode_mem, nmsetups)); }
int IMEXGARKStepGetNumMassMultSetups(void *arkode_mem, long int *nmvsetups) {
  return(arkLSGetNumMassMatvecSetups(arkode_mem, nmvsetups)); }
int IMEXGARKStepGetNumMassMult(void *arkode_mem, long int *nmvevals) {
  return(arkLSGetNumMassMult(arkode_mem, nmvevals)); }
int IMEXGARKStepGetNumMassSolves(void *arkode_mem, long int *nmsolves) {
  return(arkLSGetNumMassSolves(arkode_mem, nmsolves)); }
int IMEXGARKStepGetNumMassPrecEvals(void *arkode_mem, long int *nmpevals) {
  return(arkLSGetNumMassPrecEvals(arkode_mem, nmpevals)); }
int IMEXGARKStepGetNumMassPrecSolves(void *arkode_mem, long int *nmpsolves) {
  return(arkLSGetNumMassPrecSolves(arkode_mem, nmpsolves)); }
int IMEXGARKStepGetNumMassIters(void *arkode_mem, long int *nmiters) {
  return(arkLSGetNumMassIters(arkode_mem, nmiters)); }
int IMEXGARKStepGetNumMassConvFails(void *arkode_mem, long int *nmcfails) {
  return(arkLSGetNumMassConvFails(arkode_mem, nmcfails)); }
int IMEXGARKStepGetNumMTSetups(void *arkode_mem, long int *nmtsetups) {
  return(arkLSGetNumMTSetups(arkode_mem, nmtsetups)); }
int IMEXGARKStepGetLastMassFlag(void *arkode_mem, long int *flag) {
  return(arkLSGetLastMassFlag(arkode_mem, flag)); }

char *IMEXGARKStepGetLinReturnFlagName(long int flag) {
  return(arkLSGetReturnFlagName(flag)); }



/*===============================================================
  ARKStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  IMEXGARKStepSetDefaults:

  Resets all IMEXGARKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.  Also leaves alone any data
  structures/options related to the ARKode infrastructure itself
  (e.g. root-finding).
  ---------------------------------------------------------------*/
int IMEXGARKStepSetDefaults(void* mem)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetDefaults", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetDefaults", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* Set default values for integrator optional inputs */
  step_mem->q                = Q_DEFAULT;      /* method order */
  step_mem->p                = 0;              /* embedding order */
  step_mem->hadapt_pq        = SUNFALSE;       /* use embedding order */
  step_mem->predictor        = 0;              /* trivial predictor */
  step_mem->linear           = SUNFALSE;       /* nonlinear problem */
  step_mem->linear_timedep   = SUNTRUE;        /* dfi/dy depends on t */
  if (step_mem->hadapt_mem != NULL) {
    step_mem->hadapt_mem->etamx1      = ETAMX1;     /* max change on first step */
    step_mem->hadapt_mem->etamxf      = ETAMXF;     /* max change on error-failed step */
    step_mem->hadapt_mem->small_nef   = SMALL_NEF;  /* num error fails before ETAMXF enforced */
    step_mem->hadapt_mem->etacf       = ETACF;      /* max change on convergence failure */
    step_mem->hadapt_mem->HAdapt      = NULL;       /* step adaptivity fn */
    step_mem->hadapt_mem->HAdapt_data = NULL;       /* step adaptivity data */
    step_mem->hadapt_mem->imethod     = 0;          /* PID controller */
    step_mem->hadapt_mem->cfl         = CFLFAC;     /* explicit stability factor */
    step_mem->hadapt_mem->safety      = SAFETY;     /* step adaptivity safety factor  */
    step_mem->hadapt_mem->bias        = BIAS;       /* step adaptivity error bias */
    step_mem->hadapt_mem->growth      = GROWTH;     /* step adaptivity growth factor */
    step_mem->hadapt_mem->lbound      = HFIXED_LB;  /* step adaptivity no-change lower bound */
    step_mem->hadapt_mem->ubound      = HFIXED_UB;  /* step adaptivity no-change upper bound */
    step_mem->hadapt_mem->k1          = AD0_K1;     /* step adaptivity parameter */
    step_mem->hadapt_mem->k2          = AD0_K2;     /* step adaptivity parameter */
    step_mem->hadapt_mem->k3          = AD0_K3;     /* step adaptivity parameter */
  }
  step_mem->maxcor           = MAXCOR;         /* max nonlinear iters/stage */
  step_mem->maxnef           = MAXNEF;         /* max error test fails */
  step_mem->maxncf           = MAXNCF;         /* max convergence fails */
  step_mem->nlscoef          = NLSCOEF;        /* nonlinear tolerance coefficient */
  step_mem->crdown           = CRDOWN;         /* nonlinear convergence estimate coeff. */
  step_mem->rdiv             = RDIV;           /* nonlinear divergence tolerance */
  step_mem->dgmax            = DGMAX;          /* max step change before recomputing J or P */
  step_mem->msbp             = MSBP;           /* max steps between updates to J or P */
  step_mem->stages           = 0;              /* no stages */
  step_mem->istage           = 0;              /* current stage */
  step_mem->Bee              = NULL;           /* no Butcher tables */
  step_mem->Bii              = NULL;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetLinear:

  Specifies that the implicit portion of the problem is linear,
  and to tighten the linear solver tolerances while taking only
  one Newton iteration.  DO NOT USE IN COMBINATION WITH THE
  FIXED-POINT SOLVER.  Automatically tightens DeltaGammaMax
  to ensure that step size changes cause Jacobian recomputation.

  The argument should be 1 or 0, where 1 indicates that the
  Jacobian of fi with respect to y depends on time, and
  0 indicates that it is not time dependent.  Alternately, when
  using an iterative linear solver this flag denotes time
  dependence of the preconditioner.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetLinear(void *arkode_mem, int timedepend)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetLinear",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set parameters */
  step_mem->linear = SUNTRUE;
  step_mem->linear_timedep = (timedepend == 1);
  step_mem->dgmax = RCONST(100.0)*UNIT_ROUNDOFF;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetNonlinear:

  Specifies that the implicit portion of the problem is nonlinear.
  Used to undo a previous call to IMEXGARKStepSetLinear.  Automatically
  loosens DeltaGammaMax back to default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNonlinear(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetNonlinear",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set parameters */
  step_mem->linear = SUNFALSE;
  step_mem->linear_timedep = SUNTRUE;
  step_mem->dgmax = DGMAX;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetARKTables:

  Specifies to use customized Butcher tables for the ImEx system
  (automatically calls IMEXGARKStepSetImEx).

  If either b2e == NULL or b2i == NULL, then the method is
  automatically flagged as a fixed-step method; a user MUST also
  call either ARKodeSetFixedStep or ARKodeSetInitStep to set the
  desired time step size.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IMEXGARKStepSetTables(void *mem, int s,
                                          int q, int p,
                                          realtype *ce, realtype *ci,
                                          realtype *Aee, realtype *Aei,
                                          realtype *Aie, realtype *Aii,
                                          realtype *be, realtype *bi,
                                          realtype *de, realtype *di)
{
  int i, j;
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetButcherTables", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetButcherTables", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* check for legal inputs */
  if ((ci  == NULL) || (ce  == NULL) ||
      (Aee == NULL) || (Aei == NULL) ||
      (Aie == NULL) || (Aii == NULL) ||
      (bi  == NULL) || (be  == NULL)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetButcherTables", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  ARKodeButcherTable_Free(step_mem->Bee); step_mem->Bee = NULL;
  ARKodeButcherTable_Free(step_mem->Bei); step_mem->Bei = NULL; /* don't need a full Butcher table (only A) */
  ARKodeButcherTable_Free(step_mem->Bie); step_mem->Bie = NULL; /* don't need a full Butcher table (only A) */
  ARKodeButcherTable_Free(step_mem->Bii); step_mem->Bii = NULL;

  /* set the relevant parameters */
  step_mem->stages = s;
  step_mem->q = q;
  step_mem->p = p;
  step_mem->Bee = ARKodeButcherTable_Alloc(s, (de != NULL));
  step_mem->Bei = ARKodeButcherTable_Alloc(s, (de != NULL)); /* don't need a full Butcher table (only A) */
  step_mem->Bie = ARKodeButcherTable_Alloc(s, (de != NULL)); /* don't need a full Butcher table (only A) */
  step_mem->Bii = ARKodeButcherTable_Alloc(s, (di != NULL));
  step_mem->Bee->p = p;
  step_mem->Bii->q = q;
  step_mem->Bee->q = q;
  step_mem->Bii->p = p;
  for (i=0; i<s; i++) {
    step_mem->Bee->c[i] = ce[i];
    step_mem->Bii->c[i] = ci[i];
    step_mem->Bee->b[i] = be[i];
    step_mem->Bii->b[i] = bi[i];
    for (j=0; j<s; j++) {
      step_mem->Bee->A[i][j] = Aee[i*s + j];
      step_mem->Bei->A[i][j] = Aei[i*s + j]; /* not necessarily square */
      step_mem->Bie->A[i][j] = Aie[i*s + j]; /* not necessarily square */
      step_mem->Bii->A[i][j] = Aii[i*s + j];
    }
  }

  /* set embeddings (if applicable), otherwise set as fixed-step method */
  if ((de == NULL) || (di == NULL)) {
    arkode_mem->fixedstep = SUNTRUE;
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
    }
  } else {
    for (i=0; i<s; i++) {
      step_mem->Bee->d[i] = de[i];
      step_mem->Bii->d[i] = di[i];
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetCFLFraction:

  Specifies the safety factor to use on the maximum explicitly-
  stable step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetCFLFraction(void *arkode_mem, realtype cfl_frac)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetCFLFraction",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetCFLFraction",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if (cfl_frac >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetCFLFraction", "Illegal CFL fraction");
    return(ARK_ILL_INPUT);
  }

  /* set positive-valued parameters, otherwise set default */
  if (cfl_frac <= ZERO) {
    hadapt_mem->cfl = CFLFAC;
  } else {
    hadapt_mem->cfl = cfl_frac;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetSafetyFactor:

  Specifies the safety factor to use on the error-based predicted
  time step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetSafetyFactor(void *arkode_mem, realtype safety)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetSafetyFactor",
                                     &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetSafetyFactoy",MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if (safety >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetSafetyFactor", "Illegal safety factor");
    return(ARK_ILL_INPUT);
  }

  /* set positive-valued parameters, otherwise set default */
  if (safety <= ZERO) {
    hadapt_mem->safety = SAFETY;
  } else {
    hadapt_mem->safety = safety;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetErrorBias:

  Specifies the error bias to use when performing adaptive-step
  error control.  Allowable values must be >= 1.0.  Any illegal
  value implies a reset to the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetErrorBias(void *arkode_mem, realtype bias)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetErrorBias",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "IMEXGARKStepSetErrorBias", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (bias < 1.0) {
    hadapt_mem->bias = BIAS;
  } else {
    hadapt_mem->bias = bias;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxGrowth:

  Specifies the maximum step size growth factor to be allowed
  between successive integration steps.  Note: the first step uses
  a separate maximum growth factor.  Allowable values must be
  > 1.0.  Any illegal value implies a reset to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxGrowth(void *arkode_mem, realtype mx_growth)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetMaxGrowth",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetMaxGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (mx_growth == ZERO) {
    hadapt_mem->growth = GROWTH;
  } else {
    hadapt_mem->growth = mx_growth;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetFixedStepBounds:

  Specifies the step size growth interval within which the step
  size will remain unchanged.  Allowable values must enclose the
  value 1.0.  Any illegal interval implies a reset to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetFixedStepBounds",
                                         &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetFixedStepBounds", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* set allowable interval, otherwise set defaults */
  if ((lb <= 1.0) && (ub >= 1.0)) {
    hadapt_mem->lbound = lb;
    hadapt_mem->ubound = ub;
  } else {
    hadapt_mem->lbound = HFIXED_LB;
    hadapt_mem->ubound = HFIXED_UB;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetAdaptivityMethod:

  Specifies the built-in time step adaptivity algorithm (and
  optionally, its associated parameters) to use.  All parameters
  will be checked for validity when used by the solver.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetAdaptivityMethod(void *arkode_mem, int imethod,
                               int idefault, int pq,
                               realtype *adapt_params)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetAdaptivityMethod",
                                     &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetAdaptivityMethod", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if ((imethod > 5) || (imethod < 0)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  /* set adaptivity method */
  hadapt_mem->imethod = imethod;

  /* set flag whether to use p or q */
  step_mem->hadapt_pq = (pq != 0);

  /* set method parameters */
  if (idefault == 1) {
    switch (hadapt_mem->imethod) {
    case (0):
      hadapt_mem->k1 = AD0_K1;
      hadapt_mem->k2 = AD0_K2;
      hadapt_mem->k3 = AD0_K3; break;
    case (1):
      hadapt_mem->k1 = AD1_K1;
      hadapt_mem->k2 = AD1_K2; break;
    case (2):
      hadapt_mem->k1 = AD2_K1; break;
    case (3):
      hadapt_mem->k1 = AD3_K1;
      hadapt_mem->k2 = AD3_K2; break;
    case (4):
      hadapt_mem->k1 = AD4_K1;
      hadapt_mem->k2 = AD4_K2; break;
    case (5):
      hadapt_mem->k1 = AD5_K1;
      hadapt_mem->k2 = AD5_K2;
      hadapt_mem->k3 = AD5_K3; break;
    }
  } else {
    hadapt_mem->k1 = adapt_params[0];
    hadapt_mem->k2 = adapt_params[1];
    hadapt_mem->k3 = adapt_params[2];
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetAdaptivityFn:

  Specifies the user-provided time step adaptivity function to use.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun,
                                void *h_data)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetAdaptivityFn",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetAdaptivityFn", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* NULL hfun sets default, otherwise set inputs */
  if (hfun == NULL) {
    hadapt_mem->HAdapt      = NULL;
    hadapt_mem->HAdapt_data = NULL;
    hadapt_mem->imethod     = 0;
  } else {
    hadapt_mem->HAdapt      = hfun;
    hadapt_mem->HAdapt_data = h_data;
    hadapt_mem->imethod     = -1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxFirstGrowth:

  Specifies the user-provided time step adaptivity constant
  etamx1.  Legal values are greater than 1.0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetMaxFirstGrowth",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetMaxFirstGrowth",MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* if argument legal set it, otherwise set default */
  if (etamx1 <= ONE) {
    hadapt_mem->etamx1 = ETAMX1;
  } else {
    hadapt_mem->etamx1 = etamx1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxEFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etamxf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetMaxEFailGrowth",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "IMEXGARKStepSetMaxEFailGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* if argument legal set it, otherwise set default */
  if ((etamxf <= ZERO) || (etamxf > ONE)) {
    hadapt_mem->etamxf = ETAMXF;
  } else {
    hadapt_mem->etamxf = etamxf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetSmallNumEFails:

  Specifies the user-provided time step adaptivity constant
  small_nef.  Legal values are > 0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetSmallNumEFails(void *arkode_mem, int small_nef)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetSmallNumEFails",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetSmallNumEFails", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0) {
    hadapt_mem->small_nef = SMALL_NEF;
  } else {
    hadapt_mem->small_nef = small_nef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxCFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etacf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxCFailGrowth(void *arkode_mem, realtype etacf)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetMaxCFailGrowth",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetMaxCFailGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* if argument legal set it, otherwise set default */
  if ((etacf <= ZERO) || (etacf > ONE)) {
    hadapt_mem->etacf = ETACF;
  } else {
    hadapt_mem->etacf = etacf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNonlinCRDown(void *arkode_mem, realtype crdown)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetNonlinCRDown",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (crdown <= ZERO) {
    step_mem->crdown = CRDOWN;
  } else {
    step_mem->crdown = crdown;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNonlinRDiv(void *arkode_mem, realtype rdiv)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetNonlinRDiv",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (rdiv <= ZERO) {
    step_mem->rdiv = RDIV;
  } else {
    step_mem->rdiv = rdiv;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply
  a reset to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetDeltaGammaMax(void *arkode_mem, realtype dgmax)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetDeltaGammaMax",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (dgmax <= ZERO) {
    step_mem->dgmax = DGMAX;
  } else {
    step_mem->dgmax = dgmax;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxStepsBetweenLSet:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the number of time steps to wait
  before calling lsetup; negative values imply recomputation of
  lsetup at each nonlinear solve; a zero value implies a reset
  to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxStepsBetweenLSet(void *arkode_mem, int msbp)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetMaxStepsBetweenLSet",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (msbp == 0) {
    step_mem->msbp = MSBP;
  } else {
    step_mem->msbp = msbp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.
  Non-default choices are {1,2,3,4}, all others will use default
  (trivial) predictor.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetPredictorMethod(void *arkode_mem, int pred_method)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetPredictorMethod",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set parameters */
  step_mem->predictor = pred_method;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetStabilityFn:

  Specifies the user-provided explicit time step stability
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int IMEXGARKStepSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab,
                               void *estab_data)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IEMXGARKStepSetStabilityFn",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetStabilityFn", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* NULL argument sets default, otherwise set inputs */
  if (EStab == NULL) {
    hadapt_mem->expstab    = arkExpStab;
    hadapt_mem->estab_data = ark_mem;
  } else {
    hadapt_mem->expstab    = EStab;
    hadapt_mem->estab_data = estab_data;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxErrTestFails(void *arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetMaxErrTestFails",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    step_mem->maxnef = MAXNEF;
  } else {
    step_mem->maxnef = maxnef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxConvFails:

  Specifies the maximum number of nonlinear convergence failures
  during one step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxConvFails(void *arkode_mem, int maxncf)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetMaxConvFails",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* argument <= 0 sets default, otherwise set input */
  if (maxncf <= 0) {
    step_mem->maxncf = MAXNCF;
  } else {
    step_mem->maxncf = maxncf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxNonlinIters(void *arkode_mem, int maxcor)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetMaxNonlinIters",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return error message if no NLS module is present */
  if (step_mem->NLS == NULL) {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetMaxNonlinIters",
                    "No SUNNonlinearSolver object is present");
    return(ARK_ILL_INPUT);
  }

  /* argument <= 0 sets default, otherwise set input */
  if (maxcor <= 0) {
    step_mem->maxcor = MAXCOR;
  } else {
    step_mem->maxcor = maxcor;
  }

  /* send argument to NLS structure */
  retval = SUNNonlinSolSetMaxIters(step_mem->NLS, step_mem->maxcor);
  if (retval != SUN_NLS_SUCCESS) {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepSetMaxNonlinIters",
                    "Error setting maxcor in SUNNonlinearSolver object");
    return(ARK_NLS_OP_ERR);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNonlinConvCoef(void *arkode_mem, realtype nlscoef)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepSetNonlinConvCoef",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) {
    step_mem->nlscoef = NLSCOEF;
  } else {
    step_mem->nlscoef = nlscoef;
  }

  return(ARK_SUCCESS);
}


/*===============================================================
  IMEXGARKStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  IMEXGARKStepGetNumExpSteps:

  Returns the current number of stability-limited steps
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetNumExpSteps",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if step adaptivity structure not allocated, just return 0 */
  if (step_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = step_mem->hadapt_mem->nst_exp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumAccSteps:

  Returns the current number of accuracy-limited steps
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetNumAccSteps",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if step adaptivity structure not allocated, just return 0 */
  if (step_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = step_mem->hadapt_mem->nst_acc;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumStepAttempts:

  Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumStepAttempts(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetNumStepAttempts",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get value from step_mem */
  *nsteps = step_mem->nst_attempts;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumRhsEvals:

  Returns the current number of calls to fe and fi
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumRhsEvals(void *arkode_mem, long int *fe_evals,
                               long int *fi_evals)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetNumRhsEvals",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get values from step_mem */
  *fe_evals = step_mem->nfe;
  *fi_evals = step_mem->nfi;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumLinSolvSetups(void *arkode_mem, long int *nlinsetups)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetNumLinSolvSetups",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get value from step_mem */
  *nlinsetups = step_mem->nsetups;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumErrTestFails(void *arkode_mem, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetNumErrTestFails",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get value from step_mem */
  *netfails = step_mem->netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetCurrentButcherTables:

  Sets pointers to the explicit and implicit Butcher tables
  currently in use.
  ---------------------------------------------------------------*/
int IMEXGARKStepGetCurrentButcherTables(void *arkode_mem,
                                        ARKodeButcherTable *Bee,
                                        ARKodeButcherTable *Bei,
                                        ARKodeButcherTable *Bie,
                                        ARKodeButcherTable *Bii)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetCurrentButcherTables",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get tables from step_mem */
  *Bii = step_mem->Bii;
  *Bee = step_mem->Bee;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetEstLocalErrors: (updated to the correct vector, but
  need to verify that it is unchanged between filling the
  estimated error and the end of the time step)

  Returns an estimate of the local error
  ---------------------------------------------------------------*/
int IMEXGARKStepGetEstLocalErrors(void *arkode_mem, N_Vector ele)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetEstLocalErrors",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* copy vector to output */
  N_VScale(ONE, ark_mem->tempv1, ele);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int IMEXGARKStepGetTimestepperStats(void *arkode_mem, long int *expsteps,
                                    long int *accsteps, long int *step_attempts,
                                    long int *fe_evals, long int *fi_evals,
                                    long int *nlinsetups, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetTimestepperStats",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if step adaptivity structure not allocated,
     just set expsteps and accsteps to 0 */
  if (step_mem->hadapt_mem == NULL) {
    *expsteps = 0;
    *accsteps = 0;
  } else {
    *expsteps = step_mem->hadapt_mem->nst_exp;
    *accsteps = step_mem->hadapt_mem->nst_acc;
  }

  /* set remaining outputs from step_mem */
  *step_attempts = step_mem->nst_attempts;
  *fe_evals      = step_mem->nfe;
  *fi_evals      = step_mem->nfi;
  *nlinsetups    = step_mem->nsetups;
  *netfails      = step_mem->netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumNonlinSolvIters(void *arkode_mem, long int *nniters)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetNumNonlinSolvIters",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if a NLS object is present, set output from that; otherwise
     we took zero iterations */
  if (step_mem->NLS) {
    retval = SUNNonlinSolGetNumIters(step_mem->NLS, nniters);
    if (retval != SUN_NLS_SUCCESS) {
      arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::IMEXGARKStep",
                      "IMEXGARKStepGetNumNonlinSolvIters",
                      "Error retrieving nniters from SUNNonlinearSolver");
      return(ARK_NLS_OP_ERR);
    }
  } else {
    *nniters = 0;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumNonlinSolvConvFails(void *arkode_mem, long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetNumNonlinSolvConvFails",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set output from step_mem */
  *nncfails = step_mem->ncfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNonlinSolvStats(void *arkode_mem, long int *nniters,
                                   long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepGetNonlinSolvStats",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set outputs from NLS module and step_mem structure (if present);
     otherwise there were zero iterations and no nonlinear failures */
  if (step_mem->NLS) {
    retval = SUNNonlinSolGetNumIters(step_mem->NLS, nniters);
    if (retval != SUN_NLS_SUCCESS) {
      arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::IMEXGARKStep",
                      "IMEXGARKStepGetNonlinSolvStats",
                      "Error retrieving nniters from SUNNonlinearSolver");
      return(ARK_NLS_OP_ERR);
    }
    *nncfails = step_mem->ncfn;
  } else {
    *nniters = 0;
    *nncfails = 0;
  }

  return(ARK_SUCCESS);
}


/*===============================================================
  IMEXGARKStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  IMEXGARKStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int IMEXGARKStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int flag, retval;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepWriteParameters",
                                      &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* output ARKode infrastructure parameters first */
  flag = arkWriteParameters(ark_mem, fp);
  if (flag != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepWriteParameters",
                    "Error writing ARKode infrastructure parameters");
    return(flag);
  }

  /* print integrator parameters to file */
  fprintf(fp, "IMEXGARKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n",step_mem->q);
  if (step_mem->linear) {
    fprintf(fp, "  Linear implicit problem");
    if (step_mem->linear_timedep) {
      fprintf(fp, " (time-dependent Jacobian)\n");
    } else {
      fprintf(fp, " (time-independent Jacobian)\n");
    }
  }
  fprintf(fp, "  ImEx integrator\n");
  if (step_mem->hadapt_mem != NULL) {
    fprintf(fp, "  Maximum step increase (first step) = %"RSYM"\n",
            step_mem->hadapt_mem->etamx1);
    fprintf(fp, "  Step reduction factor on multiple error fails = %"RSYM"\n",
            step_mem->hadapt_mem->etamxf);
    fprintf(fp, "  Minimum error fails before above factor is used = %i\n",
            step_mem->hadapt_mem->small_nef);
    fprintf(fp, "  Step reduction factor on nonlinear convergence failure = %"RSYM"\n",
            step_mem->hadapt_mem->etacf);
    fprintf(fp, "  Explicit safety factor = %"RSYM"\n",
            step_mem->hadapt_mem->cfl);
    if (step_mem->hadapt_mem->HAdapt == NULL) {
      fprintf(fp, "  Time step adaptivity method %i\n", step_mem->hadapt_mem->imethod);
      fprintf(fp, "     Safety factor = %"RSYM"\n", step_mem->hadapt_mem->safety);
      fprintf(fp, "     Bias factor = %"RSYM"\n", step_mem->hadapt_mem->bias);
      fprintf(fp, "     Growth factor = %"RSYM"\n", step_mem->hadapt_mem->growth);
      fprintf(fp, "     Step growth lower bound = %"RSYM"\n", step_mem->hadapt_mem->lbound);
      fprintf(fp, "     Step growth upper bound = %"RSYM"\n", step_mem->hadapt_mem->ubound);
      fprintf(fp, "     k1 = %"RSYM"\n", step_mem->hadapt_mem->k1);
      fprintf(fp, "     k2 = %"RSYM"\n", step_mem->hadapt_mem->k2);
      fprintf(fp, "     k3 = %"RSYM"\n", step_mem->hadapt_mem->k3);
      if (step_mem->hadapt_mem->expstab == arkExpStab) {
        fprintf(fp, "  Default explicit stability function\n");
      } else {
        fprintf(fp, "  User provided explicit stability function\n");
      }
    } else {
      fprintf(fp, "  User provided time step adaptivity function\n");
    }
  }

  fprintf(fp, "  Maximum number of error test failures = %i\n",step_mem->maxnef);
  fprintf(fp, "  Maximum number of convergence test failures = %i\n",step_mem->maxncf);
  fprintf(fp, "  Implicit predictor method = %i\n",step_mem->predictor);
  fprintf(fp, "  Implicit solver tolerance coefficient = %"RSYM"\n",step_mem->nlscoef);
  fprintf(fp, "  Maximum number of nonlinear corrections = %i\n",step_mem->maxcor);
  fprintf(fp, "  Nonlinear convergence rate constant = %"RSYM"\n",step_mem->crdown);
  fprintf(fp, "  Nonlinear divergence tolerance = %"RSYM"\n",step_mem->rdiv);
  fprintf(fp, "  Gamma factor LSetup tolerance = %"RSYM"\n",step_mem->dgmax);
  fprintf(fp, "  Number of steps between LSetup calls = %i\n",step_mem->msbp);
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepWriteButcher:

  Outputs Butcher tables to the provided file pointer.
  ---------------------------------------------------------------*/
int IMEXGARKStepWriteButcher(void *arkode_mem, FILE *fp)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeIMEXGARKStepMem step_mem;
  int i, j;

  /* access ARKodeARKStepMem structure */
  retval = imexgarkStep_AccessStepMem(arkode_mem, "IMEXGARKStepWriteButcher",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  if ( (step_mem->Bee == NULL) || (step_mem->Bii == NULL) ) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::IMEXGARKStep",
                    "IMEXGARKStepWriteButcher",
                    "Butcher tables must both be non-NULL");
    return(ARK_ILL_INPUT);
  }

  /* print Butcher tables to file */
  fprintf(fp, "\nIMEXGARKStep Butcher table (stages = %i):\n\n", step_mem->stages);

  fprintf(fp, "  Explicit Butcher table:\n");
  ARKodeButcherTable_Write(step_mem->Bee, fp);

  for (i=0; i<step_mem->stages; i++) {
    fprintf(fp, "     %"RSYM"",step_mem->Bee->c[i]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++)
      fprintf(fp, " %"RSYM"",step_mem->Bee->A[i][j]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++)
      fprintf(fp, " %"RSYM"",step_mem->Bei->A[i][j]);
    fprintf(fp, " | ");
    fprintf(fp,"\n");
  }

  fprintf(fp, "      ");
  for (j=0; j < (step_mem->stages*2)+1; j++)
    fprintf(fp, "-----------");
  fprintf(fp,"\n");

  for (i=0; i<step_mem->stages; i++) {
    fprintf(fp, "     %"RSYM"",step_mem->Bii->c[i]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++)
      fprintf(fp, " %"RSYM"",step_mem->Bie->A[i][j]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++)
      fprintf(fp, " %"RSYM"",step_mem->Bii->A[i][j]);
    fprintf(fp, " | ");
    fprintf(fp,"\n");
  }

  fprintf(fp, "      ");
  for (j=0; j < (step_mem->stages*2)+1; j++)
    fprintf(fp, "-----------");
  fprintf(fp,"\n");

  fprintf(fp, "                ");
  fprintf(fp, " | ");
  for (j=0; j<step_mem->stages; j++)
    fprintf(fp, " %"RSYM"",step_mem->Bee->b[j]);
  fprintf(fp, " | ");
  for (j=0; j<step_mem->stages; j++)
    fprintf(fp, " %"RSYM"",step_mem->Bii->b[j]);
  fprintf(fp, " | ");
  fprintf(fp,"\n");

  fprintf(fp, "                ");
  fprintf(fp, " | ");
  if (step_mem->Bee->d != NULL) {
    for (j=0; j<step_mem->stages; j++)
      fprintf(fp, " %"RSYM"",step_mem->Bee->d[j]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++)
      fprintf(fp, " %"RSYM"",step_mem->Bii->d[j]);
    fprintf(fp, " | ");
    fprintf(fp,"\n");
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
