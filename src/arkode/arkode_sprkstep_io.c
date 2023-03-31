/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for the optional input and
 * output functions for the ARKODE SPRKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here,
 * with slightly different names.  The code transition will be
 * minimal, but the documentation changes will be significant.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode/arkode_sprk.h"
#include "arkode_sprkstep_impl.h"

/*===============================================================
  SPRKStep Optional input functions (wrappers for generic ARKODE
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int SPRKStepSetDenseOrder(void* arkode_mem, int dord)
{
  return (SPRKStepSetInterpolantDegree(arkode_mem, dord));
}

int SPRKStepSetInterpolantDegree(void* arkode_mem, int degree)
{
  if (degree < 0) degree = ARK_INTERP_MAX_DEGREE;
  return (arkSetInterpolantDegree(arkode_mem, degree));
}

int SPRKStepSetInterpolantType(void* arkode_mem, int itype)
{
  return (arkSetInterpolantType(arkode_mem, itype));
}

int SPRKStepSetErrHandlerFn(void* arkode_mem, ARKErrHandlerFn ehfun, void* eh_data)
{
  return (arkSetErrHandlerFn(arkode_mem, ehfun, eh_data));
}

int SPRKStepSetErrFile(void* arkode_mem, FILE* errfp)
{
  return (arkSetErrFile(arkode_mem, errfp));
}

int SPRKStepSetDiagnostics(void* arkode_mem, FILE* diagfp)
{
  return (arkSetDiagnostics(arkode_mem, diagfp));
}

int SPRKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)
{
  return (arkSetMaxNumSteps(arkode_mem, mxsteps));
}

int SPRKStepSetMaxHnilWarns(void* arkode_mem, int mxhnil)
{
  return (arkSetMaxHnilWarns(arkode_mem, mxhnil));
}

int SPRKStepSetInitStep(void* arkode_mem, realtype hin)
{
  return (arkSetInitStep(arkode_mem, hin));
}

int SPRKStepSetMinStep(void* arkode_mem, realtype hmin)
{
  return (arkSetMinStep(arkode_mem, hmin));
}

int SPRKStepSetMaxStep(void* arkode_mem, realtype hmax)
{
  return (arkSetMaxStep(arkode_mem, hmax));
}

int SPRKStepSetStopTime(void* arkode_mem, realtype tstop)
{
  return (arkSetStopTime(arkode_mem, tstop));
}

int SPRKStepSetRootDirection(void* arkode_mem, int* rootdir)
{
  return (arkSetRootDirection(arkode_mem, rootdir));
}

int SPRKStepSetNoInactiveRootWarn(void* arkode_mem)
{
  return (arkSetNoInactiveRootWarn(arkode_mem));
}

// int SPRKStepSetConstraints(void* arkode_mem, N_Vector constraints)
// {
//   return (arkSetConstraints(arkode_mem, constraints));
// }

int SPRKStepSetMaxNumConstrFails(void* arkode_mem, int maxfails)
{
  return (arkSetMaxNumConstrFails(arkode_mem, maxfails));
}

int SPRKStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep)
{
  return (arkSetPostprocessStepFn(arkode_mem, ProcessStep));
}

int SPRKStepSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage)
{
  return (arkSetPostprocessStageFn(arkode_mem, ProcessStage));
}

// int SPRKStepSetCFLFraction(void* arkode_mem, realtype cfl_frac)
// {
//   return (arkSetCFLFraction(arkode_mem, cfl_frac));
// }

// int SPRKStepSetSafetyFactor(void* arkode_mem, realtype safety)
// {
//   return (arkSetSafetyFactor(arkode_mem, safety));
// }

// int SPRKStepSetErrorBias(void* arkode_mem, realtype bias)
// {
//   return (arkSetErrorBias(arkode_mem, bias));
// }

// int SPRKStepSetMaxGrowth(void* arkode_mem, realtype mx_growth)
// {
//   return (arkSetMaxGrowth(arkode_mem, mx_growth));
// }

// int SPRKStepSetMinReduction(void* arkode_mem, realtype eta_min)
// {
//   return (arkSetMinReduction(arkode_mem, eta_min));
// }

// int SPRKStepSetFixedStepBounds(void* arkode_mem, realtype lb, realtype ub)
// {
//   return (arkSetFixedStepBounds(arkode_mem, lb, ub));
// }

int SPRKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault,
                                int pq, realtype adapt_params[3])
{
  return (arkSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, adapt_params));
}

int SPRKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data)
{
  return (arkSetAdaptivityFn(arkode_mem, hfun, h_data));
}

// int SPRKStepSetMaxFirstGrowth(void* arkode_mem, realtype etamx1)
// {
//   return (arkSetMaxFirstGrowth(arkode_mem, etamx1));
// }

// int SPRKStepSetMaxEFailGrowth(void* arkode_mem, realtype etamxf)
// {
//   return (arkSetMaxEFailGrowth(arkode_mem, etamxf));
// }

// int SPRKStepSetSmallNumEFails(void* arkode_mem, int small_nef)
// {
//   return (arkSetSmallNumEFails(arkode_mem, small_nef));
// }

// int SPRKStepSetMaxCFailGrowth(void* arkode_mem, realtype etacf)
// {
//   return (arkSetMaxCFailGrowth(arkode_mem, etacf));
// }

int SPRKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)
{
  return (arkSetStabilityFn(arkode_mem, EStab, estab_data));
}

// int SPRKStepSetMaxErrTestFails(void* arkode_mem, int maxnef)
// {
//   return (arkSetMaxErrTestFails(arkode_mem, maxnef));
// }

// int SPRKStepSetMaxConvFails(void* arkode_mem, int maxncf)
// {
//   return (arkSetMaxConvFails(arkode_mem, maxncf));
// }

int SPRKStepSetFixedStep(void* arkode_mem, realtype hfixed)
{
  return (arkSetFixedStep(arkode_mem, hfixed));
}

/*===============================================================
  SPRKStep Optional output functions (wrappers for generic ARKODE
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int SPRKStepGetNumStepAttempts(void* arkode_mem, long int* nstep_attempts)
{
  return (arkGetNumStepAttempts(arkode_mem, nstep_attempts));
}

int SPRKStepGetNumSteps(void* arkode_mem, long int* nsteps)
{
  return (arkGetNumSteps(arkode_mem, nsteps));
}

int SPRKStepGetActualInitStep(void* arkode_mem, realtype* hinused)
{
  return (arkGetActualInitStep(arkode_mem, hinused));
}

int SPRKStepGetLastStep(void* arkode_mem, realtype* hlast)
{
  return (arkGetLastStep(arkode_mem, hlast));
}

int SPRKStepGetCurrentStep(void* arkode_mem, realtype* hcur)
{
  return (arkGetCurrentStep(arkode_mem, hcur));
}

int SPRKStepGetCurrentTime(void* arkode_mem, realtype* tcur)
{
  return (arkGetCurrentTime(arkode_mem, tcur));
}

int SPRKStepGetCurrentState(void* arkode_mem, N_Vector* state)
{
  return (arkGetCurrentState(arkode_mem, state));
}

// int SPRKStepGetTolScaleFactor(void* arkode_mem, realtype* tolsfact)
// {
//   return (arkGetTolScaleFactor(arkode_mem, tolsfact));
// }

// int SPRKStepGetErrWeights(void* arkode_mem, N_Vector eweight)
// {
//   return (arkGetErrWeights(arkode_mem, eweight));
// }

// int SPRKStepGetResWeights(void* arkode_mem, N_Vector rweight)
// {
//   return (arkGetResWeights(arkode_mem, rweight));
// }

int SPRKStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)
{
  return (arkGetWorkSpace(arkode_mem, lenrw, leniw));
}

int SPRKStepGetStepStats(void* arkode_mem, long int* nsteps, realtype* hinused,
                         realtype* hlast, realtype* hcur, realtype* tcur)
{
  return (arkGetStepStats(arkode_mem, nsteps, hinused, hlast, hcur, tcur));
}

int SPRKStepGetNumConstrFails(void* arkode_mem, long int* nconstrfails)
{
  return (arkGetNumConstrFails(arkode_mem, nconstrfails));
}

int SPRKStepGetNumExpSteps(void* arkode_mem, long int* nsteps)
{
  return (arkGetNumExpSteps(arkode_mem, nsteps));
}

int SPRKStepGetNumAccSteps(void* arkode_mem, long int* nsteps)
{
  return (arkGetNumAccSteps(arkode_mem, nsteps));
}

int SPRKStepGetNumErrTestFails(void* arkode_mem, long int* netfails)
{
  return (arkGetNumErrTestFails(arkode_mem, netfails));
}

int SPRKStepGetNumStepSolveFails(void* arkode_mem, long int* nncfails)
{
  return (arkGetNumStepSolveFails(arkode_mem, nncfails));
}

int SPRKStepGetUserData(void* arkode_mem, void** user_data)
{
  return (arkGetUserData(arkode_mem, user_data));
}

char* SPRKStepGetReturnFlagName(long int flag)
{
  return (arkGetReturnFlagName(flag));
}

/*===============================================================
  SPRKStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  SPRKStepSetUserData:

  Wrapper for generic arkSetUserData and arkLSSetUserData
  routines.
  ---------------------------------------------------------------*/
int SPRKStepSetUserData(void* arkode_mem, void* user_data)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetUserData", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) return (retval);

  /* set user_data in ARKODE mem */
  retval = arkSetUserData(arkode_mem, user_data);
  if (retval != ARK_SUCCESS) return (retval);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepSetDefaults:

  Resets all SPRKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.  Also leaves alone any data
  structures/options related to the ARKODE infrastructure itself
  (e.g., root-finding and post-process step).
  ---------------------------------------------------------------*/
int SPRKStepSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetDefaults", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) return (retval);

  /* Set default ARKODE infrastructure parameters */
  retval = arkSetDefaults(ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::SPRKStep", "SPRKStepSetDefaults",
                    "Error setting ARKODE infrastructure defaults");
    return (retval);
  }

  /* Fixed step mode by default */
  SPRKStepSetFixedStep(ark_mem, 0.01);

  /* set using default method order */
  SPRKStepSetOrder(arkode_mem, 0);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepSetUseCompensatedSums:

  Turns on/off compensated summation in SPRKStep and ARKODE.
  ---------------------------------------------------------------*/
int SPRKStepSetUseCompensatedSums(void* arkode_mem, sunbooleantype onoff)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetOrder", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (onoff)
  {
    arkSetUseCompensatedSums(arkode_mem, SUNTRUE);
    ark_mem->step = sprkStep_TakeStep_Compensated;
  }
  else
  {
    arkSetUseCompensatedSums(arkode_mem, SUNFALSE);
    ark_mem->step = sprkStep_TakeStep;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepSetMethod:

  Specifies the SPRK method

  ** Note in documentation that this should not be called along
  with SPRKStepSetOrder. **
  ---------------------------------------------------------------*/
int SPRKStepSetMethod(void* arkode_mem, ARKodeSPRKMem sprk_mem)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetOrder", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (step_mem->method)
  {
    ARKodeSPRKMem_Free(step_mem->method);
    step_mem->method = NULL;
  }

  step_mem->method = sprk_mem;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepSetOrder:

  Specifies the method order

  ** Note in documentation that this should not be called along
  with SPRKStepSetMethod. **
  ---------------------------------------------------------------*/
int SPRKStepSetOrder(void* arkode_mem, int ord)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetOrder", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) return (retval);

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) { step_mem->q = 4; }
  else { step_mem->q = ord; }

  if (step_mem->method)
  {
    ARKodeSPRKMem_Free(step_mem->method);
    step_mem->method = NULL;
  }

  return (ARK_SUCCESS);
}

/*===============================================================
  SPRKStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  SPRKStepGetNumRhsEvals:

  Returns the current number of calls to f1 and f2
  ---------------------------------------------------------------*/
int SPRKStepGetNumRhsEvals(void* arkode_mem, long int* nf1, long int* nf2)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepGetNumRhsEvals",
                                  &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return (retval);

  *nf1 = step_mem->nf1;
  *nf2 = step_mem->nf2;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int SPRKStepGetTimestepperStats(void* arkode_mem, long int* expsteps,
                                long int* accsteps, long int* step_attempts,
                                long int* nf1, long int* nf2,
                                long int* nlinsetups, long int* netfails)
{
  ARKodeMem ark_mem;

  SPRKStepGetNumExpSteps(arkode_mem, expsteps);
  SPRKStepGetNumAccSteps(arkode_mem, accsteps);
  SPRKStepGetNumRhsEvals(arkode_mem, nf1, nf2);
  SPRKStepGetNumStepAttempts(arkode_mem, step_attempts);
  SPRKStepGetNumErrTestFails(arkode_mem, netfails);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepPrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
int SPRKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepPrintAllStats", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) return (retval);

  /* step and rootfinding stats */
  retval = arkPrintAllStats(arkode_mem, outfile, fmt);
  if (retval != ARK_SUCCESS) return (retval);

  switch (fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    /* function evaluations */
    fprintf(outfile, "f1 RHS fn evals              = %ld\n", step_mem->nf1);
    fprintf(outfile, "f2 RHS fn evals              = %ld\n", step_mem->nf2);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    /* function evaluations */
    fprintf(outfile, ",f1 RHS evals,%ld", step_mem->nf1);
    fprintf(outfile, ",f2 RHS fn evals,%ld", step_mem->nf2);
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "SPRKStepPrintAllStats",
                    "Invalid formatting option.");
    return (ARK_ILL_INPUT);
  }

  return (ARK_SUCCESS);
}

/*===============================================================
  SPRKStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  SPRKStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int SPRKStepWriteParameters(void* arkode_mem, FILE* fp)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int flag, retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepWriteParameters",
                                  &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return (retval);

  /* output ARKODE infrastructure parameters first */
  flag = arkWriteParameters(ark_mem, fp);
  if (flag != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::SPRKStep",
                    "SPRKStepWriteParameters",
                    "Error writing ARKODE infrastructure parameters");
    return (flag);
  }

  /* print integrator parameters to file */
  fprintf(fp, "SPRKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n", step_mem->method->q);
  fprintf(fp, "\n");

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
