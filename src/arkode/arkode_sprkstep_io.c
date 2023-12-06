/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
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

#include "arkode/arkode.h"
#include "arkode/arkode_sprk.h"
#include "arkode_sprkstep_impl.h"

/*===============================================================
  SPRKStep Optional input functions (wrappers for generic ARKODE
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int SPRKStepSetInterpolantDegree(void* arkode_mem, int degree)
{
  if (degree < 0) { degree = ARK_INTERP_MAX_DEGREE; }
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

int SPRKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)
{
  return (arkSetMaxNumSteps(arkode_mem, mxsteps));
}

int SPRKStepSetStopTime(void* arkode_mem, sunrealtype tstop)
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

int SPRKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed)
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

int SPRKStepGetLastStep(void* arkode_mem, sunrealtype* hlast)
{
  return (arkGetLastStep(arkode_mem, hlast));
}

int SPRKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur)
{
  return (arkGetCurrentStep(arkode_mem, hcur));
}

int SPRKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur)
{
  return (arkGetCurrentTime(arkode_mem, tcur));
}

int SPRKStepGetCurrentState(void* arkode_mem, N_Vector* state)
{
  return (arkGetCurrentState(arkode_mem, state));
}

int SPRKStepGetRootInfo(void* arkode_mem, int* rootsfound)
{
  return (arkGetRootInfo(arkode_mem, rootsfound));
}

int SPRKStepGetStepStats(void* arkode_mem, long int* nsteps,
                         sunrealtype* hinused, sunrealtype* hlast,
                         sunrealtype* hcur, sunrealtype* tcur)
{
  return (arkGetStepStats(arkode_mem, nsteps, hinused, hlast, hcur, tcur));
}

int SPRKStepGetNumConstrFails(void* arkode_mem, long int* nconstrfails)
{
  return (arkGetNumConstrFails(arkode_mem, nconstrfails));
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
  return (arkSetUserData(arkode_mem, user_data));
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
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetDefaults", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Set default ARKODE infrastructure parameters */
  retval = arkSetDefaults(ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Error setting ARKODE infrastructure defaults");
    return (retval);
  }

  /* use the default method order */
  SPRKStepSetOrder(arkode_mem, 0);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepSetUseCompensatedSums:

  Turns on/off compensated summation in SPRKStep and ARKODE.
  ---------------------------------------------------------------*/
int SPRKStepSetUseCompensatedSums(void* arkode_mem, sunbooleantype onoff)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetUseCompensatedSums",
                                  &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (onoff)
  {
    arkSetUseCompensatedSums(arkode_mem, SUNTRUE);
    ark_mem->step = sprkStep_TakeStep_Compensated;
    if (!step_mem->yerr)
    {
      if (!arkAllocVec(ark_mem, ark_mem->yn, &(step_mem->yerr)))
      {
        return ARK_MEM_FAIL;
      }
    }
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
int SPRKStepSetMethod(void* arkode_mem, ARKodeSPRKTable sprk_storage)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetMethod", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (step_mem->method)
  {
    ARKodeSPRKTable_Free(step_mem->method);
    step_mem->method = NULL;
  }

  step_mem->method = ARKodeSPRKTable_Copy(sprk_storage);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepSetMethodName:

  Specifies the SPRK method.
  ---------------------------------------------------------------*/
int SPRKStepSetMethodName(void* arkode_mem, const char* method)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetMethodName", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (step_mem->method)
  {
    ARKodeSPRKTable_Free(step_mem->method);
    step_mem->method = NULL;
  }

  step_mem->method = ARKodeSPRKTable_LoadByName(method);

  return step_mem->method ? ARK_SUCCESS : ARK_ILL_INPUT;
}

/*---------------------------------------------------------------
  SPRKStepSetOrder:

  Specifies the method order

  ** Note in documentation that this should not be called along
  with SPRKStepSetMethod. **
  ---------------------------------------------------------------*/
int SPRKStepSetOrder(void* arkode_mem, int ord)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepSetOrder", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Invalid orders result in the default order being used. */
  if (ord == 7 || ord == 9 || ord > 10) { ord = -1; }

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) { step_mem->q = 4; }
  else { step_mem->q = ord; }

  if (step_mem->method)
  {
    ARKodeSPRKTable_Free(step_mem->method);
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
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepGetNumRhsEvals",
                                  &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nf1 = step_mem->nf1;
  *nf2 = step_mem->nf2;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepGetCurrentMethod:

  Returns the stepper method structure.
  ---------------------------------------------------------------*/
int SPRKStepGetCurrentMethod(void* arkode_mem, ARKodeSPRKTable* sprk_storage)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepGetNumRhsEvals",
                                  &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *sprk_storage = step_mem->method;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepPrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
int SPRKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepPrintAllStats", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* step and rootfinding stats */
  retval = arkPrintAllStats(arkode_mem, outfile, fmt);
  if (retval != ARK_SUCCESS) { return (retval); }

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
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
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
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int flag                   = 0;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepWriteParameters",
                                  &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* output ARKODE infrastructure parameters first */
  flag = arkWriteParameters(ark_mem, fp);
  if (flag != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Error writing ARKODE infrastructure parameters");
    return (flag);
  }

  /* print integrator parameters to file */
  fprintf(fp, "SPRKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n", step_mem->method->q);
  fprintf(fp, "  Method stages %i\n", step_mem->method->stages);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
