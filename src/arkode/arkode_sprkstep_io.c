/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
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
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  SPRKStepSetUseCompensatedSums:

  Turns on/off compensated summation in SPRKStep and ARKODE.
  ---------------------------------------------------------------*/
int SPRKStepSetUseCompensatedSums(void* arkode_mem, sunbooleantype onoff)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeMem and ARKodeSPRKStepMem structures */
  retval = sprkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (onoff)
  {
    ark_mem->use_compensated_sums = SUNTRUE;
    ark_mem->step                 = sprkStep_TakeStep_Compensated;
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
    ark_mem->use_compensated_sums = SUNFALSE;
    ark_mem->step                 = sprkStep_TakeStep;
  }

  return (retval);
}

/*---------------------------------------------------------------
  SPRKStepSetMethod:

  Specifies the SPRK method

  ** Note in documentation that this should not be called along
  with ARKodeSetOrder. **
  ---------------------------------------------------------------*/
int SPRKStepSetMethod(void* arkode_mem, ARKodeSPRKTable sprk_storage)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeMem and ARKodeSPRKStepMem structures */
  retval = sprkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
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

  /* access ARKodeMem and ARKodeSPRKStepMem structures */
  retval = sprkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
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

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  SPRKStepGetCurrentMethod:

  Returns the stepper method structure.
  ---------------------------------------------------------------*/
int SPRKStepGetCurrentMethod(void* arkode_mem, ARKodeSPRKTable* sprk_storage)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeMem and ARKodeSPRKStepMem structures */
  retval = sprkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *sprk_storage = step_mem->method;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepGetNumRhsEvals:

  Returns the current number of calls to f1 and f2
  ---------------------------------------------------------------*/
int SPRKStepGetNumRhsEvals(void* arkode_mem, long int* nf1, long int* nf2)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeMem and ARKodeSPRKStepMem structures */
  retval = sprkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nf1 = step_mem->nf1;
  *nf2 = step_mem->nf2;

  return (ARK_SUCCESS);
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  sprkStep_SetDefaults:

  Resets all SPRKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.  Also leaves alone any data
  structures/options related to the ARKODE infrastructure itself
  (e.g., root-finding and post-process step).
  ---------------------------------------------------------------*/
int sprkStep_SetDefaults(ARKodeMem ark_mem)
{
  /* use the default method order */
  return (sprkStep_SetOrder(ark_mem, 0));
}

/*---------------------------------------------------------------
  sprkStep_SetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
int sprkStep_SetOrder(ARKodeMem ark_mem, int ord)
{
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(ark_mem, __func__, &step_mem);
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

/*---------------------------------------------------------------
  sprkStep_PrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
int sprkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt)
{
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(ark_mem, __func__, &step_mem);
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

/*---------------------------------------------------------------
  sprkStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int sprkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp)
{
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* print integrator parameters to file */
  fprintf(fp, "SPRKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n", step_mem->method->q);
  fprintf(fp, "  Method stages %i\n", step_mem->method->stages);

  return (ARK_SUCCESS);
}

/*===============================================================
  Exported-but-deprecated user-callable functions.
  ===============================================================*/

int SPRKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)
{
  return (ARKodeReset(arkode_mem, tR, yR));
}

int SPRKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)
{
  return (ARKodeRootInit(arkode_mem, nrtfn, g));
}

int SPRKStepSetRootDirection(void* arkode_mem, int* rootdir)
{
  return (ARKodeSetRootDirection(arkode_mem, rootdir));
}

int SPRKStepSetNoInactiveRootWarn(void* arkode_mem)
{
  return (ARKodeSetNoInactiveRootWarn(arkode_mem));
}

int SPRKStepSetDefaults(void* arkode_mem)
{
  return (ARKodeSetDefaults(arkode_mem));
}

int SPRKStepSetOrder(void* arkode_mem, int ord)
{
  return (ARKodeSetOrder(arkode_mem, ord));
}

int SPRKStepSetInterpolantType(void* arkode_mem, int itype)
{
  return (ARKodeSetInterpolantType(arkode_mem, itype));
}

int SPRKStepSetInterpolantDegree(void* arkode_mem, int degree)
{
  return (ARKodeSetInterpolantDegree(arkode_mem, degree));
}

int SPRKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)
{
  return (ARKodeSetMaxNumSteps(arkode_mem, mxsteps));
}

int SPRKStepSetStopTime(void* arkode_mem, sunrealtype tstop)
{
  return (ARKodeSetStopTime(arkode_mem, tstop));
}

int SPRKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed)
{
  return (ARKodeSetFixedStep(arkode_mem, hfixed));
}

int SPRKStepSetUserData(void* arkode_mem, void* user_data)
{
  return (ARKodeSetUserData(arkode_mem, user_data));
}

int SPRKStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep)
{
  return (ARKodeSetPostprocessStepFn(arkode_mem, ProcessStep));
}

int SPRKStepSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage)
{
  return (ARKodeSetPostprocessStageFn(arkode_mem, ProcessStage));
}

int SPRKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                   sunrealtype* tret, int itask)
{
  return (ARKodeEvolve(arkode_mem, tout, yout, tret, itask));
}

int SPRKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)
{
  return (ARKodeGetDky(arkode_mem, t, k, dky));
}

char* SPRKStepGetReturnFlagName(long int flag)
{
  return (ARKodeGetReturnFlagName(flag));
}

int SPRKStepGetCurrentState(void* arkode_mem, N_Vector* state)
{
  return (ARKodeGetCurrentState(arkode_mem, state));
}

int SPRKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur)
{
  return (ARKodeGetCurrentStep(arkode_mem, hcur));
}

int SPRKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur)
{
  return (ARKodeGetCurrentTime(arkode_mem, tcur));
}

int SPRKStepGetLastStep(void* arkode_mem, sunrealtype* hlast)
{
  return (ARKodeGetLastStep(arkode_mem, hlast));
}

int SPRKStepGetNumStepAttempts(void* arkode_mem, long int* nstep_attempts)
{
  return (ARKodeGetNumStepAttempts(arkode_mem, nstep_attempts));
}

int SPRKStepGetNumSteps(void* arkode_mem, long int* nsteps)
{
  return (ARKodeGetNumSteps(arkode_mem, nsteps));
}

int SPRKStepGetRootInfo(void* arkode_mem, int* rootsfound)
{
  return (ARKodeGetRootInfo(arkode_mem, rootsfound));
}

int SPRKStepGetUserData(void* arkode_mem, void** user_data)
{
  return (ARKodeGetUserData(arkode_mem, user_data));
}

int SPRKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  return (ARKodePrintAllStats(arkode_mem, outfile, fmt));
}

void SPRKStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodePrintMem(arkode_mem, outfile);
}

int SPRKStepSetMaxNumConstrFails(void* arkode_mem, int maxfails)
{
  return (ARKodeSetMaxNumConstrFails(arkode_mem, maxfails));
}

int SPRKStepWriteParameters(void* arkode_mem, FILE* fp)
{
  return (ARKodeWriteParameters(arkode_mem, fp));
}

int SPRKStepGetStepStats(void* arkode_mem, long int* nsteps,
                         sunrealtype* hinused, sunrealtype* hlast,
                         sunrealtype* hcur, sunrealtype* tcur)
{
  return (ARKodeGetStepStats(arkode_mem, nsteps, hinused, hlast, hcur, tcur));
}

int SPRKStepGetNumConstrFails(void* arkode_mem, long int* nconstrfails)
{
  return (ARKodeGetNumConstrFails(arkode_mem, nconstrfails));
}

void SPRKStepFree(void** arkode_mem) { ARKodeFree(arkode_mem); }

/*===============================================================
  EOF
  ===============================================================*/
