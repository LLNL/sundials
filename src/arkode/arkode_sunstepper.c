/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
                  Cody J. Balos @ LLNL
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
 * This is the implementation file for ARKODE's interfacing with
 * SUNStepper
 *--------------------------------------------------------------*/

#include <arkode/arkode.h>
#include <sundials/sundials_stepper.h>
#include "arkode_impl.h"
#include "sundials_macros.h"

static SUNErrCode arkSUNStepperEvolveHelper(SUNStepper stepper,
                                            SUNDIALS_MAYBE_UNUSED sunrealtype t0,
                                            sunrealtype tout, N_Vector y,
                                            sunrealtype* tret, int mode)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));

  /* evolve inner ODE */
  stepper->last_flag = ARKodeEvolve(arkode_mem, tout, y, tret, mode);
  if (stepper->last_flag < 0) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

static SUNErrCode arkSUNStepperEvolve(SUNStepper stepper, sunrealtype t0,
                                      sunrealtype tout, N_Vector y,
                                      sunrealtype* tret)
{
  return arkSUNStepperEvolveHelper(stepper, t0, tout, y, tret, ARK_NORMAL);
}

static SUNErrCode arkSUNStepperOneStep(SUNStepper stepper, sunrealtype t0,
                                       sunrealtype tout, N_Vector y,
                                       sunrealtype* tret)
{
  return arkSUNStepperEvolveHelper(stepper, t0, tout, y, tret, ARK_ONE_STEP);
}

static SUNErrCode arkSUNStepperTryStep(SUNStepper stepper, sunrealtype t0,
                                       sunrealtype tout, N_Vector y,
                                       sunrealtype* tret)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));

  /* try to evolve inner ODE */
  stepper->last_flag = arkTryStep(arkode_mem, t0, tout, y, tret, NULL);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperFullRhsFn to compute the full inner
  (fast) ODE IVP RHS.
  ----------------------------------------------------------------------------*/

static SUNErrCode arkSUNStepperFullRhs(SUNStepper stepper, sunrealtype t,
                                       N_Vector y, N_Vector f, int mode)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;

  stepper->last_flag = ark_mem->step_fullrhs(ark_mem, t, y, f, mode);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperResetFn to reset the inner (fast) stepper
  state.
  ----------------------------------------------------------------------------*/

static SUNErrCode arkSUNStepperReset(SUNStepper stepper, sunrealtype tR,
                                     N_Vector yR)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));

  stepper->last_flag = ARKodeReset(arkode_mem, tR, yR);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperStopTimeFn to set the tstop time
  ----------------------------------------------------------------------------*/

static SUNErrCode arkSUNStepperSetStopTime(SUNStepper stepper, sunrealtype tstop)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));

  stepper->last_flag = ARKodeSetStopTime(arkode_mem, tstop);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

static SUNErrCode arkSUNStepperSetForcing(SUNStepper stepper, sunrealtype tshift,
                                          sunrealtype tscale, N_Vector* forcing,
                                          int nforcing)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;

  stepper->last_flag = ark_mem->step_setforcing(ark_mem, tshift, tscale,
                                                forcing, nforcing);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

int ARKodeCreateSUNStepper(void* arkode_mem, SUNStepper* stepper)
{
  /* unpack ark_mem */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;

  SUNFunctionBegin(ark_mem->sunctx);
  SUNCheckCall(SUNStepper_Create(ark_mem->sunctx, stepper));
  SUNCheckCall(SUNStepper_SetContent(*stepper, arkode_mem));
  SUNCheckCall(SUNStepper_SetEvolveFn(*stepper, arkSUNStepperEvolve));
  SUNCheckCall(SUNStepper_SetOneStepFn(*stepper, arkSUNStepperOneStep));
  SUNCheckCall(SUNStepper_SetTryStepFn(*stepper, arkSUNStepperTryStep));
  SUNCheckCall(SUNStepper_SetFullRhsFn(*stepper, arkSUNStepperFullRhs));
  SUNCheckCall(SUNStepper_SetResetFn(*stepper, arkSUNStepperReset));
  SUNCheckCall(SUNStepper_SetStopTimeFn(*stepper, arkSUNStepperSetStopTime));
  SUNCheckCall(SUNStepper_SetForcingFn(*stepper, arkSUNStepperSetForcing));

  return ARK_SUCCESS;
}
