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

/*------------------------------------------------------------------------------
  Implementation of SUNStepperFullRhsFn to compute the full inner
  (fast) ODE IVP RHS.
  ----------------------------------------------------------------------------*/

static SUNErrCode arkSUNStepperFullRhs(SUNStepper stepper, sunrealtype t,
                                       N_Vector y, N_Vector f)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;

  stepper->last_flag = ark_mem->step_fullrhs(ark_mem, t, y, f, ARK_FULLRHS_OTHER);
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

  SUNErrCode err = SUNStepper_Create(ark_mem->sunctx, stepper);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }

  err = SUNStepper_SetContent(*stepper, arkode_mem);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }

  err = SUNStepper_SetEvolveFn(*stepper, arkSUNStepperEvolve);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }

  err = SUNStepper_SetOneStepFn(*stepper, arkSUNStepperOneStep);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }

  err = SUNStepper_SetFullRhsFn(*stepper, arkSUNStepperFullRhs);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }

  err = SUNStepper_SetResetFn(*stepper, arkSUNStepperReset);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }

  err = SUNStepper_SetStopTimeFn(*stepper, arkSUNStepperSetStopTime);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }

  if (ark_mem->step_setforcing != NULL)
  {
    err = SUNStepper_SetForcingFn(*stepper, arkSUNStepperSetForcing);
    if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }
  }

  return ARK_SUCCESS;
}
