/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sundials/sundials_errors.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"
#include "sundials_stepper_impl.h"

SUNErrCode SUNStepper_Create(SUNContext sunctx, SUNStepper* stepper_ptr)
{
  SUNFunctionBegin(sunctx);
  SUNCheck(stepper_ptr, SUN_ERR_ARG_CORRUPT);

  SUNStepper stepper = malloc(sizeof(*stepper));
  SUNAssert(stepper, SUN_ERR_MALLOC_FAIL);

  stepper->content   = NULL;
  stepper->sunctx    = sunctx;
  stepper->last_flag = SUN_SUCCESS;

  stepper->ops = malloc(sizeof(*(stepper->ops)));
  SUNAssert(stepper->ops, SUN_ERR_MALLOC_FAIL);

  stepper->ops->evolve           = NULL;
  stepper->ops->onestep          = NULL;
  stepper->ops->fullrhs          = NULL;
  stepper->ops->reset            = NULL;
  stepper->ops->setstoptime      = NULL;
  stepper->ops->setstepdirection = NULL;
  stepper->ops->setforcing       = NULL;
  stepper->ops->destroy          = NULL;

  *stepper_ptr = stepper;

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_Destroy(SUNStepper* stepper_ptr)
{
  if (stepper_ptr != NULL)
  {
    const SUNStepper_Ops ops = (*stepper_ptr)->ops;
    if (ops && ops->destroy) { ops->destroy(*stepper_ptr); }
    free(ops);
    free(*stepper_ptr);
    *stepper_ptr = NULL;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_Evolve(SUNStepper stepper, sunrealtype tout, N_Vector y,
                             sunrealtype* tret)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->evolve)
  {
    return stepper->ops->evolve(stepper, tout, y, tret);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_OneStep(SUNStepper stepper, sunrealtype tout, N_Vector y,
                              sunrealtype* tret)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->onestep)
  {
    return stepper->ops->onestep(stepper, tout, y, tret);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_FullRhs(SUNStepper stepper, sunrealtype t, N_Vector v,
                              N_Vector f, SUNFullRhsMode mode)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->fullrhs)
  {
    return stepper->ops->fullrhs(stepper, t, v, f, mode);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_Reset(SUNStepper stepper, sunrealtype tR, N_Vector yR)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->reset) { return stepper->ops->reset(stepper, tR, yR); }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_SetStopTime(SUNStepper stepper, sunrealtype tstop)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->setstoptime)
  {
    return stepper->ops->setstoptime(stepper, tstop);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_SetStepDirection(SUNStepper stepper, sunrealtype stepdir)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->setstepdirection)
  {
    return stepper->ops->setstepdirection(stepper, stepdir);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_SetForcing(SUNStepper stepper, sunrealtype tshift,
                                 sunrealtype tscale, N_Vector* forcing,
                                 int nforcing)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->setforcing)
  {
    return stepper->ops->setforcing(stepper, tshift, tscale, forcing, nforcing);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_SetContent(SUNStepper stepper, void* content)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->content = content;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_GetContent(SUNStepper stepper, void** content)
{
  SUNFunctionBegin(stepper->sunctx);
  *content = stepper->content;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetLastFlag(SUNStepper stepper, int last_flag)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->last_flag = last_flag;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_GetLastFlag(SUNStepper stepper, int* last_flag)
{
  SUNFunctionBegin(stepper->sunctx);
  *last_flag = stepper->last_flag;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetEvolveFn(SUNStepper stepper, SUNStepperEvolveFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->evolve = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetOneStepFn(SUNStepper stepper, SUNStepperOneStepFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->onestep = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetFullRhsFn(SUNStepper stepper, SUNStepperFullRhsFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->fullrhs = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetResetFn(SUNStepper stepper, SUNStepperResetFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->reset = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetStopTimeFn(SUNStepper stepper, SUNStepperSetStopTimeFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->setstoptime = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetStepDirectionFn(SUNStepper stepper,
                                         SUNStepperSetStepDirectionFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->setstepdirection = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetForcingFn(SUNStepper stepper, SUNStepperSetForcingFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->setforcing = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetDestroyFn(SUNStepper stepper, SUNStepperDestroyFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->destroy = fn;
  return SUN_SUCCESS;
}
