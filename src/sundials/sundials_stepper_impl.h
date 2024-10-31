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

#ifndef _SUNDIALS_STEPPER_IMPL_H
#define _SUNDIALS_STEPPER_IMPL_H

#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SUNStepper_Ops_* SUNStepper_Ops;

struct SUNStepper_Ops_
{
  SUNStepperEvolveFn evolve;
  SUNStepperOneStepFn onestep;
  SUNStepperFullRhsFn fullrhs;
  SUNStepperResetFn reset;
  SUNStepperSetStopTimeFn setstoptime;
  SUNStepperSetStepDirectionFn setstepdirection;
  SUNStepperSetForcingFn setforcing;
  SUNStepperDestroyFn destroy;
};

struct SUNStepper_
{
  /* stepper specific content and operations */
  void* content;
  SUNStepper_Ops ops;

  /* stepper context */
  SUNContext sunctx;

  /* last stepper return flag */
  int last_flag;
};

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_STEPPER_IMPL_H */
