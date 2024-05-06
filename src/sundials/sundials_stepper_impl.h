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

typedef struct SUNStepper_Ops_s* SUNStepper_Ops;

struct SUNStepper_Ops_s
{
  SUNStepperAdvanceFn advance;
  SUNStepperOneStepFn onestep;
  SUNStepperTryStepFn trystep;
  SUNStepperFullRhsFn fullrhs;
  SUNStepperResetFn reset;
};

struct SUNStepper_s
{
  /* stepper specific content and operations */
  void* content;
  SUNStepper_Ops ops;

  /* stepper context */
  SUNContext sunctx;

  /* base class data */
  N_Vector* forcing;      /* array of forcing vectors            */
  int nforcing;           /* number of forcing vectors active    */
  int nforcing_allocated; /* number of forcing vectors allocated */
  int last_flag;          /* last stepper return flag            */
  sunrealtype tshift;     /* time normalization shift            */
  sunrealtype tscale;     /* time normalization scaling          */

  /* fused op workspace */
  sunrealtype* fused_scalars;
  N_Vector* fused_vectors;
};

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_STEPPER_IMPL_H */
