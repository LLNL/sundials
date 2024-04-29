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

#ifndef _SUNDIALS_STEPPER_H
#define _SUNDIALS_STEPPER_H

#include <sundials/sundials_core.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef _SUNDIALS_STRUCT_ SUNStepper_s* SUNStepper;

typedef int (*SUNStepperEvolveFn)(SUNStepper stepper, sunrealtype t0,
                                  sunrealtype tout, N_Vector y);

typedef int (*SUNStepperFullRhsFn)(SUNStepper stepper, sunrealtype t,
                                   N_Vector y, N_Vector f, int mode);

typedef int (*SUNStepperResetFn)(SUNStepper stepper, sunrealtype tR, N_Vector yR);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Create(SUNContext sunctx, SUNStepper* stepper);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Free(SUNStepper* stepper);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetContent(SUNStepper stepper, void* content);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_GetContent(SUNStepper stepper, void** content);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetEvolveFn(SUNStepper stepper, SUNStepperEvolveFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetFullRhsFn(SUNStepper stepper, SUNStepperFullRhsFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetResetFn(SUNStepper stepper, SUNStepperResetFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_AddForcing(SUNStepper stepper, sunrealtype t, N_Vector f);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_GetForcingData(SUNStepper stepper, sunrealtype* tshift,
                                     sunrealtype* tscale, N_Vector** forcing,
                                     int* nforcing);

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_STEPPER_H */
