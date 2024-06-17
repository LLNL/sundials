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

typedef int (*SUNJacFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                        void* user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3);

typedef int (*SUNJacTimesFn)(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                             N_Vector fy, void* user_data, N_Vector tmp);

typedef _SUNDIALS_STRUCT_ SUNStepper_* SUNStepper;

typedef int (*SUNStepperAdvanceFn)(SUNStepper stepper, sunrealtype t0,
                                   sunrealtype tout, N_Vector y,
                                   sunrealtype* tret, int* stop_reason);

typedef int (*SUNStepperOneStepFn)(SUNStepper stepper, sunrealtype t0,
                                   sunrealtype tout, N_Vector y,
                                   sunrealtype* tret, int* stop_reason);

typedef int (*SUNStepperTryStepFn)(SUNStepper stepper, sunrealtype t0,
                                   sunrealtype tout, N_Vector y,
                                   sunrealtype* tret, int* stop_reason);

typedef int (*SUNStepperFullRhsFn)(SUNStepper stepper, sunrealtype t,
                                   N_Vector y, N_Vector f, int mode);

typedef int (*SUNStepperResetFn)(SUNStepper stepper, sunrealtype tR, N_Vector yR);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Create(SUNContext sunctx, SUNStepper* stepper);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Destroy(SUNStepper* stepper);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetContent(SUNStepper stepper, void* content);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_GetContent(SUNStepper stepper, void** content);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetAdvanceFn(SUNStepper stepper, SUNStepperAdvanceFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetOneStepFn(SUNStepper stepper, SUNStepperOneStepFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetTryStepFn(SUNStepper stepper, SUNStepperTryStepFn fn);

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

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Advance(SUNStepper stepper, sunrealtype t0,
                              sunrealtype tout, N_Vector y, sunrealtype* tret,
                              int* stop_reason);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Step(SUNStepper stepper, sunrealtype t0, sunrealtype tout,
                           N_Vector y, sunrealtype* tret, int* stop_reason);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_TryStep(SUNStepper stepper, sunrealtype t0,
                              sunrealtype tout, N_Vector y, sunrealtype* tret,
                              int* stop_reason);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Reset(SUNStepper stepper, sunrealtype tR, N_Vector yR);

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_STEPPER_H */
