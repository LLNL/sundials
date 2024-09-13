/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
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
  * This is the header file for the ARKODE SplittingStep module.
 *--------------------------------------------------------------*/

#ifndef ARKODE_SPLITTINGSTEP_H_
#define ARKODE_SPLITTINGSTEP_H_

#include <arkode/arkode_execution_policy.h>
#include <arkode/arkode_splittingstep_coefficients.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_stepper.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* TODO(SBR): would (t0, y0, steppers, partitions, sunctx) be a better arg
   order? That would be more consistent with MRIStep but less with others */
SUNDIALS_EXPORT void* SplittingStepCreate(SUNStepper* steppers, int partitions,
                                          sunrealtype t0, N_Vector y0,
                                          SUNContext sunctx);

SUNDIALS_EXPORT int SplittingStepCreateForcing(SUNStepper stepper1,
                                               SUNStepper stepper2,
                                               sunrealtype t0, N_Vector y0,
                                               SUNContext sunctx);

SUNDIALS_EXPORT int SplittingStep_SetCoefficients(
  void* arkode_mem, SplittingStepCoefficients coefficients);

SUNDIALS_EXPORT int SplittingStep_SetExecutionPolicy(
  void* arkode_mem, ARKodeSplittingExecutionPolicy policy);

SUNDIALS_EXPORT int SplittingStep_GetNumEvolves(void* arkode_mem, int partition,
                                                long int* evolves);

#ifdef __cplusplus
}
#endif

#endif
