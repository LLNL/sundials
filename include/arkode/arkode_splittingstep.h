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
 * TODO
 *--------------------------------------------------------------*/

#ifndef ARKODE_SPLITTINGSTEP_H_
#define ARKODE_SPLITTINGSTEP_H_

#include <arkode/arkode_execution_policy.h>
#include <arkode/arkode_splitting_coefficients.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* TODO: remove once SUNStepper is ready */
typedef MRIStepInnerStepper SUNStepper;

/* ARKODE functions */
/* TODO: can we use `const SUNStepper* steppers` since we don't make copies? */
/* TODO: would (t0, y0, steppers, partitions, sunctx) be a better arg order? Seems slightly more consistent with MRIStepCreate */
SUNDIALS_EXPORT void* SplittingStepCreate(SUNStepper* steppers, int partitions,
                                          sunrealtype t0, N_Vector y0,
                                          SUNContext sunctx);

SUNDIALS_EXPORT int SplittingStep_SetCoefficients(
  void* arkode_mem, ARKodeSplittingCoefficients coefficients);

SUNDIALS_EXPORT int SplittingStep_SetExecutionPolicy(
  void* arkode_mem, ARKodeSplittingExecutionPolicy policy);

#ifdef __cplusplus
}
#endif

#endif
