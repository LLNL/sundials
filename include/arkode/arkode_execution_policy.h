/* -----------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This header defines execution policy data structures and functions.
 * ---------------------------------------------------------------------------*/

#ifndef ARKODE_EXECUTION_POLICY_H_
#define ARKODE_EXECUTION_POLICY_H_

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Parallelization Policy */
typedef int (*ARKExecutionPolicyFn)(int i, N_Vector y, void* user_data);

typedef _SUNDIALS_STRUCT_ ARKodeSplittingExecutionPolicyMem* ARKodeSplittingExecutionPolicy;

struct ARKodeSplittingExecutionPolicyMem
{
  int (*setup)(ARKodeSplittingExecutionPolicy policy, N_Vector y,
               int sequential_methods);

  int (*execute)(ARKodeSplittingExecutionPolicy policy, ARKExecutionPolicyFn fn,
                 N_Vector yn, N_Vector ycur, N_Vector tmp, sunrealtype* alpha,
                 int sequential_methods, void* user_data);

  void (*free)(ARKodeSplittingExecutionPolicy policy);

  void* data;
};

// TODO: Return error code? Accept pointer?
SUNDIALS_EXPORT void ARKodeSplittingExecutionPolicy_Free(
  ARKodeSplittingExecutionPolicy* policy);

#ifdef __cplusplus
}
#endif

#endif
