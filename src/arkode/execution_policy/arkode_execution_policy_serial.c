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
 * This is the implementation file for the serial execution
 * policy.
 *--------------------------------------------------------------*/

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <arkode/execution_policy/arkode_execution_policy_serial.h>

#include "sundials_macros.h"

/*---------------------------------------------------------------
  This routine does the setup for serial execution for which
  nothing is needed.
  ---------------------------------------------------------------*/
static int setup_serial(
  SUNDIALS_MAYBE_UNUSED const ARKodeSplittingExecutionPolicy policy,
  SUNDIALS_MAYBE_UNUSED const N_Vector y,
  SUNDIALS_MAYBE_UNUSED const int sequential_methods)
{
  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  This routine evaluates the policy function fn in serial over
  the range of sequential methods
  ---------------------------------------------------------------*/
static int execute_serial(
  SUNDIALS_MAYBE_UNUSED const ARKodeSplittingExecutionPolicy policy,
  const ARKExecutionPolicyFn fn, const N_Vector yn, const N_Vector ycur,
  const N_Vector tmp, sunrealtype* const alpha, const int sequential_methods,
  void* const user_data)
{
  N_VScale(1, yn, ycur);
  int retval = fn(0, ycur, user_data);
  if (retval != ARK_SUCCESS) { return retval; }

  if (alpha[0] != SUN_RCONST(1.0)) { N_VScale(alpha[0], ycur, ycur); }

  for (int i = 1; i < sequential_methods; i++)
  {
    N_VScale(SUN_RCONST(1.0), yn, tmp);
    retval = fn(i, tmp, user_data);
    if (retval != ARK_SUCCESS) { return retval; }
    N_VLinearSum(SUN_RCONST(1.0), ycur, alpha[i], tmp, ycur);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  This routine does the freeing for serial execution which happens
  to be nothing
  ---------------------------------------------------------------*/
static void free_serial(
  SUNDIALS_MAYBE_UNUSED const ARKodeSplittingExecutionPolicy policy)
{
  // Nothing needed
}

/*---------------------------------------------------------------
  This routine creates a serial execution policy
  ---------------------------------------------------------------*/
ARKodeSplittingExecutionPolicy ARKodeSplittingExecutionPolicy_New_Serial()
{
  ARKodeSplittingExecutionPolicy policy = malloc(sizeof(*policy));
  if (policy == NULL) { return NULL; }
  policy->setup   = setup_serial;
  policy->execute = execute_serial;
  policy->free    = free_serial;
  policy->data    = NULL;
  return policy;
}
