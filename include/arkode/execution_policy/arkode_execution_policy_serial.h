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
 * The header file for the serial execution policy
 * ---------------------------------------------------------------------------*/

#ifndef ARKODE_EXECUTION_POLICY_SERIAL_H_
#define ARKODE_EXECUTION_POLICY_SERIAL_H_

#include <arkode/arkode_execution_policy.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT ARKodeSplittingExecutionPolicy
ARKodeSplittingExecutionPolicy_New_Serial();

#ifdef __cplusplus
}
#endif

#endif
