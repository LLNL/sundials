/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This header files defines internal utility functions and macros
 * for working with CUDA.
 * -----------------------------------------------------------------
 */

#include <assert.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cusolverSp.h>
#include <stdio.h>
#include <sundials/sundials_types.h>

#ifndef _SUNDIALS_CUSOLVER_H
#define _SUNDIALS_CUSOLVER_H

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * Utility macros
 * ---------------------------------------------------------------------------*/

#define SUNDIALS_CUSOLVER_VERIFY(cuerr) \
  SUNDIALS_CUSOLVER_Assert(cuerr, __FILE__, __LINE__)

/* ---------------------------------------------------------------------------
 * Utility functions
 * ---------------------------------------------------------------------------*/

inline sunbooleantype SUNDIALS_CUSOLVER_Assert(cusolverStatus_t status,
                                               const char* file, int line)
{
  if (status != CUSOLVER_STATUS_SUCCESS)
  {
#ifdef SUNDIALS_DEBUG
    fprintf(stderr, "ERROR in cuSOLVER runtime operation: cusolverStatus_t = %d %s:%d\n",
            status, file, line);
#ifdef SUNDIALS_DEBUG_ASSERT
    assert(false);
#endif
#endif
    return SUNFALSE; /*  Assert failed */
  }
  return SUNTRUE; /* Assert OK */
}

#ifdef __cplusplus /* wrapper to enable C++ usage */
}
#endif

#endif /* _SUNDIALS_CUSOLVER_H */
