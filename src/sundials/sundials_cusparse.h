/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
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
#include <cusparse.h>
#include <stdio.h>
#include <sundials/sundials_types.h>

#ifndef _SUNDIALS_CUSPARSE_H
#define _SUNDIALS_CUSPARSE_H

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * Utility macros
 * ---------------------------------------------------------------------------*/

#define SUNDIALS_CUSPARSE_VERIFY(cuerr) \
  SUNDIALS_CUSPARSE_Assert(cuerr, __FILE__, __LINE__)

/* ---------------------------------------------------------------------------
 * Utility functions
 * ---------------------------------------------------------------------------*/

inline sunbooleantype SUNDIALS_CUSPARSE_Assert(cusparseStatus_t status,
                                               const char* file, int line)
{
  if (status != CUSPARSE_STATUS_SUCCESS)
  {
#ifdef SUNDIALS_DEBUG
    fprintf(stderr, "ERROR in cuSPARSE runtime operation: cusparseStatus_t = %d %s:%d\n",
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

#endif /* _SUNDIALS_CUSPARSE_H */
