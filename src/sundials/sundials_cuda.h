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
#include <stdio.h>
#include <sundials/sundials_types.h>

#ifndef _SUNDIALS_CUDA_H
#define _SUNDIALS_CUDA_H

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * Utility macros
 * ---------------------------------------------------------------------------*/

#define SUNDIALS_CUDA_VERIFY(cuerr) \
  SUNDIALS_CUDA_Assert(cuerr, __FILE__, __LINE__)

#define SUNDIALS_KERNEL_NAME(...) __VA_ARGS__
#ifndef SUNDIALS_DEBUG_CUDA_LASTERROR
#define SUNDIALS_LAUNCH_KERNEL(kernel, gridDim, blockDim, shMem, stream, ...) \
  {                                                                           \
    kernel<<<gridDim, blockDim, shMem, stream>>>(__VA_ARGS__);                \
  }
#else
#define SUNDIALS_LAUNCH_KERNEL(kernel, gridDim, blockDim, shMem, stream, ...) \
  {                                                                           \
    kernel<<<gridDim, blockDim, shMem, stream>>>(__VA_ARGS__);                \
    cudaDeviceSynchronize();                                                  \
    SUNDIALS_CUDA_VERIFY(cudaGetLastError());                                 \
  }
#endif

/* ---------------------------------------------------------------------------
 * Utility functions
 * ---------------------------------------------------------------------------*/

inline sunbooleantype SUNDIALS_CUDA_Assert(cudaError_t cuerr, const char* file,
                                           int line)
{
  if (cuerr != cudaSuccess)
  {
#ifdef SUNDIALS_DEBUG
    fprintf(stderr, "ERROR in CUDA runtime operation: %s %s:%d\n",
            cudaGetErrorString(cuerr), file, line);
#ifdef SUNDIALS_DEBUG_ASSERT
    assert(false);
#endif
#endif
    return SUNFALSE; /* Assert failed */
  }
  return SUNTRUE; /* Assert OK */
}

#ifdef __cplusplus /* wrapper to enable C++ usage */
}
#endif

#endif /* _SUNDIALS_CUDA_H */
