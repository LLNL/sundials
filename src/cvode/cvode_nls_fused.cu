 /*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file implements fused CUDA kernels for cvode_nls.c.
 * -----------------------------------------------------------------
 */

#include <cuda_runtime.h>

#include "cvode_impl.h"
#include <nvector/nvector_cuda.h>
#include "sundials_cuda_kernels.cuh"


/*
 * -----------------------------------------------------------------
 * Compute the nonlinear residual.
 * -----------------------------------------------------------------
 */


__global__
void cvNlsResid_cukernel(const sunindextype length,
                         const realtype rl1,
                         const realtype ngamma,
                         const realtype* zn1,
                         const realtype* ycor,
                         const realtype* ftemp,
                         realtype* res)
{
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VLinearSum(cv_mem->cv_rl1, cv_mem->cv_zn[1], ONE, ycor, res);
    // N_VLinearSum(-cv_mem->cv_gamma, cv_mem->cv_ftemp, ONE, res, res);
    realtype tmp = rl1*zn1[i] + ycor[i];
    res[i] = ngamma*ftemp[i] + tmp;
  }
}

extern "C"
int cvNlsResid_fused(const realtype rl1,
                     const realtype ngamma,
                     const N_Vector zn1,
                     const N_Vector ycor,
                     const N_Vector ftemp,
                     N_Vector res)
{
  const SUNCudaExecPolicy* exec_policy = ((N_VectorContent_Cuda)res->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(res);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvNlsResid_cukernel<<<grid, block, 0, exec_policy->stream()>>>
  (
    N,
    rl1,
    ngamma,
    N_VGetDeviceArrayPointer_Cuda(zn1),
    N_VGetDeviceArrayPointer_Cuda(ycor),
    N_VGetDeviceArrayPointer_Cuda(ftemp),
    N_VGetDeviceArrayPointer_Cuda(res)
  );

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return -1;
#endif

  return 0;
}
