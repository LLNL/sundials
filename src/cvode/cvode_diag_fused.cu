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
 * This file implements fused CUDA kernels for CVDiag.
 * -----------------------------------------------------------------
 */

#include <cuda_runtime.h>

#include <nvector/nvector_cuda.h>
#include "sundials_cuda_kernels.cuh"


/* 
 * -----------------------------------------------------------------
 * Form y with perturbation = FRACT*(func. iter. correction)
 * -----------------------------------------------------------------
 */

__global__
void cvDiagSetup_formY_kernel(const sunindextype length,
                              const realtype h,
                              const realtype r,
                              const realtype* fpred,
                              const realtype* zn1,
                              const realtype* ypred,
                              realtype* ftemp,
                              realtype* y)
{
  // N_VLinearSum(h, fpred, -ONE, zn[1], ftemp);
  // N_VLinearSum(r, ftemp, ONE, ypred, y);
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    ftemp[i] = h*fpred[i] - zn1[i];
    y[i] = r*ftemp[i] + ypred[i];
  }
}

extern "C"
int cvDiagSetup_formY(const realtype h,
                      const realtype r,
                      const N_Vector fpred,
                      const N_Vector zn1,
                      const N_Vector ypred,
                      N_Vector ftemp,
                      N_Vector y)
{
  const SUNCudaExecPolicy* exec_policy = ((N_VectorContent_Cuda)y->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(y);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvDiagSetup_formY_kernel<<<grid, block, 0, exec_policy->stream()>>>
  (
    N,
    h,
    r,
    N_VGetDeviceArrayPointer_Cuda(fpred),
    N_VGetDeviceArrayPointer_Cuda(zn1),
    N_VGetDeviceArrayPointer_Cuda(ypred),
    N_VGetDeviceArrayPointer_Cuda(ftemp),
    N_VGetDeviceArrayPointer_Cuda(y)
  );

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return -1;
#endif

  return 0;
}

/* 
 * -----------------------------------------------------------------
 * Construct M = I - gamma*J with J = diag(deltaf_i/deltay_i)
 * protecting against deltay_i being at roundoff level.
 * -----------------------------------------------------------------
 */

__global__
void cvDiagSetup_buildM_kernel(const sunindextype length,
                               const realtype fract,
                               const realtype uround,
                               const realtype h,
                               const realtype* ftemp,
                               const realtype* fpred,
                               const realtype* ewt,
                               realtype* bit,
                               realtype* bitcomp,
                               realtype* y,
                               realtype* M)
{
  static const realtype zero = 0.0;
  static const realtype one = 1.0;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VLinearSum(ONE, M, -ONE, fpred, M);
    // N_VLinearSum(FRACT, ftemp, -h, M, M);
    // N_VProd(ftemp, ewt, y);
    M[i] = fract*ftemp[i] - h*(M[i] - fpred[i]);
    y[i] = ftemp[i] * ewt[i];
        
    // N_VCompare(uround, y, bit);
    // N_VAddConst(bit, -ONE, bitcomp);
    bool test = (abs(y[i]) > uround);
    bit[i] = test ? one : zero;
    bitcomp[i] = test ? zero : -one;

    // N_VProd(ftemp, bit, y);
    // N_VLinearSum(FRACT, y, -ONE, bitcomp, y);
    // N_VDiv(M, y, M);
    // N_VProd(M, bit, M);
    // N_VLinearSum(ONE, M, -ONE, bitcomp, M);
    y[i] = fract*ftemp[i]*bit[i] - bitcomp[i];
    M[i] = M[i]/y[i] * bit[i] - bitcomp[i];
  }
}

extern "C"
int cvDiagSetup_buildM(const realtype fract,
                       const realtype uround,
                       const realtype h,
                       const N_Vector ftemp,
                       const N_Vector fpred,
                       const N_Vector ewt,
                       N_Vector bit,
                       N_Vector bitcomp,
                       N_Vector y,
                       N_Vector M)
{
  const SUNCudaExecPolicy* exec_policy = ((N_VectorContent_Cuda)M->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(M);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvDiagSetup_buildM_kernel<<<grid, block, 0, exec_policy->stream()>>>
  (
    N,
    fract,
    uround,
    h,
    N_VGetDeviceArrayPointer_Cuda(ftemp),
    N_VGetDeviceArrayPointer_Cuda(fpred),
    N_VGetDeviceArrayPointer_Cuda(ewt),
    N_VGetDeviceArrayPointer_Cuda(bit),
    N_VGetDeviceArrayPointer_Cuda(bitcomp),
    N_VGetDeviceArrayPointer_Cuda(y),
    N_VGetDeviceArrayPointer_Cuda(M)
  );

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return -1;
#endif

  return 0;
}


/*
 * -----------------------------------------------------------------
 *  Update M with changed gamma so that M = I - gamma*J.
 * -----------------------------------------------------------------
 */

 __global__
void cvDiagSolve_updateM_kernel(const sunindextype length, const realtype r, realtype* M)
{
  static const realtype one = 1.0;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VInv(M, M);
    // N_VAddConst(M, -ONE, M);
    // N_VScale(r, M, M);
    // N_VAddConst(M, ONE, M);
    realtype a = one/M[i] - one;
    M[i] = r*a + one;
  }
}


extern "C"
int cvDiagSolve_updateM(const realtype r, N_Vector M)
{
  const SUNCudaExecPolicy* exec_policy = ((N_VectorContent_Cuda)M->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(M);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvDiagSolve_updateM_kernel<<<grid, block, 0, exec_policy->stream()>>>
  (
    N,
    r,
    N_VGetDeviceArrayPointer_Cuda(M)
  );

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return -1;
#endif

  return 0;
}