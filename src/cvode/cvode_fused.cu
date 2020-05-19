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
 * This file implements fused CUDA kernels for CVODE.
 * -----------------------------------------------------------------
 */

#include <cuda_runtime.h>

#include "cvode_impl.h"
#include <nvector/nvector_cuda.h>
#include "sundials_cuda_kernels.cuh"


/*
 * -----------------------------------------------------------------
 * Compute the ewt vector when the tol type is CV_SS.
 * -----------------------------------------------------------------
 */


__global__
void cvEwtSetSS_cukernel(const sunindextype length,
                         const realtype reltol,
                         const realtype Sabstol,
                         const realtype* ycur,
                         realtype* tempv,
                         realtype* weight)
{
  const realtype one = 1.0;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VAbs(ycur, cv_mem->cv_tempv);
    // N_VScale(cv_mem->cv_reltol, cv_mem->cv_tempv, cv_mem->cv_tempv);
    // N_VAddConst(cv_mem->cv_tempv, cv_mem->cv_Sabstol, cv_mem->cv_tempv);
    // N_VInv(cv_mem->cv_tempv, weight);
    realtype tmp = abs(ycur[i]);
    tempv[i] = reltol*tmp + Sabstol;
    weight[i] = one/tempv[i];
  }
}

extern "C"
int cvEwtSetSS_fused(const realtype reltol,
                     const realtype Sabstol,
                     const N_Vector ycur,
                     N_Vector tempv,
                     N_Vector weight)
{
  const SUNCudaExecPolicy* exec_policy = ((N_VectorContent_Cuda)weight->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(weight);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvEwtSetSS_cukernel<<<grid, block, 0, exec_policy->stream()>>>
  (
    N,
    reltol,
    Sabstol,
    N_VGetDeviceArrayPointer_Cuda(ycur),
    N_VGetDeviceArrayPointer_Cuda(tempv),
    N_VGetDeviceArrayPointer_Cuda(weight)
  );

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return -1;
#endif

  return 0;
}


/*
 * -----------------------------------------------------------------
 * Compute the ewt vector when the tol type is CV_SV.
 * -----------------------------------------------------------------
 */


__global__
void cvEwtSetSV_cukernel(const sunindextype length,
                        const realtype reltol,
                        const realtype* Vabstol,
                        const realtype* ycur,
                        realtype* tempv,
                        realtype* weight)
{
  const realtype one = 1.0;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VAbs(ycur, cv_mem->cv_tempv);
    // N_VLinearSum(cv_mem->cv_reltol, cv_mem->cv_tempv, ONE,
    //             cv_mem->cv_Vabstol, cv_mem->cv_tempv);
    // N_VInv(cv_mem->cv_tempv, weight);
    realtype tmp = abs(ycur[i]);
    tempv[i] = reltol*tmp + Vabstol[i];
    weight[i] = one/tempv[i];
  }
}

extern "C"
int cvEwtSetSV_fused(const realtype reltol,
                     const N_Vector Vabstol,
                     const N_Vector ycur,
                     N_Vector tempv,
                     N_Vector weight)
{
  const SUNCudaExecPolicy* exec_policy = ((N_VectorContent_Cuda)weight->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(weight);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvEwtSetSV_cukernel<<<grid, block, 0, exec_policy->stream()>>>
  (
    N,
    reltol,
    N_VGetDeviceArrayPointer_Cuda(Vabstol),
    N_VGetDeviceArrayPointer_Cuda(ycur),
    N_VGetDeviceArrayPointer_Cuda(tempv),
    N_VGetDeviceArrayPointer_Cuda(weight)
  );

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return -1;
#endif

  return 0;
}


/*
 * -----------------------------------------------------------------
 * Determine if the constraints of the problem are satisfied by
 * the proposed step.
 * -----------------------------------------------------------------
 */


__global__
void cvCheckConstraints_cukernel(const sunindextype length,
                                 const realtype* c,
                                 const realtype* ewt,
                                 const realtype* y,
                                 const realtype* mm,
                                 realtype* tempv)
{
  static const realtype zero = 0.0;
  static const realtype pt1 = 0.1;
  static const realtype one = 1.0;
  static const realtype onept5 = 1.5;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VCompare(ONEPT5, cv_mem->cv_constraints, tmp); /* a[i]=1 when |c[i]|=2  */
    // N_VProd(tmp, cv_mem->cv_constraints, tmp);       /* a * c                 */
    // N_VDiv(tmp, cv_mem->cv_ewt, tmp);                /* a * c * wt            */
    // N_VLinearSum(ONE, cv_mem->cv_y, -PT1, tmp, tmp); /* y - 0.1 * a * c * wt  */
    // N_VProd(tmp, mm, tmp);                           /* v = mm*(y-0.1*a*c*wt) */
    realtype tmp = (abs(c[i]) >= onept5) ? one : zero;
    tmp = tmp*c[i];
    tmp = tmp/ewt[i];
    tmp = y[i] - pt1*tmp;
    tempv[i] = tmp*mm[i];
  }
}

extern "C"
int cvCheckConstraints_fused(const N_Vector c,
                             const N_Vector ewt,
                             const N_Vector y,
                             const N_Vector mm,
                             N_Vector tempv)
{
  const SUNCudaExecPolicy* exec_policy = ((N_VectorContent_Cuda)c->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(c);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvCheckConstraints_cukernel<<<grid, block, 0, exec_policy->stream()>>>
  (
    N,
    N_VGetDeviceArrayPointer_Cuda(c),
    N_VGetDeviceArrayPointer_Cuda(ewt),
    N_VGetDeviceArrayPointer_Cuda(y),
    N_VGetDeviceArrayPointer_Cuda(mm),
    N_VGetDeviceArrayPointer_Cuda(tempv)
  );

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return -1;
#endif

  return 0;
}

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

/*
 * -----------------------------------------------------------------
 * Attach CUDA integrator fused kernels
 * -----------------------------------------------------------------
 */

extern "C"
int CVodeSetUseIntegratorFusedKernels_CUDA(void *cvode_mem, booleantype onoff)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE",
                   "CVodeSetUseIntegratorFusedKernels", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

#ifdef SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS
  if (!cv_mem->cv_MallocDone ||
      N_VGetVectorID(cv_mem->cv_ewt) != SUNDIALS_NVEC_CUDA) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE",
                   "CVodeSetUseIntegratorFusedKernels", MSGCV_BAD_NVECTOR);
    return(CV_MEM_NULL);
  }
  cv_mem->cv_usefused = onoff;
  if (onoff)
  {
    cv_mem->EwtSS_fused = cvEwtSetSS_fused;
    cv_mem->EwtSV_fused = cvEwtSetSV_fused;
    cv_mem->CheckConstr_fused = cvCheckConstraints_fused;
    cv_mem->Resid_fused = cvNlsResid_fused;
  }
  else
  {
    cv_mem->EwtSS_fused = NULL;
    cv_mem->EwtSV_fused = NULL;
    cv_mem->CheckConstr_fused = NULL;
    cv_mem->Resid_fused = NULL;
  }
  return(CV_SUCCESS);
#else
  cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE",
                 "CVodeSetUseIntegratorFusedKernels",
                 "CVODE was not built with fused integrator kernels enabled");
  return(CV_ILL_INPUT);
#endif
}
