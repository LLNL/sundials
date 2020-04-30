/*
 * -----------------------------------------------------------------
 * Programmer(s): David Gardner @ LLNL
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
 */


#ifndef _VECTOR_ARRAY_KERNELS_CUH_
#define _VECTOR_ARRAY_KERNELS_CUH_

#include <limits>
#include <cuda_runtime.h>

#include "sundials_cuda.h"

namespace sundials
{
namespace nvector_cuda
{

/* -----------------------------------------------------------------
 * The namespace for CUDA kernels
 *
 * Reduction CUDA kernels in nvector are based in part on "reduction"
 * example in NVIDIA Corporation CUDA Samples, and parallel reduction
 * examples in textbook by J. Cheng at al. "CUDA C Programming".
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------------------
 * fused vector operation kernels
 * -----------------------------------------------------------------------------
 */

/*
 * Computes the linear combination of nv vectors
 */
template <typename T, typename I>
__global__ void
linearCombinationKernel(int nv, T* c, T** xd, T* zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n) {
    zd[i] = c[0]*xd[0][i];
    for (int j=1; j<nv; j++)
      zd[i] += c[j]*xd[j][i];
  }
}

/*
 * Computes the scaled sum of one vector with nv other vectors
 */
template <typename T, typename I>
__global__ void
scaleAddMultiKernel(int nv, T* c, T* xd, T** yd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int j=0; j<nv; j++)
      zd[j][i] = c[j] * xd[i] + yd[j][i];
}


/*
 * Dot product of one vector with nv other vectors.
 *
 */
template <typename T, typename I>
__global__ void
dotProdMultiKernel(int nv, T* xd, T** yd, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // Initialize shared memory to zero
  for (int k=0; k<nv; k++)
    shmem[tid + k*blockDim.x] = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] = xd[i] * yd[k][i];
  if (i + blockDim.x < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] += (xd[i + blockDim.x] * yd[k][i + blockDim.x]);

  __syncthreads();

  // Perform blockwise reduction in shared memory
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j)
      for (int k=0; k<nv; k++)
        shmem[tid + k*blockDim.x] += shmem[tid + j + k*blockDim.x];

    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    for (int k=0; k<nv; k++)
      out[blockIdx.x + k*gridDim.x] = shmem[k*blockDim.x];
}


/*
 * Sums all elements of the vector.
 *
 */
template <typename T, typename I>
__global__ void
sumReduceVectorKernel(int nv, T* x, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // First reduction step before storing data in shared memory.
  if (i < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] = x[i];
  if (i + blockDim.x < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] += x[i+blockDim.x];

  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j)
      for (int k=0; k<nv; k++)
        shmem[tid + k*blockDim.x] += shmem[tid + j + k*blockDim.x];

    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    for (int k=0; k<nv; k++)
      out[blockIdx.x + k*gridDim.x] = shmem[k*blockDim.x];
}



/*
 * -----------------------------------------------------------------------------
 * vector array operation kernels
 * -----------------------------------------------------------------------------
 */

/*
 * Computes the linear sum of multiple vectors
 */
template <typename T, typename I>
__global__ void
linearSumVectorArrayKernel(int nv, T a, T** xd, T b, T** yd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int j=0; j<nv; j++)
      zd[j][i] = a * xd[j][i] + b * yd[j][i];
}


/*
 * Scales multiple vectors
 */
template <typename T, typename I>
__global__ void
scaleVectorArrayKernel(int nv, T* c, T** xd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int j=0; j<nv; j++)
      zd[j][i] = c[j] * xd[j][i];
}


/*
 * Sets multiple vectors equal to a constant
 */
template <typename T, typename I>
__global__ void
constVectorArrayKernel(int nv, T c, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int j=0; j<nv; j++)
      zd[j][i] = c;
}


/*
 * WRMS norm of nv vectors.
 *
 */
template <typename T, typename I>
__global__ void
wL2NormSquareVectorArrayKernel(int nv, T** xd, T** wd, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // Initialize shared memory to zero
  for (int k=0; k<nv; k++)
    shmem[tid + k*blockDim.x] = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] = xd[k][i] * wd[k][i] * xd[k][i] * wd[k][i];
  if (i + blockDim.x < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] += (xd[k][i + blockDim.x] * wd[k][i + blockDim.x]
                                    * xd[k][i + blockDim.x] * wd[k][i + blockDim.x]);

  __syncthreads();

  // Perform blockwise reduction in shared memory
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j)
      for (int k=0; k<nv; k++)
        shmem[tid + k*blockDim.x] += shmem[tid + j + k*blockDim.x];

    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    for (int k=0; k<nv; k++)
      out[blockIdx.x + k*gridDim.x] = shmem[k*blockDim.x];
}


/*
 * Masked WRMS norm of nv vectors.
 *
 */
template <typename T, typename I>
__global__ void
wL2NormSquareMaskVectorArrayKernel(int nv, T** xd, T** wd, T* id, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // Initialize shared memory to zero
  for (int k=0; k<nv; k++)
    shmem[tid + k*blockDim.x] = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n && id[i] > 0.0)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] = xd[k][i] * wd[k][i] * xd[k][i] * wd[k][i];
  if (i + blockDim.x < n && id[i + blockDim.x] > 0.0)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] += (xd[k][i + blockDim.x] * wd[k][i + blockDim.x]
                                    * xd[k][i + blockDim.x] * wd[k][i + blockDim.x]);

  __syncthreads();

  // Perform blockwise reduction in shared memory
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j)
      for (int k=0; k<nv; k++)
        shmem[tid + k*blockDim.x] += shmem[tid + j + k*blockDim.x];

    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    for (int k=0; k<nv; k++)
      out[blockIdx.x + k*gridDim.x] = shmem[k*blockDim.x];
}


/*
 * Computes the scaled sum of a vector array with multiple other vector arrays
 */
template <typename T, typename I>
__global__ void
scaleAddMultiVectorArrayKernel(int nv, int ns, T* c, T** xd, T** yd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int k=0; k<nv; k++)
      for (int j=0; j<ns; j++)
        zd[k*ns+j][i] = c[j] * xd[k][i] + yd[k*ns+j][i];
}


/*
 * Computes the scaled sum of a vector array with multiple other vector arrays
 */
template <typename T, typename I>
__global__ void
linearCombinationVectorArrayKernel(int nv, int ns, T* c, T** xd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n) {
    for (int k=0; k<nv; k++) {
      zd[k][i] = c[0]*xd[k*ns][i];
      for (int j=1; j<ns; j++) {
        zd[k][i] += c[j]*xd[k*ns+j][i];
      }
    }
  }
}

} // namespace nvector_cuda
} // namespace sundials

#endif // _VECTOR_ARRAY_KERNELS_CUH_
