/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, Cody J. Balos @ LLNL
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


#ifndef _VECTOR_KERNELS_CUH_
#define _VECTOR_KERNELS_CUH_

#include <limits>
#include <cuda_runtime.h>

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
 * Sets all elements of the vector X to constant value a.
 *
 */

template <typename T, typename I>
__global__ void
setConstKernel(T a, T *X, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
  {
    X[i] = a;
  }
}


/*
 * Computes linear sum (combination) of two vectors.
 *
 */

template <typename T, typename I>
__global__ void
linearSumKernel(T a, const T *X, T b, const T *Y, T *Z, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
  {
    Z[i] = a*X[i] + b*Y[i];
  }
}


/*
 * Elementwise product of two vectors.
 *
 */

template <typename T, typename I>
__global__ void
prodKernel(const T *X, const T *Y, T *Z, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
  {
    Z[i] = X[i]*Y[i];
  }
}


/*
 * Elementwise division of two vectors.
 *
 */

template <typename T, typename I>
__global__ void
divKernel(const T *X, const T *Y, T *Z, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
  {
    Z[i] = X[i]/Y[i];
  }
}


/*
 * Scale vector with scalar value 'a'.
 *
 */

template <typename T, typename I>
__global__ void
scaleKernel(T a, const T *X, T *Z, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
  {
    Z[i] = a*X[i];
  }
}


/*
 * Stores absolute values of vector X elements into vector Z.
 *
 */

template <typename T, typename I>
__global__ void
absKernel(const T *X, T *Z, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
  {
    Z[i] = abs(X[i]);
  }
}


/*
 * Elementwise inversion.
 *
 */

template <typename T, typename I>
__global__ void
invKernel(const T *X, T *Z, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
  {
    Z[i] = 1.0/(X[i]);
  }
}


/*
 * Add constant 'c' to each vector element.
 *
 */

template <typename T, typename I>
__global__ void
addConstKernel(T a, const T *X, T *Z, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
  {
    Z[i] = a + X[i];
  }
}


/*
 * Compare absolute values of vector 'X' with constant 'c'.
 *
 */

template <typename T, typename I>
__global__ void
compareKernel(T c, const T *X, T *Z, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
  {
    Z[i] = (abs(X[i]) >= c) ? 1.0 : 0.0;
  }
}


/*
 * Sums all elements of the vector.
 *
 */
template <typename T, typename I>
__global__ void
sumReduceKernel(const T *x, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T sum = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n)
    sum = x[i];
  if (i + blockDim.x < n)
    sum += x[i+blockDim.x];
  shmem[tid] = sum;
  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j) {
      sum += shmem[tid + j];
      shmem[tid] = sum;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    out[blockIdx.x] = sum;
}


/*
 * Dot product of two vectors.
 *
 */
template <typename T, typename I>
__global__ void
dotProdKernel(const T *x, const T *y, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T sum = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n)
    sum = x[i] * y[i];
  if (i + blockDim.x < n)
    sum += ( x[i+blockDim.x] * y[i+blockDim.x]);
  shmem[tid] = sum;
  __syncthreads();

  // Perform blockwise reduction in shared memory
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j) {
      sum += shmem[tid + j];
      shmem[tid] = sum;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    out[blockIdx.x] = sum;
}


/*
 * Finds max norm the vector.
 *
 */
template <typename T, typename I>
__global__ void
maxNormKernel(const T *x, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T maximum = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n)
    maximum = abs(x[i]);
  if (i + blockDim.x < n)
    maximum = max(abs(x[i+blockDim.x]), maximum);
  shmem[tid] = maximum;
  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j) {
      maximum = max(shmem[tid + j], maximum);
      shmem[tid] = maximum;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    out[blockIdx.x] = maximum;
}


/*
 * Weighted L2 norm squared.
 *
 */
template <typename T, typename I>
__global__ void
wL2NormSquareKernel(const T *x, const T *w, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T sum = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n)
    sum = x[i] * w[i] * x[i] * w[i];
  if (i + blockDim.x < n)
    sum += ( x[i+blockDim.x] * w[i+blockDim.x] * x[i+blockDim.x] * w[i+blockDim.x] );

  shmem[tid] = sum;
  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1)
  {
    if (tid < j) {
      sum += shmem[tid + j];
      shmem[tid] = sum;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    out[blockIdx.x] = sum;
}

/*
 * Weighted L2 norm squared with mask. Vector id specifies the mask.
 *
 */
template <typename T, typename I>
__global__ void
wL2NormSquareMaskKernel(const T *x, const T *w, const T *id, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T sum = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n && id[i] > 0.0)
    sum = x[i] * w[i] * x[i] * w[i];
  if ((i + blockDim.x < n) && (id[i+blockDim.x] > 0.0))
    sum += ( x[i+blockDim.x] * w[i+blockDim.x] * x[i+blockDim.x] * w[i+blockDim.x]);
  shmem[tid] = sum;
  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j) {
      sum += shmem[tid + j];
      shmem[tid] = sum;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    out[blockIdx.x] = sum;
}

/*
 * Finds min value in the vector.
 *
 */
template <typename T, typename I>
__global__ void
findMinKernel(T MAX_VAL, const T *x, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T minimum = MAX_VAL;

  // First reduction step before storing data in shared memory.
  if (i < n)
    minimum = x[i];
  if (i + blockDim.x < n)
    minimum = min((x[i+blockDim.x]), minimum);
  shmem[tid] = minimum;
  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j) {
      minimum = min(shmem[tid + j], minimum);
      shmem[tid] = minimum;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    out[blockIdx.x] = minimum;
}


/*
 * Computes L1 norm of vector
 *
 */
template <typename T, typename I>
__global__ void
L1NormKernel(const T *x, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T sum = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n)
    sum = abs(x[i]);
  if (i + blockDim.x < n)
    sum += abs(x[i+blockDim.x]);
  shmem[tid] = sum;
  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j) {
      sum += shmem[tid + j];
      shmem[tid] = sum;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    out[blockIdx.x] = sum;
}

/*
 * Vector inverse  z[i] = 1/x[i] with check for zeros. Reduction is performed
 * to flag the result if any x[i] = 0.
 *
 */
template <typename T, typename I>
__global__ void
invTestKernel(const T *x, T *z, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T flag = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n) {
    if (x[i] == 0.0) {
      flag = 1.0;
    } else {
      flag = 0.0;
      z[i] = 1.0/x[i];
    }
  }

  if (i + blockDim.x < n) {
    if (x[i + blockDim.x] == 0.0) {
      flag += 1.0;
    } else {
      z[i + blockDim.x] = 1.0/x[i + blockDim.x];
    }
  }

  shmem[tid] = flag;
  __syncthreads();

  // Inverse calculation is done. Perform reduction block-wise in shared
  // to find if any x[i] = 0.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j) {
      flag += shmem[tid + j];
      shmem[tid] = flag;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    out[blockIdx.x] = flag;
}

/*
 * Checks if inequality constraints are satisfied. Constraint check
 * results are stored in vector 'm'. A sum reduction over all elements
 * of 'm' is performed to find if any of the constraints is violated.
 * If all constraints are satisfied sum == 0.
 *
 */
template <typename T, typename I>
__global__ void
constrMaskKernel(const T *c, const T *x, T *m, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T sum = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n){
    // test1 = true if constraints violated
    bool test1 = (std::abs(c[i]) > 1.5 && c[i]*x[i] <= 0.0) ||
                 (std::abs(c[i]) > 0.5 && c[i]*x[i] <  0.0);
    m[i] = test1 ? 1.0 : 0.0;
    sum = m[i];
  }

  if (i + blockDim.x < n) {
    // test2 = true if constraints violated
    bool test2 = (std::abs(c[i + blockDim.x]) > 1.5 && c[i + blockDim.x]*x[i + blockDim.x] <= 0.0) ||
                 (std::abs(c[i + blockDim.x]) > 0.5 && c[i + blockDim.x]*x[i + blockDim.x] <  0.0);
    m[i+blockDim.x] = test2 ? 1.0 : 0.0;
    sum += m[i+blockDim.x];
  }

  shmem[tid] = sum;
  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j) {
      sum += shmem[tid + j];
      shmem[tid] = sum;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    out[blockIdx.x] = sum;
}

/*
 * Finds minimum component-wise quotient.
 *
 */
template <typename T, typename I>
__global__ void
minQuotientKernel(const T MAX_VAL, const T *num, const T *den, T *min_quotient, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // Initialize "minimum" to maximum floating point value.
  T minimum = MAX_VAL;
  const T zero = static_cast<T>(0.0);

  // Load vector quotient in the shared memory. Skip if the denominator
  // value is zero.
  if (i < n && den[i] != zero)
    minimum = num[i]/den[i];

  // First level of reduction is upon storing values to shared memory.
  if (i + blockDim.x < n && den[i + blockDim.x] != zero)
    minimum = min(num[i+blockDim.x]/den[i+blockDim.x], minimum);

  shmem[tid] = minimum;
  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j) {
      minimum = min(shmem[tid + j], minimum);
      shmem[tid] = minimum;
    }
    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    min_quotient[blockIdx.x] = minimum;
}

} // namespace nvector_cuda
} // namespace sundials

#endif // _VECTOR_KERNELS_CUH_
