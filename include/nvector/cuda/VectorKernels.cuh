/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 */


#ifndef _VECTOR_KERNELS_CUH_
#define _VECTOR_KERNELS_CUH_

#include <limits>
#include <cuda_runtime.h>


namespace suncudavec
{

/* -----------------------------------------------------------------
 * The namespace for CUDA kernels
 *
 * Reduction CUDA kernels in nvector are based in part on "reduction"
 * example in NVIDIA Corporation CUDA Samples, and parallel reduction
 * examples in textbook by J. Cheng at al. "CUDA C Programming".
 * -----------------------------------------------------------------
 */
namespace math_kernels
{


/**
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


/**
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


/**
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


/**
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


/**
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


/**
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


/**
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


/**
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


/**
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
 * Weighted root mean square norm of a vector.
 *    
 */
template <typename T, typename I>
__global__ void
wrmsNormKernel(const T *x, const T *w, T *out, I n)
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
 * Weighted root mean square norm of a vector values selected by id.
 *
 */
template <typename T, typename I>
__global__ void
wrmsNormMaskKernel(const T *x, const T *w, const T *id, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  T sum = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n && id[i] > 0.0)
    sum = x[i] * w[i] * x[i] * w[i];
  if ((i + blockDim.x < n) && (id[i + blockDim.x] > 0.0))
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
 * Weighted root mean square notm of a vector.
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

  T flag;

  // First reduction step before storing data in shared memory.
  if (i < n && x[i] == 0.0) {
    flag = 1.0;
  } else {
    flag = 0.0;
    z[i] = 1.0/x[i];
  }

  if (i + blockDim.x < n && x[i + blockDim.x] == 0.0)
  {
    flag += 1.0;
  }
  else
  {
    z[i + blockDim.x] = 1.0/x[i + blockDim.x];
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
 * results are stored in vector 'm'. A reduction is performed to set a
 * flag > 0 if any of the constraints is violated.
 *
 */
template <typename T, typename I>
__global__ void
constrMaskKernel(const T *c, const T *x, T *m, T *out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // First reduction step before storing data in shared memory.

  // test1 = true if test failed
  bool test1 = (abs(c[i]) > 1.5 && c[i]*x[i] <= 0.0) ||
      (abs(c[i]) > 0.5 && c[i]*x[i] <  0.0);
  T sum = m[i] = (i < n && test1) ? 1.0 : 0.0;

  // test2 = true if test failed
  bool test2 = (abs(c[i + blockDim.x]) > 1.5 && c[i + blockDim.x]*x[i + blockDim.x] <= 0.0) ||
      (abs(c[i + blockDim.x]) > 0.5 && c[i + blockDim.x]*x[i + blockDim.x] <  0.0);
  m[i+blockDim.x] = (i+blockDim.x < n && test2) ? 1.0 : 0.0;
  sum += m[i+blockDim.x];

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
wrmsNormVectorArrayKernel(int nv, T** xd, T** wd, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    
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
wrmsNormMaskVectorArrayKernel(int nv, T** xd, T** wd, T* id, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    
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

} // namespace math_kernels




template <typename T, typename I>
inline cudaError_t setConst(T a, Vector<T,I>& X)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::setConstKernel<<<grid, block>>>(a, X.device(), X.size());
  return cudaGetLastError();
}

template <typename T, typename I>
inline cudaError_t linearSum(T a, const Vector<T,I>& X, T b, const Vector<T,I>& Y, Vector<T,I>& Z)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::linearSumKernel<<<grid, block>>>(a, X.device(), b, Y.device(), Z.device(), X.size());
  return cudaGetLastError();
}

template <typename T, typename I>
inline cudaError_t prod(const Vector<T,I>& X, const Vector<T,I>& Y, Vector<T,I>& Z)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::prodKernel<<<grid, block>>>(X.device(), Y.device(), Z.device(), X.size());
  return cudaGetLastError();
}

template <typename T, typename I>
inline cudaError_t div(const Vector<T,I>& X, const Vector<T,I>& Y, Vector<T,I>& Z)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::divKernel<<<grid, block>>>(X.device(), Y.device(), Z.device(), X.size());
  return cudaGetLastError();
}

template <typename T, typename I>
inline cudaError_t scale(T const a, const Vector<T,I>& X, Vector<T,I>& Z)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::scaleKernel<<<grid, block>>>(a, X.device(), Z.device(), X.size());
  return cudaGetLastError();
}

template <typename T, typename I>
inline cudaError_t absVal(const Vector<T,I>& X, Vector<T,I>& Z)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::absKernel<<<grid, block>>>(X.device(), Z.device(), X.size());
  return cudaGetLastError();
}

template <typename T, typename I>
inline cudaError_t inv(const Vector<T,I>& X, Vector<T,I>& Z)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::invKernel<<<grid, block>>>(X.device(), Z.device(), X.size());
  return cudaGetLastError();
}

template <typename T, typename I>
inline cudaError_t addConst(T const a, const Vector<T,I>& X, Vector<T,I>& Z)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::addConstKernel<<<grid, block>>>(a, X.device(), Z.device(), X.size());
  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t compare(T const c, const Vector<T,I>& X, Vector<T,I>& Z)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::compareKernel<<<grid, block>>>(c, X.device(), Z.device(), X.size());
  return cudaGetLastError();
}


template <typename T, typename I>
inline T dotProd(const Vector<T,I>& x, const Vector<T,I>& y)
{
  // Set partitioning
  ReducePartitioning<T, I>& p = x.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::dotProdKernel<T,I><<< grid, block, shMemSize >>>(x.device(), y.device(), p.devBuffer(), x.size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // Rerun reduction kernel
    math_kernels::sumReduceKernel<T,I><<< grid, block, shMemSize >>>(p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result += p.hostBuffer()[i];
  }
  return gpu_result;
}

template <typename T, typename I>
inline T maxNorm(const Vector<T,I>& x)
{
  // Set partitioning
  ReducePartitioning<T, I>& p = x.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::maxNormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), p.devBuffer(), x.size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // (Re)run reduction kernel
    math_kernels::maxNormKernel<T,I><<< grid, block, shMemSize >>>(p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    if (p.hostBuffer()[i] > gpu_result)
      gpu_result = p.hostBuffer()[i];
  }
  return gpu_result;
}

template <typename T, typename I>
inline T wrmsNorm(const Vector<T,I>& x, const Vector<T,I>& w)
{
  // Set partitioning
  ReducePartitioning<T, I>& p = x.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::wrmsNormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), w.device(), p.devBuffer(), x.size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // (Re)run reduction kernel
    math_kernels::sumReduceKernel<T,I><<< grid, block, shMemSize >>>(p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result += p.hostBuffer()[i];
  }
  return sqrt(gpu_result/x.size());
}

template <typename T, typename I>
inline T wrmsNormMask(const Vector<T,I>& x, const Vector<T,I>& w, const Vector<T,I>& id)
{
  // Set partitioning
  ReducePartitioning<T, I>& p = x.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::wrmsNormMaskKernel<T,I><<< grid, block, shMemSize >>>(x.device(), w.device(), id.device(), p.devBuffer(), x.size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // (Re)run reduction kernel
    math_kernels::sumReduceKernel<T,I><<< grid, block, shMemSize >>>(p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result += p.hostBuffer()[i];
  }
  return sqrt(gpu_result/x.size());
}

template <typename T, typename I>
inline T findMin(const Vector<T,I>& x)
{
  T maxVal = std::numeric_limits<T>::max();

  // Set partitioning
  ReducePartitioning<T, I>& p = x.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::findMinKernel<T,I><<< grid, block, shMemSize >>>(maxVal, x.device(), p.devBuffer(), x.size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // Rerun reduction kernel
    math_kernels::findMinKernel<T,I><<< grid, block, shMemSize >>>(maxVal, p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    if (p.hostBuffer()[i] < gpu_result)
      gpu_result = p.hostBuffer()[i];
  }
  return gpu_result;
}


template <typename T, typename I>
inline T wL2Norm(const Vector<T,I>& x, const Vector<T,I>& y)
{
  // Set partitioning
  ReducePartitioning<T, I>& p = x.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::wrmsNormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), y.device(), p.devBuffer(), x.size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // Rerun reduction kernel
    math_kernels::sumReduceKernel<T,I><<< grid, block, shMemSize >>>(p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result += p.hostBuffer()[i];
  }
  return sqrt(gpu_result);
}


template <typename T, typename I>
inline T L1Norm(const Vector<T,I>& x)
{
  // Set partitioning
  ReducePartitioning<T, I>& p = x.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::L1NormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), p.devBuffer(), x.size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // Rerun reduction kernel
    math_kernels::sumReduceKernel<T,I><<< grid, block, shMemSize >>>(p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result += p.hostBuffer()[i];
  }
  return gpu_result;
}


template <typename T, typename I>
inline bool invTest(const Vector<T,I>& x, Vector<T,I>& z)
{
  // Set partitioning
  ReducePartitioning<T, I>& p = x.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::invTestKernel<T,I><<< grid, block, shMemSize >>>(x.device(), z.device(), p.devBuffer(), x.size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // Rerun reduction kernel
    math_kernels::sumReduceKernel<T,I><<< grid, block, shMemSize >>>(p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result += p.hostBuffer()[i];
  }
  return !(gpu_result > 0.0);
}


template <typename T, typename I>
inline bool constrMask(const Vector<T,I>& c, const Vector<T,I>& x, Vector<T,I>& m)
{
  // Set partitioning
  ReducePartitioning<T, I>& p = x.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::constrMaskKernel<T,I><<< grid, block, shMemSize >>>(c.device(), x.device(), m.device(), p.devBuffer(), x.size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // Rerun reduction kernel
    math_kernels::sumReduceKernel<T,I><<< grid, block, shMemSize >>>(p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result += p.hostBuffer()[i];
  }
  return (gpu_result < 0.5);
}


template <typename T, typename I>
inline T minQuotient(const Vector<T,I>& num, const Vector<T,I>& den)
{
  // Starting value for min reduction
  T maxVal = std::numeric_limits<T>::max();

  // Set partitioning
  ReducePartitioning<T, I>& p = num.partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = p.shmem();

  math_kernels::minQuotientKernel<T,I><<< grid, block, shMemSize >>>(maxVal, num.device(), den.device(), p.devBuffer(), num.size());

  // All quotients are computed by now. Find the minimum.
  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    p.setPartitioning(n, grid, block, shMemSize);

    // Rerun reduction kernel
    math_kernels::findMinKernel<T,I><<< grid, block, shMemSize >>>(maxVal, p.devBuffer(), p.devBuffer(), n);
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  p.copyFromDevBuffer(n);

  T gpu_result = p.hostBuffer()[0];
  for (unsigned int i=1; i<n; i++)
  {
    if (p.hostBuffer()[i] < gpu_result)
      gpu_result = p.hostBuffer()[i];
  }
  return gpu_result;
}


/*
 * -----------------------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------------------
 */

template <typename T, typename I>
inline cudaError_t linearCombination(int nvec, T* c, Vector<T,I>** X, Vector<T,I>* Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  StreamPartitioning<T, I>& p = X[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::linearCombinationKernel<<<grid, block>>>(nvec, d_c, d_Xd,
                                                         Z->device(), Z->size());

  // Free host array
  delete[] h_Xd;

  // Free device arrays
  err = cudaFree(d_c);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t scaleAddMulti(int nvec, T* c, Vector<T,I>* X,
                                 Vector<T,I>** Y, Vector<T,I>** Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Yd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = Y[i]->device();

  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  StreamPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::scaleAddMultiKernel<<<grid, block>>>(nvec, d_c, X->device(),
                                                     d_Yd, d_Zd, X->size());

  // Free host array
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_c);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t dotProdMulti(int nvec, Vector<T,I>* x, Vector<T,I>** Y,
                                T* dots)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Yd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = Y[i]->device();

  // Copy array of device pointers to device from host
  T** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ReducePartitioning<T, I>& p = x->partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = nvec*block*sizeof(T);

  // Allocate reduction buffer on device
  T* d_buff;
  err = cudaMalloc((void**) &d_buff, nvec*grid*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  
  math_kernels::dotProdMultiKernel<T,I><<< grid, block, shMemSize >>>(nvec,
                                                                      x->device(),
                                                                      d_Yd,
                                                                      d_buff,
                                                                      x->size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax) {

    // Recompute partitioning
    grid = (n + block - 1)/block;

    // Rerun reduction kernel
    math_kernels::sumReduceVectorKernel<T,I><<< grid, block, shMemSize >>>(nvec,
                                                                           d_buff,
                                                                           d_buff,
                                                                           n);

    // update buffer array working length
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  T* h_buff = new T[nvec*n*sizeof(T)];
  err = cudaMemcpy(h_buff, d_buff, nvec*n*sizeof(T), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) return cudaGetLastError();

  for (int k=0; k<nvec; k++) {
    dots[k] = h_buff[k*n];
    for (int i=1; i<n; i++){
      dots[k] += h_buff[i + k*n];
    }
  }

  // Free host array
  delete[] h_Yd;
  delete[] h_buff;

  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_buff);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

template <typename T, typename I>
inline cudaError_t linearSumVectorArray(int nvec, T a, Vector<T,I>** X, T b,
                                        Vector<T,I>** Y, Vector<T,I>** Z)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Yd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = Y[i]->device();

  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  StreamPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::linearSumVectorArrayKernel<<<grid, block>>>(nvec, a, d_Xd, b,
                                                            d_Yd, d_Zd, Z[0]->size());

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t scaleVectorArray(int nvec, T* c, Vector<T,I>** X,
                                    Vector<T,I>** Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  StreamPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::scaleVectorArrayKernel<<<grid, block>>>(nvec, d_c, d_Xd, d_Zd,
                                                        Z[0]->size());

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t constVectorArray(int nvec, T c, Vector<T,I>** Z)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  StreamPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::constVectorArrayKernel<<<grid, block>>>(nvec, c, d_Zd,
                                                        Z[0]->size());

  // Free host array
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t wrmsNormVectorArray(int nvec, Vector<T,I>** X,
                                       Vector<T,I>** W, T* nrm)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Wd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = W[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Wd;
  err = cudaMalloc((void**) &d_Wd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Wd, h_Wd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ReducePartitioning<T, I>& p = X[0]->partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = nvec*block*sizeof(T);

  // Allocate reduction buffer on device
  T* d_buff;
  err = cudaMalloc((void**) &d_buff, nvec*grid*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  
  math_kernels::wrmsNormVectorArrayKernel<<< grid, block, shMemSize >>>(nvec,
                                                                        d_Xd,
                                                                        d_Wd,
                                                                        d_buff,
                                                                        X[0]->size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax) {

    // Recompute partitioning
    grid = (n + block - 1)/block;

    // Rerun reduction kernel
    math_kernels::sumReduceVectorKernel<T,I><<< grid, block, shMemSize >>>(nvec,
                                                                           d_buff,
                                                                           d_buff,
                                                                           n);

    // update buffer array working length
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  T* h_buff = new T[nvec*n*sizeof(T)];
  err = cudaMemcpy(h_buff, d_buff, nvec*n*sizeof(T), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) return cudaGetLastError();

  for (int k=0; k<nvec; k++) {
    nrm[k] = h_buff[k*n];
    for (int i=1; i<n; i++){
      nrm[k] += h_buff[i + k*n];
    }
    nrm[k] = sqrt(nrm[k]/X[0]->size());
  }

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;
  delete[] h_buff;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Wd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_buff);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t wrmsNormMaskVectorArray(int nvec, Vector<T,I>** X,
                                           Vector<T,I>** W, Vector<T,I>* ID,
                                           T* nrm)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Wd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = W[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Wd;
  err = cudaMalloc((void**) &d_Wd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Wd, h_Wd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ReducePartitioning<T, I>& p = X[0]->partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = nvec*block*sizeof(T);

  // Allocate reduction buffer on device
  T* d_buff;
  err = cudaMalloc((void**) &d_buff, nvec*grid*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  
  math_kernels::wrmsNormMaskVectorArrayKernel<<< grid, block, shMemSize >>>(nvec,
                                                                            d_Xd,
                                                                            d_Wd,
                                                                            ID->device(),
                                                                            d_buff,
                                                                            X[0]->size());

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax) {

    // Recompute partitioning
    grid = (n + block - 1)/block;

    // Rerun reduction kernel
    math_kernels::sumReduceVectorKernel<T,I><<< grid, block, shMemSize >>>(nvec,
                                                                           d_buff,
                                                                           d_buff,
                                                                           n);

    // update buffer array working length
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  T* h_buff = new T[nvec*n*sizeof(T)];
  err = cudaMemcpy(h_buff, d_buff, nvec*n*sizeof(T), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) return cudaGetLastError();

  for (int k=0; k<nvec; k++) {
    nrm[k] = h_buff[k*n];
    for (int i=1; i<n; i++){
      nrm[k] += h_buff[i + k*n];
    }
    nrm[k] = sqrt(nrm[k]/X[0]->size());
  }

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;
  delete[] h_buff;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Wd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_buff);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t scaleAddMultiVectorArray(int nvec, int nsum, T* c,
                                            Vector<T,I>** X, Vector<T,I>** Y,
                                            Vector<T,I>** Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nsum*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Yd = new T*[nsum*nvec];
  for (int i=0; i<nsum*nvec; i++)
    h_Yd[i] = Y[i]->device();

  T** h_Zd = new T*[nsum*nvec];
  for (int i=0; i<nsum*nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nsum*nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nsum*nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nsum*nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nsum*nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  StreamPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::scaleAddMultiVectorArrayKernel<<<grid, block>>>(nvec, nsum, d_c,
                                                                d_Xd, d_Yd, d_Zd,
                                                                Z[0]->size());

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t linearCombinationVectorArray(int nvec, int nsum, T* c,
                                                Vector<T,I>** X, Vector<T,I>** Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nsum*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Xd = new T*[nsum*nvec];
  for (int i=0; i<nsum*nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nsum*nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nsum*nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  StreamPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::linearCombinationVectorArrayKernel<<<grid, block>>>(nvec, nsum,
                                                                    d_c, d_Xd,
                                                                    d_Zd,
                                                                    Z[0]->size());

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


} // namespace nvec



#endif // _VECTOR_KERNELS_CUH_
