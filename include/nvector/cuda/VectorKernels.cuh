/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
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
 * This software contains source code provided by NVIDIA Corporation.
 * Reduction CUDA kernels in nvector are based in part on "reduction"
 * example in NVIDIA Corporation CUDA Samples. The NVIDIA CUDA Samples
 * License Agreement is available in Chapter 2 in NVIDIA-EULA.txt
 * document.
 * -----------------------------------------------------------------
 */


#ifndef _VECTOR_KERNELS_CUH_
#define _VECTOR_KERNELS_CUH_

#include <limits>

// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>

/// Forward declarations of Vector and ThreadPartitioning classes
namespace suncudavec
{
    template <typename T, typename I>
    class Vector;

    template <typename T, typename I>
    class ReducePartitioning;

    template <typename T, typename I>
    class StreamPartitioning;
}

//#define abs(x) ((x)<0 ? -(x) : (x))

namespace suncudavec
{
namespace math_kernels
{
// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
    __device__ inline operator T *()
    {
        extern __shared__ float sundials_shmem_ptr[];
        return (T *) sundials_shmem_ptr;
    }
    
    __device__ inline operator const T *() const
    {
        extern __shared__ float sundials_shmem_ptr[];
        return (T *) sundials_shmem_ptr;
    }
};

// specialize for double to avoid unaligned memory
// access compile errors
template<>
struct SharedMemory<double>
{
    __device__ inline operator double *()
    {
        extern __shared__ double sundials_shmem_ptr_double[];
        return (double *) sundials_shmem_ptr_double;
    }

    __device__ inline operator const double *() const
    {
        extern __shared__ double sundials_shmem_ptr_double[];
        return (double *) sundials_shmem_ptr_double;
    }
};
        
// specialize for long double to avoid unaligned memory
// access compile errors
template<>
struct SharedMemory<long double>
{
    __device__ inline operator long double *()
    {
        extern __shared__ long double sundials_shmem_ptr_long_double[];
        return (long double *) sundials_shmem_ptr_long_double;
    }

    __device__ inline operator const long double *() const
    {
        extern __shared__ long double sundials_shmem_ptr_long_double[];
        return (long double *) sundials_shmem_ptr_long_double;
    }
};


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


template <typename T, typename I>
__global__ void
axpyKernel(T a, const T *X, T *Y, I n)
{
    I i = blockDim.x * blockIdx.x + threadIdx.x;
    
    if (i < n)
    {
        Y[i] += a*X[i];
    }
}


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
 * It is based off reduce3 kernel from CUDA examples. Uses n/2 threads.
 * Performs the first level of reduction when reading from global memory.
 *
 */
template <typename T, typename I>
__global__ void
sumReduceKernel(const T *g_idata, T *g_odata, I n)
{
    T *sdata = SharedMemory<T>();
    
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    I tid = threadIdx.x;
    I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    
    T mySum = (i < n) ? g_idata[i] : 0;
    
    if (i + blockDim.x < n)
        mySum += g_idata[i+blockDim.x];
    
    sdata[tid] = mySum;
    __syncthreads();
    
    // do reduction in shared mem (blockDim.x is a power of 2)
    for (I s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }
        
        __syncthreads();
    }
    
    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = mySum;
}


/*
 * Dot product of two vectors.
 *    
 * It is based off reduce3 kernel from CUDA examples. Uses n/2 threads.
 * Performs the first level of reduction when reading from global memory.
 *
 */
template <typename T, typename I>
__global__ void
dotProdKernel(const T *g_idata1, const T *g_idata2, T *g_odata, I n)
{
    T *sdata = SharedMemory<T>();
    
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    I tid = threadIdx.x;
    I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    
    T mySum = (i < n) ? g_idata1[i] * g_idata2[i] : 0;
    
    if (i + blockDim.x < n)
        mySum += ( g_idata1[i+blockDim.x] * g_idata2[i+blockDim.x]);
    
    sdata[tid] = mySum;
    __syncthreads();
    
    // do reduction in shared mem
    for (I s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }
        
        __syncthreads();
    }
    
    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = mySum;
}


/*
 * Finds max norm the vector.
 *    
 * It is based off reduce3 kernel from CUDA examples. Uses n/2 threads.
 * Performs the first level of reduction when reading from global memory.
 *
 */
template <typename T, typename I>
__global__ void
maxNormKernel(const T *g_idata, T *g_odata, I n)
{
    T *sdata = SharedMemory<T>();
    
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    I tid = threadIdx.x;
    I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    
    T maximum = (i < n) ? abs(g_idata[i]) : 0;
    
    if (i + blockDim.x < n)
        maximum = (abs(g_idata[i+blockDim.x]) > maximum) ? abs(g_idata[i+blockDim.x]) : maximum;
    
    sdata[tid] = maximum;
    __syncthreads();
    
    // do reduction in shared mem (blockDim.x is a power of 2)
    for (I s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] = maximum = sdata[tid + s] > maximum ? sdata[tid + s] : maximum;
        }
        
        __syncthreads();
    }
    
    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = maximum;
}


/*
 * Weighted root mean square norm of a vector.
 *    
 * It is based off reduce3 kernel from CUDA examples. Uses n/2 threads.
 * Performs the first level of reduction when reading from global memory.
 *
 */
template <typename T, typename I>
__global__ void
wrmsNormKernel(const T *g_idata1, const T *g_idata2, T *g_odata, I n)
{
    T *sdata = SharedMemory<T>();
    
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    I tid = threadIdx.x;
    I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    
    T mySum = (i < n) ? g_idata1[i] * g_idata2[i] * g_idata1[i] * g_idata2[i] : 0;
    
    if (i + blockDim.x < n)
        mySum += ( g_idata1[i+blockDim.x] * g_idata2[i+blockDim.x] * g_idata1[i+blockDim.x] * g_idata2[i+blockDim.x] );
    
    sdata[tid] = mySum;
    __syncthreads();
    
    // do reduction in shared mem
    for (I s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }
        
        __syncthreads();
    }
    
    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = mySum;
}

/*
 * Weighted root mean square notm of a vector.
 *
 * It is based off reduce3 kernel from CUDA examples. Uses n/2 threads.
 * Performs the first level of reduction when reading from global memory.
 *
 */
template <typename T, typename I>
__global__ void
wrmsNormMaskKernel(const T *g_idata1, const T *g_idata2, const T *g_id, T *g_odata, I n)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    I tid = threadIdx.x;
    I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

    T mySum = (i < n && g_id[i] > 0.0) ? g_idata1[i] * g_idata2[i] * g_idata1[i] * g_idata2[i] : 0.0;

    if ((i + blockDim.x < n) && (g_id[i] > 0.0))
        mySum += ( g_idata1[i+blockDim.x] * g_idata2[i+blockDim.x] * g_idata1[i+blockDim.x] * g_idata2[i+blockDim.x]);

    sdata[tid] = mySum;
    __syncthreads();

    // do reduction in shared mem
    for (I s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }

        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = mySum;
}

/*
 * Finds min value in the vector.
 *    
 * It is based off reduce3 kernel from CUDA examples. Uses n/2 threads.
 * Performs the first level of reduction when reading from global memory.
 *
 */
template <typename T, typename I>
__global__ void
findMinKernel(T MAX, const T *g_idata, T *g_odata, I n)
{
    T *sdata = SharedMemory<T>();
    
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    I tid = threadIdx.x;
    I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    
    T minimum = (i < n) ? g_idata[i] : MAX;
    
    if (i + blockDim.x < n)
        minimum = (g_idata[i+blockDim.x]) < minimum ? g_idata[i+blockDim.x] : minimum;
    
    sdata[tid] = minimum;
    __syncthreads();
    
    // do reduction in shared mem (blockDim.x is a power of 2)
    for (I s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] = minimum = sdata[tid + s] < minimum ? sdata[tid + s] : minimum;
        }
        
        __syncthreads();
    }
    
    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = minimum;
}


/*
 * Weighted root mean square notm of a vector.
 *
 * It is based off reduce3 kernel from CUDA examples. Uses n/2 threads.
 * Performs the first level of reduction when reading from global memory.
 *
 */
template <typename T, typename I>
__global__ void
L1NormKernel(const T *g_idata1, T *g_odata, I n)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    I tid = threadIdx.x;
    I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

    T mySum = (i < n) ? abs(g_idata1[i]) : 0.0;

    if (i + blockDim.x < n)
        mySum += abs(g_idata1[i+blockDim.x]);

    sdata[tid] = mySum;
    __syncthreads();

    // do reduction in shared mem
    for (I s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }

        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = mySum;
}


template <typename T, typename I>
__global__ void
invTestKernel(const T *g_idata, T *g_odata1, T *g_buffer, I n)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    I tid = threadIdx.x;
    I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

    T mySum; // = (i < n) ? g_idata[i] : 0;

    if (i < n && g_idata[i] == 0.0)
    {
      mySum = 1.0;
    }
    else
    {
      mySum = 0.0;
      g_odata1[i] = 1.0/g_idata[i];
    }

    if (i + blockDim.x < n && g_idata[i + blockDim.x] == 0.0)
    {
      mySum += 1.0;
    }
    else
    {
      g_odata1[i + blockDim.x] = 1.0/g_idata[i + blockDim.x];
    }

    sdata[tid] = mySum;
    __syncthreads();

    // do reduction in shared mem (blockDim.x is a power of 2)
    for (I s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }

        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_buffer[blockIdx.x] = mySum;
}

/*
 * Checks if inequality constraints are satisfied.
 *
 * It is based off reduce3 kernel from CUDA examples. Uses n/2 threads.
 * Performs the first level of reduction when reading from global memory.
 *
 */
template <typename T, typename I>
__global__ void
constrMaskKernel(const T *g_c, const T *g_x, T *g_m, T *g_odata, I n)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    // computing mask for failed constarints g_m on the fly
    I tid = threadIdx.x;
    I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

    // test1 = true if test failed
    bool test1 = (abs(g_c[i]) > 1.5 && g_c[i]*g_x[i] <= 0.0) ||
                 (abs(g_c[i]) > 0.5 && g_c[i]*g_x[i] <  0.0);
    T mySum = g_m[i] = (i < n && test1) ? 1.0 : 0.0;

    // test2 = true if test failed
    bool test2 = (abs(g_c[i + blockDim.x]) > 1.5 && g_c[i + blockDim.x]*g_x[i + blockDim.x] <= 0.0) ||
                 (abs(g_c[i + blockDim.x]) > 0.5 && g_c[i + blockDim.x]*g_x[i + blockDim.x] <  0.0);
    g_m[i+blockDim.x] = (i+blockDim.x < n && test2) ? 1.0 : 0.0;
    mySum += g_m[i+blockDim.x];

    sdata[tid] = mySum;
    __syncthreads();

    // do reduction in shared mem (blockDim.x is a power of 2)
    for (I s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }

        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = mySum;
}



/*
 * Finds minimum component-wise quotient.
 *
 * TODO: Replace '0.0' with 'ZERO'
 *
 * It is based off reduce3 kernel from CUDA examples. Uses n/2 threads.
 * Performs the first level of reduction when reading from global memory.
 *
 */
template <typename T, typename I>
__global__ void
minQuotientKernel(const T MAX_VAL, const T *num, const T *den, T *min_quotient, I n)
{
    T *sdata = SharedMemory<T>();

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

    sdata[tid] = minimum;
    __syncthreads();

    // Perform reduction block-wise in shared memory.
    for (I s = blockDim.x/2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
          minimum = min(sdata[tid + s], minimum);
          sdata[tid] = minimum; // = sdata[tid + s] < minimum ? sdata[tid + s] : minimum;
        }

        __syncthreads();
    }

    // Copy reduction result for each block to global memory
    if (tid == 0)
      min_quotient[blockIdx.x] = minimum;
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
inline cudaError_t axpy(T a, const Vector<T,I>& X, Vector<T,I>& Y)
{
  // Set partitioning
  StreamPartitioning<T, I>& p = X.partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();

  math_kernels::axpyKernel<<<grid, block>>>(a, X.device(), Y.device(), X.size());
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();
    
    math_kernels::dotProdKernel<T,I><<< grid, block, shMemSize >>>(x.device(), y.device(), p.devBuffer(), x.size());
    
    I n = grid;
    I nmax = 2*block;
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();
    
    math_kernels::maxNormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), p.devBuffer(), x.size());
    
    I n = grid;
    I nmax = 2*block;
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();
    
    math_kernels::wrmsNormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), w.device(), p.devBuffer(), x.size());
    
    I n = grid;
    I nmax = 2*block;
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();

    math_kernels::wrmsNormMaskKernel<T,I><<< grid, block, shMemSize >>>(x.device(), w.device(), id.device(), p.devBuffer(), x.size());

    I n = grid;
    I nmax = 2*block;
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();

    math_kernels::findMinKernel<T,I><<< grid, block, shMemSize >>>(maxVal, x.device(), p.devBuffer(), x.size());

    I n = grid;
    I nmax = 2*block;
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();

    math_kernels::wrmsNormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), y.device(), p.devBuffer(), x.size());

    I n = grid;
    I nmax = 2*block;
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();

    math_kernels::L1NormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), p.devBuffer(), x.size());

    I n = grid;
    I nmax = 2*block;
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();

    math_kernels::invTestKernel<T,I><<< grid, block, shMemSize >>>(x.device(), z.device(), p.devBuffer(), x.size());

    I n = grid;
    I nmax = 2*block;
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();

    math_kernels::constrMaskKernel<T,I><<< grid, block, shMemSize >>>(c.device(), x.device(), m.device(), p.devBuffer(), x.size());

    I n = grid;
    I nmax = 2*block;
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
    I grid                      = p.grid();
    unsigned block              = p.block();
    unsigned shMemSize          = p.shmem();

    math_kernels::minQuotientKernel<T,I><<< grid, block, shMemSize >>>(maxVal, num.device(), den.device(), p.devBuffer(), num.size());

    // All quotients are computed by now. Find the minimum.
    I n = grid;
    I nmax = 2*block;
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



} // namespace nvec



#endif // _VECTOR_KERNELS_CUH_
