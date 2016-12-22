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


#ifndef _VECTOR_KERNELS_RAJA_CUH_
#define _VECTOR_KERNELS_RAJA_CUH_

#include <limits>

// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
#include <RAJA/RAJA.hxx>

/// Forward declarations of Vector and ThreadPartitioning classes
namespace rvec
{
    template <typename T, typename I>
    class Vector;

    template <typename T, typename I>
    class ThreadPartitioning;
}


namespace rvec
{
namespace math_kernels
{
// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};

// specialize for double to avoid unaligned memory
// access compile errors
template<>
struct SharedMemory<double>
{
    __device__ inline operator       double *()
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }

    __device__ inline operator const double *() const
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }
};

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



} // namespace math_kernels


template <typename T, typename I>
inline cudaError_t linearSum(T a, const Vector<T,I>& X, T b, const Vector<T,I>& Y, Vector<T,I>& Z)
{
//  /* Thread block size is chosen at compile time in current RAJA, */
//  /* used to be a runtime choice.  Easy to switch back to runtime. */
//  const T* xdata = X.device();
//  const T* ydata = Y.device();
//  T* zdata = Z.device();
//  I N = X.size();
//  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(I i) {
//     zdata[i] = a*xdata[i] + b*ydata[i];
//  });
//  return cudaGetLastError();
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::linearSumKernel<<<grid, block>>>(a, X.device(), b, Y.device(), Z.device(), X.size());
    return cudaGetLastError();
}


    
template <typename T, typename I>
inline T dotProd(const rvec::Vector<T,I>& x, const rvec::Vector<T,I>& y)
{
    // Reduction result storage on CPU
    T gpu_result = 0;
    
    // Set partitioning
    ThreadPartitioning<T, I>* p = x.partReduce();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();
    unsigned int shMemSize = p->getShMemSize();
    
    math_kernels::dotProdKernel<T,I><<< grid, block, shMemSize >>>(x.device(), y.device(), p->devBuffer(), x.size());
    
    unsigned int n = grid.x * grid.y * grid.z;
    while (n>2048)
    {
        // Recompute partitioning
        unsigned int grid2;
        unsigned int block2;
        unsigned int shMemSize2;
        p->setPartitioningReduce(n, grid2, block2, shMemSize2);
        
        // Rerun reduction kernel
        math_kernels::sumReduceKernel<T,I><<< grid2, block2, shMemSize2 >>>(p->devBuffer(), p->devBuffer(), n);
        n = (n + 2*block2 - 1) / (2*block2);
    }
    
    // sum partial sums from each block on CPU
    // copy result from device to host
    p->copyFromDevBuffer(n);
    
    for (unsigned int i=0; i<n; i++)
    {
        gpu_result += p->hostBuffer()[i];
    }
    return gpu_result;
}


} // namespace rvec



#endif // _VECTOR_KERNELS_CUH_
