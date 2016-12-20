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
namespace nvec 
{
    template <typename T, typename I>
    class Vector;

    template <typename T, typename I>
    class ThreadPartitioning;
}


namespace nvec
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
 * Weighted root mean square notm of a vector.
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




} // namespace math_kernels

template <typename T, typename I>
inline cudaError_t setConst(T a, Vector<T,I>& X)
{
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::setConstKernel<<<grid, block>>>(a, X.device(), X.size());
    return cudaGetLastError();    
}

template <typename T, typename I>
inline cudaError_t linearSum(T a, const Vector<T,I>& X, T b, const Vector<T,I>& Y, Vector<T,I>& Z)
{
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::linearSumKernel<<<grid, block>>>(a, X.device(), b, Y.device(), Z.device(), X.size());
    return cudaGetLastError();    
}

template <typename T, typename I>
inline cudaError_t axpy(T a, const Vector<T,I>& X, Vector<T,I>& Y)
{
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::axpyKernel<<<grid, block>>>(a, X.device(), Y.device(), X.size());
    return cudaGetLastError();    
}

template <typename T, typename I>
inline cudaError_t prod(const Vector<T,I>& X, const Vector<T,I>& Y, Vector<T,I>& Z)
{
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::prodKernel<<<grid, block>>>(X.device(), Y.device(), Z.device(), X.size());
    return cudaGetLastError();    
}

template <typename T, typename I>
inline cudaError_t div(const Vector<T,I>& X, const Vector<T,I>& Y, Vector<T,I>& Z)
{
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::divKernel<<<grid, block>>>(X.device(), Y.device(), Z.device(), X.size());
    return cudaGetLastError();    
}

template <typename T, typename I>
inline cudaError_t scale(T const a, const Vector<T,I>& X, Vector<T,I>& Z)
{
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::scaleKernel<<<grid, block>>>(a, X.device(), Z.device(), X.size());
    return cudaGetLastError();    
}

template <typename T, typename I>
inline cudaError_t absVal(const Vector<T,I>& X, Vector<T,I>& Z)
{
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::absKernel<<<grid, block>>>(X.device(), Z.device(), X.size());
    return cudaGetLastError();    
}

template <typename T, typename I>
inline cudaError_t inv(const Vector<T,I>& X, Vector<T,I>& Z)
{
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::invKernel<<<grid, block>>>(X.device(), Z.device(), X.size());
    return cudaGetLastError();    
}

template <typename T, typename I>
inline cudaError_t addConst(T const a, const Vector<T,I>& X, Vector<T,I>& Z)
{
    // Set partitioning
    ThreadPartitioning<T, I>* p = X.partStream();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();

    math_kernels::addConstKernel<<<grid, block>>>(a, X.device(), Z.device(), X.size());
    return cudaGetLastError();    
}

    
template <typename T, typename I>
inline T dotProd(const nvec::Vector<T,I>& x, const nvec::Vector<T,I>& y)
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
    
template <typename T, typename I>
inline T maxNorm(const nvec::Vector<T,I>& x)
{
    // Reduction result storage on CPU
    T gpu_result = 0;
    
    // Set partitioning
    ThreadPartitioning<T, I>* p = x.partReduce();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();
    unsigned int shMemSize = p->getShMemSize();
    
    math_kernels::maxNormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), p->devBuffer(), x.size());
    
    unsigned int n = grid.x * grid.y * grid.z;
    while (n>2048)
    {
        // Recompute partitioning
        unsigned int grid2;
        unsigned int block2;
        unsigned int shMemSize2;
        p->setPartitioningReduce(n, grid2, block2, shMemSize2);
        
        // Rerun reduction kernel
        math_kernels::maxNormKernel<T,I><<< grid2, block2, shMemSize2 >>>(p->devBuffer(), p->devBuffer(), n);
        n = (n + 2*block2 - 1) / (2*block2);
    }
    
    // sum partial sums from each block on CPU
    // copy result from device to host
    p->copyFromDevBuffer(n);
    
    gpu_result = p->hostBuffer()[0];
    for (unsigned int i=1; i<n; i++)
    {
        if (p->hostBuffer()[i])
            gpu_result = p->hostBuffer()[i];
    }
    return gpu_result;
}
    
template <typename T, typename I>
inline T wrmsNorm(const nvec::Vector<T,I>& x, const nvec::Vector<T,I>& y)
{
    // Reduction result storage on CPU
    T gpu_result = 0;
    
    // Set partitioning
    ThreadPartitioning<T, I>* p = x.partReduce();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();
    unsigned int shMemSize = p->getShMemSize();
    
    math_kernels::wrmsNormKernel<T,I><<< grid, block, shMemSize >>>(x.device(), y.device(), p->devBuffer(), x.size());
    
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
    return sqrt(gpu_result/x.size());
}
    
template <typename T, typename I>
inline T findMin(const nvec::Vector<T,I>& x)
{
    // Reduction result storage on CPU
    T gpu_result = 0;
    T maxVal = std::numeric_limits<T>::max();
    
    // Set partitioning
    ThreadPartitioning<T, I>* p = x.partReduce();
    const dim3& grid       = p->getGrid();
    const dim3& block      = p->getBlock();
    unsigned int shMemSize = p->getShMemSize();
    
    math_kernels::findMinKernel<T,I><<< grid, block, shMemSize >>>(maxVal, x.device(), p->devBuffer(), x.size());
    
    unsigned int n = grid.x * grid.y * grid.z;
    while (n>2048)
    {
        // Recompute partitioning
        unsigned int grid2;
        unsigned int block2;
        unsigned int shMemSize2;
        p->setPartitioningReduce(n, grid2, block2, shMemSize2);
        
        // Rerun reduction kernel
        math_kernels::findMinKernel<T,I><<< grid2, block2, shMemSize2 >>>(maxVal, p->devBuffer(), p->devBuffer(), n);
        n = (n + 2*block2 - 1) / (2*block2);
    }
    
    // sum partial sums from each block on CPU
    // copy result from device to host
    p->copyFromDevBuffer(n);
    
    gpu_result = p->hostBuffer()[0];
    for (unsigned int i=1; i<n; i++)
    {
        if (p->hostBuffer()[i])
            gpu_result = p->hostBuffer()[i];
    }
    return gpu_result;
}
    

    

} // namespace nvec



#endif // _VECTOR_KERNELS_CUH_
