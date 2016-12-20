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
 */



#ifndef _THREAD_PARTITIONING_HPP_
#define _THREAD_PARTITIONING_HPP_

#include <iostream>
#include <cuda_runtime.h>

namespace nvec
{

/// Finds next power of two larger than value
template<class T>
T nextPow2(T value)
{
    //std::cout << "Integer size is " << sizeof(T) << " bytes.\n"; 
    --value;
    for(size_t i = 1; i < sizeof(T) * CHAR_BIT; i*=2)
        value |= value >> i;
    return ++value;
}

 
template<class T, class I=int>
class ThreadPartitioning
{
public:
    enum KernelTypeEnumeration {STREAM=0, REDUCTION};

    /// Constructor
    ThreadPartitioning(KernelTypeEnumeration t, I N) 
    : grid_(1,1,1), 
      block_(1,1,1), 
      shMemSize_(0), 
      maxThreads_(256), 
      d_buffer_(NULL), 
      h_buffer_(NULL),
      bufferSize_(0),
      kernelType_(t)
    {
        setPartitioning(N);
    }
    
    /// Copy constructor
    explicit ThreadPartitioning(ThreadPartitioning<T,I>& p) 
    : grid_(p.grid_), 
      block_(p.block_), 
      shMemSize_(p.shMemSize_), 
      maxThreads_(p.maxThreads_), 
      d_buffer_(NULL), 
      h_buffer_(NULL),
      bufferSize_(0),
      kernelType_(p.kernelType_)
    {
        allocateBuffer();
    }
    


    ~ThreadPartitioning()
    {    
        cudaError_t err;
        if (h_buffer_ != NULL) // this if statement is a bad idea, think of something better...
            free(h_buffer_);
        if (d_buffer_ != NULL)
        {
            err = cudaFree(d_buffer_);
            if(err != cudaSuccess)
                std::cout << "Failed to free device vector (error code " << err << ")!\n";
        }
    }

    const dim3& getGrid() const
    {
        return grid_;
    }
    
    const dim3& getBlock() const
    {
        return block_;
    }
    
    unsigned int getShMemSize()
    {
        return shMemSize_;
    }
        
    T* devBuffer()
    {
        return d_buffer_;   
    }
    
    const T* devBuffer() const
    {
        return d_buffer_;   
    }
    
    T* hostBuffer()
    {
        return h_buffer_;   
    }
    
    const T* hostBuffer() const
    {
        return h_buffer_;   
    }
    
    void copyFromDevBuffer(unsigned int n) const
    {
        cudaError_t err = cudaMemcpy(h_buffer_, d_buffer_, n*sizeof(T), cudaMemcpyDeviceToHost);
        if(err != cudaSuccess)
            std::cout << "Failed to copy vector from device to host (error code " << err << ")!\n";
    }

    int setPartitioning(I N)
    {
        switch (kernelType_)
        {
        case STREAM:
            setPartitioningStream(N);
            break;
        case REDUCTION:
            setPartitioningReduce(N);
            break;
        }
        return 0;
    }
    
    int setPartitioningStream(I N)
    {
       block_.x = maxThreads_;
       grid_.x = (N + maxThreads_ - 1) / maxThreads_;
       
       return 0;
    }
    
    int setPartitioningReduce(I N)
    {
        block_.x = (N < maxThreads_*2) ? nextPow2((N + 1)/ 2) : maxThreads_;
        grid_.x  = (N + (block_.x * 2 - 1)) / (block_.x * 2);
        
        // when there is only one warp per block, we need to allocate two warps
        // worth of shared memory so that we don't index shared memory out of bounds
        shMemSize_ = (block_.x <= 32) ? 2 * block_.x * sizeof(T) : block_.x * sizeof(T);
        
        // allocate reduction buffer
        allocateBuffer();
        return 0;
    }
    
    int setPartitioningReduce(I N, unsigned int& grid, unsigned int& block, unsigned int& shMemSize)
    {
        block = (N < maxThreads_*2) ? nextPow2((N + 1)/ 2) : maxThreads_;
        grid  = (N + (block * 2 - 1)) / (block * 2);
        
        // when there is only one warp per block, we need to allocate two warps
        // worth of shared memory so that we don't index shared memory out of bounds
        shMemSize = (block <= 32) ? 2 * block * sizeof(T) : block * sizeof(T);

        return 0;
    }
    
    int allocateBuffer()
    {
        // Streaming kernel does not need a buffer
        if (kernelType_ == STREAM)
            return 0;
        
        bufferSize_ = grid_.x * sizeof(T);
        h_buffer_ = static_cast<T*>(malloc(bufferSize_));
        if(h_buffer_ == NULL)
            std::cout << "Failed to allocate host vector!\n";

        cudaError_t err;
        err = cudaMalloc((void**) &d_buffer_, bufferSize_); 
        if(err != cudaSuccess)
            std::cout << "Failed to allocate device vector (error code " << err << ")!\n";
        
        return 0;
    }
    
    unsigned int buffSize()
    {
        return bufferSize_;
    }
   
private:
    dim3 grid_;
    dim3 block_;
    unsigned int shMemSize_;  ///< Shared memory size
    unsigned int maxThreads_; ///< Number of threads per block
    T* d_buffer_;
    T* h_buffer_;
    unsigned int bufferSize_;
    KernelTypeEnumeration kernelType_;
};


} // namespace nvec

#endif // _THREAD_PARTITIONING_HPP_
