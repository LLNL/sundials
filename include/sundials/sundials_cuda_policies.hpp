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
 * This header files defines the CudaExecPolicy classes which
 * are utilized to determine CUDA kernel launch paramaters.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_CUDAEXECPOLICIES_HPP
#define _SUNDIALS_CUDAEXECPOLICIES_HPP

#include <stdio.h>
#include <cuda_runtime.h>

namespace sundials
{

class CudaExecPolicy
{
public:
  virtual size_t gridSize(size_t numWorkElements = 0, size_t blockDim = 0) const = 0;
  virtual size_t blockSize(size_t numWorkElements = 0, size_t gridDim = 0) const = 0; 
  virtual cudaStream_t stream() const = 0;
  virtual CudaExecPolicy* clone() const = 0;
};


/* 
 * A simple execution policy for streaming kernels.
 * Will map each work element to a thread.
 * The number of theads per block is fixed.
 */
class CudaStreamingExecPolicy : public CudaExecPolicy
{
public:
  CudaStreamingExecPolicy(const size_t blockDim, const cudaStream_t stream = 0)
    : blockDim_(blockDim), stream_(stream)
  {}

  CudaStreamingExecPolicy(const CudaStreamingExecPolicy& ex)
    : blockDim_(ex.blockDim_), stream_(ex.stream_)
  {} 

  virtual size_t gridSize(size_t numWorkElements = 0, size_t blockDim = 0) const
  {
    if (numWorkElements == 0)
    {
      return 0;
    }
    return (numWorkElements + blockSize() - 1) / blockSize();
  }

  virtual size_t blockSize(size_t numWorkElements = 0, size_t gridDim = 0) const
  {
    return blockDim_;
  }

  virtual cudaStream_t stream() const
  {
    return stream_;
  }

  virtual CudaExecPolicy* clone() const
  {
    return static_cast<CudaExecPolicy*>(new CudaStreamingExecPolicy(*this));
  }

private:
  const cudaStream_t stream_;
  const size_t blockDim_;
};


/* 
 * A simple execution policy for kernels that do reductions across thread blocks.
 * The number of theads per block is fixed.
 */
class CudaReductionExecPolicy : public CudaExecPolicy
{
public:
  CudaReductionExecPolicy(const size_t blockDim, const cudaStream_t stream = 0)
    : blockDim_(blockDim), stream_(stream)
  {}

  CudaReductionExecPolicy(const CudaReductionExecPolicy& ex)
    : blockDim_(ex.blockDim_), stream_(ex.stream_)
  {} 

  virtual size_t gridSize(size_t numWorkElements = 0, size_t blockDim = 0) const
  {
    if (numWorkElements == 0)
    {
      return 0;
    }
    return (numWorkElements + (blockSize() * 2 - 1)) / (blockSize() * 2);
  }

  virtual size_t blockSize(size_t numWorkElements = 0, size_t gridDim = 0) const
  {
    return blockDim_;
  }
  
  virtual cudaStream_t stream() const
  {
    return stream_;
  }

  virtual CudaExecPolicy* clone() const
  {
    return static_cast<CudaExecPolicy*>(new CudaReductionExecPolicy(*this));
  }

private:
  const cudaStream_t stream_;
  const size_t blockDim_;
};

} // namespace sundials

typedef sundials::CudaExecPolicy SUNCudaExecPolicy;
typedef sundials::CudaStreamingExecPolicy SUNCudaStreamingExecPolicy;
typedef sundials::CudaReductionExecPolicy SUNCudaReductionExecPolicy;
  
#endif