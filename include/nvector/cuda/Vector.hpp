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


/**
 * Vector class
 */

#ifndef _NVECTOR_HPP_
#define _NVECTOR_HPP_

#include <cstdlib>
#include <iostream>

#include <cuda_runtime.h>
#include "ThreadPartitioning.hpp"

#include <nvector/nvector_cuda.h>

namespace nvec
{

template <typename T, typename I=int>
class Vector : public _N_VectorContent_Cuda
{
public:
    Vector(I N)
    : size_(N),
      mem_size_(N*sizeof(T)),
      ownPartitioning_(true),
      isClone_(false),
      ownHostData_(true),
      ownDevData_(true)
    {
        // Set partitioning
        partStream_ = new ThreadPartitioning<T, I>(ThreadPartitioning<T,I>::STREAM, N);
        partReduce_ = new ThreadPartitioning<T, I>(ThreadPartitioning<T,I>::REDUCTION, N);

        allocate();
    }

    /// This is temporary solution will be removed.
    /// Use only if you know what you are doing.
    Vector(I N, T* data)
    : size_(N),
      mem_size_(N*sizeof(T)),
      ownPartitioning_(true),
      isClone_(false),
      ownHostData_(false),
      ownDevData_(false)
    {
        // Set partitioning
        partStream_ = new ThreadPartitioning<T, I>(ThreadPartitioning<T,I>::STREAM, N);
        partReduce_ = new ThreadPartitioning<T, I>(ThreadPartitioning<T,I>::REDUCTION, N);

        if (data == NULL)
        {
          h_vec_ = NULL;
          d_vec_ = NULL;
        }
        else
        {
          // Sets pointer to device memory and allocates host memory of size N
          allocate(data);
          ownHostData_ = true;
        }
    }

    /// Copy constructor does not copy values
    explicit Vector(const Vector& v)
    : size_(v.size()),
      mem_size_(size_*sizeof(T)),
      partStream_(v.partStream_),
      partReduce_(v.partReduce_),
      ownPartitioning_(false),
      isClone_(true), ///< temporary, will be removed!
      ownHostData_(true),
      ownDevData_(true)
    {
        allocate();
    }

    ~Vector()
    {
        if (ownPartitioning_)
        {
            delete partReduce_;
            delete partStream_;
        }
        clear();
    }


    void allocate()
    {
        cudaError_t err;
        h_vec_ = static_cast<T*>(malloc(mem_size_));
        if(h_vec_ == NULL)
            std::cout << "Failed to allocate host vector!\n";
        err = cudaMalloc((void**) &d_vec_, mem_size_);
        if(err != cudaSuccess)
            std::cout << "Failed to allocate device vector (error code " << err << ")!\n";
    }

    /// This is temporary solution! Use only if you know what you are doing.
    void allocate(T* data)
    {
      // Allocate host data
      h_vec_ = static_cast<T*>(malloc(mem_size_));
      if(h_vec_ == NULL)
        std::cout << "Failed to allocate host vector!\n";

      // Set pointer to device data
      d_vec_ = data;
    }

    void clear()
    {
      if(ownHostData_)
      {
        free(h_vec_);
      }
      if (ownDevData_)
      {
        cudaError_t err = cudaFree(d_vec_);
        if(err != cudaSuccess)
          std::cout << "Failed to free device vector (error code " << err << ")!\n";
      }
    }

    int size() const
    {
        return size_;
    }

    bool isClone()
    {
        return isClone_;
    }

    T* host()
    {
        return h_vec_;
    }

    const T* host() const
    {
        return h_vec_;
    }

    T* device()
    {
        return d_vec_;
    }

    const T* device() const
    {
        return d_vec_;
    }

    void copyToDev()
    {
        cudaError_t err = cudaMemcpy(d_vec_, h_vec_, mem_size_, cudaMemcpyHostToDevice);
        if(err != cudaSuccess)
            std::cout << "Failed to copy vector from host to device (error code " << err << ")!\n";
    }

    void copyFromDev()
    {
        cudaError_t err = cudaMemcpy(h_vec_, d_vec_, mem_size_, cudaMemcpyDeviceToHost);
        if(err != cudaSuccess)
            std::cout << "Failed to copy vector from device to host (error code " << err << ")!\n";
    }

    /// This is dangerous function and is here only temporary. It will be removed.
    void setFromHost(T* h_vec)
    {
      if (ownHostData_)
      {
        free(h_vec_);
        ownHostData_ = false;
      }
      h_vec_ = h_vec;
      copyToDev();
    }

    /// This is dangerous function and is here only temporary. It will be removed.
    void setFromDevice(T* d_vec)
    {
      if(ownDevData_)
      {
        cudaError_t err = cudaFree(d_vec_);
        if(err != cudaSuccess)
          std::cout << "Failed to free device vector (error code " << err << ")!\n";
        ownDevData_ = false;
      }
      d_vec_ = d_vec;
      // Do not copy to host
    }

    ThreadPartitioning<T, I>* partStream()
    {
        return partStream_;
    }

    ThreadPartitioning<T, I>* partStream() const
    {
        return partStream_;
    }

    ThreadPartitioning<T, I>* partReduce()
    {
        return partReduce_;
    }

    ThreadPartitioning<T, I>* partReduce() const
    {
        return partReduce_;
    }

private:
    I size_;
    I mem_size_;
    T* h_vec_;
    T* d_vec_;
    ThreadPartitioning<T, I>* partStream_;
    ThreadPartitioning<T, I>* partReduce_;
    bool ownPartitioning_;
    bool isClone_;    ///< temporary, will be removed!
    bool ownHostData_;
    bool ownDevData_;
};





// Vector extractor
inline nvec::Vector<double, long int>* extract(N_Vector v) 
{ 
    return static_cast<nvec::Vector<double, long int>*>(v->content); 
}

} // namespace nvec




#endif // _NVECTOR_HPP_
