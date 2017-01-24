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

#ifndef _NVECTOR_RAJA_HPP_
#define _NVECTOR_RAJA_HPP_

#include <cstdlib>
#include <iostream>

#include <nvector/nvector_raja.h>

namespace rvec
{

template <typename T, typename I=int>
class Vector : public _N_VectorContent_Raja
{
public:
    Vector(I N) : size_(N), mem_size_(N*sizeof(T)), isClone_(false)
    {
        allocate();
    }
    
    /// Copy constructor does not copy values
    explicit Vector(const Vector& v) 
    : size_(v.size()), 
      mem_size_(size_*sizeof(T)), 
      isClone_(true) ///< temporary, will be removed!
    {
        allocate();
    }
    
    ~Vector()
    {
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
    
    void clear()
    {
        cudaError_t err;
        free(h_vec_);
        err = cudaFree(d_vec_);
        if(err != cudaSuccess)
            std::cout << "Failed to free device vector (error code " << err << ")!\n";        
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
    
private:
    I size_;
    I mem_size_;
    T* h_vec_;
    T* d_vec_;
    bool isClone_;    ///< temporary, will be removed! 
};



    

} // namespace rvec

// Vector extractor
inline rvec::Vector<double, long int>* extract_raja(N_Vector v)
{
    return static_cast<rvec::Vector<double, long int>*>(v->content);
}




#endif // _NVECTOR_RAJA_HPP_
