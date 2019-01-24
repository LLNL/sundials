..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2019, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3


.. _NVectors.CUDA:

The NVECTOR_CUDA Module
======================================

The NVECTOR_CUDA module is an experimental NVECTOR implementation
in the CUDA language. This module allows for SUNDIALS vector kernels
to run on GPU devices. It is intended for users who are already
familiar with CUDA and GPU programming.  Building this vector  
module requires a CUDA compiler and, by extension, C++ compiler.
The class ``Vector``  in the namespace ``suncudavec`` manages
the vector data layout. 

.. code-block:: c++

  template <class T, class I>
  class Vector {
    I size_;
    I mem_size_;
    I global_size_;
    T* h_vec_;
    T* d_vec_;
    ThreadPartitioning<T, I>* partStream_;
    ThreadPartitioning<T, I>* partReduce_;
    bool ownPartitioning_;
    bool ownData_;
    bool managed_mem_;
    SUNMPI_Comm comm_;
    ...
  };

The class members are vector size (length), size of the vector data memory block, pointers
to vector data on the host and the device, pointers to classes ``StreamPartitioning``
and ``ReducePartitioning``, which handle thread partitioning for streaming and 
reduction vector kernels, respectively, a boolean flag that signals if the
vector owns the thread partitioning, a boolean flag that signals if the vector
owns the data, a boolean flag that signals if managed memory is used for the
data arrays, and the MPI communicator. he class ``Vector`` inherits from empty
structure

.. code-block:: c++

   struct _N_VectorContent_Cuda {
   };

to interface the C++ class with ``N_Vector`` C code.
Due to rapid progress in of CUDA development, we expect
that ``suncudavec::Vector`` class will change frequently in the future
SUNDIALS releases. The code is structured so that it can tolerate
significant changes in the ``suncudavec::Vector`` class without
requiring changes to user API. 

When instantiated, the class ``Vector`` will allocate memory on both, host
and device by default. Optionally, managed memory can be allocated instead
(see ``N_VNewManaged_Cuda``), or a user can provide data arrays
(see ``N_VMake_Cuda`` and ``N_VMakeManaged_Cuda``).

The NVECTOR_CUDA module can be utilized for single-node parallelism or in
a distributed context with MPI. The header file to include when using this
module for single-node parallelism is ``nvector_cuda.h``. The header file
to include when using this module in the distributed case is
``nvector_mpicuda.h``. The installed module libraries to link to are
``libsundials_nveccuda.lib`` in the single-node case, or
``libsundials_nvecmpicuda.lib`` in the distributed case. Only one one of
these libraries may be linked to when creating an executable or library.
SUNDIALS must be built with MPI support if the distributed library is desired.

Unlike other native SUNDIALS vector types, the NVECTOR_CUDA module does not
provide macros to access its member variables. Instead, user should use the
accessor functions:



.. c:function:: sunindextype N_VGetLength_Cuda(N_Vector v)

   This function returns the global length of the vector.


.. c:function:: sunindextype N_VGetLocalLength_Cuda(N_Vector v)

   This function returns the local length of the vector.

   Note: This function is for use in a *distributed* context
   and is defined in the header ``nvector_mpicuda.h`` and the
   library to link to is ``libsundials_nvecmpicuda.lib``.


.. c:function:: realtype* N_VGetHostArrayPointer_Cuda(N_Vector v)

   This function returns pointer to the vector data on the host.


.. c:function:: realtype* N_VGetDeviceArrayPointer_Cuda(N_Vector v)

   This function returns pointer to the vector data on the device.


.. c:function:: MPI_Comm N_VGetMPIComm_Cuda(N_Vector v)

   This function returns the MPI communicator for the vector.

   Note: This function is for use in a *distributed* context
   and is defined in the header ``nvector_mpicuda.h`` and the
   library to link to is ``libsundials_nvecmpicuda.lib``.


.. c:function:: booleantype N_VIsManagedMemory_Cuda(N_Vector v)

   This function returns a boolean flag indiciating if the vector
   data array is in managed memory or not.


The NVECTOR_CUDA module defines implementations of all standard vector
operations defined in the sections :ref:`NVectors.Ops`,
:ref:`NVectors.FusedOps`, and  :ref:`NVectors.ArrayOps`, except for
``N_VGetArrayPointer`` and ``N_VSetArrayPointer``.  As such, this
vector cannot be used with SUNDIALS Fortran interfaces, nor with
SUNDIALS direct solvers and preconditioners. This support will be
added in subsequent SUNDIALS releases.  The NVECTOR_CUDA module
provides separate functions to access data on the host and on the
device. It also provides methods for copying from the host to the
device and vice versa. Usage examples of NVECTOR_CUDA are provided in  
example programs for CVODE [HSR2017]_.

The names of vector operations are obtained from those in the sections
:ref:`NVectors.Ops`, :ref:`NVectors.FusedOps` and
:ref:`NVectors.ArrayOps` by appending the suffix ``_Cuda``
(e.g. ``N_VDestroy_Cuda``).  The module NVECTOR_CUDA provides the
following additional user-callable routines:



.. c:function:: N_Vector N_VNew_Cuda(sunindextype length)
                N_Vector N_VNew_Cuda(MPI_Comm comm, sunindextype local_length, sunindextype global_length)

   This function creates and allocates memory for a CUDA ``N_Vector``.
   The vector data array is allocated on both the host and device.

   In the *single-node* setting, the only input is the vector length.
   This constructor is defined in the header ``nvector_cuda.h`` and
   the library to link to is is ``libsundials_nveccuda.lib``.
 
   When used in a *distributed* context with MPI, the arguments are the
   MPI communicator, the local vector length, and the global vector length.
   This constructor is defined in the header ``nvector_mpicuda.h`` and
   the library to link to is ``libsundials_nvecmpicuda.lib``.


.. c:function:: N_Vector N_VNewManaged_Cuda(sunindextype vec_length)
                N_Vector N_VNewManaged_Cuda(MPI_Comm comm, sunindextype local_length, sunindextype global_length)

   This function creates and allocates memory for a CUDA
   ``N_Vector``. The vector data array is allocated in managed memory.
   
   When used in the *single-node* setting, the only input is the vector length.
   this constructor is defined in the header ``nvector_cuda.h`` and
   the library to link to is is ``libsundials_nveccuda.lib``.
   
   When used in a *distributed* context with MPI, the arguments are the
   MPI communicator, the local vector length, and the global vector length.
   This constructor is defined in the header ``nvector_mpicuda.h`` and
   the library to link to is ``libsundials_nvecmpicuda.lib``.


.. c:function:: N_Vector N_VNewEmpty_Cuda(sunindextype vec_length)

   This function creates a new ``N_Vector`` wrapper with the pointer
   to the wrapped CUDA vector set to ``NULL``.  It is used by
   :c:func:`N_VNew_Cuda()`, :c:func:`N_VMake_Cuda()`, and
   :c:func:`N_VClone_Cuda()` implementations. 
   

.. c:function:: N_Vector N_VMake_Cuda(sunindextype vec_length, realtype *h_vdata, realtype *d_vdata)
                N_Vector N_VMake_Cuda(MPI_Comm comm, sunindextype global_length, sunindextype local_length, realtype *h_vdata, realtype *d_vdata)

   
   This function creates a CUDA ``N_Vector`` with user-supplied vector data arrays
   for the host and the device.
   
   When used in the *single-node* setting, the arguments are the
   the vector length, the host data array, and the device data array.
   This constructor is defined in the header ``nvector_cuda.h`` and
   the library to link to is is ``libsundials_nveccuda.lib``.

   When used in a *distributed* context with MPI, the arguments are the
   MPI communicator, the global vector length, the local vector length,
   the host data array, the device data array.
   This constructor is defined in the header ``nvector_mpicuda.h`` and
   the library to link to is ``libsundials_nvecmpicuda.lib``.


.. c:function:: N_Vector N_VMakeManaged_Cuda(sunindextype vec_length, realtype *vdata)
                N_Vector N_VMakeManaged_Cuda(MPI_Comm comm, sunindextype global_length, sunindextype local_length, realtype *vdata)

   This function creates a CUDA ``N_Vector`` with a user-supplied
   managed memory data array.

   When used in the *single-node* setting, the arguments are the
   the vector length, and the managed data array. This constructor
   is defined in the header ``nvector_cuda.h`` and
   the library to link to is is ``libsundials_nveccuda.lib``.

   When used in a *distributed* context with MPI, the arguments are the
   MPI communicator, the global vector length, the local vector length,
   the managed data array. This constructor is defined in the header
   ``nvector_mpicuda.h`` and the library to link to is
   ``libsundials_nvecmpicuda.lib``.


.. c:function:: N_Vector* N_VCloneVectorArray_Cuda(int count, N_Vector w)

   This function creates (by cloning) an array of *count* NVECTOR_CUDA
   vectors. 


.. c:function:: N_Vector* N_VCloneVectorArrayEmpty_Cuda(int count, N_Vector w)

   This function creates (by cloning) an array of *count* NVECTOR_CUDA
   vectors, each with pointers to CUDA vectors set to ``NULL``. 


.. c:function:: void N_VDestroyVectorArray_Cuda(N_Vector* vs, int count)
  
   This function frees memory allocated for the array of *count*
   variables of type ``N_Vector`` created with
   :c:func:`N_VCloneVectorArray_Cuda()` or with
   :c:func:`N_VCloneVectorArrayEmpty_Cuda()`. 


.. c:function:: realtype* N_VCopyToDevice_Cuda(N_Vector v)

   This function copies host vector data to the device.


.. c:function:: realtype* N_VCopyFromDevice_Cuda(N_Vector v)

   This function copies vector data from the device to the host.


.. c:function:: void N_VPrint_Cuda(N_Vector v)

   This function prints the content of a CUDA vector to ``stdout``.


.. c:function:: void N_VPrintFile_Cuda(N_Vector v, FILE *outfile)

   This function prints the content of a CUDA vector to ``outfile``.


By default all fused and vector array operations are disabled in the NVECTOR_CUDA
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with :c:func:`N_VNew_Cuda`, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using :c:func:`N_VClone`. This guarantees the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned from while vectors created with
:c:func:`N_VNew_Cuda` will have the default settings for the NVECTOR_CUDA module.

.. c:function:: void N_VEnableFusedOps_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.
   
.. c:function:: void N_VEnableLinearCombination_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableScaleAddMulti_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the CUDA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableDotProdMulti_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the CUDA vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: void N_VEnableLinearSumVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableScaleVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableConstVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the CUDA vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: void N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the CUDA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the CUDA vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When there is a need to access components of an ``N_Vector_Cuda``, ``v``, 
  it is recommeded to use functions :c:func:`N_VGetDeviceArrayPointer_Cuda()` or 
  :c:func:`N_VGetHostArrayPointer_Cuda()`.        

* To maximize efficiency, vector operations in the NVECTOR_CUDA implementation
  that have more than one ``N_Vector`` argument do not check for
  consistent internal representations of these vectors. It is the user's 
  responsibility to ensure that such routines are called with ``N_Vector``
  arguments that were all created with the same internal representations.
