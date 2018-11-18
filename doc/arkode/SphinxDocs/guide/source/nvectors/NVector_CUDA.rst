..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _NVectors.CUDA:

The NVECTOR_CUDA Module
======================================

The NVECTOR_CUDA module is an experimental implementation of
``N_Vector`` in CUDA language.  It allows for SUNDIALS vector kernels
to run on GPU devices. It is intended for users who are already
familiar with CUDA and GPU programming.  Building this vector  
module requires CUDA compiler and, by extension, C++ compiler. Class ``Vector`` 
in namespace ``suncudavec`` manages vector data layout. 

.. code-block:: c++

   template <class T, class I>
   class Vector {
     I size_;
     I mem_size_;
     T* h_vec_;
     T* d_vec_;
     StreamPartitioning<T, I>* partStream_;
     ReducePartitioning<T, I>* partReduce_;
     bool ownPartitioning_;
  
     ...
   };

The class members are vector size (length), size of the vector data memory block, pointers
to vector data on the host and the device, pointers to classes ``StreamPartitioning``
and ``ReducePartitioning``, which handle thread partitioning for streaming and 
reduction vector kernels, respectively, and the boolean flag that signals if the
vector owns thread partitioning. The class ``Vector`` inherits from empty structure


.. code-block:: c++

   struct _N_VectorContent_Cuda {
   };

to interface the C++ class with ``N_Vector`` C code. When
instantiated, the class ``Vector`` will allocate memory on both, host
and device. Due to rapid progress in of CUDA development, we expect
that ``suncudavec::Vector`` class will change frequently in the future
SUNDIALS releases. The code is structured so that it can tolerate
significant changes in the ``suncudavec::Vector`` class without
requiring changes to user API. 

The header file to be included when using this module is
``nvector/nvector_cuda.h``.  Unlike other native SUNDIALS vector
types, NVECTOR_CUDA does not provide macros to access its member
variables. Note that NVECTOR_CUDA requires SUNDIALS to be built with
MPI support. 


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
     


.. c:function:: N_Vector N_VNew_Cuda(sunindextype vec_length)

   This function creates and allocates memory for a CUDA
   ``N_Vector``. The memory is allocated on both, host and device. Its
   only argument is the vector length. 


.. c:function:: N_Vector N_VNewEmpty_Cuda(sunindextype vec_length)

   This function creates a new ``N_Vector`` wrapper with the pointer
   to the wrapped CUDA vector set to ``NULL``.  It is used by
   :c:func:`N_VNew_Cuda()`, :c:func:`N_VMake_Cuda()`, and
   :c:func:`N_VClone_Cuda()` implementations. 

      
.. c:function:: N_Vector N_VMake_Cuda(N_VectorContent_Cuda c)

   This function creates and allocates memory for an NVECTOR_CUDA
   wrapper around a user-provided ``suncudavec::Vector`` class.  
   Its only argument is of type ``N_VectorContent_Cuda``, which
   is the pointer to the class.

 
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


.. c:function:: sunindextype N_VGetLength_Cuda(N_Vector v)

   This function returns the length of the vector.


.. c:function:: realtype* N_VGetHostArrayPointer_Cuda(N_Vector v)

   This function returns pointer to the vector data on the host.


.. c:function:: realtype* N_VGetDeviceArrayPointer_Cuda(N_Vector v)

   This function returns pointer to the vector data on the device.


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
