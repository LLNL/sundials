..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2020, Lawrence Livermore National Security
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

The NVECTOR_CUDA module is an experimental NVECTOR implementation in the CUDA language.
The module allows for SUNDIALS vector kernels to run on GPU devices. It is intended for users
who are already familiar with CUDA and GPU programming. Building this vector
module requires a CUDA compiler and, by extension, a C++ compiler. The class ``Vector``
in the namespace ``suncudavec`` manages the vector data layout:

.. code-block:: c++

  template <class T, class I>
  class Vector {
    I size_;
    I mem_size_;
    T* h_vec_;
    T* d_vec_;
    ThreadPartitioning<T, I>* partStream_;
    ThreadPartitioning<T, I>* partReduce_;
    bool ownPartitioning_;
    bool ownData_;
    bool managed_mem_;
    ...
  };

The class members are vector size (length), size of the vector data memory
block, pointers to vector data on the host and the device, pointers
to ``ThreadPartitioning`` implementations that handle thread partitioning for
streaming and reduction vector kernels, a boolean flag that signals if the
vector owns the thread partitioning, a boolean flag that signals if the vector
owns the data, and a boolean flag that signals if managed memory is used for the
data arrays. The class ``Vector`` inherits from the empty structure

.. code-block:: c++

   struct _N_VectorContent_Cuda {};

to interface the C++ class with the NVECTOR C code. Due to the rapid progress
of CUDA development, we expect that the ``suncudavec::Vector`` class will
change frequently in future SUNDIALS releases. The code is structured so that
it can tolerate significant changes in the ``suncudavec::Vector`` class without
requiring changes to the user API.

When instantiated with ``N_VNew_Cuda``, the class ``Vector`` will allocate
memory on both the host and the device. Alternatively, a user can provide host
and device data arrays by using the ``N_VMake_Cuda`` constructor. To use CUDA
managed memory, the constructors ``N_VNewManaged_Cuda`` and
``N_VMakeManaged_Cuda`` are provided. Details on each of these constructors
are provided below.

To use the NVECTOR_CUDA module, include ``nvector_cuda.h`` and link to
the library ``libsundials_nveccuda.lib``. The extension, ``.lib``, is
typically ``.so`` for shared libraries and ``.a`` for static libraries.


NVECTOR_CUDA functions
-----------------------------------

Unlike other native SUNDIALS vector types, the NVECTOR_CUDA module does not
provide macros to access its member variables. Instead, user should use the
accessor functions:


.. c:function:: realtype* N_VGetHostArrayPointer_Cuda(N_Vector v)

   This function returns pointer to the vector data on the host.


.. c:function:: realtype* N_VGetDeviceArrayPointer_Cuda(N_Vector v)

   This function returns pointer to the vector data on the device.


.. c:function:: booleantype N_VIsManagedMemory_Cuda(N_Vector v)

   This function returns a boolean flag indiciating if the vector
   data array is in managed memory or not.


The NVECTOR_CUDA module defines implementations of all standard vector
operations defined in the sections :ref:`NVectors.Ops`, :ref:`NVectors.FusedOps`,
:ref:`NVectors.ArrayOps`, and :ref:`NVectors.LocalOps`, except for
``N_VSetArrayPointer``, and, if using unmanaged memory, ``N_VGetArrayPointer``.
As such, this vector can only be used with SUNDIALS Fortran interfaces, and the
SUNDIALS direct solvers and preconditioners when using managed memory.
The NVECTOR_CUDA module provides separate functions to access data on the host
and on the device for the unmanaged memory use case. It also provides methods for
copying from the host to the device and vice versa. Usage examples of NVECTOR_CUDA
are provided in example programs for CVODE [HSR2017]_.

The names of vector operations are obtained from those in the sections
:ref:`NVectors.Ops`, :ref:`NVectors.FusedOps`, :ref:`NVectors.ArrayOps`, and
:ref:`NVectors.LocalOps` by appending the suffix ``_Cuda``
(e.g. ``N_VDestroy_Cuda``).  The module NVECTOR_CUDA provides the
following additional user-callable routines:



.. c:function:: N_Vector N_VNew_Cuda(sunindextype length)

   This function creates and allocates memory for a CUDA ``N_Vector``.
   The vector data array is allocated on both the host and device.


.. c:function:: N_Vector N_VNewManaged_Cuda(sunindextype vec_length)

   This function creates and allocates memory for a CUDA
   ``N_Vector``. The vector data array is allocated in managed memory.


.. c:function:: N_Vector N_VNewEmpty_Cuda(sunindextype vec_length)

   This function creates a new ``N_Vector`` wrapper with the pointer
   to the wrapped CUDA vector set to ``NULL``.  It is used by
   :c:func:`N_VNew_Cuda()`, :c:func:`N_VMake_Cuda()`, and
   :c:func:`N_VClone_Cuda()` implementations.


.. c:function:: N_Vector N_VMake_Cuda(sunindextype vec_length, realtype \*h_vdata, realtype \*d_vdata)


   This function creates a CUDA ``N_Vector`` with user-supplied vector data arrays
   for the host and the device.


.. c:function:: N_Vector N_VMakeManaged_Cuda(sunindextype vec_length, realtype \*vdata)

   This function creates a CUDA ``N_Vector`` with a user-supplied
   managed memory data array.


.. c:function:: N_Vector N_VMakeWithManagedAllocator_Cuda(sunindextype length, void* (\*allocfn)(size_t size), void (\*freefn)(void* ptr))

   This function creates a CUDA ``N_Vector`` with a user-supplied memory allocator.
   It requires the user to provide a corresponding free function as well.
   The memory allocated by the allocator function must behave like CUDA managed memory.
   


The module NVECTOR_CUDA also provides the following user-callable routines:

.. c:function:: void N_VSetCudaStream_Cuda(N_Vector v, cudaStream_t \*stream)

   This function sets the CUDA stream that all vector kernels will be launched on.
   By default an NVECTOR_CUDA uses the default CUDA stream.

   *Note: All vectors used in a single instance of a SUNDIALS solver must
   use the same CUDA stream, and the CUDA stream must be set prior to
   solver initialization. Additionally, if manually instantiating the stream and
   reduce ``ThreadPartitioning`` of a ``suncudavec::Vector``, ensure that they
   use the same CUDA stream.*


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

.. c:function:: int N_VEnableFusedOps_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the CUDA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the CUDA vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the CUDA vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the CUDA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the CUDA vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When there is a need to access components of an ``N_Vector_Cuda``, ``v``,
  it is recommeded to use functions :c:func:`N_VGetDeviceArrayPointer_Cuda()` or
  :c:func:`N_VGetHostArrayPointer_Cuda()`. However, when using managed memory,
  the function :c:func:`N_VGetArrayPointer` may also be used.

* To maximize efficiency, vector operations in the NVECTOR_CUDA implementation
  that have more than one ``N_Vector`` argument do not check for
  consistent internal representations of these vectors. It is the user's
  responsibility to ensure that such routines are called with ``N_Vector``
  arguments that were all created with the same internal representations.
