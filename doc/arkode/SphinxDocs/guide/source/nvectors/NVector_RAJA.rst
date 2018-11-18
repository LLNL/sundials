..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _NVectors.RAJA:

The NVECTOR_RAJA Module
======================================

The NVECTOR_RAJA module is an experimental implementation of
``N_Vector`` using the RAJA hardware abstraction layer
`https://software.llnl.gov/RAJA/ <https://software.llnl.gov/RAJA/>`_.
In this implementation, RAJA allows for SUNDIALS vector kernels to run
on GPU devices. The module is intended for users who are already
familiar with RAJA and GPU programming. Building this vector module
requires a C++11 compliant compiler and a CUDA software development
toolkit.  Besides the CUDA backend, RAJA has other backends such as
serial, OpenMP and OpenAC. These backends are not used in this SUNDIALS release. 
Class ``Vector`` in namespace ``sunrajavec`` manages the vector data layout:

.. code-block:: c++

   template <class T, class I>
   class Vector {
     I size_;
     I mem_size_;
     T* h_vec_;
     T* d_vec_;
  
     ...
   };

The class members are: vector size (length), size of the vector data
memory block, and pointers to vector data on the host and on the
device. The class ``Vector`` inherits from an empty structure 

.. code-block:: c++

   struct _N_VectorContent_Raja {
   };

to interface the C++ class with the ``N_Vector`` C code. When
instantiated, the class ``Vector`` will allocate memory on both the host
and the device. Due to the rapid progress of RAJA development, we expect
that the ``sunrajavec::Vector`` class will change frequently in the future
SUNDIALS releases. The code is structured so that it can tolerate
significant changes in the ``sunrajavec::Vector`` class without
requiring changes to the user API. 

The header file to be included when using this module is
``nvector/nvector_raja.h``.  Unlike other native SUNDIALS vector
types, NVECTOR_RAJA does not provide macros to access its member
variables. Note that NVECTOR_RAJA requires SUNDIALS to be built with
MPI support. 


The NVECTOR_RAJA module defines the implementations of all vector
operations listed in the sections :ref:`NVectors.Ops`,
:ref:`NVectors.FusedOps` and :ref:`NVectors.ArrayOps`, except for 
``N_VDotProdMulti``, ``N_VWrmsNormVectorArray``,
``N_VWrmsNormMaskVectorArray`` as support for arrays of reduction
vectors is not yet supported in RAJA.  These functions will be added
to the NVECTOR_RAJA implementation in the future.  Additionally, the
operations ``N_VGetArrayPointer`` and ``N_VSetArrayPointer`` are not
implemented by the RAJA vector.  As such, this
vector cannot be used with SUNDIALS Fortran interfaces, nor with
SUNDIALS direct solvers and preconditioners. The NVECTOR_RAJA module
provides separate functions to access data on the host and on the
device. It also provides methods for copying from the host to the
device and vice versa. Usage examples of NVECTOR_RAJA are provided in  
some example programs for CVODE [HSR2017]_.

The names of vector operations are obtained from those in the sections
:ref:`NVectors.Ops`, :ref:`NVectors.FusedOps` and
:ref:`NVectors.ArrayOps` by appending the suffix ``_Raja`` 
(e.g. ``N_VDestroy_Raja``).  The module NVECTOR_RAJA 
provides the following additional user-callable routines:


.. c:function:: N_Vector N_VNew_Raja(sunindextype vec_length)

   This function creates and allocates memory for a RAJA
   ``N_Vector``. The memory is allocated on both the host and the
   device. Its only argument is the vector length. 


.. c:function:: N_Vector N_VNewEmpty_Raja(sunindextype vec_length)

   This function creates a new ``N_Vector`` wrapper with the pointer
   to the wrapped RAJA vector set to ``NULL``.  It is used by
   :c:func:`N_VNew_Raja()`, :c:func:`N_VMake_Raja()`, and
   :c:func:`N_VClone_Raja()` implementations. 

      
.. c:function:: N_Vector N_VMake_Raja(N_VectorContent_Raja c)

   This function creates and allocates memory for an NVECTOR_RAJA
   wrapper around a user-provided ``sunrajavec::Vector`` class.  
   Its only argument is of type ``N_VectorContent_Raja``, which
   is the pointer to the class.

 
.. c:function:: N_Vector* N_VCloneVectorArray_Raja(int count, N_Vector w)

   This function creates (by cloning) an array of *count* NVECTOR_RAJA
   vectors. 


.. c:function:: N_Vector* N_VCloneVectorArrayEmpty_Raja(int count, N_Vector w)

   This function creates (by cloning) an array of *count* NVECTOR_RAJA
   vectors, each with pointers to RAJA vectors set to ``NULL``. 


.. c:function:: void N_VDestroyVectorArray_Raja(N_Vector* vs, int count)
  
   This function frees memory allocated for the array of *count*
   variables of type ``N_Vector`` created with
   :c:func:`N_VCloneVectorArray_Raja()` or with
   :c:func:`N_VCloneVectorArrayEmpty_Raja()`. 


.. c:function:: sunindextype N_VGetLength_Raja(N_Vector v)

   This function returns the length of the vector.


.. c:function:: realtype* N_VGetHostArrayPointer_Raja(N_Vector v)

   This function returns a pointer to the vector data on the host.


.. c:function:: realtype* N_VGetDeviceArrayPointer_Raja(N_Vector v)

   This function returns a pointer to the vector data on the device.


.. c:function:: realtype* N_VCopyToDevice_Raja(N_Vector v)

   This function copies host vector data to the device.


.. c:function:: realtype* N_VCopyFromDevice_Raja(N_Vector v)

   This function copies vector data from the device to the host.


.. c:function:: void N_VPrint_Raja(N_Vector v)

   This function prints the content of a RAJA vector to ``stdout``.


.. c:function:: void N_VPrintFile_Raja(N_Vector v, FILE *outfile)

   This function prints the content of a RAJA vector to ``outfile``.

    
By default all fused and vector array operations are disabled in the NVECTOR_RAJA
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with :c:func:`N_VNew_Raja`, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using :c:func:`N_VClone`. This guarantees the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned from while vectors created with
:c:func:`N_VNew_Raja` will have the default settings for the NVECTOR_RAJA module.

.. c:function:: void N_VEnableFusedOps_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.
   
.. c:function:: void N_VEnableLinearCombination_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableScaleAddMulti_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the RAJA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

..
   .. c:function:: void N_VEnableDotProdMulti_Raja(N_Vector v, booleantype tf)

      This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
      dot products fused operation in the RAJA vector. The return value is ``0``
      for success and ``-1`` if the input vector or its ``ops`` structure are
      ``NULL``.

.. c:function:: void N_VEnableLinearSumVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableScaleVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableConstVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

..
   .. c:function:: void N_VEnableWrmsNormVectorArray_Raja(N_Vector v, booleantype tf)

      This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
      operation for vector arrays in the RAJA vector. The return value is ``0`` for
      success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

   .. c:function:: void N_VEnableWrmsNormMaskVectorArray_Raja(N_Vector v, booleantype tf)

      This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
      norm operation for vector arrays in the RAJA vector. The return value is
      ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
      ``NULL``.

.. c:function:: void N_VEnableScaleAddMultiVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the RAJA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: void N_VEnableLinearCombinationVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the RAJA vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When there is a need to access components of an ``N_Vector_Raja``, ``v``, 
  it is recommeded to use functions :c:func:`N_VGetDeviceArrayPointer_Raja()` or 
  :c:func:`N_VGetHostArrayPointer_Raja()`.        

* To maximize efficiency, vector operations in the NVECTOR_RAJA implementation
  that have more than one ``N_Vector`` argument do not check for
  consistent internal representations of these vectors. It is the user's 
  responsibility to ensure that such routines are called with ``N_Vector``
  arguments that were all created with the same internal representations.
