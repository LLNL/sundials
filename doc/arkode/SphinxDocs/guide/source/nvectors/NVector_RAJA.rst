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


.. _NVectors.RAJA:

The NVECTOR_RAJA Module
======================================

The NVECTOR_RAJA module is an experimental {\nvector} implementation using the
`RAJA <https://software.llnl.gov/RAJA/>`_ hardware abstraction
layer. In this implementation, RAJA
allows for SUNDIALS vector kernels to run on GPU devices. The module is intended for users
who are already familiar with RAJA and GPU programming. Building this vector
module requires a C++11 compliant compiler and a CUDA software development toolkit.
Besides the CUDA backend, RAJA has other backends such as serial, OpenMP,
and OpenACC. These backends are not used in this SUNDIALS release.
The vector content layout is as follows:

.. code-block:: c++

   struct _N_VectorContent_Raja
   {
      sunindextype length;
      booleantype  own_data;
      realtype*    host_data;
      realtype*    device_data;
      void*        priv; /* 'private' data */
   };


The content members are the vector length (size), a boolean flag that signals if
the vector owns the data (i.e., it is in charge of freeing the data), pointers to
vector data on the host and the device, and a private data structure which holds
the memory management type, which should not be accessed directly.

When instantiated with ``N_VNew_Raja``, the underlying data will be allocated
on both the host and the device. Alternatively, a user can provide host
and device data arrays by using the ``N_VMake_Raja`` constructor. To use CUDA
managed memory, the constructors ``N_VNewManaged_Raja`` and
``N_VMakeManaged_Raja`` are provided. Details on each of these constructors
are provided below.

The header file to include when using this is ``nvector_raja.h``.
The installed module library to link to is ``libsundials_nveccudaraja.lib``.
The extension ``.lib`` is typically ``.so`` for shared libraries ``.a`` for
static libraries.


NVECTOR_RAJA functions
-----------------------------------

Unlike other native SUNDIALS vector types, the NVECTOR_RAJA module does not
provide macros to access its member variables. Instead, user should use the
accessor functions:



.. c:function:: realtype* N_VGetHostArrayPointer_Raja(N_Vector v)

   This function returns pointer to the vector data on the host.


.. c:function:: realtype* N_VGetDeviceArrayPointer_Raja(N_Vector v)

   This function returns pointer to the vector data on the device.

.. c:function:: booleantype N_VIsManagedMemory_Raja(N_Vector v)

   This function returns a boolean flag indicating if the vector
   data is allocated in managed memory or not.


The NVECTOR_RAJA module defines the implementations of all vector
operations listed in the sections :ref:`NVectors.Ops`,
:ref:`NVectors.FusedOps`, :ref:`NVectors.ArrayOps`, and
:ref:`NVectors.LocalOps`, except for
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
:ref:`NVectors.Ops`, :ref:`NVectors.FusedOps`, :ref:`NVectors.ArrayOps`, and
:ref:`NVectors.LocalOps` by appending the suffix ``_Raja``
(e.g. ``N_VDestroy_Raja``).  The module NVECTOR_RAJA
provides the following additional user-callable routines:


.. c:function:: N_Vector N_VNew_Raja(sunindextype vec_length)

   This function creates and allocates memory for a RAJA
   ``N_Vector``. The memory is allocated on both the host and the
   device. Its only argument is the vector length.


.. c:function:: N_Vector N_VNewManaged_Raja(sunindextype vec_length)

   This function creates and allocates memory for a RAJA ``N_Vector``.
   The vector data array is allocated in managed memory.


.. c:function:: N_Vector N_VNewEmpty_Raja(sunindextype vec_length)

   This function creates a new ``N_Vector`` wrapper with the pointer
   to the wrapped RAJA vector set to ``NULL``.  It is used by
   :c:func:`N_VNew_Raja()`, :c:func:`N_VMake_Raja()`, and
   :c:func:`N_VClone_Raja()` implementations.


.. c:function:: N_Vector N_VMake_Raja(sunindextype length, realtype *vdata)

   This function creates an NVECTOR_RAJA with a user-supplied managed
   memory data array. This function does not allocate memory for data itself.


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

.. c:function:: int N_VEnableFusedOps_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the RAJA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

..
   .. c:function:: int N_VEnableDotProdMulti_Raja(N_Vector v, booleantype tf)

      This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
      dot products fused operation in the RAJA vector. The return value is ``0``
      for success and ``-1`` if the input vector or its ``ops`` structure are
      ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the RAJA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

..
   .. c:function:: int N_VEnableWrmsNormVectorArray_Raja(N_Vector v, booleantype tf)

      This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
      operation for vector arrays in the RAJA vector. The return value is ``0`` for
      success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

   .. c:function:: int N_VEnableWrmsNormMaskVectorArray_Raja(N_Vector v, booleantype tf)

      This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
      norm operation for vector arrays in the RAJA vector. The return value is
      ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
      ``NULL``.

.. c:function:: int N_VEnableScaleAddMultiVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the RAJA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombinationVectorArray_Raja(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the RAJA vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When there is a need to access components of an ``N_Vector_Raja``, ``v``,
  it is recommended to use functions :c:func:`N_VGetDeviceArrayPointer_Raja()` or
  :c:func:`N_VGetHostArrayPointer_Raja()`. However, when using managed memory,
  the function ``N_VGetArrayPointer`` may also be used.

* To maximize efficiency, vector operations in the NVECTOR_RAJA implementation
  that have more than one ``N_Vector`` argument do not check for
  consistent internal representations of these vectors. It is the user's
  responsibility to ensure that such routines are called with ``N_Vector``
  arguments that were all created with the same internal representations.
