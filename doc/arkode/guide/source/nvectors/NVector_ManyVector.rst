..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3

.. _NVectors.ManyVector:

The NVECTOR_MANYVECTOR Module
================================

The NVECTOR_MANYVECTOR implementation of the NVECTOR module provided with
SUNDIALS is designed to facilitate problems with an inherent
data partitioning for the solution vector within a computational node.
These data partitions are entirely user-defined, through construction
of distinct NVECTOR modules for each component, that are then combined
together to form the NVECTOR_MANYVECTOR.  We envision two generic use
cases for this implementation:

A. *Heterogenous computational architectures*:
   for users who wish to partition data on a node between different
   computing resources, they may create architecture-specific
   subvectors for each partition.  For example, a user could create
   one serial component based on NVECTOR_SERIAL, another component for
   GPU accelerators based on NVECTOR_CUDA, and another threaded
   component based on NVECTOR_OPENMP.

B. *Structure of arrays (SOA) data layouts*:
   for users who wish to create separate subvectors for each
   solution component, e.g., in a Navier-Stokes simulation they
   could have separate subvectors for density, velocities and
   pressure, which are combined together into a single
   NVECTOR_MANYVECTOR for the overall "solution".

We note that the above use cases are not mutually exclusive, and the
NVECTOR_MANYVECTOR implementation should support arbitrary combinations
of these cases.

The NVECTOR_MANYVECTOR implementation is designed to work with any
NVECTOR subvectors that implement the minimum *required* set
of operations.  Additionally, NVECTOR_MANYVECTOR sets no limit on the
number of subvectors that may be attached (aside from the limitations
of using ``sunindextype`` for indexing, and standard per-node memory
limitations).  However, while this ostensibly supports subvectors
with one entry each (i.e., one subvector for each solution entry), we
anticipate that this extreme situation will hinder performance due to
non-stride-one memory accesses and increased function call overhead.
We therefore recommend a relatively coarse partitioning of the
problem, although actual performance will likely be
problem-dependent.

As a final note, in the coming years we plan to introduce additional
algebraic solvers and time integration modules that will leverage the
problem partitioning enabled by NVECTOR_MANYVECTOR.  However, even at
present we anticipate that users will be able to leverage such data
partitioning in their problem-defining ODE right-hand side, DAE
residual, or nonlinear solver residual functions.



NVECTOR_MANYVECTOR structure
-------------------------------

The NVECTOR_MANYVECTOR implementation defines the *content* field
of ``N_Vector`` to be a structure containing the number of
subvectors comprising the ManyVector, the global length of the
ManyVector (including all subvectors), a pointer to
the beginning of the array of subvectors, and a boolean flag
``own_data`` indicating ownership of the subvectors that populate
``subvec_array``.

.. code-block:: c

   struct _N_VectorContent_ManyVector {
     sunindextype  num_subvectors;  /* number of vectors attached      */
     sunindextype  global_length;   /* overall manyvector length       */
     N_Vector*     subvec_array;    /* pointer to N_Vector array       */
     booleantype   own_data;        /* flag indicating data ownership  */
   };

The header file to include when using this module is
``nvector_manyvector.h``. The installed module library to link against is
``libsundials_nvecmanyvector.lib`` where ``.lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.



NVECTOR_MANYVECTOR functions
-------------------------------

The NVECTOR_MANYVECTOR module implements all vector operations listed
in the sections :ref:`NVectors.Ops`, :ref:`NVectors.FusedOps`,
:ref:`NVectors.ArrayOps`, and :ref:`NVectors.LocalOps`, except for
:c:func:`N_VGetArrayPointer()`, :c:func:`N_VSetArrayPointer()`,
:c:func:`N_VScaleAddMultiVectorArray()`, and
:c:func:`N_VLinearCombinationVectorArray()`.  As such, this vector
cannot be used with the SUNDIALS Fortran-77 interfaces, nor with the
SUNDIALS direct solvers and preconditioners. Instead, the
NVECTOR_MANYVECTOR module provides functions to access subvectors,
whose data may in turn be accessed according to their NVECTOR
implementations.

The names of vector operations are obtained from those in the sections
:ref:`NVectors.Ops`, :ref:`NVectors.FusedOps`,
:ref:`NVectors.ArrayOps`, and :ref:`NVectors.LocalOps` by
appending the suffix ``_ManyVector`` (e.g. ``N_VDestroy_ManyVector``).
The module NVECTOR_MANYVECTOR provides the following additional
user-callable routines:

.. c:function:: N_Vector N_VNew_ManyVector(sunindextype num_subvectors, N_Vector *vec_array)

   This function creates a ManyVector from a set of existing
   NVECTOR objects.

   This routine will copy all ``N_Vector`` pointers from the input
   ``vec_array``, so the user may modify/free that pointer array
   after calling this function.  However, this routine does *not*
   allocate any new subvectors, so the underlying NVECTOR objects
   themselves should not be destroyed before the ManyVector that
   contains them.

   Upon successful completion, the new ManyVector is returned;
   otherwise this routine returns ``NULL`` (e.g., a memory allocation
   failure occurred).

   Users of the Fortran 2003 interface to this function will first need to use
   the generic ``N_Vector`` utility functions ``N_VNewVectorArray``, and
   ``N_VSetVecAtIndexVectorArray`` to create the ``N_Vector*`` argument.  This is
   further explained in Chapter :ref:`Fortran2003.Differences.NVectorArrays`,
   and the functions are documented in Chapter :ref:`NVectors.utilities`.


.. c:function:: N_Vector N_VGetSubvector_ManyVector(N_Vector v, sunindextype vec_num)

   This function returns the *vec_num* subvector from the NVECTOR array.


.. c:function:: realtype *N_VGetSubvectorArrayPointer_ManyVector(N_Vector v, sunindextype vec_num)

   This function returns the data array pointer for the *vec_num*
   subvector from the NVECTOR array.

   If the input *vec_num* is invalid, or if the subvector does not
   support the ``N_VGetArrayPointer`` operation, then ``NULL`` is
   returned.


.. c:function:: int N_VSetSubvectorArrayPointer_ManyVector(realtype *v_data, N_Vector v, sunindextype vec_num)

   This function sets the data array pointer for the *vec_num*
   subvector from the NVECTOR array.

   If the input *vec_num* is invalid, or if the subvector does not
   support the ``N_VSetArrayPointer`` operation, then ``-1`` is
   returned; otherwise it returns ``0``.


.. c:function:: sunindextype N_VGetNumSubvectors_ManyVector(N_Vector v)

   This function returns the overall number of subvectors in the ManyVector object.


By default all fused and vector array operations are disabled in the
NVECTOR_MANYVECTOR module, except for :c:func:`N_VWrmsNormVectorArray()`
and :c:func:`N_VWrmsNormMaskVectorArray()`, that are enabled by
default. The following additional user-callable routines are provided
to enable or disable fused and vector array operations for a specific
vector. To ensure consistency across vectors it is recommended to
first create a vector with :c:func:`N_VNew_ManyVector()`,
enable/disable the desired operations
for that vector with the functions below, and create any additional
vectors from that vector using :c:func:`N_VClone()`. This guarantees
that the new vectors will have the same operations enabled/disabled,
since cloned vectors inherit those configuration options from the
vector they are cloned from, while vectors created with
:c:func:`N_VNew_ManyVector()` will
have the default settings for the NVECTOR_MANYVECTOR module.  We note
that these routines *do not* call the corresponding routines on
subvectors, so those should be set up as desired *before* attaching
them to the ManyVector in :c:func:`N_VNew_ManyVector()`.

.. c:function:: int N_VEnableFusedOps_ManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the manyvector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_ManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the manyvector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_ManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the manyvector vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_ManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the manyvector vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_ManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the manyvector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_ManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the manyvector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_ManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the manyvector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_ManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the manyvector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_ManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the manyvector vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.


**Notes**

* :c:func:`N_VNew_ManyVector()` sets
  the field ``own_data = SUNFALSE``.
  :c:func:`N_VDestroy_ManyVector()` will not attempt to call
  :c:func:`N_VDestroy()` on any subvectors contained in the
  subvector array for any ``N_Vector`` with ``own_data`` set to
  ``SUNFALSE``. In such a case, it is the user's responsibility to
  deallocate the subvectors.

* To maximize efficiency, arithmetic vector operations in the
  NVECTOR_MANYVECTOR implementation that have more than one
  ``N_Vector`` argument do not check for consistent internal
  representation of these vectors. It is the user's responsibility to
  ensure that such routines are called with ``N_Vector`` arguments
  that were all created with the same subvector representations.
