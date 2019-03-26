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

.. _NVectors.ManyVector:

The NVECTOR_MANYVECTOR Module
================================

The NVECTOR_MANYVECTOR implementation of the NVECTOR module provided with
SUNDIALS is designed to facilitate problems with an inherent
data partitioning for the solution vector.  These data partitions are
entirely user-defined, through construction of distinct NVECTOR
modules for each component, that are then combined together to form
the NVECTOR_MANYVECTOR.  We envision three generic use cases for this
implementation:

A. *Heterogenous computational architectures (serial or parallel)*:
   for users who wish to partition data on a node between different
   computing resources, they may create architecture-specific
   subvectors for each partition.  For example, a user could create
   one MPI-parallel component based on NVECTOR_PARALLEL, another
   single-node component for GPU accelerators based on NVECTOR_CUDA,
   and another threaded single-node component based on
   NVECTOR_OPENMP.

B. *Process-based multiphysics decompositions (parallel)*: for users
   who wish to combine separate simulations together, e.g., where
   one subvector resides on one subset of MPI processes, while
   another subvector resides on a different subset of MPI processes,
   and where the user has created a MPI *intercommunicator* to
   connect these distinct process sets together.

C. *Structure of arrays (SOA) data layouts (serial or parallel)*:
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
of operations, however significant performance benefits may be
obtained when subvectors additionally implement the optional local
reduction operations listed in the section :ref:`NVectors.LocalOps`.

Additionally, NVECTOR_MANYVECTOR sets no limit on the number of
subvectors that may be attached (aside from the limitations of using
``sunindextype`` for indexing, and standard per-node memory
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
of ``N_Vector`` to be a structure containing the MPI communicator
(or ``SUNMPI_COMM_NULL`` if running in serial), the number of
subvectors comprising the ManyVector, the global length of the
ManyVector (including all subvectors on all MPI tasks), a pointer to
the beginning of the array of subvectors, and a boolean flag
``own_data`` indicating ownership of the subvectors that populate
``subvec_array``.

.. code-block:: c

   struct _N_VectorContent_ManyVector {
     SUNMPI_Comm   comm;            /* overall MPI communicator        */
     sunindextype  num_subvectors;  /* number of vectors attached      */
     sunindextype  global_length;   /* overall manyvector length       */
     N_Vector*     subvec_array;    /* pointer to N_Vector array       */
     booleantype   own_data;        /* flag indicating data ownership  */
   };

The header file to include when using this module is
``nvector_manyvector.h``. The installed module library to link against is 
``libsundials_nvecmanyvector.lib`` where ``.lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.

*Note:* If SUNDIALS is configured with MPI enabled, then the
ManyVector library will be built for the parallel use case and ``SUNMPI_Comm``
is set to ``MPI_Comm``. As such an MPI-aware compiler will become necessary
even in single node uses of the ManyVector library. If SUNDIALS is configured
with MPI disabled, then the ManyVector library is built for the single-node
(serial) use case and ``SUNMPI_Comm`` is set to ``int``. As such, users need
not include ``mpi.h``, nor must executables be built with an MPI-aware
compiler. See the :ref:`_Installation` section for details.


NVECTOR_MANYVECTOR functions
-------------------------------

The NVECTOR_MANYVECTOR module implements all vector operations listed
in the sections :ref:`NVectors.Ops`, :ref:`NVectors.FusedOps`,
:ref:`NVectors.ArrayOps`, and :ref:`NVectors.LocalOps`, except for
:c:func:`N_VGetArrayPointer()`, :c:func:`N_VSetArrayPointer()`,
:c:func:`N_VScaleAddMultiVectorArray()`, and
:c:func:`N_VLinearCombinationVectorArray()`.  As such, this vector
cannot be used with the SUNDIALS Fortran interfaces, nor with the
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
   NVECTOR objects, under the requirement that all MPI-aware
   subvectors use the same MPI communicator (this is checked
   internally).  If none of the subvectors are MPI-aware, then this
   may equivalently be used to describe data partitioning within a
   single node.  We note that this routine is designed to support use
   cases A and C above.

   This routine will copy all ``N_Vector`` pointers from the input
   ``vec_array``, so the user may modify/free that pointer array
   after calling this function.  However, this routine does *not*
   allocate any new subvectors, so the underlying NVECTOR objects
   themselves should not be destroyed before the ManyVector that
   contains them.

   Upon successful completion, the new ManyVector is returned;
   otherwise this routine returns ``NULL`` (e.g., if two MPI-aware
   subvectors use different MPI communicators).


.. c:function:: N_Vector N_VMake_ManyVector(void *comm, sunindextype num_subvectors, N_Vector *vec_array)

   This function creates a ManyVector from a set of existing NVECTOR
   objects, and a user-created MPI communicator that "connects" these
   subvectors.  Any MPI-aware subvectors may use different MPI
   communicators than the input *comm*.  We note that this routine
   is designed to support any combination of the use cases above.

   The input *comm* should be the memory reference to this
   user-created MPI communicator.  We note that since many MPI
   implementations ``#define``  ``MPI_COMM_WORLD`` to be a specific
   integer *value* (that has no memory reference), users who wish
   to supply ``MPI_COMM_WORLD`` to this routine should first
   duplicate this to a specific ``MPI_Comm`` variable before passing
   in the reference, e.g.

   .. code-block:: c

      MPI_Comm comm;
      N_Vector x;
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);
      x = N_VMake_ManyVector(&comm, ...);

   This routine will internally call ``MPI_Comm_dup`` to create a
   copy of the input ``comm``, so the user-supplied ``comm`` argument
   need not be retained after the call to
   :c:func:`N_VMake_ManyVector()`.

   If all subvectors are MPI-unaware, then the input *comm* argument
   should be ``NULL``, although in this case, it would be simpler to
   call :c:func:`N_VNew_ManyVector()` instead.

   This routine will copy all ``N_Vector`` pointers from the input
   *vec_array*, so the user may modify/free that pointer array
   after calling this function.  However, this routine does *not*
   allocate any new subvectors, so the underlying NVECTOR objects
   themselves should not be destroyed before the ManyVector that
   contains them.

   Upon successful completion, the new ManyVector is returned;
   otherwise this routine returns ``NULL`` (e.g., if the input
   *vec_array* is ``NULL``).


.. c:function:: N_Vector N_VGetSubvector_ManyVector(N_Vector v, sunindextype vec_num)

   This function returns the *vec_num* subvector from the NVECTOR array.


.. c:function:: sunindextype N_VGetNumSubvectors_ManyVector(N_Vector v)

   This function returns the overall number of subvectors in the ManyVector object.


By default all fused and vector array operations are disabled in the NVECTOR_MANYVECTOR
module, except for :c:func:`N_VWrmsNormVectorArray()` and
:c:func:`N_VWrmsNormMaskVectorArray()`, that are enabled by default. The
following additional user-callable routines are provided to enable or
disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a
vector with :c:func:`N_VNew_ManyVector()` or
:c:func:`N_VMake_ManyVector()`, enable/disable the desired operations
for that vector with the functions below, and create any additional
vectors from that vector using :c:func:`N_VClone()`. This guarantees
that the new vectors will have the same operations enabled/disabled,
since cloned vectors inherit those configuration options from the
vector they are cloned from, while vectors created with
:c:func:`N_VNew_ManyVector()` and :c:func:`N_VMake_ManyVector()` will
have the default settings for the NVECTOR_MANYVECTOR module.  We note
that these routines *do not* call the corresponding routines on
subvectors, so those should be set up as desired *before* attaching
them to the ManyVector in :c:func:`N_VNew_ManyVector()` or
:c:func:`N_VMake_ManyVector()`.

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

* :c:func:`N_VNew_ManyVector()` and :c:func:`N_VMake_ManyVector()` set
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
