..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2015, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _NVectors.ParHyp:

The NVECTOR_PARHYP Module
======================================

The NVECTOR_PARHYP implementation of the NVECTOR  module provided with
SUNDIALS is a wrapper around HYPRE's ParVector class. 
Most of the vector kernels simply call HYPRE vector operations. 
The implementation defines the *content* field of ``N_Vector`` to 
be a structure containing the global and local lengths of the vector, a 
pointer to an object of type ``hypre_ParVector``, an MPI communicator, 
and a boolean flag *own_parvector* indicating ownership of the
HYPRE parallel vector object *x*.


.. code-block:: c

   struct _N_VectorContent_ParHyp {
     sunindextype local_length;
     sunindextype global_length;
     booleantype own_data;
     booleantype own_parvector;
     realtype *data;
     MPI_Comm comm;
     hypre_ParVector *x;
   };

The header file to be included when using this module is ``nvector_parhyp.h``.
Unlike native SUNDIALS vector types, NVECTOR_PARHYP does not provide macros 
to access its member variables.


The NVECTOR_PARHYP module defines implementations of all vector
operations listed in the section :ref:`NVectors.Ops`, except for
``N_VSetArrayPointer`` and ``N_VGetArrayPointer``, because accessing
raw vector data is handled by low-level HYPRE functions.  As such,
this vector is not available for use with SUNDIALS Fortran
interfaces.  When access to raw vector data is needed, one should
extract the HYPRE HYPRE vector first, and then use HYPRE methods to
access the data.  Usage examples of NVECTOR_PARHYP are provided in
the ``cvAdvDiff_non_ph.c`` example programs for CVODE and the
``ark_diurnal_kry_ph.c`` example program for ARKode.

The names of parhyp methods are obtained from those in the section
:ref:`NVectors.Ops` by appending the suffix ``_ParHyp``
(e.g. ``N_VDestroy_ParHyp``).  The module {\nvecph} provides the
following additional user-callable routines:


.. c:function:: N_Vector N_VNewEmpty_ParHyp(MPI_Comm comm, sunindextype local_length, sunindextype global_length)

   This function creates a new parhyp ``N_Vector`` with the pointer to the
   HYPRE vector set to ``NULL``. 


.. c:function:: N_Vector N_VMake_ParHyp(hypre_ParVector *x)

   This function creates an ``N_Vector`` wrapper around an existing
   HYPRE parallel vector.  It does *not* allocate memory for ``x`` itself.


.. c:function:: hypre_ParVector *N_VGetVector_ParHyp(N_Vector v)
  
   This function returns a pointer to the underlying HYPRE vector.


.. c:function:: N_Vector* N_VCloneVectorArray_ParHyp(int count, N_Vector w)

   This function creates (by cloning) an array of *count* parhyp
   vectors. 


.. c:function:: N_Vector* N_VCloneVectorArrayEmpty_ParHyp(int count, N_Vector w)

   This function creates (by cloning) an array of *count* parhyp
   vectors, each with an empty (```NULL``) data array.


.. c:function:: void N_VDestroyVectorArray_ParHyp(N_Vector* vs, int count)
  
   This function frees memory allocated for the array of *count*
   variables of type ``N_Vector`` created with
   :c:func:`N_VCloneVectorArray_ParHyp()` or with
   :c:func:`N_VCloneVectorArrayEmpty_ParHyp()`. 


.. c:function:: void N_VPrint_ParHyp(N_Vector v)

   This function prints the local content of a parhyp vector to ``stdout``.


.. c:function:: void N_VPrintFile_ParHyp(N_Vector v, FILE *outfile)

   This function prints the local content of a parhyp vector to ``outfile``.

    

**Notes**

* When there is a need to access components of an ``N_Vector_ParHyp v``, 
  it is recommended to extract the HYPRE vector via 
  ``x_vec = N_VGetVector_ParHyp(v)`` and then access components using 
  appropriate HYPRE functions.

* :c:func:`N_VNewEmpty_ParHyp()`, :c:func:`N_VMake_ParHyp()`, and
  :c:func:`N_VCloneVectorArrayEmpty_ParHyp()` set the field *own_parvector*
  to ``SUNFALSE``.  The functions :c:func:`N_VDestroy_ParHyp()` and
  :c:func:`N_VDestroyVectorArray_ParHyp()` will not attempt to delete an
  underlying HYPRE vector for any ``N_Vector`` with *own_parvector*
  set to ``SUNFALSE``.  In such a case, it is the user's responsibility
  to delete the underlying vector.

* To maximize efficiency, vector operations in the NVECTOR_PARHYP
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representations of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  internal representations.


