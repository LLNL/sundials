..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _NVectors.NVPETSc:

The NVECTOR_PETSC Module
================================

The NVECTOR_PETSC module is an NVECTOR wrapper around the PETSc vector. It
defines the *content* field of a ``N_Vector`` to be a structure
containing the global and local lengths of the vector, a pointer to
the PETSc vector, an MPI communicator, and a boolean flag  *own_data*
indicating ownership of the wrapped PETSc vector.

.. code-block:: c

   struct _N_VectorContent_Petsc { 
      sunindextype local_length; 
      sunindextype global_length; 
      booleantype own_data;
      Vec *pvec;
      MPI_Comm comm; 
   };

The header file to be included when using this module is
``nvector_petsc.h``.  Unlike native SUNDIALS vector types,
NVECTOR_PETSC does not provide macros to access its member variables.
Note that NVECTOR_PETSC requires SUNDIALS to be built with MPI support.



The NVECTOR_PETSC module defines implementations of all vector
operations listed in the section :ref:`NVectors.Ops`, except for
``N_VGetArrayPointer`` and ``N_VSetArrayPointer``.  As such, this
vector cannot be used with SUNDIALS Fortran interfaces.  When access
to raw vector data is needed, it is recommended to extract the PETSc
vector first, and then use PETSc methods to access the data.  Usage
examples of NVECTOR_PETSC is provided in example programs for IDA.

The names of vector operations are obtained from those in the section
:ref:`NVectors.Ops` by appending the suffice ``_Petsc``
(e.g. ``N_VDestroy_Petsc``).  The module NVECTOR_PETSC provides the
following additional user-callable routines:


.. c:function:: N_Vector N_VNewEmpty_Petsc(MPI_Comm comm, sunindextype local_length, sunindextype global_length)

   This function creates a new PETSC ``N_Vector`` with the pointer to
   the wrapped PETSc vector set to ``NULL``. It is used by the
   ``N_VMake_Petsc`` and ``N_VClone_Petsc`` implementations.  It
   should be used only with great caution.
 

.. c:function:: N_Vector N_VMake_Petsc(Vec* pvec)

   This function creates and allocates memory for an NVECTOR_PETSC
   wrapper with a user-provided PETSc vector.  It does *not* allocate
   memory for the vector ``pvec`` itself.


.. c:function:: Vec *N_VGetVector_Petsc(N_Vector v)

   This function returns a pointer to the underlying PETSc vector.


.. c:function:: N_Vector* N_VCloneVectorArray_Petsc(int count, N_Vector w)

   This function creates (by cloning) an array of *count*
   NVECTOR_PETSC vectors.


.. c:function:: N_Vector* N_VCloneVectorArrayEmpty_Petsc(int count, N_Vector w)

   This function creates (by cloning) an array of *count*
   NVECTOR_PETSC vectors, each with pointers to PETSc vectors set to ``NULL``. 


.. c:function:: void N_VDestroyVectorArray_Petsc(N_Vector* vs, int count)

   This function frees memory allocated for the array of *count*
   variables of type ``N_Vector`` created with
   :c:func:`N_VCloneVectorArray_Petsc()` or with
   :c:func:`N_VCloneVectorArrayEmpty_Petsc()`. 


.. c:function:: void N_VPrint_Petsc(N_Vector v)

   This function prints the global content of a wrapped PETSc vector to ``stdout``. 


.. c:function:: void N_VPrintFile_Petsc(N_Vector v, const char fname[])

   This function prints the global content of a wrapped PETSc vector to ``fname``. 




**Notes**

* When there is a need to access components of an ``N_Vector_Petsc v``, it
  is recommeded to extract the PETSc vector via 

  ``x_vec = N_VGetVector_Petsc(v);`` 

  and then access components using appropriate PETSc functions.

* The functions :c:func:`N_VNewEmpty_Petsc()`, :c:func:`N_VMake_Petsc()`, 
  and :c:func:`N_VCloneVectorArrayEmpty_Petsc()` set the field
  *own_data* to ``SUNFALSE``. The routines :c:func:`N_VDestroy_Petsc()` and
  :c:func:`N_VDestroyVectorArray_Petsc()` will not attempt to free the
  pointer ``pvec`` for any ``N_Vector`` with *own_data* set to
  ``SUNFALSE``. In such a case, it is the user's responsibility to
  deallocate the ``pvec`` pointer. 

* To maximize efficiency, vector operations in the NVECTOR_PETSC
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representations of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  internal representations.
