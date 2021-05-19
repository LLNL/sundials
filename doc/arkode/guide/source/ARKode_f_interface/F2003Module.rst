..
   Programmer(s): Cody J. Balos @ LLNL
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

.. _Fortran2003:

=====================================
ARKode Fortran 2003 Interface Modules
=====================================

The ARKode Fortran 2003 modules define interfaces to most of the ARKode C
API using the intrinsic ``iso_c_binding`` module which provides a standardized
mechanism for interoperating with C. AKRode provides four Fortran 2003 modules:

* ``farkode_arkstep_mod``, ``farkode_erkstep_mod``, ``farkode_mristep_mod`` provide
  interfaces to the ARKStep, ERKStep, and MRIStep time-stepping modules respectively
* ``farkode_mod`` which interfaces to the components of ARKode which are shared by the
  time-stepping modules

All interfaced functions are named after the corresponding C function, but with a
leading 'F'. For example. the ARKStep function ``ARKStepCreate`` is interfaced as
``FARKStepCreate``. Thus, the steps to use an ARKode time-stepping module from Fortran
are identical (ignoring language differences) to using it from C/C++.

The Fortran 2003 ARKode interface modules can be accessed by the ``use`` statement,
i.e. ``use farkode_mod``, and linking to the library ``libsundials_farkode_mod.lib``
in addition to ``libsundials_farkode.lib``. Further information on the location of
installed modules is provided in the Chapter :ref:`Installation`.

The Fortran 2003 interface modules were generated with SWIG Fortran, a
fork of SWIG [JPE2019]_. Users who are interested in the SWIG code used
in the generation process should contact the SUNDIALS development team.

.. _Fortran2003.SundialsModules:

SUNDIALS Fortran 2003 Interface Modules
----------------------------------------

All of the generic SUNDIALS modules provide Fortran 2003 interface modules.
Many of the generic module implementations provide Fortran 2003 interfaces
(a complete list of modules with Fortran 2003 interfaces is given in
:ref:`Fortran2003.InterfacesTable`. A module can be accessed with the ``use``
statement, e.g. ``use fnvector_openmp_mod``, and linking to the Fortran
2003 library in addition to the C library, e.g.
``libsundials_fnvecpenmp_mod.lib`` and ``libsundials_nvecopenmp.lib``.

The Fortran 2003 interfaces leverage the ``iso_c_binding`` module and the
``bind(C)`` attribute to closely follow the SUNDIALS C API (ignoring
language differences). The generic SUNDIALS structures, e.g. ``N_Vector``,
are interfaced as Fortran derived types, and function signatures are matched
but with an ``F`` prepending the name, e.g. ``FN_VConst`` instead of
``N_VConst``. Constants are named exactly as they are in the C API.
Accordingly, using SUNDIALS via the Fortran 2003 interfaces looks just like
using it in C. Some caveats stemming from the language differences are
discussed in the section :ref:`Fortran2003.Differences`. A discussion on the
topic of equivalent data types in C and Fortran 2003 is presented in
section :ref:`Fortran2003.DataTypes`.

Further information on the Fortran 2003 interfaces specific to modules is given
in the NVECTOR, SUNMatrix, SUNLinearSolver, and SUNNonlinearSolver sections
alongside the C documentation (chapters :ref:`NVectors`, :ref:`SUNMatrix`,
:ref:`SUNLinSol`, and :ref:`SUNNonlinSol` respectively). For details on where
the Fortran 2003 module (``.mod``) files and libraries are installed see Appendix
:ref:`Installation`.

.. _Fortran2003.InterfacesTable:

Table: SUNDIALS Fortran 2003 Interface Modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=======================  ====================================
   **Module**               **Fortran 2003 Module Name**
=======================  ====================================
NVECTOR                  ``fsundials_nvector_mod``
NVECTOR_SERIAL           ``fnvector_serial_mod``
NVECTOR_OPENMP           ``fnvector_openmp_mod``
NVECTOR_PTHREADS         ``fnvector_pthreads_mod``
NVECTOR_PARALLEL         ``fnvector_parallel_mod``
NVECTOR_PARHYP           Not interfaced
NVECTOR_PETSC            Not interfaced
NVECTOR_CUDA             Not interfaced
NVECTOR_RAJA             Not interfaced
NVECTOR_MANVECTOR        ``fnvector_manyvector_mod``
NVECTOR_MPIMANVECTOR     ``fnvector_mpimanyvector_mod``
NVECTOR_MPIPLUSX         ``fnvector_mpiplusx_mod``
SUNMATRIX                ``fsundials_matrix_mod``
SUNMATRIX_BAND           ``fsunmatrix_band_mod``
SUNMATRIX_DENSE          ``fsunmatrix_dense_mod``
SUNMATRIX_SPARSE         ``fsunmatrix_sparse_mod``
SUNLINSOL                ``fsundials_linearsolver_mod``
SUNLINSOL_BAND           ``fsunlinsol_band_mod``
SUNLINSOL_DENSE          ``fsunlinsol_dense_mod``
SUNLINSOL_LAPACKBAND     Not interfaced
SUNLINSOL_LAPACKDENSE    Not interfaced
SUNLINSOL_KLU            ``fsunlinsol_klu_mod``
SUNLINSOL_SLUMT          Not interfaced
SUNLINSOL_SLUDIST        Not interfaced
SUNLINSOL_SPGMR          ``fsunlinsol_spgmr_mod``
SUNLINSOL_SPFGMR         ``fsunlinsol_spfgmr_mod``
SUNLINSOL_SPBCGS         ``fsunlinsol_spbcgs_mod``
SUNLINSOL_SPTFQMR        ``fsunlinsol_sptfqmr_mod``
SUNLINSOL_PCG            ``fsunlinsol_pcg_mof``
SUNNONLINSOL             ``fsundials_nonlinearsolver_mod``
SUNNONLINSOL_NEWTON      ``fsunnonlinsol_newton_mod``
SUNNONLINSOL_FIXEDPOINT  ``fsunnonlinsol_fixedpoint_mod``
=======================  ====================================


.. _Fortran2003.DataTypes:

Data Types
----------

Generally, the Fortran 2003 type that is equivalent to the C type is what one
would expect. Primitive types map to the ``iso_c_binding`` type equivalent.
SUNDIALS generic types map to a Fortran derived type. However, the handling
of pointer types is not always clear as they can depend on the parameter direction.
ref:`Fortran2003.DataTypesTable` presents a summary of the type equivalencies
with the parameter direction in mind.

*NOTE*: Currently, the Fortran 2003 interfaces are only compatible with SUNDIALS builds
where the ``realtype`` is double-precision the ``sunindextype`` size is 64-bits.

.. _Fortran2003.DataTypesTable:

Table: C/Fortran-2003 Equivalent Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-------------------------+-------------------------------+-------------------------------------------+
| **C Type**              | **Parameter Direction**       | **Fortran 2003 type**                     |
+=========================+===============================+===========================================+
|``double``               | in, inout, out, return        | ``real(c_double)``                        |
+-------------------------+-------------------------------+-------------------------------------------+
|``int``                  | in, inout, out, return        | ``integer(c_int)``                        |
+-------------------------+-------------------------------+-------------------------------------------+
|``long``                 | in, inout, out, return        | ``integer(c_long)``                       |
+-------------------------+-------------------------------+-------------------------------------------+
|``booleantype``          | in, inout, out, return        | ``integer(c_int)``                        |
+-------------------------+-------------------------------+-------------------------------------------+
|``realtype``             | in, inout, out, return        | ``real(c_double)``                        |
+-------------------------+-------------------------------+-------------------------------------------+
|``sunindextype``         | in, inout, out, return        | ``integer(c_long)``                       |
+-------------------------+-------------------------------+-------------------------------------------+
|``double*``              | in, inout, out                | ``real(c_double), dimension(*)``          |
+-------------------------+-------------------------------+-------------------------------------------+
|``double*``              | return                        | ``real(c_double), pointer, dimension(:)`` |
+-------------------------+-------------------------------+-------------------------------------------+
|``int*``                 | in, inout, out                | ``real(c_int), dimension(*)``             |
+-------------------------+-------------------------------+-------------------------------------------+
|``int*``                 | return                        | ``real(c_int), pointer, dimension(:)``    |
+-------------------------+-------------------------------+-------------------------------------------+
|``long*``                | in, inout, out                | ``real(c_long), dimension(*)``            |
+-------------------------+-------------------------------+-------------------------------------------+
|``long*``                | return                        | ``real(c_long), pointer, dimension(:)``   |
+-------------------------+-------------------------------+-------------------------------------------+
|``realtype*``            | in, inout, out                | ``real(c_double), dimension(*)``          |
+-------------------------+-------------------------------+-------------------------------------------+
|``realtype*``            | return                        | ``real(c_double), pointer, dimension(:)`` |
+-------------------------+-------------------------------+-------------------------------------------+
|``sunindextype*``        | in, inout, out                | ``real(c_long), dimension(*)``            |
+-------------------------+-------------------------------+-------------------------------------------+
|``sunindextype*``        | return                        | ``real(c_long), pointer, dimension(:)``   |
+-------------------------+-------------------------------+-------------------------------------------+
|``realtype[]``           | in, inout, out                | ``real(c_double), dimension(*)``          |
+-------------------------+-------------------------------+-------------------------------------------+
|``sunindextype[]``       | in, inout, out                | ``integer(c_long), dimension(*)``         |
+-------------------------+-------------------------------+-------------------------------------------+
|``N_Vector``             | in, inout, out                | ``type(N_Vector)``                        |
+-------------------------+-------------------------------+-------------------------------------------+
|``N_Vector``             | return                        | ``type(N_Vector), pointer``               |
+-------------------------+-------------------------------+-------------------------------------------+
|``SUNMatrix``            | in, inout, out                | ``type(SUNMatrix)``                       |
+-------------------------+-------------------------------+-------------------------------------------+
|``SUNMatrix``            | return                        | ``type(SUNMatrix), pointer``              |
+-------------------------+-------------------------------+-------------------------------------------+
|``SUNLinearSolver``      | in, inout, out                | ``type(SUNLinearSolver)``                 |
+-------------------------+-------------------------------+-------------------------------------------+
|``SUNLinearSolver``      | return                        | ``type(SUNLinearSolver), pointer``        |
+-------------------------+-------------------------------+-------------------------------------------+
|``SUNNonlinearSolver``   | in, inout, out                | ``type(SUNNonlinearSolver)``              |
+-------------------------+-------------------------------+-------------------------------------------+
|``SUNNonlinearSolver``   | return                        | ``type(SUNNonlinearSolver), pointer``     |
+-------------------------+-------------------------------+-------------------------------------------+
|``FILE*``                | in, inout, out, return        | ``type(c_ptr)``                           |
+-------------------------+-------------------------------+-------------------------------------------+
|``void*``                | in, inout, out, return        | ``type(c_ptr)``                           |
+-------------------------+-------------------------------+-------------------------------------------+
|``T**``                  | in, inout, out, return        | ``type(c_ptr)``                           |
+-------------------------+-------------------------------+-------------------------------------------+
|``T***``                 | in, inout, out, return        | ``type(c_ptr)``                           |
+-------------------------+-------------------------------+-------------------------------------------+
|``T****``                | in, inout, out, return        | ``type(c_ptr)``                           |
+-------------------------+-------------------------------+-------------------------------------------+


.. _Fortran2003.Differences:

Notable Fortran/C usage differences
-----------------------------------

While the Fortran 2003 interface to SUNDIALS closely follows the C API,
some differences are inevitable due to the differences between Fortran and C.
In this section, we note the most critical differences. Additionally, section
:ref:`Fortran2003.DataTypes` discusses equivalencies of data types in the
two languages.


.. _Fortran2003.Differences.CreatingObjects:

Creating generic SUNDIALS objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the C API a generic SUNDIALS object, such as an ``N_Vector``, is actually
a pointer to an underlying C struct. However, in the Fortran 2003 interface,
the derived type is bound to the C struct, not the pointer to the struct. E.g.,
``type(N_Vector)`` is bound to the C struct ``_generic_N_Vector`` not the
``N_Vector`` type. The consequence of this is that creating and declaring SUNDIALS
objects in Fortran is nuanced. This is illustrated in the code snippets below:

C code:

.. sourcecode:: c

   N_Vector x;
   x = N_VNew_Serial(N);

Fortran code:

.. sourcecode:: Fortran

   type(N_Vector), pointer :: x
   x => FN_VNew_Serial(N)

Note that in the Fortran declaration, the vector is a ``type(N_Vector), pointer``,
and that the pointer assignment operator is then used.


.. _Fortran2003.Differences.ArraysAndPointers:

Arrays and pointers
~~~~~~~~~~~~~~~~~~~

Unlike in the C API, in the Fortran 2003 interface, arrays and pointers are
treated differently when they are return values versus arguments to a function.
Additionally, pointers which are meant to be out parameters, not arrays,
in the C API must still be declared as a rank-1 array in Fortran.
The reason for this is partially due to the Fortran 2003 standard for C bindings,
and partially due to the tool used to generate the interfaces. Regardless, the
code snippets below illustrate the differences.

C code:

.. sourcecode:: c

   N_Vector x
   realtype* xdata;
   long int leniw, lenrw;

   x = N_VNew_Serial(N);

   /* capturing a returned array/pointer */
   xdata = N_VGetArrayPointer(x)

   /* passing array/pointer to a function */
   N_VSetArrayPointer(xdata, x)

   /* pointers that are out-parameters */
   N_VSpace(x, &leniw, &lenrw);


Fortran code:

.. sourcecode:: Fortran

   type(N_Vector), pointer :: x
   real(c_double), pointer :: xdataptr(:)
   real(c_double)          :: xdata(N)
   integer(c_long)         :: leniw(1), lenrw(1)

   x => FN_VNew_Serial(x)

   ! capturing a returned array/pointer
   xdataptr => FN_VGetArrayPointer(x)

   ! passing array/pointer to a function
   call FN_VSetArrayPointer(xdata, x)

   ! pointers that are out-parameters
   call FN_VSpace(x, leniw, lenrw)


.. _Fortran2003.Differences.ProcedurePointers:

Passing procedure pointers and user data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since functions/subroutines passed to SUNDIALS will be called from within
C code, the Fortran procedure must have the attribute ``bind(C)``
Additionally, when providing them as arguments to a Fortran 2003 interface
routine, it is required to convert a procedure's Fortran address to C with
the Fortran intrinsic ``c_funloc``.

Typically when passing user data to a SUNDIALS function, a user may
simply cast some custom data structure as a ``void*``. When using the
Fortran 2003 interfaces, the same thing can be achieved. Note, the
custom data structure *does not* have to be ``bind(C)`` since
it is never accessed on the C side.

C code:

.. sourcecode:: c

   MyUserData* udata;
   void *cvode_mem;

   ierr = CVodeSetUserData(cvode_mem, udata);

Fortran code:

.. sourcecode:: Fortran

   type(MyUserData) :: udata
   type(c_ptr)      :: arkode_mem

   ierr = FARKStepSetUserData(arkode_mem, c_loc(udata))

On the other hand, Fortran users may instead choose to store problem-specific data, e.g.
problem parameters, within modules, and thus do not need the SUNDIALS-provided ``user_data``
pointers to pass such data back to user-supplied functions. These users should supply the
``c_null_ptr`` input for user_data arguments to the relevant SUNDIALS functions.

.. _Fortran2003.Differences.OptionalParameters:

Passing NULL to optional parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the SUNDIALS C API some functions have optional parameters that a
caller can pass ``NULL`` to. If the optional parameter is of a type that is
equivalent to a Fortran ``type(c_ptr)`` (see section :ref:`Fortran2003.DataTypes`),
then a Fortran user can pass the intrinsic ``c_null_ptr``. However, if the
optional parameter is of a type that is not equivalent to ``type(c_ptr)``,
then a caller must provide a Fortran pointer that is dissociated. This
is demonstrated in the code example below.

C code:

.. sourcecode:: c

   SUNLinearSolver LS;
   N_Vector x, b;

   ! SUNLinSolSolve expects a SUNMatrix or NULL
   ! as the second parameter.
   ierr = SUNLinSolSolve(LS, NULL, x, b);

Fortran code:

.. sourcecode:: Fortran

   type(SUNLinearSolver), pointer :: LS
   type(SUNMatrix), pointer :: A
   type(N_Vector), pointer :: x, b

   A => null()

   ! SUNLinSolSolve expects a type(SUNMatrix), pointer
   ! as the second parameter. Therefore, we cannot
   ! pass a c_null_ptr, rather we pass a disassociated A.
   ierr = FSUNLinSolSolve(LS, A, x, b)

.. _Fortran2003.Differences.NVectorArrays:

Working with ``N_Vector`` arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Arrays of ``N_Vector`` objects are interfaced to Fortran 2003 as opaque
``type(c_ptr)``.  As such, it is not possible to directly index an array of
``N_Vector`` objects returned by the ``N_Vector`` "VectorArray" operations, or
packages with sensitivity capabilities.  Instead, SUNDIALS provides a utility
function ``FN_VGetVecAtIndexVectorArray`` that can be called for accessing a
vector in a vector array. The example below demonstrates this:

C code:

.. sourcecode:: c

   N_Vector x;
   N_Vector* vecs;

   vecs = N_VCloneVectorArray(count, x);
   for (int i=0; i < count; ++i)
   N_VConst(vecs[i]);

Fortran code:

.. sourcecode:: Fortran

   type(N_Vector), pointer :: x, xi
   type(c_ptr)             :: vecs

   vecs = FN_VCloneVectorArray(count, x)
   do index, count
   xi => FN_VGetVecAtIndexVectorArray(vecs, index)
   call FN_VConst(xi)
   enddo

SUNDIALS also provides the functions ``FN_VSetVecAtIndexVectorArray`` and
``FN_VNewVectorArray`` for working with ``N_Vector`` arrays. These functions are
particularly useful for users of the Fortran interface to the ``NVECTOR_MANYVECTOR``
or ``NVECTOR_MPIMANYVECTOR`` when creating the subvector array.  Both of these
functions along with ``FN_VGetVecAtIndexVectorArray`` are further described in Chapter
:ref:`NVectors.utilities`.


Providing file pointers
~~~~~~~~~~~~~~~~~~~~~~~

Expert SUNDIALS users may notice that there are a few advanced functions in the SUNDIALS C
API which take a ``FILE*`` argument. Since there is no portable way to convert between a
Fortran file descriptor and a C file pointer, SUNDIALS provides two utility functions
for creating a ``FILE*`` and destroying it. These functions are defined in the module
``fsundials_futils_mod``. 

.. f:function:: FSUNDIALSFileOpen(filename, mode)

  The function allocates a ``FILE*`` by calling the C function
  ``fopen`` with the provided filename and I/O mode. 

  The function argument ``filename`` is the full path to the file and has the type
  ``character(kind=C_CHAR, len=*)``.

  The function argument ``mode`` has the type ``character(kind=C_CHAR, len=*)``.
  The string begins with one of the following characters: 
      
  * "r"  - open text file for reading
  * "r+" - open text file for reading and writing
  * "w"  - truncate text file to zero length or create it for writing
  * "w+" - open text file for reading or writing, create it if it does not exist
  * "a"  - open for appending, see documentation of ``fopen`` for your system/compiler
  * "a+ - open for reading and appending, see documentation for ``fopen`` for your system/compiler
  
  The function returns a ``type(C_PTR)`` which holds a C ``FILE*``.

.. f:subroutine:: FSUNDIALSFileClose(fp)

  The function deallocates a C ``FILE*`` by calling the C function ``fclose``
  with the provided pointer.

  The function argument ``fp`` has the type ``type(c_ptr)`` and should be
  the C ``FILE*`` obtained from ``fopen``.
  

.. _Fortran2003.Portability:

Important notes on portability
------------------------------

The SUNDIALS Fortran 2003 interface *should* be compatible with any compiler
supporting the Fortran 2003 ISO standard. However, it has only been tested
and confirmed to be working with GNU Fortran 4.9+ and Intel Fortran 18.0.1+.

Upon compilation of SUNDIALS, Fortran module (``.mod``) files are generated
for each Fortran 2003 interface. These files are highly compiler specific, and
thus it is almost always necessary to compile a consuming application with the
same compiler used to generate the modules.
