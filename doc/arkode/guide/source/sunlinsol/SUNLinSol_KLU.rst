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


.. _SUNLinSol_KLU:

The SUNLinSol_KLU Module
======================================

The KLU implementation of the ``SUNLinearSolver`` module provided with
SUNDIALS, SUNLinSol_KLU, is designed to be used with the
corresponding SUNMATRIX_SPARSE matrix type, and one of the serial or
shared-memory ``N_Vector`` implementations (NVECTOR_SERIAL, NVECTOR_OPENMP, or
NVECTOR_PTHREADS).

.. _SUNLinSol_KLU.Usage:

SUNLinSol_KLU Usage
------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_klu.h``.  The installed module
library to link to is ``libsundials_sunlinsolklu`` *.lib*
where *.lib* is typically ``.so`` for shared libraries and
``.a`` for static libraries.

The module SUNLinSol_KLU provides the following additional
user-callable routines:


.. c:function:: SUNLinearSolver SUNLinSol_KLU(N_Vector y, SUNMatrix A)

   This constructor function creates and allocates memory for a SUNLinSol_KLU
   object.  Its arguments are an ``N_Vector`` and ``SUNMatrix``, that it
   uses to determine the linear system size and to assess compatibility
   with the linear solver implementation.

   This routine will perform consistency checks to ensure that it is
   called with consistent ``N_Vector`` and ``SUNMatrix`` implementations.
   These are currently limited to the SUNMATRIX_SPARSE matrix type
   (using either CSR or CSC storage formats) and the NVECTOR_SERIAL,
   NVECTOR_OPENMP, and NVECTOR_PTHREADS vector types.  As additional
   compatible matrix and vector implementations are added to
   SUNDIALS, these will be included within this compatibility
   check.

   If either ``A`` or ``y`` are incompatible then this routine will
   return ``NULL``.


.. c:function:: int SUNLinSol_KLUReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz, int reinit_type)

   This function reinitializes memory and flags for a new factorization
   (symbolic and numeric) to be conducted at the next solver setup
   call.  This routine is useful in the cases where the number of
   nonzeroes has changed or if the structure of the linear system has
   changed which would require a new symbolic (and numeric
   factorization).

   The ``reinit_type`` argument governs the level of
   reinitialization.  The allowed values are:

   1. The Jacobian matrix will be destroyed and a new one will be
      allocated based on the ``nnz`` value passed to this call.  New
      symbolic and numeric factorizations will be completed at the next
      solver setup.

   2. Only symbolic and numeric factorizations will be completed.
      It is assumed that the Jacobian size has not exceeded the size of
      ``nnz`` given in the sparse matrix provided to the original
      constructor routine (or the previous ``SUNKLUReInit`` call).

   This routine assumes no other changes to solver use are necessary.

   The return values from this function are ``SUNLS_MEM_NULL``
   (either ``S`` or ``A`` are ``NULL``), ``SUNLS_ILL_INPUT``
   (``A`` does not have type ``SUNMATRIX_SPARSE`` or
   ``reinit_type`` is invalid), ``SUNLS_MEM_FAIL`` (reallocation
   of the sparse matrix failed) or ``SUNLS_SUCCESS``.


.. c:function:: int SUNLinSol_KLUSetOrdering(SUNLinearSolver S, int ordering_choice)

   This function sets the ordering used by KLU for reducing fill in
   the linear solve.  Options for ``ordering_choice`` are:

   0. AMD,

   1. COLAMD, and

   2. the natural ordering.

   The default is 1 for COLAMD.

   The return values from this function are ``SUNLS_MEM_NULL``
   (``S`` is ``NULL``), ``SUNLS_ILL_INPUT``
   (invalid ``ordering_choice``), or ``SUNLS_SUCCESS``.


.. c:function:: sun_klu_symbolic* SUNLinSol_KLUGetSymbolic(SUNLinearSolver S)

   This function returns a pointer to the KLU symbolic factorization
   stored in the SUNLinSol_KLU ``content`` structure.

   When SUNDIALS is compiled with 32-bit indices (``SUNDIALS_INDEX_SIZE=32``),
   ``sun_klu_symbolic`` is mapped to the KLU type ``klu_symbolic``; when
   SUNDIALS compiled with 64-bit indices (``SUNDIALS_INDEX_SIZE=64``) this is
   mapped to the KLU type ``klu_l_symbolic``.


.. c:function:: sun_klu_numeric* SUNLinSol_KLUGetNumeric(SUNLinearSolver S)

   This function returns a pointer to the KLU numeric factorization
   stored in the SUNLinSol_KLU ``content`` structure.

   When SUNDIALS is compiled with 32-bit indices (``SUNDIALS_INDEX_SIZE=32``),
   ``sun_klu_numeric`` is mapped to the KLU type ``klu_numeric``; when
   SUNDIALS is compiled with 64-bit indices (``SUNDIALS_INDEX_SIZE=64``) this is
   mapped to the KLU type ``klu_l_numeric``.


.. c:function:: sun_klu_common* SUNLinSol_KLUGetCommon(SUNLinearSolver S)

   This function returns a pointer to the KLU common structure
   stored in the SUNLinSol_KLU ``content`` structure.

   When SUNDIALS is compiled with 32-bit indices (``SUNDIALS_INDEX_SIZE=32``),
   ``sun_klu_common`` is mapped to the KLU type ``klu_common``; when
   SUNDIALS is compiled with 64-bit indices  (``SUNDIALS_INDEX_SIZE=64``) this is
   mapped to the KLU type ``klu_l_common``.


For backwards compatibility, we also provide the wrapper functions,
each with identical input and output arguments to the routines that
they wrap:

.. c:function:: SUNLinearSolver SUNKLU(N_Vector y, SUNMatrix A)

   Wrapper function for :c:func:`SUNLinSol_KLU()`

.. c:function:: int SUNKLUReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz, int reinit_type)

   Wrapper function for :c:func:`SUNLinSol_KLUReInit()`

.. c:function:: int SUNKLUSetOrdering(SUNLinearSolver S, int ordering_choice)

   Wrapper function for :c:func:`SUNLinSol_KLUSetOrdering()`



For solvers that include a Fortran interface module, the
SUNLinSol_KLU module also includes the Fortran-callable
function :f:func:`FSUNKLUInit()` to initialize this SUNLinSol_KLU
module for a given SUNDIALS solver.

.. f:subroutine:: FSUNKLUInit(CODE, IER)

   Initializes a KLU sparse ``SUNLinearSolver`` structure for
   use in a SUNDIALS package.

   This routine must be called *after* both the ``N_Vector`` and
   ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *CODE* (``int``, input) -- flag denoting the SUNDIALS solver
        this matrix will be used for: CVODE=1, IDA=2, KINSOL=3, ARKode=4.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).


Additionally, when using ARKode with a non-identity mass matrix, the
Fortran-callable function :f:func:`FSUNMassKLUInit()` initializes this
SUNLinSol_KLU module for solving mass matrix linear systems.

.. f:subroutine:: FSUNMassKLUInit(IER)

   Initializes a KLU sparse ``SUNLinearSolver`` structure for
   use in solving mass matrix systems in ARKode.

   This routine must be called *after* both the ``N_Vector`` and
   ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).

The :c:func:`SUNLinSol_KLUReInit()` and :c:func:`SUNLinSol_KLUSetOrdering()`
routines also support Fortran interfaces for the system and mass
matrix solvers:

.. f:subroutine:: FSUNKLUReInit(CODE, NNZ, REINIT_TYPE, IER)

   Fortran interface to :c:func:`SUNLinSol_KLUReInit()` for system
   linear solvers.

   This routine must be called *after*
   :f:func:`FSUNKLUInit()` has been called.

   **Arguments:** *NNZ* should have type ``long int``, all others
   should have type ``int``; all arguments have meanings identical to
   those listed above.


.. f:subroutine:: FSUNMassKLUReInit(NNZ, REINIT_TYPE, IER)

   Fortran interface to :c:func:`SUNLinSol_KLUReInit()` for mass matrix
   linear solvers in ARKode.

   This routine must be called *after*
   :f:func:`FSUNMassKLUInit()` has been called.

   **Arguments:** *NNZ* should have type ``long int``, all others
   should have type ``int``; all arguments have meanings identical to
   those listed above.

.. f:subroutine:: FSUNKLUSetOrdering(CODE, ORDERING, IER)

   Fortran interface to :c:func:`SUNLinSol_KLUSetOrdering()` for system
   linear solvers.

   This routine must be called *after* :f:func:`FSUNKLUInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.

.. f:subroutine:: FSUNMassKLUSetOrdering(ORDERING, IER)

   Fortran interface to :c:func:`SUNLinSol_KLUSetOrdering()` for mass matrix
   linear solvers in ARKode.

   This routine must be called *after* :f:func:`FSUNMassKLUInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.





.. _SUNLinSol_KLU.Description:

SUNLinSol_KLU Description
--------------------------


The SUNLinSol_KLU module defines the *content*
field of a ``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_KLU {
     int              last_flag;
     int              first_factorize;
     sun_klu_symbolic *symbolic;
     sun_klu_numeric  *numeric;
     sun_klu_common   common;
     sunindextype     (*klu_solver)(sun_klu_symbolic*, sun_klu_numeric*,
                                    sunindextype, sunindextype,
                                    double*, sun_klu_common*);
   };

These entries of the *content* field contain the following
information:

* ``last_flag`` - last error return flag from internal function
  evaluations,

* ``first_factorize`` - flag indicating whether the factorization
  has ever been performed,

* ``Symbolic`` - KLU storage structure for symbolic
  factorization components, with underlying type ``klu_symbolic``
  or ``klu_l_symbolic``, depending on whether SUNDIALS was
  installed with 32-bit versus 64-bit indices, respectively,

* ``Numeric`` - KLU storage structure for numeric factorization
  components, with underlying type ``klu_numeric``
  or ``klu_l_numeric``, depending on whether SUNDIALS was
  installed with 32-bit versus 64-bit indices, respectively,

* ``Common`` - storage structure for common KLU solver
  components, with underlying type ``klu_common``
  or ``klu_l_common``, depending on whether SUNDIALS was
  installed with 32-bit versus 64-bit indices, respectively,

* ``klu_solver`` -- pointer to the appropriate KLU solver function
  (depending on whether it is using a CSR or CSC sparse matrix, and
  on whether SUNDIALS was installed with 32-bit or 64-bit indices).


The SUNLinSol_KLU module is a ``SUNLinearSolver`` wrapper for
the KLU sparse matrix factorization and solver library written by Tim
Davis ([KLU]_, [DP2010]_).  In order to use the
SUNLinSol_KLU interface to KLU, it is assumed that KLU has
been installed on the system prior to installation of SUNDIALS, and
that SUNDIALS has been configured appropriately to link with KLU
(see section :ref:`Installation.CMake.ExternalLibraries` for details).
Additionally, this wrapper only supports double-precision
calculations, and therefore cannot be compiled if SUNDIALS is
configured to have ``realtype`` set to either ``extended`` or
``single`` (see section :ref:`ARKStep_CInterface.DataTypes` for
details). Since the KLU library supports both 32-bit and 64-bit
integers, this interface will be compiled for either of the available
``sunindextype`` options.

The KLU library has a symbolic factorization routine that computes
the permutation of the linear system matrix to block triangular form
and the permutations that will pre-order the diagonal blocks (the only
ones that need to be factored) to reduce fill-in (using AMD, COLAMD,
CHOLAMD, natural, or an ordering given by the user).  Of these
ordering choices, the default value in the SUNLinSol_KLU
module is the COLAMD ordering.

KLU breaks the factorization into two separate parts.  The first is
a symbolic factorization and the second is a numeric factorization
that returns the factored matrix along with final pivot information.
KLU also has a refactor routine that can be called instead of the numeric
factorization.  This routine will reuse the pivot information.  This routine
also returns diagnostic information that a user can examine to determine if
numerical stability is being lost and a full numerical factorization should
be done instead of the refactor.

Since the linear systems that arise within the context of SUNDIALS
calculations will typically have identical sparsity patterns, the
SUNLinSol_KLU module is constructed to perform the
following operations:

* The first time that the "setup" routine is called, it
  performs the symbolic factorization, followed by an initial
  numerical factorization.

* On subsequent calls to the "setup" routine, it calls the
  appropriate KLU "refactor" routine, followed by estimates of
  the numerical conditioning using the relevant "rcond", and if
  necessary "condest", routine(s).  If these estimates of the
  condition number are larger than :math:`\varepsilon^{-2/3}` (where
  :math:`\varepsilon` is the double-precision unit roundoff), then a new
  factorization is performed.

* The module includes the routine ``SUNKLUReInit``, that
  can be called by the user to force a full refactorization at the
  next "setup" call.

* The "solve" call performs pivoting and forward and
  backward substitution using the stored KLU data structures.  We
  note that in this solve KLU operates on the native data arrays
  for the right-hand side and solution vectors, without requiring
  costly data copies.


The SUNLinSol_KLU module defines implementations of all
"direct" linear solver operations listed in the section
:ref:`SUNLinSol.API`:

* ``SUNLinSolGetType_KLU``

* ``SUNLinSolInitialize_KLU`` -- this sets the
  ``first_factorize`` flag to 1, forcing both symbolic and numerical
  factorizations on the subsequent "setup" call.

* ``SUNLinSolSetup_KLU`` -- this performs either a :math:`LU`
  factorization or refactorization of the input matrix.

* ``SUNLinSolSolve_KLU`` -- this calls the appropriate KLU
  solve routine to utilize the :math:`LU` factors to solve the linear
  system.

* ``SUNLinSolLastFlag_KLU``

* ``SUNLinSolSpace_KLU`` -- this only returns information for
  the storage within the solver *interface*, i.e. storage for the
  integers ``last_flag`` and ``first_factorize``.  For additional
  space requirements, see the KLU documentation.

* ``SUNLinSolFree_KLU``
