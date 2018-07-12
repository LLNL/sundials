..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _SUNLinSol_LapackDense:

The SUNLINSOL_LAPACKDENSE Module
======================================

The LAPACK dense implementation of the ``SUNLinearSolver`` module provided
with SUNDIALS, SUNLINSOL_LAPACKDENSE, is designed to be used with the
corresponding SUNMATRIX_DENSE matrix type, and one of the serial or
shared-memory ``N_Vector`` implementations (NVECTOR_SERIAL, NVECTOR_OPENMP, or
NVECTOR_PTHREADS).  The SUNLINSOL_LAPACKDENSE module defines the
*content* field of a ``SUNLinearSolver`` to be the following
structure:

.. code-block:: c

   struct _SUNLinearSolverContent_Dense {
     sunindextype N;
     sunindextype *pivots;
     long int last_flag;
   };

These entries of the *content* field contain the following
information:

* ``N`` - size of the linear system,

* ``pivots`` - index array for partial pivoting in LU
  factorization,

* ``last_flag`` - last error return flag from internal function
  evaluations.


The SUNLINSOL_LAPACKDENSE module is a ``SUNLinearSolver`` wrapper for
the LAPACK dense matrix factorization and solve routines, ``*GETRF``
and ``*GETRS``, where ``*`` is either ``D`` or ``S``, depending on
whether SUNDIALS was configured to have ``realtype`` set to
``double`` or ``single``, respectively (see section
:ref:`ARKStep_CInterface.DataTypes` for details).  In order to use the
SUNLINSOL_LAPACKDENSE module it is assumed that LAPACK has been
installed on the system prior to installation of
SUNDIALS, and that SUNDIALS has been configured appropriately to
link with LAPACK (see section
:ref:`Installation.CMake.ExternalLibraries` for details).
We note that since there do not exist 128-bit floating-point
factorization and solve routines in LAPACK, this interface cannot be
compiled when using ``extended`` precision for ``realtype``.
Similarly, since there do not exist 64-bit integer LAPACK routines,
the SUNLINSOL_LAPACKDENSE module also cannot be compiled when using
``int64_t`` for the ``sunindextype``.

This solver is constructed to perform the following operations:

* The "setup" call performs a :math:`LU` factorization with
  partial (row) pivoting (:math:`\mathcal O(N^3)` cost),
  :math:`PA=LU`, where :math:`P` is a permutation matrix, :math:`L` is
  a lower triangular matrix with 1's on the diagonal, and :math:`U` is
  an upper triangular matrix.  This factorization is stored in-place
  on the input SUNMATRIX_DENSE object :math:`A`, with pivoting
  information encoding :math:`P` stored in the ``pivots`` array.

* The "solve" call performs pivoting and forward and
  backward substitution using the stored ``pivots`` array and the
  :math:`LU` factors held in the SUNMATRIX_DENSE object
  (:math:`\mathcal O(N^2)` cost).

The header file to be included when using this module
is ``sunlinsol/sunlinsol_lapackdense.h``.

The SUNLINSOL_LAPACKDENSE module defines dense implementations of all
"direct" linear solver operations listed in the section
:ref:`SUNLinSol.Ops`:

* ``SUNLinSolGetType_LapackDense``

* ``SUNLinSolInitialize_LapackDense`` -- this does nothing, since all
  consistency checks are performed at solver creation.

* ``SUNLinSolSetup_LapackDense`` -- this calls either
  ``DGETRF`` or ``SGETRF`` to perform the :math:`LU` factorization.

* ``SUNLinSolSolve_LapackDense`` -- this calls either
  ``DGETRS`` or ``SGETRS`` to use the :math:`LU` factors and
  ``pivots`` array to perform the solve.

* ``SUNLinSolLastFlag_LapackDense``

* ``SUNLinSolSpace_LapackDense`` -- this only returns information for
  the storage *within* the solver object, i.e. storage
  for ``N``, ``last_flag``, and ``pivots``.

* ``SUNLinSolFree_LapackDense``

The module SUNLINSOL_LAPACKDENSE provides the following additional
user-callable constructor routine:


.. c:function:: SUNLinearSolver SUNLapackDense(N_Vector y, SUNMatrix A)

   This function creates and allocates memory for a LAPACK dense
   ``SUNLinearSolver``.  Its arguments are an ``N_Vector`` and
   ``SUNMatrix``, that it uses to determine the linear system size and
   to assess compatibility with the linear solver implementation.

   This routine will perform consistency checks to ensure that it is
   called with consistent ``N_Vector`` and ``SUNMatrix`` implementations.
   These are currently limited to the SUNMATRIX_DENSE matrix type and
   the NVECTOR_SERIAL, NVECTOR_OPENMP, and NVECTOR_PTHREADS vector
   types.  As additional compatible matrix and vector implementations
   are added to SUNDIALS, these will be included within this
   compatibility check.

   If either ``A`` or ``y`` are incompatible then this routine will
   return ``NULL``.

For solvers that include a Fortran interface module, the
SUNLINSOL_LAPACKDENSE module also includes the Fortran-callable
function :f:func:`FSUNLapackDenseInit()` to initialize
this SUNLINSOL_LAPACKDENSE module for a given SUNDIALS solver.

.. f:subroutine:: FSUNLapackDenseInit(CODE, IER)

   Initializes a dense LAPACK ``SUNLinearSolver`` structure for
   use in a SUNDIALS package.

   This routine must be called *after* both the ``N_Vector`` and
   ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *CODE* (``int``, input) -- flag denoting the SUNDIALS solver
        this matrix will be used for: CVODE=1, IDA=2, KINSOL=3, ARKode=4.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).

Additionally, when using ARKode with a non-identity mass matrix, the
Fortran-callable function :f:func:`FSUNMassLapackDenseInit()`
initializes this SUNLINSOL_LAPACKDENSE module for solving mass matrix
linear systems.

.. f:subroutine:: FSUNMassLapackDenseInit(IER)

   Initializes a dense LAPACK ``SUNLinearSolver`` structure for
   use in solving mass matrix systems in ARKode.

   This routine must be called *after* both the ``N_Vector`` and
   ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).
