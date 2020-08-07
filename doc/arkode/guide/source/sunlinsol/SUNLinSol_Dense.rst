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


.. _SUNLinSol_Dense:

The SUNLinSol_Dense Module
======================================

The dense implementation of the ``SUNLinearSolver`` module provided with
SUNDIALS, SUNLinSol_Dense, is designed to be used with the
corresponding SUNMATRIX_DENSE matrix type, and one of the serial or
shared-memory ``N_Vector`` implementations (NVECTOR_SERIAL, NVECTOR_OPENMP or
NVECTOR_PTHREADS).

.. _SUNLinSol_Dense.Usage:

SUNLinSol_Dense Usage
------------------------

The header file to be included when using this module is
``sunlinsol/sunlinsol_dense.h``.  The SUNLinSol_Dense module is
accessible from all SUNDIALS solvers *without*
linking to the ``libsundials_sunlinsoldense`` module library.


The module SUNLinSol_Dense provides the following user-callable constructor routine:


.. c:function:: SUNLinearSolver SUNLinSol_Dense(N_Vector y, SUNMatrix A)

   This function creates and allocates memory for a dense ``SUNLinearSolver``.
   Its arguments are an ``N_Vector`` and ``SUNMatrix``, that it uses to
   determine the linear system size and to assess compatibility with
   the linear solver implementation.

   This routine will perform consistency checks to ensure that it is
   called with consistent ``N_Vector`` and ``SUNMatrix`` implementations.
   These are currently limited to the SUNMATRIX_DENSE matrix type and
   the NVECTOR_SERIAL, NVECTOR_OPENMP, and NVECTOR_PTHREADS vector types.  As
   additional compatible matrix and vector implementations are added to
   SUNDIALS, these will be included within this compatibility check.

   If either ``A`` or ``y`` are incompatible then this routine will
   return ``NULL``.


For backwards compatibility, we also provide the wrapper function,

.. c:function:: SUNLinearSolver SUNDenseLinearSolver(N_Vector y, SUNMatrix A)

   Wrapper function for :c:func:`SUNLinSol_Dense()`, with identical input and
   output arguments


For solvers that include a Fortran interface module, the
SUNLinSol_Dense module also includes the Fortran-callable
function :f:func:`FSUNDenseLinSolInit()` to initialize
this SUNLinSol_Dense module for a given SUNDIALS solver.

.. f:subroutine:: FSUNDenseLinSolInit(CODE, IER)

   Initializes a dense ``SUNLinearSolver`` structure for
   use in a SUNDIALS package.

   This routine must be called *after* both the ``N_Vector`` and
   ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *CODE* (``int``, input) -- flag denoting the SUNDIALS solver
        this matrix will be used for: CVODE=1, IDA=2, KINSOL=3, ARKode=4.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).

Additionally, when using ARKode with a non-identity mass matrix, the
Fortran-callable function :f:func:`FSUNMassDenseLinSolInit()`
initializes this SUNLinSol_Dense module for solving mass matrix linear
systems.

.. f:subroutine:: FSUNMassDenseLinSolInit(IER)

   Initializes a dense ``SUNLinearSolver`` structure for
   use in solving mass matrix systems in ARKode.

   This routine must be called *after* both the ``N_Vector`` and
   ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).



.. _SUNLinSol_Dense.Description:

SUNLinSol_Dense Description
-----------------------------


The SUNLinSol_Dense module defines the *content*
field of a ``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_Dense {
     sunindextype N;
     sunindextype *pivots;
     sunindextype last_flag;
   };

These entries of the *content* field contain the following
information:

* ``N`` - size of the linear system,

* ``pivots`` - index array for partial pivoting in LU factorization,

* ``last_flag`` - last error return flag from internal function evaluations.


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


The SUNLinSol_Dense module defines dense implementations of all
"direct" linear solver operations listed in the section
:ref:`SUNLinSol.API`:

* ``SUNLinSolGetType_Dense``

* ``SUNLinSolInitialize_Dense`` -- this does nothing, since all
  consistency checks are performed at solver creation.

* ``SUNLinSolSetup_Dense`` -- this performs the :math:`LU` factorization.

* ``SUNLinSolSolve_Dense`` -- this uses the :math:`LU` factors
  and ``pivots`` array to perform the solve.

* ``SUNLinSolLastFlag_Dense``

* ``SUNLinSolSpace_Dense`` -- this only returns information for
  the storage *within* the solver object, i.e. storage
  for ``N``, ``last_flag``, and ``pivots``.

* ``SUNLinSolFree_Dense``
