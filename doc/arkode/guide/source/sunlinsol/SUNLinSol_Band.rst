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


.. _SUNLinSol_Band:

The SUNLinSol_Band Module
======================================

The band implementation of the ``SUNLinearSolver`` module provided with
SUNDIALS, SUNLinSol_Band, is designed to be used with the
corresponding SUNMATRIX_BAND matrix type, and one of the serial or
shared-memory ``N_Vector`` implementations (NVECTOR_SERIAL, NVECTOR_OPENMP or
NVECTOR_PTHREADS).

.. _SUNLinSol_Band.Usage:

SUNLinSol_Band Usage
---------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_band.h``.  The SUNLinSol_Band module
is accessible from all SUNDIALS solvers *without*
linking to the
``libsundials_sunlinsolband`` module library.

The module SUNLinSol_Band provides the following user-callable constructor routine:

.. c:function:: SUNLinearSolver SUNLinSol_Band(N_Vector y, SUNMatrix A)

   This function creates and allocates memory for a band ``SUNLinearSolver``.
   Its arguments are an ``N_Vector`` and ``SUNMatrix``, that it uses to
   determine the linear system size and to assess compatibility with
   the linear solver implementation.

   This routine will perform consistency checks to ensure that it is
   called with consistent ``N_Vector`` and ``SUNMatrix`` implementations.
   These are currently limited to the SUNMATRIX_BAND matrix type and
   the NVECTOR_SERIAL, NVECTOR_OPENMP, and NVECTOR_PTHREADS vector types.  As
   additional compatible matrix and vector implementations are added to
   SUNDIALS, these will be included within this compatibility check.

   Additionally, this routine will verify that the input matrix ``A``
   is allocated with appropriate upper bandwidth storage for the :math:`LU`
   factorization.

   If either ``A`` or ``y`` are incompatible then this routine will
   return ``NULL``.


For backwards compatibility, we also provide the wrapper function,

.. c:function:: SUNLinearSolver SUNBandLinearSolver(N_Vector y, SUNMatrix A)

   Wrapper function for :c:func:`SUNLinSol_Band()`, with identical input and
   output arguments.


For solvers that include a Fortran interface module, the
SUNLinSol_Band module also includes the Fortran-callable
function :f:func:`FSUNBandLinSolInit()` to initialize
this SUNLinSol_Band module for a given SUNDIALS solver.

.. f:subroutine:: FSUNBandLinSolInit(CODE, IER)

   Initializes a banded ``SUNLinearSolver`` structure for
   use in a SUNDIALS package.

   This routine must be called *after* both the ``N_Vector`` and
   ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *CODE* (``int``, input) -- flag denoting the SUNDIALS solver
        this matrix will be used for: CVODE=1, IDA=2, KINSOL=3, ARKode=4.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).

Additionally, when using ARKode with a non-identity mass matrix, the
Fortran-callable function :f:func:`FSUNMassBandLinSolInit()`
initializes this SUNLinSol_Band module for solving mass matrix linear
systems.

.. f:subroutine:: FSUNMassBandLinSolInit(IER)

   Initializes a banded ``SUNLinearSolver`` structure for
   use in solving mass matrix systems in ARKode.

   This routine must be called *after* both the ``N_Vector`` and
   ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).



.. _SUNLinSol_Band.Description:

SUNLinSol_Band Description
---------------------------


The SUNLinSol_Band module defines the *content*
field of a ``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_Band {
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
  partial (row) pivoting, :math:`PA=LU`, where :math:`P` is a permutation matrix,
  :math:`L` is a lower triangular matrix with 1's on the diagonal, and :math:`U`
  is an upper triangular matrix.  This factorization is stored
  in-place on the input SUNMATRIX_BAND object :math:`A`, with pivoting
  information encoding :math:`P` stored in the ``pivots`` array.

* The "solve" call performs pivoting and forward and
  backward substitution using the stored ``pivots`` array and the
  :math:`LU` factors held in the SUNMATRIX_BAND object.

* :math:`A` must be allocated to accommodate the increase in upper
  bandwidth that occurs during factorization.  More precisely, if :math:`A`
  is a band matrix with upper bandwidth ``mu`` and lower bandwidth
  ``ml``, then the upper triangular factor :math:`U` can have upper
  bandwidth as big as ``smu = MIN(N-1,mu+ml)``. The lower triangular
  factor :math:`L` has lower bandwidth ``ml``.


The SUNLinSol_Band module defines band implementations of all
"direct" linear solver operations listed in the section
:ref:`SUNLinSol.API`:

* ``SUNLinSolGetType_Band``

* ``SUNLinSolInitialize_Band`` -- this does nothing, since all
  consistency checks are performed at solver creation.

* ``SUNLinSolSetup_Band`` -- this performs the :math:`LU` factorization.

* ``SUNLinSolSolve_Band`` -- this uses the :math:`LU` factors
  and ``pivots`` array to perform the solve.

* ``SUNLinSolLastFlag_Band``

* ``SUNLinSolSpace_Band`` -- this only returns information for
  the storage *within* the solver object, i.e. storage
  for ``N``, ``last_flag``, and ``pivots``.

* ``SUNLinSolFree_Band``
