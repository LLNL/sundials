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


.. _SUNLinSol_SuperLUMT:

The SUNLinSol_SuperLUMT Module
======================================

The SuperLU_MT implementation of the ``SUNLinearSolver`` module
provided with SUNDIALS, SUNLinSol_SuperLUMT, is designed to be used
with the corresponding SUNMATRIX_SPARSE matrix type, and one of the
serial or shared-memory ``N_Vector`` implementations (NVECTOR_SERIAL,
NVECTOR_OPENMP, or NVECTOR_PTHREADS).  While these are compatible, it
is not recommended to use a threaded vector module with
SUNLinSol_SuperLUMT unless it is the NVECTOR_OPENMP module and the
SuperLU_MT library has also been compiled with OpenMP.


.. _SUNLinSol_SuperLUMT.Usage:

SUNLinSol_SuperLUMT Usage
-----------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_superlumt.h``.  The installed module
library to link to is ``libsundials_sunlinsolsuperlumt`` *.lib*
where *.lib* is typically ``.so`` for shared libraries and
``.a`` for static libraries.

The module SUNLinSol_SuperLUMT provides the following user-callable routines:

.. c:function:: SUNLinearSolver SUNLinSol_SuperLUMT(N_Vector y, SUNMatrix A, int num_threads)

   This constructor function creates and allocates memory for a
   SUNLinSol_SuperLUMT object.  Its arguments are an ``N_Vector``, a
   ``SUNMatrix``, and a desired number of threads (OpenMP or Pthreads,
   depending on how SuperLU_MT was installed) to use during the
   factorization steps. This routine analyzes the input matrix and
   vector to determine the linear system size and to assess
   compatibility with the SuperLU_MT library.

   This routine will perform consistency checks to ensure that it is
   called with consistent ``N_Vector`` and ``SUNMatrix``
   implementations.  These are currently limited to the
   SUNMATRIX_SPARSE matrix type (using either CSR or CSC storage
   formats) and the NVECTOR_SERIAL, NVECTOR_OPENMP, and
   NVECTOR_PTHREADS vector types.  As additional compatible matrix and
   vector implementations are added to SUNDIALS, these will be
   included within this compatibility check.

   If either ``A`` or ``y`` are incompatible then this routine will
   return ``NULL``.  The ``num_threads`` argument is not checked
   and is passed directly to SuperLU_MT routines.


.. c:function:: int SUNLinSol_SuperLUMTSetOrdering(SUNLinearSolver S, int ordering_choice)

   This function sets the ordering used by SuperLU_MT for reducing fill in
   the linear solve.  Options for ``ordering_choice`` are:

   0. natural ordering

   1. minimal degree ordering on :math:`A^TA`

   2. minimal degree ordering on :math:`A^T+A`

   3. COLAMD ordering for unsymmetric matrices

   The default is 3 for COLAMD.

   The return values from this function are ``SUNLS_MEM_NULL``
   (``S`` is ``NULL``), ``SUNLS_ILL_INPUT``
   (invalid ``ordering_choice``), or ``SUNLS_SUCCESS``.


For backwards compatibility, we also provide the wrapper functions,
each with identical input and output arguments to the routines that
they wrap:

.. c:function:: SUNLinearSolver SUNSuperLUMT(N_Vector y, SUNMatrix A, int num_threads)

   Wrapper for :c:func:`SUNLinSol_SuperLUMT()`.

and

.. c:function:: int SUNSuperLUMTSetOrdering(SUNLinearSolver S, int ordering_choice)

   Wrapper for :c:func:`SUNLinSol_SuperLUMTSetOrdering()`.


For solvers that include a Fortran interface module, the
SUNLinSol_SuperLUMT module also includes the Fortran-callable
function :f:func:`FSUNSuperLUMTInit()` to initialize this
SUNLinSol_SuperLUMT module for a given SUNDIALS solver.

.. f:subroutine:: FSUNSuperLUMTInit(CODE, NUM_THREADS, IER)

   Initializes a SuperLU_MT sparse ``SUNLinearSolver`` structure for
   use in a SUNDIALS package.

   This routine must be called *after* both the ``N_Vector`` and
   ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *CODE* (``int``, input) -- flag denoting the SUNDIALS solver
        this matrix will be used for: CVODE=1, IDA=2, KINSOL=3, ARKode=4.
      * *NUM_THREADS* (``int``, input) -- desired number of
        OpenMP/Pthreads threads to use in the factorization.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).

Additionally, when using ARKode with a non-identity mass matrix, the
Fortran-callable function :f:func:`FSUNMassSuperLUMTInit()`
initializes this SUNLinSol_SuperLUMT module for solving mass matrix
linear systems.

.. f:subroutine:: FSUNMassSuperLUMTInit(NUM_THREADS, IER)

   Initializes a SuperLU_MT sparse ``SUNLinearSolver`` structure for
   use in solving mass matrix systems in ARKode.

   This routine must be called *after* both the ``N_Vector`` and
   the mass ``SUNMatrix`` objects have been initialized.

   **Arguments:**
      * *NUM_THREADS* (``int``, input) -- desired number of
        OpenMP/Pthreads threads to use in the factorization.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).

The :c:func:`SUNLinSol_SuperLUMTSetOrdering()` routine also supports Fortran
interfaces for the system and mass matrix solvers:

.. f:subroutine:: FSUNSuperLUMTSetOrdering(CODE, ORDERING, IER)

   Fortran interface to :c:func:`SUNLinSol_SuperLUMTSetOrdering()` for system
   linear solvers.

   This routine must be called *after*
   :f:func:`FSUNSuperLUMTInit()` has been called

   **Arguments:** all should have type ``int`` and have meanings
   identical to those listed above

.. f:subroutine:: FSUNMassSuperLUMTSetOrdering(ORDERING, IER)

   Fortran interface to :c:func:`SUNLinSol_SuperLUMTSetOrdering()` for mass
   matrix linear solves in ARKode.

   This routine must be called *after*
   :f:func:`FSUNMassSuperLUMTInit()` has been called

   **Arguments:** all should have type ``int`` and have meanings
   identical to those listed above




.. _SUNLinSol_SuperLUMT.Description:

SUNLinSol_SuperLUMT Description
----------------------------------

The SUNLinSol_SuperLUMT module defines the *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_SuperLUMT {
     int          last_flag;
     int          first_factorize;
     SuperMatrix  *A, *AC, *L, *U, *B;
     Gstat_t      *Gstat;
     sunindextype *perm_r, *perm_c;
     sunindextype N;
     int          num_threads;
     realtype     diag_pivot_thresh;
     int          ordering;
     superlumt_options_t *options;
   };

These entries of the *content* field contain the following
information:

* ``last_flag`` - last error return flag from internal function
  evaluations,

* ``first_factorize`` - flag indicating whether the factorization
  has ever been performed,

* ``A, AC, L, U, B`` - ``SuperMatrix`` pointers used in solve,

* ``Gstat`` - ``GStat_t`` object used in solve,

* ``perm_r, perm_c`` - permutation arrays used in solve,

* ``N`` - size of the linear system,

* ``num_threads`` - number of OpenMP/Pthreads threads to use,

* ``diag_pivot_thresh`` - threshold on diagonal pivoting,

* ``ordering`` - flag for which reordering algorithm to use,

* ``options`` - pointer to SuperLU_MT options structure.

The SUNLinSol_SuperLUMT module is a ``SUNLinearSolver`` wrapper for
the SuperLU_MT sparse matrix factorization and solver library
written by X. Sherry Li ([SuperLUMT]_, [L2005]_, [DGL1999]_).  The
package performs matrix factorization using threads to enhance
efficiency in shared memory parallel environments.  It should be noted
that threads are only used in the factorization step.  In
order to use the SUNLinSol_SuperLUMT interface to SuperLU_MT, it is
assumed that SuperLU_MT has been installed on the system prior to
installation of SUNDIALS, and that SUNDIALS has been configured
appropriately to link with SuperLU_MT (see section
:ref:`Installation.CMake.ExternalLibraries` for details).
Additionally, this wrapper only supports single- and
double-precision calculations, and therefore cannot be compiled if
SUNDIALS is configured to have ``realtype`` set to ``extended``
(see section :ref:`ARKStep_CInterface.DataTypes` for details).  Moreover,
since the SuperLU_MT library may be installed to support either 32-bit
or 64-bit integers, it is assumed that the SuperLU_MT library is
installed using the same integer precision as the SUNDIALS
``sunindextype`` option.

The SuperLU_MT library has a symbolic factorization routine that
computes the permutation of the linear system matrix to reduce fill-in
on subsequent :math:`LU` factorizations (using COLAMD, minimal degree
ordering on :math:`A^T*A`, minimal degree ordering on :math:`A^T+A`,
or natural ordering).  Of these ordering choices, the default value in
the SUNLinSol_SuperLUMT module is the COLAMD ordering.

Since the linear systems that arise within the context of SUNDIALS
calculations will typically have identical sparsity patterns, the
SUNLinSol_SuperLUMT module is constructed to perform the
following operations:

* The first time that the "setup" routine is called, it
  performs the symbolic factorization, followed by an initial
  numerical factorization.

* On subsequent calls to the "setup" routine, it skips the
  symbolic factorization, and only refactors the input matrix.

* The "solve" call performs pivoting and forward and
  backward substitution using the stored SuperLU_MT data
  structures.  We note that in this solve SuperLU_MT operates on the
  native data arrays for the right-hand side and solution vectors,
  without requiring costly data copies.


The SUNLinSol_SuperLUMT module defines implementations of all
"direct" linear solver operations listed in the section
:ref:`SUNLinSol.API`:


* ``SUNLinSolGetType_SuperLUMT``

* ``SUNLinSolInitialize_SuperLUMT`` -- this sets the
  ``first_factorize`` flag to 1 and resets the internal SuperLU_MT
  statistics variables.

* ``SUNLinSolSetup_SuperLUMT`` -- this performs either a :math:`LU`
  factorization or refactorization of the input matrix.

* ``SUNLinSolSolve_SuperLUMT`` -- this calls the appropriate
  SuperLU_MT solve routine to utilize the :math:`LU` factors to solve the
  linear system.

* ``SUNLinSolLastFlag_SuperLUMT``

* ``SUNLinSolSpace_SuperLUMT`` -- this only returns information for
  the storage within the solver *interface*, i.e. storage for the
  integers ``last_flag`` and ``first_factorize``.  For additional
  space requirements, see the SuperLU_MT documentation.

* ``SUNLinSolFree_SuperLUMT``
