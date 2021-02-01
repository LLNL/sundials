..
   Daniel R. Reynolds @ SMU
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


.. _SUNLinSol.API:

The SUNLinearSolver API
=============================

The SUNLinSol API defines several linear solver operations that enable
SUNDIALS packages to utilize any SUNLinSol implementation that
provides the required functions. These functions can be divided into
three categories. The first are the core linear solver functions. The
second group of functions consists of set routines to supply the
linear solver with functions provided by the SUNDIALS time integrators
and to modify solver parameters. The final group consists of get
routines for retrieving linear solver statistics. All of these
functions are defined in the header file
``sundials/sundials_linearsolver.h``.



.. _SUNLinSol.CoreFn:

SUNLinearSolver core functions
-----------------------------------------------------

The core linear solver functions consist of two required functions to get the
linear solver type (:c:func:`SUNLinSolGetType`) and solve the linear system
:math:`Ax=b` (:c:func:`SUNLinSolSolve`). The remaining functions are for
getting the solver ID (:c:func:`SUNLinSolGetID`), initializing the linear solver
object once all solver-specific options have been set
(:c:func:`SUNLinSolInitialize`), setting up the linear solver object to utilize
an updated matrix :math:`A` (:c:func:`SUNLinSolSetup`), and for destroying the
linear solver object (:c:func:`SUNLinSolFree`) are optional.


.. c:function:: SUNLinearSolver_Type SUNLinSolGetType(SUNLinearSolver LS)

   Returns the type identifier for the linear solver *LS*. It is used to
   determine the solver type (direct, iterative, or matrix-iterative) from
   the abstract ``SUNLinearSolver`` interface.  Returned values are
   one of the following:

   * ``SUNLINEARSOLVER_DIRECT`` -- ``0``, the SUNLinSol module
     requires a matrix, and computes an 'exact' solution to the linear
     system defined by that matrix.

   * ``SUNLINEARSOLVER_ITERATIVE`` -- ``1``, the SUNLinSol module does
     not require a matrix (though one may be provided), and computes
     an inexact solution to the linear system using a matrix-free
     iterative algorithm. That is it solves the linear system defined
     by the package-supplied ``ATimes`` routine (see
     :c:func:`SUNLinSolSetATimes()` below), even if that linear system
     differs from the one encoded in the matrix object (if one is
     provided). As the solver computes the solution only inexactly (or
     may diverge), the linear solver should check for solution
     convergence/accuracy as appropriate.

   * ``SUNLINEARSOLVER_MATRIX_ITERATIVE`` -- ``2``, the SUNLinSol
     module requires a matrix, and computes an inexact solution to the
     linear system defined by that matrix using an iterative
     algorithm. That is it solves the linear system defined by the
     matrix object even if that linear system differs from that
     encoded by the package-supplied ``ATimes`` routine. As the solver
     computes the solution only inexactly (or may diverge), the linear
     solver should check for solution convergence/accuracy as
     appropriate.

   Usage:

   .. code-block:: c

      type = SUNLinSolGetType(LS);

   Notes: See section :ref:`SUNLinSol.Intended` for more information
   on intended use cases corresponding to the linear solver type.

.. c:function:: SUNLinearSolver_ID SUNLinSolGetID(SUNLinearSolver LS)

   Returns the identifier for the linear solver *LS*. It is recommended that a
   user-supplied ``SUNLinearSolver`` implementation return the
   ``SUNLINEARSOLVER_CUSTOM`` identifier.

   Usage:

   .. code-block:: c

      id = SUNLinSolGetID(LS);


.. c:function:: int SUNLinSolInitialize(SUNLinearSolver LS)

   Performs linear solver initialization (assuming that all
   solver-specific options have been set).  This should return zero for a
   successful call, and a negative value for a failure, ideally
   returning one of the generic error codes listed in section
   :ref:`SUNLinSol.ErrorCodes`.

   Usage:

   .. code-block:: c

      retval = SUNLinSolInitialize(LS);


.. c:function:: int SUNLinSolSetup(SUNLinearSolver LS, SUNMatrix A)

   Performs any linear solver setup needed, based on an updated system
   ``SUNMatrix`` *A*.  This may be called frequently (e.g., with a full
   Newton method) or infrequently (for a modified Newton method), based
   on the type of integrator and/or nonlinear solver requesting the
   solves.  This should return zero for a successful call, a positive
   value for a recoverable failure and a negative value for an
   unrecoverable failure, ideally returning one of the generic error
   codes listed in section :ref:`SUNLinSol.ErrorCodes`.

   Usage:

   .. code-block:: c

      retval = SUNLinSolSetup(LS, A);


.. c:function:: int SUNLinSolSolve(SUNLinearSolver LS, SUNMatrix A, N_Vector x, N_Vector b, realtype tol)

   This *required* function Solves a linear system :math:`Ax = b`.

   **Arguments:**
      * *LS* -- a SUNLinSol object.
      * *A* -- a ``SUNMatrix`` object.
      * *x* -- a ``N_Vector`` object containing the initial guess for
        the solution of the linear system, and the solution to the
        linear system upon return.
      * *b* -- a ``N_Vector`` object containing the linear system
        right-hand side.
      * *tol* -- the desired linear solver tolerance.

   **Return value:**  This should return zero for a successful call, a
   positive value for a recoverable failure and a negative value for
   an unrecoverable failure, ideally returning one of the generic
   error codes listed in section :ref:`SUNLinSol.ErrorCodes`.

   **Direct solvers:** can ignore the *tol* argument.

   **Matrix-free solvers:** (those that identify as
   ``SUNLINEARSOLVER_ITERATIVE``) can ignore the ``SUNMatrix`` input
   *A*, and should rely on the matrix-vector product function supplied
   through the routine :c:func:`SUNLinSolSetATimes()`.

   **Iterative solvers:** (those that identify as
   ``SUNLINEARSOLVER_ITERATIVE`` or
   ``SUNLINEARSOLVER_MATRIX_ITERATIVE``) should attempt to solve to
   the specified tolerance *tol* in a weighted 2-norm. If the solver
   does not support scaling then it should just use a 2-norm.

   Usage:

   .. code-block:: c

      retval = SUNLinSolSolve(LS, A, x, b, tol);


.. c:function:: int SUNLinSolFree(SUNLinearSolver LS)

   Frees memory allocated by the linear solver.  This should return
   zero for a successful call, and a negative value for a failure.

   Usage:

   .. code-block:: c

      retval = SUNLinSolFree(LS);




.. _SUNLinSol.SetFn:

SUNLinearSolver set functions
-------------------------------------

The following set functions are used to supply linear solver modules with
functions defined by the SUNDIALS packages and to modify solver
parameters.  Only the routine for setting the matrix-vector product
routine is required, and that is only for matrix-free linear solver
modules.  Otherwise, all other set functions are optional.  SUNLinSol
implementations that do not provide the functionality for any optional
routine should leave the corresponding function pointer ``NULL``
instead of supplying a dummy routine.


.. c:function:: int SUNLinSolSetATimes(SUNLinearSolver LS, void* A_data, ATimesFn ATimes)

   This function is *required for matrix-free linear solvers*;
   otherwise it is optional.

   Provides a :c:type:`ATimesFn` function pointer, as well as a ``void*``
   pointer to a data structure used by this routine, to a linear
   solver object.  SUNDIALS packages will call this function to set the
   matrix-vector product function to either a solver-provided
   difference-quotient via vector operations or a user-supplied
   solver-specific routine.  This routine should return zero for a
   successful call, and a negative value for a failure, ideally
   returning one of the generic error codes listed in section
   :ref:`SUNLinSol.ErrorCodes`.

   Usage:

   .. code-block:: c

      retval = SUNLinSolSetATimes(LS, A_data, ATimes);


.. c:function:: int SUNLinSolSetPreconditioner(SUNLinearSolver LS, void* P_data, PSetupFn Pset, PSolveFn Psol)

   This *optional* routine provides :c:type:`PSetupFn` and
   :c:type:`PSolveFn` function pointers that implement the
   preconditioner solves :math:`P_1^{-1}` and :math:`P_2^{-1}`. This
   routine will be called by a SUNDIALS package, which will provide
   translation between the generic *Pset* and *Psol* calls and the
   package- or user-supplied routines.
   This routine should return zero for a successful call, and a
   negative value for a failure, ideally returning one of the generic
   error codes listed in section :ref:`SUNLinSol.ErrorCodes`.

   Usage:

   .. code-block:: c

      retval = SUNLinSolSetPreconditioner(LS, Pdata, Pset, Psol);


.. c:function:: int SUNLinSolSetScalingVectors(SUNLinearSolver LS, N_Vector s1, N_Vector s2)

   This *optional* routine provides left/right scaling vectors for the
   linear system solve.  Here, *s1* and *s2* are ``N_Vectors`` of positive
   scale factors containing the diagonal of the matrices :math:`S_1`
   and :math:`S_2`, respectively.  Neither of these vectors need
   to be tested for positivity, and a ``NULL`` argument for either
   indicates that the corresponding scaling matrix is the
   identity. This routine should return zero for a successful call,
   and a negative value for a failure, ideally returning one of the
   generic error codes listed in section :ref:`SUNLinSol.ErrorCodes`.

   Usage:

   .. code-block:: c

      retval = SUNLinSolSetScalingVectors(LS, s1, s2);






.. _SUNLinSol.GetFn:

SUNLinearSolver get functions
----------------------------------

The following get functions allow SUNDIALS packages to retrieve
results from a linear solve.  All routines are optional.


.. c:function:: int SUNLinSolNumIters(SUNLinearSolver LS)

   This *optional* routine should return the number of linear
   iterations performed in the last "solve" call.

   Usage:

   .. code-block:: c

      its = SUNLinSolNumIters(LS);


.. c:function:: realtype SUNLinSolResNorm(SUNLinearSolver LS)

   This *optional* routine should return the final residual norm from
   the last "solve" call.

   Usage:

   .. code-block:: c

      rnorm = SUNLinSolResNorm(LS);


.. c:function:: N_Vector SUNLinSolResid(SUNLinearSolver LS)

   If an iterative method computes the preconditioned initial residual
   and returns with a successful solve without performing any
   iterations (i.e., either the initial guess or the preconditioner is
   sufficiently accurate), then this *optional* routine may be called
   by the SUNDIALS package.  This routine should return the ``N_Vector``
   containing the preconditioned initial residual vector.

   Usage:

   .. code-block:: c

      rvec = SUNLinSolResid(LS);

   Note: since ``N_Vector`` is actually a pointer, and the results are
   not modified, this routine should *not* require additional memory
   allocation.  If the SUNLinSol object does not retain a vector for
   this purpose, then this function pointer should be set to ``NULL``
   in the implementation.

.. c:function:: sunindextype SUNLinSolLastFlag(SUNLinearSolver LS)

   This *optional* routine should return the last error flag
   encountered within the linear solver. This is not called by the
   SUNDIALS packages directly; it allows the user to investigate
   linear solver issues after a failed solve.

   Usage:

   .. code-block:: c

      lflag = SUNLinLastFlag(LS);


.. c:function:: int SUNLinSolSpace(SUNLinearSolver LS, long int *lenrwLS, long int *leniwLS)

   This *optional* routine should return the storage requirements for
   the linear solver *LS*.  *lrw* is a ``long int`` containing the
   number of realtype words and *liw* is a ``long int`` containing the
   number of integer words.  The return value is an integer flag
   denoting success/failure of the operation.

   This function is advisory only, for use by users to help determine
   their total space requirements.

   Usage:

   .. code-block:: c

      retval = SUNLinSolSpace(LS, &lrw, &liw);





.. _SUNLinSol.SUNSuppliedFn:

Functions provided by SUNDIALS packages
---------------------------------------------

To interface with SUNLinSol modules, the SUNDIALS packages supply a
variety of routines for evaluating the matrix-vector product, and
setting up and applying the preconditioniner.  These package-provided
routines translate between the user-supplied ODE, DAE, or nonlinear
systems and the generic interfaces to the linear systems of equations
that result in their solution. The types for functions provided to a
SUNLinSol module are defined in the header file
``sundials/sundials_iterative.h``, and are described below.


.. c:type:: typedef int (*ATimesFn)(void *A_data, N_Vector v, N_Vector z)

   These functions compute the action of a matrix on a vector,
   performing the operation :math:`z = Av`.  Memory for *z* will
   already be allocated prior to calling this function.  The parameter
   *A_data* is a pointer to any information about :math:`A` which the
   function needs in order to do its job. The vector :math:`v` should
   be left unchanged.  This routine should return 0 if successful and a
   non-zero value if unsuccessful.


.. c:type:: typedef int (*PSetupFn)(void *P_data)

   These functions set up any requisite problem data in preparation
   for calls to the corresponding :c:type:`PSolveFn`. This routine
   should return 0 if successful and a non-zero value if
   unsuccessful.


.. c:type:: typedef int (*PSolveFn)(void *P_data, N_Vector r, N_Vector z, realtype tol, int lr)

   These functions solve the preconditioner equation :math:`Pz = r`
   for the vector :math:`z`.  Memory for *z* will already be allocated
   prior to calling this function.  The parameter *P_data* is a
   pointer to any information about :math:`P` which the function needs
   in order to do its job (set up by the corresponding
   :c:type:`PSetupFn`). The parameter *lr* is input, and indicates
   whether :math:`P` is to be taken as the left or right
   preconditioner: *lr* = 1 for left and *lr* = 2 for right.  If
   preconditioning is on one side only, *lr* can be ignored.  If the
   preconditioner is iterative, then it should strive to solve the
   preconditioner equation so that

   .. math::

      \| Pz - r \|_{\text{wrms}} < tol

   where the error weight vector for the WRMS norm may be accessed
   from the main package memory structure.  The vector *r* should not
   be modified by the *PSolveFn*.  This routine should return 0 if
   successful and a non-zero value if unsuccessful.  On a failure, a
   negative return value indicates an unrecoverable condition, while a
   positive value indicates a recoverable one, in which the calling
   routine may reattempt the solution after updating preconditioner
   data.


.. _SUNLinSol.ErrorCodes:

SUNLinearSolver return codes
------------------------------------

The functions provided to SUNLinSol modules by each SUNDIALS package,
and functions within the SUNDIALS-provided SUNLinSol implementations
utilize a common set of return codes, listed below.  These adhere to a
common pattern: 0 indicates success, a postitive value corresponds to
a recoverable failure, and a negative value indicates a
non-recoverable failure.  Aside from this pattern, the actual values
of each error code are primarily to provide additional information to
the user in case of a linear solver failure.

* ``SUNLS_SUCCESS`` (0) -- successful call or converged solve

* ``SUNLS_MEM_NULL`` (-801) -- the memory argument to the function is ``NULL``

* ``SUNLS_ILL_INPUT`` (-802) -- an illegal input has been provided to the function

* ``SUNLS_MEM_FAIL`` (-803) -- failed memory access or allocation

* ``SUNLS_ATIMES_NULL`` (-804) -- the Atimes function is ``NULL``

* ``SUNLS_ATIMES_FAIL_UNREC`` (-805) -- an unrecoverable failure occurred in the ``ATimes`` routine

* ``SUNLS_PSET_FAIL_UNREC`` (-806) -- an unrecoverable failure occurred in the ``Pset`` routine

* ``SUNLS_PSOLVE_NULL`` (-807) -- the preconditioner solve function is ``NULL``

* ``SUNLS_PSOLVE_FAIL_UNREC`` (-808) -- an unrecoverable failure occurred in the ``Psolve`` routine

* ``SUNLS_PACKAGE_FAIL_UNREC`` (-809) -- an unrecoverable failure occurred in an external linear solver package

* ``SUNLS_GS_FAIL`` (-810) -- a failure occurred during Gram-Schmidt orthogonalization (SPGMR/SPFGMR)

* ``SUNLS_QRSOL_FAIL`` (-811) -- a singular $R$ matrix was encountered in a QR factorization (SPGMR/SPFGMR)

* ``SUNLS_VECTOROP_ERR`` (-812) -- a vector operation error occurred

* ``SUNLS_RES_REDUCED`` (801) -- an iterative solver reduced the residual, but did not converge to the desired tolerance

* ``SUNLS_CONV_FAIL`` (802) -- an iterative solver did not converge (80and the residual was not reduced)

* ``SUNLS_ATIMES_FAIL_REC`` (803) -- a recoverable failure occurred in the ``ATimes`` routine

* ``SUNLS_PSET_FAIL_REC`` (804) -- a recoverable failure occurred in the ``Pset`` routine

* ``SUNLS_PSOLVE_FAIL_REC`` (805) -- a recoverable failure occurred in the ``Psolve`` routine

* ``SUNLS_PACKAGE_FAIL_REC`` (806) -- a recoverable failure occurred in an external linear solver package

* ``SUNLS_QRFACT_FAIL`` (807) -- a singular matrix was encountered during a QR factorization (SPGMR/SPFGMR)

* ``SUNLS_LUFACT_FAIL`` (808) -- a singular matrix was encountered during a LU factorization









.. _SUNLininSol.Generic:

The generic SUNLinearSolver module
-----------------------------------------

SUNDIALS packages interact with specific SUNLinSol implementations
through the generic SUNLinSol module on which all other SUNLinSol
iplementations are built.  The ``SUNLinearSolver`` type is a pointer
to a structure containing an implementation-dependent *content* field,
and an *ops* field.  The type ``SUNLinearSolver`` is defined as

.. code-block:: c

   typedef struct _generic_SUNLinearSolver *SUNLinearSolver;

   struct _generic_SUNLinearSolver {
     void *content;
     struct _generic_SUNLinearSolver_Ops *ops;
   };

where the ``_generic_SUNLinearSolver_Ops`` structure is a list of
pointers to the various actual linear solver operations provided by a
specific implementation.  The ``_generic_SUNLinearSolver_Ops``
structure is defined as

.. code-block:: c

   struct _generic_SUNLinearSolver_Ops {
     SUNLinearSolver_Type (*gettype)(SUNLinearSolver);
     SUNLinearSolver_ID   (*getid)(SUNLinearSolver);
     int                  (*setatimes)(SUNLinearSolver, void*, ATimesFn);
     int                  (*setpreconditioner)(SUNLinearSolver, void*,
                                               PSetupFn, PSolveFn);
     int                  (*setscalingvectors)(SUNLinearSolver,
                                               N_Vector, N_Vector);
     int                  (*initialize)(SUNLinearSolver);
     int                  (*setup)(SUNLinearSolver, SUNMatrix);
     int                  (*solve)(SUNLinearSolver, SUNMatrix, N_Vector,
                                   N_Vector, realtype);
     int                  (*numiters)(SUNLinearSolver);
     realtype             (*resnorm)(SUNLinearSolver);
     sunindextype         (*lastflag)(SUNLinearSolver);
     int                  (*space)(SUNLinearSolver, long int*, long int*);
     N_Vector             (*resid)(SUNLinearSolver);
     int                  (*free)(SUNLinearSolver);
   };


The generic SUNLinSol module defines and implements the linear solver
operations defined in Sections :ref:`SUNLinSol.CoreFn` through
:ref:`SUNLinSol.GetFn`.  These routines are in fact only wrappers to
the linear solver operations defined by a particular SUNLinSol
implementation, which are accessed through the *ops* field of the
``SUNLinearSolver`` structure.  To illustrate this point we show below
the implementation of a typical linear solver operation from the
generic ``SUNLinearSolver`` module, namely ``SUNLinSolInitialize``,
which initializes a ``SUNLinearSolver`` object for use after it has
been created and configured, and returns a flag denoting a
successful or failed operation:

.. code-block:: c

   int SUNLinSolInitialize(SUNLinearSolver S)
   {
     return ((int) S->ops->initialize(S));
   }



.. _SUNLinSol.Compatibility:

Compatibility of SUNLinearSolver modules
---------------------------------------------

We note that not all ``SUNLinearSolver`` types are compatible with all
``SUNMatrix`` and ``N_Vector`` types provided with SUNDIALS.  In Table
:ref:`SUNLinSol.linsol-matrix` we show the matrix-based linear solvers
available as ``SUNLinearSolver`` modules, and the compatible matrix
implementations.  Recall that Table :ref:`ARKStep_CInterface.solver-vector`
shows the compatibility between all ``SUNLinearSolver`` modules and vector
implementations.


.. _SUNLinSol.linsol-matrix:

Compatible SUNLinearSolver and SUNMatrix implementations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

================ ===== ====== ====== =============
Linear Solver    Dense Banded Sparse User Supplied
================ ===== ====== ====== =============
Dense            X                   X
LapackDense      X                   X
Band                   X             X
LapackBand             X             X
KLU                           X      X
SuperLU_MT                    X      X
User supplied    X     X      X      X
================ ===== ====== ====== =============





.. _SUNLinSol.Custom:

Implementing a custom SUNLinearSolver module
--------------------------------------------------

A particular implementation of the ``SUNLinearSolver`` module must:

* Specify the *content* field of the SUNLinSol module.

* Define and implement the required linear solver operations.  See the
  section :ref:`SUNLinSol.ARKode` to determine which SUNLinSol
  operations are required for this SUNDIALS package.

  Note that the names of these routines should be unique to that
  implementation in order to permit using more than one
  SUNLinSol module (each with different ``SUNLinearSolver``
  internal data representations) in the same code.

* Define and implement user-callable constructor and destructor
  routines to create and free a ``SUNLinearSolver`` with
  the new *content* field and with *ops* pointing to the
  new linear solver operations.

We note that the function pointers for all unsupported optional
routines should be set to ``NULL`` in the *ops* structure.  This
allows the SUNDIALS package that is using the SUNLinSol object
to know that the associated functionality is not supported.

To aid in the creation of custom ``SUNLinearSolver`` modules the generic
``SUNLinearSolver`` module provides the utility function
:c:func:`SUNLinSolNewEmpty`. When used in custom ``SUNLinearSolver``
constructors this function will ease the introduction of any new optional linear
solver operations to the ``SUNLinearSolver`` API by ensuring only required
operations need to be set.

.. c:function:: SUNLinearSolver SUNLinSolNewEmpty()

  This function allocates a new generic ``SUNLinearSolver`` object and
  initializes its content pointer and the function pointers in the operations
  structure to ``NULL``.

  **Return value:** If successful, this function returns a ``SUNLinearSolver``
  object. If an error occurs when allocating the object, then this routine will
  return ``NULL``.

.. c:function:: void SUNLinSolFreeEmpty(SUNLinearSolver LS)

  This routine frees the generic ``SUNLinearSolver`` object, under the assumption that any
  implementation-specific data that was allocated within the underlying content structure
  has already been freed. It will additionally test whether the ops pointer is ``NULL``,
  and, if it is not, it will free it as well.

   **Arguments:**
      * *LS* -- a SUNLinearSolver object


Additionally, a ``SUNLinearSolver`` implementation *may* do the following:

* Define and implement additional user-callable "set" routines
  acting on the ``SUNLinearSolver``, e.g., for setting various
  configuration options to tune the linear solver to a particular
  problem.

* Provide additional user-callable "get" routines acting on the
  ``SUNLinearSolver`` object, e.g., for returning various solve
  statistics.


.. _SUNLinSol.Intended:


Intended use cases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The SUNLinSol (and SUNMATRIX) APIs are designed to require a minimal set
of routines to ease interfacing with custom or third-party linear solver
libraries. External solvers provide similar routines with
the necessary functionality and thus will require minimal effort to wrap within
custom SUNMATRIX and SUNLinSol implementations. Sections
:ref:`SUNMatrix.ARKode` and :ref:`SUNLinSol.ARKode` include a list of
the required set of routines that compatible SUNMATRIX and SUNLinSol
implementations must provide. As SUNDIALS packages utilize generic
SUNLinSol modules allowing for user-supplied ``SUNLinearSolver``
implementations, there exists a wide range of possible linear solver
combinations. Some intended use cases for both the SUNDIALS-provided
and user-supplied SUNLinSol modules are discussd in the following
sections.


Direct linear solvers
""""""""""""""""""""""""""""""""

Direct linear solver modules require a matrix and compute an 'exact' solution to
the linear system *defined by the matrix*. Multiple matrix formats and
associated direct linear solvers are supplied with SUNDIALS through different
SUNMATRIX and SUNLinSol implementations. SUNDIALS packages strive to
amortize the high cost of matrix construction by reusing matrix information for
multiple nonlinear iterations. As a result, each package's linear solver
interface recomputes Jacobian information as infrequently as possible.

Alternative matrix storage formats and compatible linear solvers that are not
currently provided by or interfaced with SUNDIALS can leverage this
infrastructure with minimal effort. To do so, a user must implement custom
SUNMATRIX and SUNLinSol wrappers for the desired matrix format and/or linear
solver following the APIs described in the sections :ref:`SUNMatrix`
and :ref:`SUNLinSol`.  *This user-supplied SUNLinSol module must then
self-identify as having* ``SUNLINEARSOLVER_DIRECT`` *type*.


Matrix-free iterative linear solvers
""""""""""""""""""""""""""""""""""""""

Matrix-free iterative linear solver modules do not require a matrix and compute
an inexact solution to the linear system *defined by the package-supplied*
``ATimes`` *routine*. SUNDIALS supplies multiple scaled, preconditioned
iterative linear solver (spils) SUNLinSol modules that support scaling to
allow users to handle non-dimensionalization (as best as possible) within each
SUNDIALS package and retain variables and define equations as desired in
their applications. For linear solvers that do not support left/right scaling,
the tolerance supplied to the linear solver is adjusted to compensate (see
section :ref:`SUNLinSol.Iterative_Tolerance` for more details); however, this
use case may be non-optimal and cannot handle situations where the magnitudes of
different solution components or equations vary dramatically within a single
problem.

To utilize alternative linear solvers that are not currently provided by or
interfaced with SUNDIALS a user must implement a custom SUNLinSol wrapper
for the linear solver following the API described in the section
:ref:`SUNLinSol`.  *This user-supplied SUNLinSol module must then
self-identify as having* ``SUNLINEARSOLVER_ITERATIVE`` *type*.


Matrix-based iterative linear solvers (reusing :math:`A`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Matrix-based iterative linear solver modules require a matrix and compute an
inexact solution to the linear system *defined by the matrix*.  This
matrix will be updated infrequently and resued across multiple solves
to amortize cost of matrix construction. As in the direct linear
solver case, only wrappers for the matrix and linear solver in
SUNMATRIX and SUNLinSol implementations need to be created to utilize
a new linear solver. *This user-supplied SUNLinSol module must then
self-identify as having* ``SUNLINEARSOLVER_MATRIX_ITERATIVE`` *type*.

At present, SUNDIALS has one example problem that uses this approach for
wrapping a structured-grid matrix, linear solver, and preconditioner from the
*hypre* library that may be used as a template for other customized
implementations
(see ``examples/arkode/CXX_parhyp/ark_heat2D_hypre.cpp``).


Matrix-based iterative linear solvers (current :math:`A`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For users who wish to utilize a matrix-based iterative linear solver module
where the matrix is *purely for preconditioning* and the linear system is
*defined by the package-supplied* ``ATimes`` *routine*, we envision two
current possibilities.

The preferred approach is for users to employ one of the SUNDIALS
scaled, preconditioned iterative linear solver (spils) implementations
(:c:func:`SUNLinSol_SPGMR()`, :c:func:`SUNLinSol_SPFGMR()`,
:c:func:`SUNLinSol_SPBCGS()`, :c:func:`SUNLinSol_SPTFQMR()`, or
:c:func:`SUNLinSol_PCG()`) as the outer solver. The creation and storage of the
preconditioner matrix, and interfacing with the corresponding linear solver, can
be handled through a package's preconditioner 'setup' and 'solve' functionality
(see the sections :ref:`ARKStep_CInterface.PrecSetupFn` and
:ref:`ARKStep_CInterface.PrecSolveFn`, respectively) without creating
SUNMATRIX and SUNLinSol implementations. This usage mode is
recommended primarily because the SUNDIALS-provided spils modules
support the scaling as described above.

A second approach supported by the linear solver APIs is as follows. If the
SUNLinSol implementation is matrix-based, *self-identifies
as having* ``SUNLINEARSOLVER_ITERATIVE`` *type*, and *also provides a non-NULL
:c:func:`SUNLinSolSetATimes()` routine*, then each SUNDIALS package
will call that routine to attach its package-specific matrix-vector
product routine to the SUNLinSol object. The SUNDIALS package will
then call the SUNLinSol-provided :c:func:`SUNLinSolSetup()` routine
(infrequently) to update matrix information, but will provide current
matrix-vector products to the SUNLinSol implementation through the
package-supplied ``ATimesFn`` routine.
