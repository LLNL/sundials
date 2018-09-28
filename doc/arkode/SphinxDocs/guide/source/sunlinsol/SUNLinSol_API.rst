..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _SUNLinSol.API:

The SUNLinearSolver API
=============================

For problems that require the solution of linear systems of equations,
the SUNDIALS packages operate using generic linear solver modules
defined through the SUNLinSol API.  This allows SUNDIALS
packages to utilize any valid SUNLinSol implementation that provides
a set of required functions.  These functions can be divided into
three categories.  The first are the core linear solver functions.  The
second group consists of "set" routines to supply the linear solver object
with functions provided by the SUNDIALS package, or for modification
of solver parameters.  The last group consists of "get" routines for
retrieving artifacts (statistics, residual vectors, etc.) from the
linear solver.  All of these functions are defined in the header file
``sundials/sundials_linearsolver.h``.

The implementations provided with SUNDIALS work in coordination
with the SUNDIALS generic ``N_Vector`` and ``SUNMatrix`` modules to
provide a set of compatible data structures and solvers for the
solution of linear systems using direct or matrix-free iterative
methods. Moreover, advanced users can provide a customized
``SUNLinearSolver`` implementation to any SUNDIALS package,
particularly in cases where they provide their own ``N_Vector`` and/or
``SUNMatrix`` modules.  

Historically, the SUNDIALS packages have been designed to specifically
leverage the use of either *direct linear solvers* or matrix-free,
*scaled, preconditioned, iterative linear solvers*.  However,
matrix-based iterative linear solvers are supported.

The iterative linear solvers packaged with SUNDIALS leverage scaling
and preconditioning, as applicable, to balance error between solution
components and to accelerate convergence of the linear solver.  To
this end, instead of solving the linear system :math:`Ax = b`
directly, these apply the underlying iterative algorithm to the
transformed system   

.. math::
   :label: eq:transformed_linear_system
   
   \tilde{A} \tilde{x} = \tilde{b}

where

.. math::
  :label: eq:transformed_linear_system_components
   
  \tilde{A} &= S_1 P_1^{-1} A P_2^{-1} S_2^{-1},\\
  \tilde{b} &= S_1 P_1^{-1} b,\\
  \tilde{x} &= S_2 P_2 x,

and where

* :math:`P_1` is the left preconditioner,
  
* :math:`P_2` is the right preconditioner,
    
* :math:`S_1` is a diagonal matrix of scale factors for
  :math:`P_1^{-1} b`,
      
* :math:`S_2` is a diagonal matrix of scale factors for :math:`P_2 x`.

SUNDIALS solvers request that iterative linear solvers stop
based on the 2-norm of the scaled preconditioned residual meeting a
prescribed tolerance

.. math::
   
   \left\| \tilde{b} - \tilde{A} \tilde{x} \right\|_2  <  \text{tol}.


When provided an iterative SUNLinSol implementation that does not
support the scaling matrices :math:`S_1` and :math:`S_2`, SUNDIALS'
packages will adjust the value of :math:`\text{tol}` accordingly.  In
this case, they instead request that iterative linear solvers stop
based on the criteria

.. math::
   
   \left\| P_1^{-1} b - P_1^{-1} A x \right\|_2  <  \text{tol}.

We note that the corresponding adjustments to :math:`\text{tol}` in
this case are non-optimal, in that they cannot balance error between
specific entries of the solution :math:`x`, only the aggregate error
in the overall solution vector.
   
We further note that not all of the SUNDIALS-provided iterative linear
solvers support the full range of the above options (e.g., separate
left/right preconditioning), and that some of the SUNDIALS packages
only utilize a subset of these options.  Further details on these
exceptions are described in the documentation for each
``SUNLinearSolver`` implementation, or for each SUNDIALS package.




.. _SUNLinSol.CoreFn:

SUNLinearSolver core functions
-----------------------------------------------------

The core linear solver functions consist of four required routines to get
the linear solver type (:c:func:`SUNLinSolGetType()`), initialize
the linear solver object once all solver-specific options have been
set (:c:func:`SUNLinSolInitialize()`), set up the linear solver object
to utilize an updated matrix :math:`A` (:c:func:`SUNLinSolSetup()`),
and solve the linear system :math:`Ax=b` (:c:func:`SUNLinSolSolve()`).
The remaining routine for destruction of the linear solver object
(:c:func:`SUNLinSolFree()`) is optional.


.. c:function:: SUNLinearSolver_Type SUNLinSolGetType(SUNLinearSolver LS)

   Returns the type identifier for the linear solver *LS*. It is used
   to determine the solver type (direct or iterative) from
   the abstract ``SUNLinearSolver`` interface.  This is used to
   assess compatibility with SUNDIALS-provided linear solver
   interfaces.  Returned values are one of the following:

   * ``SUNLINEARSOLVER_DIRECT`` -- ``0``, the SUNLinSol module uses
     direct methods to solve the linear system.

   * ``SUNLINEARSOLVER_ITERATIVE`` -- ``1``, the SUNLinSol module
     iteratively solves the linear system, stopping when the linear
     residual is within a prescribed tolerance.

   Usage:

   .. code-block:: c

      type = SUNLinSolGetType(LS);


.. c:function:: int SUNLinSolInitialize(SUNLinearSolver LS)

   Performs linear solver initialization (assumes that all
   solver-specific options have been set).  This should return zero for a
   successful call, and a negative value for a failure, ideally
   returning one of the generic error codes listed in section
   :ref:`SUNLinSol.ErrorCodes`. 
   
   Usage:

   .. code-block:: c

      ier = SUNLinSolInitialize(LS);


.. c:function:: int SUNLinSolSetup(SUNLinearSolver LS, SUNMatrix A)

   Performs any linear solver setup needed, based on an updated system
   ``SUNMatrix`` *A*.  This may be called frequently (e.g. with a full
   Newton method) or infrequently (for a modified Newton method), based
   on the type of integrator and/or nonlinear solver requesting the
   solves.  This should return zero for a successful call, a positive
   value for a recoverable failure and a negative value for an
   unrecoverable failure, ideally returning one of the generic error
   codes listed in section :ref:`SUNLinSol.ErrorCodes`. 

   Usage:

   .. code-block:: c

      ier = SUNLinSolSetup(LS, A);


.. c:function:: int SUNLinSolSolve(SUNLinearSolver LS, SUNMatrix A, N_Vector x, N_Vector b, realtype tol)
 
   Solves a linear system :math:`Ax = b`.  This should return zero for
   a successful call, a positive value for a recoverable failure and a
   negative value for an unrecoverable failure, ideally returning one
   of the generic error codes listed in section :ref:`SUNLinSol.ErrorCodes`. 

   **Direct solvers:** can ignore the  *tol* argument.

   **Matrix-free solvers:** can ignore the ``SUNMatrix`` input *A*
   since a ``NULL`` argument will be passed (these should instead rely on
   the matrix-vector product function supplied through the routine
   :c:func:`SUNLinSolSetATimes()`).

   **Iterative solvers:** These should attempt to solve to
   the specified ``realtype`` tolerance *tol* in a weighted 2-norm. 
   If the solver does not support scaling then it should just use a
   2-norm.
  
   Usage:

   .. code-block:: c

      ier = SUNLinSolSolve(LS, A, x, b, tol);


.. c:function:: int SUNLinSolFree(SUNLinearSolver LS)

   Frees memory allocated by the linear solver.  This should return
   zero for a successful call, and a negative value for a failure.
 
   Usage:

   .. code-block:: c
 
      ier = SUNLinSolFree(LS);




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

   (Required for matrix-free linear solvers; otherwise optional)

   Provides a :c:type:`ATimesFn` function pointer, as well as a ``void *``
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

      ier = SUNLinSolSetATimes(LS, A_data, ATimes);


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

      ier = SUNLinSolSetPreconditioner(LS, Pdata, Pset, Psol);


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

      ier = SUNLinSolSetScalingVectors(LS, s1, s2);






.. _SUNLinSol.GetFn:

SUNLinearSolver get functions
----------------------------------

The following get functions allow SUNDIALS packages to retrieve
results from the linear solve.  All routines are optional.
      

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
   iterations (i.e. either the initial guess or the preconditioner is
   sufficiently accurate), then this *optional* routine may be called
   by the SUNDIALS package.  This routine should return the ``N_Vector``
   containing the preconditioned initial residual vector (note: since
   ``N_Vector`` is actually a pointer, and the results are not
   modified, this should *not* require additional memory allocation).

   Usage: 

   .. code-block:: c

      rvec = SUNLinSolResid(LS);


.. c:function:: long int SUNLinSolLastFlag(SUNLinearSolver LS)

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

      ier = SUNLinSolSpace(LS, &lrw, &liw);





.. _SUNLinSol.SUNSuppliedFn:
      
Functions provided by SUNDIALS packages
---------------------------------------------

To interface with SUNLinSol modules, the SUNDIALS packages supply a
variety of routines for evaluating the matrix-vector product, and
setting up and applying the preconditioniner.  These package-provided
routines translate between the user-supplied ODE, DAE or nonlinear
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

* ``SUNLS_MEM_NULL`` (-1) -- the memory argument to the function is ``NULL``

* ``SUNLS_ILL_INPUT`` (-2) -- an illegal input has been provided to the function 

* ``SUNLS_MEM_FAIL`` (-3) -- failed memory access or allocation

* ``SUNLS_ATIMES_FAIL_UNREC`` (-4) -- an unrecoverable failure occurred in the ``ATimes`` routine

* ``SUNLS_PSET_FAIL_UNREC`` (-5) -- an unrecoverable failure occurred in the ``Pset`` routine

* ``SUNLS_PSOLVE_FAIL_UNREC`` (-6) -- an unrecoverable failure occurred in the ``Psolve`` routine

* ``SUNLS_PACKAGE_FAIL_UNREC`` (-7) -- an unrecoverable failure occurred in an external linear solver package

* ``SUNLS_GS_FAIL`` (-8) -- a failure occurred during Gram-Schmidt orthogonalization (SPGMR/SPFGMR)

* ``SUNLS_QRSOL_FAIL`` (-9) -- a singular $R$ matrix was encountered in a QR factorization (SPGMR/SPFGMR)

* ``SUNLS_RES_REDUCED`` (1) -- an iterative solver reduced the residual, but did not converge to the desired tolerance

* ``SUNLS_CONV_FAIL`` (2) -- an iterative solver did not converge (and the residual was not reduced)

* ``SUNLS_ATIMES_FAIL_REC`` (3) -- a recoverable failure occurred in the ``ATimes`` routine

* ``SUNLS_PSET_FAIL_REC`` (4) -- a recoverable failure occurred in the ``Pset`` routine

* ``SUNLS_PSOLVE_FAIL_REC`` (5) -- a recoverable failure occurred in the ``Psolve`` routine

* ``SUNLS_PACKAGE_FAIL_REC`` (6) -- a recoverable failure occurred in an external linear solver package

* ``SUNLS_QRFACT_FAIL`` (7) -- a singular matrix was encountered during a QR factorization (SPGMR/SPFGMR)

* ``SUNLS_LUFACT_FAIL`` (8) -- a singular matrix was encountered during a LU factorization


   



   


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
     long int             (*lastflag)(SUNLinearSolver);
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
  
* Define and implement the required linear solver operations defined
  in Sections :ref:`SUNLinSol.CoreFn` through :ref:`SUNLinSol.SetFn`.
  Note that the names of these routines should be unique to that
  implementation in order to permit using more than one
  SUNLinSol module (each with different ``SUNLinearSolver``
  internal data representations) in the same code.

* Define and implement a user-callable constructor to create the
  ``SUNLinearSolver`` object.

We note that the function pointers for all unsupported optional
routines should be set to ``NULL`` in the *ops* structure, to indicate
to the SUNDIALS package that they are not provided.

  
Additionally, a ``SUNLinearSolver`` implementation *may* do the following:
  
* Define and implement additional user-callable "set" routines
  acting on the ``SUNLinearSolver``, e.g., for setting various
  configuration options to tune the linear solver to a particular
  problem. 

* Provide additional user-callable "get" routines acting on the
  ``SUNLinearSolver`` object, e.g., for returning various solve
  statistics.


