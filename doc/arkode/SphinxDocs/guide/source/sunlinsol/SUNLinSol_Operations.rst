..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _SUNLinSol.Ops:

Description of the SUNLinearSolver operations
===============================================

For each of the ``SUNLinearSolver`` operations, we give the name,
usage of the function, and a description of its mathematical
operations below.


.. c:function:: SUNLinearSolver_Type SUNLinSolGetType(SUNLinearSolver LS)

   Returns the type identifier for the linear solver *LS*. It is used
   to determine the solver type (direct, iterative, or custom) from
   the abstract ``SUNLinearSolver`` interface.  This is used to
   assess compatibility with SUNDIALS-provided linear solver
   interfaces.  Returned values are given in the Table
   :ref:`SUNLinSol.linsolIDs`.
   
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

   **Iterative solvers:** can ignore the ``SUNMatrix`` input *A* since
   a ``NULL`` argument will be passed (these should instead rely on
   the matrix-vector product function supplied through the routine
   :c:func:`SUNLinSolSetATimes()`).  These should attempt to solve to
   the specified ``realtype`` tolerance *tol* in a weighted 2-norm. 
   If the solver does not support scaling then it should just use a
   2-norm.
      
   **Custom solvers:** all arguments will be supplied, and if the
   solver is approximate then it should attempt to solve to the
   specified ``realtype`` tolerance *tol* in a weighted 2-norm. If the
   solver does not support scaling then it should just use a 2-norm.
  
   Usage:

   .. code-block:: c

      ier = SUNLinSolSolve(LS, A, x, b, tol);


.. c:function:: int SUNLinSolFree(SUNLinearSolver LS)

   Frees memory allocated by the linear solver.  This should return
   zero for a successful call, and a negative value for a failure.
 
   Usage:

   .. code-block:: c
 
      ier = SUNLinSolFree(LS);


.. c:function:: int SUNLinSolSetATimes(SUNLinearSolver LS, void* A_data, ATimesFn ATimes)

   Provides :c:type:`ATimesFn` function pointer, as well as a ``void *``
   pointer to a data structure used by this routine, to a linear
   solver object.  SUNDIALS solvers will call this function to set the
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

   (Optional; Iterative/Custom linear solvers only)

   Provides :c:type:`PSetupFn` and :c:type:`PSolveFn` function
   pointers that implement the preconditioner solves :math:`P_1^{-1}`
   and :math:`P_2^{-1}`. This routine will be called by a SUNDIALS
   solver, which will provide translation between the generic
   *Pset* and *Psol* calls and the integrator-specific and
   integrator- or user-supplied routines.  This routine should return
   zero for a successful call, and a negative value for a failure,
   ideally returning one of the generic error codes listed in section
   :ref:`SUNLinSol.ErrorCodes`.  

   Usage:

   .. code-block:: c

      ier = SUNLinSolSetPreconditioner(LS, Pdata, Pset, Psol);


.. c:function:: int SUNLinSolSetScalingVectors(SUNLinearSolver LS, N_Vector s1, N_Vector s2)

   (Optional; Iterative/Custom linear solvers only)

   Sets pointers to left/right scaling vectors for the linear system
   solve.  Here, *s1* is an ``N_Vector`` of positive scale factors 
   containing the diagonal of the :math:`S_1` scaling matrix.
   Similarly, *s2* is an ``N_Vector`` containing the diagonal of
   the :math:`S_2` scaling matrix.  Neither of these vectors are
   tested for positivity, and a ``NULL`` argument for either indicates
   that the corresponding scaling matrix is the identity. This routine
   should return zero for a successful call, and a negative value for
   a failure, ideally returning one of the generic error codes listed
   in section :ref:`SUNLinSol.ErrorCodes`. 

   Usage:

   .. code-block:: c

      ier = SUNLinSolSetScalingVectors(LS, s1, s2);


.. c:function:: int SUNLinSolNumIters(SUNLinearSolver LS)

   (Optional; Iterative/Custom linear solvers only)

   Should return the number of linear iterations performed in the last
   "solve" call.  

   Usage:

   .. code-block:: c

      its = SUNLinSolNumIters(LS);


.. c:function:: realtype SUNLinSolResNorm(SUNLinearSolver LS)

   (Optional; Iterative/Custom linear solvers only)

   Should return the final residual norm from the last "solve" call.
 
   Usage:
 
   .. code-block:: c

      rnorm = SUNLinSolResNorm(LS);


.. c:function:: N_Vector SUNLinSolResid(SUNLinearSolver LS)

   (Optional; Iterative/Custom linear solvers only)

   If an iterative method computes the preconditioned initial residual
   and returns with a successful solve without performing any
   iterations (i.e. either the initial guess or the preconditioner is
   sufficiently accurate), then this function may be called by the
   SUNDIALS solver.  This routine should return the ``N_Vector``
   containing the preconditioned initial residual vector. 

   Usage: 

   .. code-block:: c

      rvec = SUNLinSolResid(LS);


.. c:function:: long int SUNLinSolLastFlag(SUNLinearSolver LS)

   (Optional)

   Should return the last error flag encountered within the linear
   solver. This is not called by the SUNDIALS solvers directly; it
   allows the user to investigate linear solver issues after a failed
   solve. 

   Usage:

   .. code-block:: c

      lflag = SUNLinLastFlag(LS);


.. c:function:: int SUNLinSolSpace(SUNLinearSolver LS, long int *lenrwLS, long int *leniwLS)
 
   (Optional)

   Returns the storage requirements for the linear solver *LS*.
   *lrw* is a ``long int`` containing the number of realtype words
   and *liw* is a ``long int`` containing the number of integer
   words.  The return value is an integer flag denoting
   success/failure of the operation.  

   This function is advisory only, for use in determining a user's
   total space requirements.

   Usage:

   .. code-block:: c

      ier = SUNLinSolSpace(LS, &lrw, &liw);



.. _SUNLinSol.Supplied:
      
Description of the client-supplied SUNLinearSolver routines
==============================================================

The SUNDIALS packages provide the *ATimes*, *Pset* and *Psol* routines
utilized by the ``SUNLinearSolver`` modules.  These function types are
defined in the header file ``sundials/sundials_iterative.h``, and
are described here in case a user wishes to interact directly with an
iterative ``SUNLinearSolver`` object.

.. c:type:: typedef int (*ATimesFn)(void *A_data, N_Vector v, N_Vector z)

   These functions compute the action of a matrix on a vector,
   performing the operation :math:`z = Av`.  Memory for *z* should already be
   allocted prior to calling this function.  The parameter *A_data* is
   a pointer to any information about :math:`A` which the function
   needs in order to do its job. The vector :math:`v` should be left
   unchanged.  This routine should return 0 if successful and a         
   non-zero value if unsuccessful.


.. c:type:: typedef int (*PSetupFn)(void *P_data)

   These functions set up any requisite problem data in preparation
   for calls to the corresponding :c:type:`PSolveFn`. This routine
   should return 0 if successful and a non-zero value if
   unsuccessful. 
   

.. c:type:: typedef int (*PSolveFn)(void *P_data, N_Vector r, N_Vector z, realtype tol, int lr)

   These functions solve the preconditioner equation :math:`Pz = r`
   for the vector :math:`z`.  Memory for *z* should already be
   allocted prior to calling this function.  The parameter *P_data* is
   a pointer to any information about :math:`P` which the function
   needs in order to do its job (set up by the corresponding
   :c:type:`PSetupFn`. The parameter *lr* is input, and indicates
   whether :math:`P` is to be taken as the left preconditioner or the
   right preconditioner: *lr* = 1 for left and *lr* = 2 for right.  If
   preconditioning is on one side only, *lr* can be ignored.  If the
   preconditioner is iterative, then it should strive to solve the
   preconditioner equation so that 

   .. math::
         
      \| Pz - r \|_{\text{wrms}} < tol

   where the weight vector for the WRMS norm may be accessed from the
   main package memory structure.  The vector *r* should not be
   modified by the *PSolveFn*.  This routine should return 0 if
   successful and a non-zero value if unsuccessful.  On a failure, a
   negative return value indicates an unrecoverable condition, while a
   positive value indicates a recoverable one, in which the calling
   routine may reattempt the solution after updating preconditioner
   data.
