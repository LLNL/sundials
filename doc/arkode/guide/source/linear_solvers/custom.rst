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


.. _LinearSolvers.custom:

Providing Alternate Linear Solver Modules
==================================================


Newton system linear solver
------------------------------

The central ARKode module interfaces with the Newton system linear
solver module using calls to one of four routines. These are denoted
here by :c:func:`linit()`, :c:func:`lsetup()`, :c:func:`lsolve()`, and
:c:func:`lfree()`. Briefly, their purposes are as follows:

* :c:func:`linit()`: initializes memory specific to the linear solver;
* :c:func:`lsetup()`: evaluates and preprocesses the Jacobian or
  preconditioner in preparation for solves;
* :c:func:`lsolve()`: solves the linear system;
* :c:func:`lfree()`: frees the linear solver memory.

A linear solver module must also provide a user-callable **specification
function** (like those described in the section
:ref:`ARKStep_CInterface.LinearSolvers`) which will attach the above four
routines to the main ARKode memory block. The ARKode memory block is a
structure defined in the header file ``arkode_impl.h``. A pointer to
such a structure is defined as the type ``ARKodeMem``. The four
fields in the ``ARKodeMem`` structure that refer to the Newton system
linear solver's functions are ``ark_linit``, ``ark_lsetup``,
``ark_lsolve``, and ``ark_lfree``, respectively.  Note that of these
interface functions , only the :c:func:`lsolve()` function is
required. The :c:func:`lfree()` routine must be provided only if the
solver specification routine makes any memory allocation.  For any of
the functions that are *not* provided, the corresponding field should
be set to ``NULL``. The linear
solver specification function must also set the value of the field
``ark_setupNonNull`` in the ARKode memory block -- to ``SUNTRUE`` if
:c:func:`lsetup()` is used, or ``SUNFALSE`` otherwise.

Typically, the linear solver will require a block of memory specific
to the solver, and a principal function of the specification function
is to allocate that memory block, and initialize it.  Then the field
``ark_lmem`` in the ARKode memory block ``ARKodeMem`` is available to
attach a pointer to that linear solver memory.  This block can then be
used to facilitate the exchange of data between the four interface
functions.

If the linear solver involves adjustable parameters, the specification
function should set the default values of those.  User-callable
functions may be defined that could, optionally, override the default
parameter values.

We encourage the use of performance counters in connection with the various
operations involved with the linear solver.  Such counters would be
members of the linear solver memory block, would be initialized in the
:c:func:`linit()` function, and would be incremented by the
:c:func:`lsetup()` and :c:func:`lsolve()` functions.  Then
user-callable functions would be needed to obtain the values of these
counters.

For consistency with the existing ARKode linear solver modules, we
recommend that the return value of the specification function be 0 for
a successful return, and a negative value if an error occurs.
Possible error conditions include: the pointer to the main ARKode
memory block is ``NULL``, an input is illegal, the NVECTOR
implementation is not compatible, or a memory allocation fails.





Mass matrix linear solver
------------------------------

Similarly, for problems involving a non-identity mass matrix
:math:`M\ne I`, the main ARKode module interfaces with the mass matrix
linear solver module using calls to one of four routines:
:c:func:`minit()`, :c:func:`msetup()`, :c:func:`msolve()`, and
:c:func:`mfree()`. Briefly, their purposes are as follows:

* :c:func:`minit()`: initializes memory specific to the mass matrix
  linear solver;
* :c:func:`msetup()`: evaluates and preprocesses the mass matrix or
  associated preconditioner in preparation for solves;
* :c:func:`msolve()`: solves the mass matrix system;
* :c:func:`mfree()`: frees the mass matrix linear solver memory.

As with the Newton system linear solver, a mass matrix linear solver
module must also provide a user-callable **specification function** (like
those described in the section :ref:`ARKStep_CInterface.LinearSolvers`) which
will attach the above four functions to the main ARKode memory
block.  The four fields in the ``ARKodeMem`` structure that refer to
the mass matrix system linear solver's functions are ``ark_minit``,
``ark_msetup``, ``ark_msolve``, and ``ark_mfree``, respectively.  As
with the Newton system solver, only :c:func:`msolve()` is required,
and :c:func:`mfree()` must be provided only if the solver
specification function makes any memory allocation.  For any of the
functions that are *not* provided, the corresponding field should be
set to ``NULL``.  The mass matrix linear solver specification function
must also set the value of the field ``ark_MassSetupNonNull`` in the
ARKode memory block -- to ``SUNTRUE`` if :c:func:`msetup()` is used, or
``SUNFALSE`` otherwise.

As with the Newton system linear solver, the mass matrix linear solver
will require a block of memory specific to the solver, so a principal
function of the specification function is to allocate that memory
block, and initialize it.  Then the field ``ark_mass_mem`` in the
ARKode memory block ``ARKodeMem`` is available to attach a pointer to
that mass matrix solver memory.  This block can then be used to
facilitate the exchange of data between the various interface functions.

If the linear solver involves adjustable parameters, the specification
function should set the default values of those.  User-callable
functions may be defined that could, optionally, override the default
parameter values.

We encourage the use of performance counters in connection with the various
operations involved with the linear solver.  Such counters would be
members of the linear solver memory block, would be initialized in the
:c:func:`minit()` function, and would be incremented by the
:c:func:`msetup()` and :c:func:`msolve()` functions.  Then
user-callable functions would be needed to obtain the values of these
counters.

For consistency with the existing ARKode linear solver modules, we
recommend that the return value of the specification function be 0 for
a successful return, and a negative value if an error occurs.
Possible error conditions include: the pointer to the main ARKode
memory block is ``NULL``, an input is illegal, the NVECTOR
implementation is not compatible, or a memory allocation fails.



These above functions, which interface between ARKode and the Newton
system or mass matrix linear solver module necessarily have fixed call
sequences.  Thus, a user wishing to implement another linear solver
within the ARKode package must adhere to this set of interfaces.  The
following is a complete description of the call list for each of these
functions.  Note that the call list of each function includes a pointer
to the main ARKode memory block, by which the function can access
various data related to the ARKode solution. The contents of this
memory block are given in the file ``arkode_impl.h`` (but not
reproduced here, for the sake of space).





Initialization function
-----------------------------------

The type definition of :c:func:`linit()` is

.. c:function:: typedef int (*linit)(ARKodeMem ark_mem)

   Completes initializations for the specific linear solver, such as
   counters and statistics.  It should also set pointers to data
   blocks that will later be passed to functions associated with the
   linear solver.  The :c:func:`linit()` function is called once only,
   at the start of the problem, during the first call to ARKode.

   **Arguments:**
      * *ark_mem* -- pointer to the ARKode memory block.

   **Return value:**  Should return 0 if it has successfully
   initialized the ARKode linear solver and a negative value
   otherwise.


Similarly, the type definition of :c:func:`minit()` is

.. c:function:: typedef int (*minit)(ARKodeMem ark_mem)

   Completes initializations for the specific mass matrix linear
   solver, such as counters and statistics.  It should also set
   pointers to data blocks that will later be passed to functions
   associated with the linear solver.  The :c:func:`minit()` function
   is called once only, at the start of the problem, during the first
   call to ARKode.

   **Arguments:**
      * *ark_mem* -- pointer to the ARKode memory block.

   **Return value:**  Should return 0 if it has successfully
   initialized the ARKode linear solver and a negative value
   otherwise.



Setup function
-----------------------------------


The type definition of :c:func:`lsetup()` is

.. c:function:: typedef int (*lsetup)(ARKodeMem ark_mem, int convfail, N_Vector ypred, N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)

   Prepares the linear solver for subsequent calls to
   :c:func:`lsolve()`, in the solution of systems :math:`A x = b`,
   where :math:`A` is some approximation to the Newton matrix,
   :math:`M-\gamma \frac{\partial f}{\partial y}`.  Here,
   :math:`\gamma` is available as ``ark_mem->ark_gamma``.

   The :c:func:`lsetup()` function may call a user-supplied function,
   or a function within the linear solver module, to compute needed
   data related to the Jacobian matrix :math:`\frac{\partial
   f}{\partial y}`.  Alterntively, it may choose to retrieve and use
   stored values of this data.

   In either case, :c:func:`lsetup()` may also preprocess that data as
   needed for :c:func:`lsolve()`, which may involve calling a generic
   function (such as for LU factorization).  This data may be intended
   either for direct use (in a direct linear solver) or for use in a
   preconditioner (in a preconditioned iterative linear solver).

   The :c:func:`lsetup()` function is not called at every stage solve
   (or even every time step), but only as frequently as the solver
   determines that it is appropriate to perform the setup task.  In
   this way, Jacobian-related data generated by :c:func:`lsetup()` is
   expected to be used over a number of time steps.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *convfail* -- an input flag used to indicate any problem that
	occurred during the solution of the nonlinear equation on the
	current time step for which the linear solver is being
	used. This flag can be used to help decide whether the
	Jacobian data kept by a linear solver needs to be
	updated or not. Its possible values are:

        - *ARK_NO_FAILURES*: this value is passed if either this is the
	  first call for this step, or the local error test failed on
	  the previous attempt at this step (but the Newton iteration
	  converged).
        - *ARK_FAIL_BAD_J*: this value is passed if (a) the previous
	  Newton corrector iteration did not converge and the linear
	  solver's setup function indicated that its Jacobian-related
	  data is not current, or (b) during the previous Newton
	  corrector iteration, the linear solver's solve function
	  failed in a recoverable manner and the linear solver's setup
	  function indicated that its Jacobian-related data is not
	  current.
        - *ARK_FAIL_OTHER*: this value is passed if during the current
	  internal step try, the previous Newton iteration failed to
	  converge even though the linear solver was using current
	  Jacobian-related data.

      * *ypred* -- is the predicted :math:`y` vector for the current
	ARKode internal step.
      * *fpred* -- is the value of the implicit right-hand side at
	*ypred*, :math:`f^I(t_n,ypred)`.
      * *jcurPtr* -- is a pointer to a boolean to be filled in by
	:c:func:`lsetup()`. The function should set ``*jcurPtr = SUNTRUE``
        if its Jacobian data is current after the call, and should set
	``*jcurPtr = SUNFALSE`` if its Jacobian data is not current. If
	:c:func:`lsetup()` calls for re-evaluation of Jacobian data
	(based on *convfail* and ARKode state data), it should return
	``*jcurPtr = SUNTRUE`` unconditionally; otherwise an infinite
	loop can result.
      * *vtemp1*, *vtemp2*, *vtemp3* -- are temporary variables of
	type ``N_Vector`` provided for use by :c:func:`lsetup()`.

   **Return value:**
   Should return 0 if successful, a positive value
   for a recoverable error, and a negative value for an unrecoverable
   error.  On a recoverable error return, the solver will attempt to
   recover by reducing the step size.


Similarly, the type definition of :c:func:`msetup()` is

.. c:function:: typedef int (*msetup)(ARKodeMem ark_mem, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)

   Prepares the mass matrix linear solver for subsequent calls to
   :c:func:`msolve()`, in the solution of systems :math:`M x = b`,
   where :math:`M` is the system mass matrix.

   The :c:func:`msetup()` function may call a user-supplied function,
   or a function within the linear solver module, to compute needed
   data related to the mass matrix.  Alterntively, it may choose to
   retrieve and use stored values of this data.

   In either case, :c:func:`msetup()` may also preprocess that data as
   needed for :c:func:`msolve()`, which may involve calling a generic
   function (such as for LU factorization).  This data may be intended
   either for direct use (in a direct linear solver) or for use in a
   preconditioner (in a preconditioned iterative linear solver).

   The :c:func:`msetup()` function is called at every time step, as
   discussed in section :ref:`Mathematics.MassSolve`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *vtemp1*, *vtemp2*, *vtemp3* -- are temporary variables of
	type ``N_Vector`` provided for use by :c:func:`msetup()`.

   **Return value:**
   Should return 0 if successful, a positive value
   for a recoverable error, and a negative value for an unrecoverable
   error.  On a recoverable error return, the solver will attempt to
   recover by reducing the step size.





Solve function
-----------------------------------

The type definition of :c:func:`lsolve()` is

.. c:function:: typedef int (*lsolve)(ARKodeMem ark_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector fcur)

   Solves the linear equation :math:`{\mathcal A} x = b`, where
   :math:`{\mathcal A}` arises  in the Newton iteration (see the
   section :ref:`Mathematics.Linear`) and gives some approximation to
   the Newton matrix :math:`M - \gamma J`, :math:`J =
   \frac{\partial}{\partial y} f^I(t_n, ycur)`. Note, the right-hand
   side vector  :math:`b` is input, and :math:`\gamma` is available as
   ``ark_mem->ark_gamma``.

   :c:func:`lsolve()` is called once per Newton iteration, hence possibly
   several times per time step.

   If there is an :c:func:`lsetup()` function, this :c:func:`lsolve()`
   function should make use of any Jacobian data that was computed and
   preprocessed by :c:func:`lsetup()`, either for direct use, or for
   use in a preconditioner.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *b* -- is the right-hand side vector :math:`b`. The solution
	is also to be returned in the vector :math:`b`.
      * *weight* -- is a vector that contains the residual weights. These
	are the :math:`rwt_i` of :ref:`ARKStep_CInterface.ResidualWeight`.
	This weight vector is included here to enable the computation
	of weighted norms needed to test for the convergence of
	iterative methods (if any) within the linear solver.
      * *ycur* -- is a vector that contains the solver's current
	approximation to :math:`y(t_n)`.
      * *fcur* -- is a vector that contains the current right-hand
         side, :math:`f^I(t_n, ycur)`.

   **Return value:**  Should return 0 if successful, a positive value
   for a recoverable error, and a negative value for an unrecoverable
   error.  On a recoverable error return, the solver will attempt to
   recover, such as by calling the :c:func:`lsetup()` function with
   the current arguments.


Similarly, the type definition of :c:func:`msolve()` is

.. c:function:: typedef int (*msolve)(ARKodeMem ark_mem, N_Vector b, N_Vector weight)

   Solves the linear equation :math:`M x = b`, where :math:`M` is the
   system mass matrix.  Note, the right-hand side vector :math:`b` is
   input, and holds the solution :math:`x` on output.

   :c:func:`msolve()` is called at least once per time step (if
   :math:`M\ne I`), as discussed in section :ref:`Mathematics.MassSolve`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *b* -- is the right-hand side vector :math:`b`. The solution
	is also to be returned in the vector :math:`b`.
      * *weight* -- is a vector that contains the error weights. These
	are the :math:`rwt_i` of :ref:`ARKStep_CInterface.ResidualWeight`.
	This weight vector is included here to enable the computation
	of weighted norms needed to test for the convergence of
	iterative methods (if any) within the linear solver.

   **Return value:**  Should return 0 if successful, and a nonzero
   value for an unrecoverable error.



Memory deallocation function
-----------------------------------

The type definition of :c:func:`lfree()` is

.. c:function:: typedef int (*lfree)(ARKodeMem ark_mem)

   free up any memory allocated by the linear solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.

   **Return value:**  This function should return 0 if successful, or
   a nonzero if not.

   **Notes:**  This function is called once a problem has been
   completed and the linear solver is no longer needed.


Similarly, the type definition of :c:func:`mfree()` is

.. c:function:: typedef int (*mfree)(ARKodeMem ark_mem)

   free up any memory allocated by the mass matrix linear solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.

   **Return value:**  This function should return 0 if successful, or
   a nonzero if not.

   **Notes:**  This function is called once a problem has been
   completed and the mass matrix solver is no longer needed.
