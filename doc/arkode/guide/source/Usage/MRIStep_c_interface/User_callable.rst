.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.MRIStep.UserCallable:

MRIStep User-callable functions
==================================

This section describes the functions that are called by the
user to setup and then solve an IVP using the MRIStep time-stepping
module. Some of these are required; however, starting with
:numref:`ARKODE.Usage.MRIStep.OptionalInputs`, the functions listed involve
optional inputs/outputs or restarting, and those paragraphs may be
skipped for a casual use of ARKODE's MRIStep module. In any case,
refer to the preceding section, :numref:`ARKODE.Usage.MRIStep.Skeleton`,
for the correct order of these calls.

On an error, each user-callable function returns a negative value  (or
``NULL`` if the function returns a pointer) and sends an error message
to the error handler routine, which prints the message to ``stderr``
by default. However, the user can set a file as error output or can
provide their own error handler function (see
:numref:`ARKODE.Usage.MRIStep.OptionalInputs` for details).



.. _ARKODE.Usage.MRIStep.Initialization:

MRIStep initialization and deallocation functions
------------------------------------------------------


.. c:function:: void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, realtype t0, N_Vector y0, MRIStepInnerStepper stepper, SUNContext sunctx)

   This function allocates and initializes memory for a problem to
   be solved using the MRIStep time-stepping module in ARKODE.

   **Arguments:**
      * *fse* -- the name of the function (of type :c:func:`ARKRhsFn()`)
        defining the explicit slow portion of the right-hand side function in
        :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.
      * *fsi* -- the name of the function (of type :c:func:`ARKRhsFn()`)
        defining the implicit slow portion of the right-hand side function in
        :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
      * *stepper* -- an :c:type:`MRIStepInnerStepper` for integrating the fast
        time scale.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a pointer to initialized problem memory of type ``void*``, to
      be passed to all user-facing MRIStep routines listed below.  If unsuccessful,
      a ``NULL`` pointer will be returned, and an error message will be printed to
      ``stderr``.

   **Example usage:**

      .. code-block:: C

         /* fast (inner) and slow (outer) ARKODE objects */
         void *inner_arkode_mem = NULL;
         void *outer_arkode_mem = NULL;

         /* MRIStepInnerStepper to wrap the inner (fast) ARKStep object */
         MRIStepInnerStepper stepper = NULL;

         /* create an ARKStep object, setting fast (inner) right-hand side
            functions and the initial condition */
         inner_arkode_mem = ARKStepCreate(ffe, ffi, t0, y0, sunctx);

         /* setup ARKStep */
         . . .

         /* create MRIStepInnerStepper wrapper for the ARKStep memory block */
         flag = ARKStepCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);

         /* create an MRIStep object, setting the slow (outer) right-hand side
            functions and the initial condition */
         outer_arkode_mem = MRIStepCreate(fse, fsi, t0, y0, stepper, sunctx)

   **Example codes:**
      * ``examples/arkode/C_serial/ark_brusselator_mri.c``
      * ``examples/arkode/C_serial/ark_twowaycouple_mri.c``
      * ``examples/arkode/C_serial/ark_brusselator_1D_mri.c``
      * ``examples/arkode/C_serial/ark_onewaycouple_mri.c``
      * ``examples/arkode/C_serial/ark_reaction_diffusion_mri.c``
      * ``examples/arkode/C_serial/ark_kpr_mri.c``
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``


.. c:function:: void MRIStepFree(void** arkode_mem)

   This function frees the problem memory *arkode_mem* created by
   :c:func:`MRIStepCreate`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.

   **Return value:**  None



.. _ARKODE.Usage.MRIStep.Tolerances:

MRIStep tolerance specification functions
------------------------------------------------------

These functions specify the integration tolerances. One of them
**should** be called before the first call to
:c:func:`MRIStepEvolve()`; otherwise default values of ``reltol =
1e-4`` and ``abstol = 1e-9`` will be used, which may be entirely
incorrect for a specific problem.

The integration tolerances ``reltol`` and ``abstol`` define a vector
of error weights, ``ewt``.  In the case of
:c:func:`MRIStepSStolerances()`, this vector has components

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol);

whereas in the case of :c:func:`MRIStepSVtolerances()` the vector components
are given by

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol[i]);

This vector is used in all error tests, which use a weighted RMS norm
on all error-like vectors :math:`v`:

.. math::
    \|v\|_{WRMS} = \left( \frac{1}{N} \sum_{i=1}^N (v_i\; ewt_i)^2 \right)^{1/2},

where :math:`N` is the problem dimension.

Alternatively, the user may supply a custom function to supply the
``ewt`` vector, through a call to :c:func:`MRIStepWFtolerances()`.



.. c:function:: int MRIStepSStolerances(void* arkode_mem, realtype reltol, realtype abstol)

   This function specifies scalar relative and absolute tolerances.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *reltol* -- scalar relative tolerance.
      * *abstol* -- scalar absolute tolerance.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the MRIStep memory was ``NULL``
      * *ARK_NO_MALLOC*  if the MRIStep memory was not allocated by the time-stepping module
      * *ARK_ILL_INPUT* if an argument has an illegal value (e.g. a negative tolerance).



.. c:function:: int MRIStepSVtolerances(void* arkode_mem, realtype reltol, N_Vector abstol)

   This function specifies a scalar relative tolerance and a vector
   absolute tolerance (a potentially different absolute tolerance for
   each vector component).

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *reltol* -- scalar relative tolerance.
      * *abstol* -- vector containing the absolute tolerances for each
        solution component.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the MRIStep memory was ``NULL``
      * *ARK_NO_MALLOC*  if the MRIStep memory was not allocated by the time-stepping module
      * *ARK_ILL_INPUT* if an argument has an illegal value (e.g. a negative tolerance).



.. c:function:: int MRIStepWFtolerances(void* arkode_mem, ARKEwtFn efun)

   This function specifies a user-supplied function *efun* to compute
   the error weight vector ``ewt``.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *efun* -- the name of the function (of type :c:func:`ARKEwtFn()`)
        that implements the error weight vector computation.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the MRIStep memory was ``NULL``
      * *ARK_NO_MALLOC*  if the MRIStep memory was not allocated by the time-stepping module




General advice on the choice of tolerances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For many users, the appropriate choices for tolerance values in
``reltol`` and ``abstol`` are a concern. The following pieces
of advice are relevant.

(1) The scalar relative tolerance ``reltol`` is to be set to control
    relative errors. So a value of :math:`10^{-4}` means that errors
    are controlled to .01%. We do not recommend using ``reltol`` larger
    than :math:`10^{-3}`. On the other hand, ``reltol`` should not be so
    small that it is comparable to the unit roundoff of the machine
    arithmetic (generally around :math:`10^{-15}` for double-precision).

(2) The absolute tolerances ``abstol`` (whether scalar or vector) need
    to be set to control absolute errors when any components of the
    solution vector :math:`y` may be so small that pure relative error
    control is meaningless.  For example, if :math:`y_i` starts at some
    nonzero value, but in time decays to zero, then pure relative
    error control on :math:`y_i` makes no sense (and is overly costly)
    after :math:`y_i` is below some noise level. Then ``abstol`` (if
    scalar) or ``abstol[i]`` (if a vector) needs to be set to that
    noise level. If the different components have different noise
    levels, then ``abstol`` should be a vector.  For example, see the
    example problem ``ark_robertson.c``, and the discussion
    of it in the ARKODE Examples Documentation :cite:p:`arkode_ex`.  In that
    problem, the three components vary between 0 and 1, and have
    different noise levels; hence the ``atols`` vector therein. It is
    impossible to give any general advice on ``abstol`` values,
    because the appropriate noise levels are completely
    problem-dependent. The user or modeler hopefully has some idea as
    to what those noise levels are.

(3) Finally, it is important to pick all the tolerance values
    conservatively, because they control the error committed on each
    individual step. The final (global) errors are an accumulation of
    those per-step errors, where that accumulation factor is
    problem-dependent.  A general rule of thumb is to reduce the
    tolerances by a factor of 10 from the actual desired limits on
    errors.  So if you want .01% relative accuracy (globally), a good
    choice for ``reltol`` is :math:`10^{-5}`.  In any case, it is
    a good idea to do a few experiments with the tolerances to see how
    the computed solution values vary as tolerances are reduced.



Advice on controlling nonphysical negative values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many applications, some components in the true solution are always
positive or non-negative, though at times very small.  In the
numerical solution, however, small negative (nonphysical) values
can then occur. In most cases, these values are harmless, and simply
need to be controlled, not eliminated, but in other cases any value
that violates a constraint may cause a simulation to halt. For both of
these scenarios the following pieces of advice are relevant.

(1) The best way to control the size of unwanted negative computed
    values is with tighter absolute tolerances.  Again this requires
    some knowledge of the noise level of these components, which may
    or may not be different for different components. Some
    experimentation may be needed.

(2) If output plots or tables are being generated, and it is important
    to avoid having negative numbers appear there (for the sake of
    avoiding a long explanation of them, if nothing else), then
    eliminate them, but only in the context of the output medium. Then
    the internal values carried by the solver are unaffected. Remember
    that a small negative value in :math:`y` returned by MRIStep, with
    magnitude comparable to ``abstol`` or less, is equivalent to zero
    as far as the computation is concerned.

(3) The user's right-hand side routine :math:`f^I`
    should never change a negative value in the solution vector :math:`y`
    to a non-negative value in attempt to "fix" this problem,
    since this can lead to numerical instability.  If the :math:`f^I`
    routine cannot tolerate a zero or negative value (e.g. because
    there is a square root or log), then the offending value should be
    changed to zero or a tiny positive number in a temporary variable
    (not in the input :math:`y` vector) for the purposes of computing
    :math:`f^I(t, y)`.

..
   (4) Positivity and non-negativity constraints on components can be
       enforced by use of the recoverable error return feature in the
       user-supplied right-hand side function, :math:`f^I`. When a
       recoverable error is encountered, MRIStep will retry the step with
       a smaller step size, which typically alleviates the problem.
       However, because this option involves some additional overhead
       cost, it should only be exercised if the use of absolute
       tolerances to control the computed values is unsuccessful.



.. _ARKODE.Usage.MRIStep.LinearSolvers:

Linear solver interface functions
-------------------------------------------

As previously explained, the Newton iterations used in solving
implicit systems within MRIStep require the solution of linear
systems of the form

.. math::
   \mathcal{A}\left(z_i^{(m)}\right) \delta^{(m+1)} = -G\left(z_i^{(m)}\right)

where

.. math::
   \mathcal{A} \approx I - \gamma J, \qquad J = \frac{\partial f^I}{\partial y}.

ARKODE's ARKLS linear solver interface supports all valid
``SUNLinearSolver`` modules for this task.

Matrix-based ``SUNLinearSolver`` modules utilize ``SUNMatrix`` objects
to store the approximate Jacobian matrix :math:`J`, the Newton matrix
:math:`\mathcal{A}`, and, when using direct solvers, the factorizations
used throughout the solution process.

Matrix-free ``SUNLinearSolver`` modules instead use iterative methods
to solve the Newton systems of equations, and only require the
*action* of the matrix on a vector, :math:`\mathcal{A}v`.  With most
of these methods, preconditioning can be done on the left only, on the
right only, on both the left and the right, or not at all.  The
exceptions to this rule are SPFGMR that supports right preconditioning
only and PCG that performs symmetric preconditioning.  For the
specification of a preconditioner, see the iterative linear solver
portions of :numref:`ARKODE.Usage.MRIStep.OptionalInputs` and
:numref:`ARKODE.Usage.UserSupplied`.

If preconditioning is done, user-supplied functions should be used to
define left and right preconditioner matrices :math:`P_1` and
:math:`P_2` (either of which could be the identity matrix), such that
the product :math:`P_{1}P_{2}` approximates the Newton matrix
:math:`\mathcal{A} = I - \gamma J`.

To specify a generic linear solver for MRIStep to use for the Newton
systems, after the call to :c:func:`MRIStepCreate()` but before any
calls to :c:func:`MRIStepEvolve()`, the user's program must create the
appropriate ``SUNLinearSolver`` object and call the function
:c:func:`MRIStepSetLinearSolver()`, as documented below.  To create
the ``SUNLinearSolver`` object, the user may call one of the
SUNDIALS-packaged SUNLinSol module constructor routines via a call of
the form

.. code:: c

   SUNLinearSolver LS = SUNLinSol_*(...);

The current list of SUNDIALS-packaged SUNLinSol modules, and their
constructor routines, may be found in chapter :numref:`SUNLinSol`.
Alternately, a user-supplied ``SUNLinearSolver`` module may be created
and used.  Specific information on how to create such user-provided
modules may be found in :numref:`SUNLinSol.API.Custom`.

Once this solver object has been constructed, the user should attach
it to MRIStep via a call to :c:func:`MRIStepSetLinearSolver()`. The
first argument passed to this function is the MRIStep memory pointer
returned by :c:func:`MRIStepCreate()`; the second argument is the
``SUNLinearSolver`` object created above.  The third argument is an
optional ``SUNMatrix`` object to accompany matrix-based
``SUNLinearSolver`` inputs (for matrix-free linear solvers, the third
argument should be ``NULL``).  A call to this function initializes the
ARKLS linear solver interface, linking it to the MRIStep integrator,
and allows the user to specify additional parameters and routines
pertinent to their choice of linear solver.

.. c:function:: int MRIStepSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix J)

   This function specifies the ``SUNLinearSolver`` object that MRIStep
   should use, as well as a template Jacobian ``SUNMatrix`` object (if
   applicable).

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *LS* -- the ``SUNLinearSolver`` object to use.
      * *J* -- the template Jacobian ``SUNMatrix`` object to use (or
        ``NULL`` if not applicable).

   **Return value:**
      * *ARKLS_SUCCESS*   if successful
      * *ARKLS_MEM_NULL*  if the MRIStep memory was ``NULL``
      * *ARKLS_MEM_FAIL*  if there was a memory allocation failure
      * *ARKLS_ILL_INPUT* if ARKLS is incompatible with the
        provided *LS* or *J* input objects, or the current
        ``N_Vector`` module.

   **Notes:**  If *LS* is a matrix-free linear solver, then the *J*
   argument should be ``NULL``.

   If *LS* is a matrix-based linear solver, then the template Jacobian
   matrix *J* will be used in the solve process, so if additional
   storage is required within the ``SUNMatrix`` object (e.g. for
   factorization of a banded matrix), ensure that the input object is
   allocated with sufficient size (see the documentation of
   the particular SUNMATRIX type in :numref:`SUNMatrix` for
   further information).

   When using sparse linear solvers, it is typically much more
   efficient to supply *J* so that it includes the full sparsity
   pattern of the Newton system matrices :math:`\mathcal{A} =
   I-\gamma J`, even if *J* itself has zeros in nonzero
   locations of :math:`I`.  The reasoning for this is
   that :math:`\mathcal{A}` is constructed in-place, on top of the
   user-specified values of *J*, so if the sparsity pattern in *J* is
   insufficient to store :math:`\mathcal{A}` then it will need to be
   resized internally by MRIStep.



.. _ARKODE.Usage.MRIStep.NonlinearSolvers:

Nonlinear solver interface functions
-------------------------------------------

When changing the nonlinear solver in MRIStep, after the
call to :c:func:`MRIStepCreate()` but before any calls to
:c:func:`MRIStepEvolve()`, the user's program must create the
appropriate SUNNonlinSol object and call
:c:func:`MRIStepSetNonlinearSolver()`, as documented below.  If any
calls to :c:func:`MRIStepEvolve()` have been made, then MRIStep will
need to be reinitialized by calling :c:func:`MRIStepReInit()` to
ensure that the nonlinear solver is initialized correctly before any
subsequent calls to :c:func:`MRIStepEvolve()`.

The first argument passed to the routine
:c:func:`MRIStepSetNonlinearSolver()` is the MRIStep memory pointer
returned by :c:func:`MRIStepCreate()`; the second argument passed
to this function is the desired ``SUNNonlinearSolver`` object to use for
solving the nonlinear system for each implicit stage. A call to this
function attaches the nonlinear solver to the main MRIStep integrator.


.. c:function:: int MRIStepSetNonlinearSolver(void* arkode_mem, SUNNonlinearSolver NLS)

   This function specifies the ``SUNNonlinearSolver`` object
   that MRIStep should use for implicit stage solves.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *NLS* -- the ``SUNNonlinearSolver`` object to use.

   **Return value:**
      * *ARK_SUCCESS*   if successful
      * *ARK_MEM_NULL*  if the MRIStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if there was a memory allocation failure
      * *ARK_ILL_INPUT* if MRIStep is incompatible with the
        provided *NLS* input object.

   **Notes:**  MRIStep will use the Newton ``SUNNonlinearSolver`` module by
   default; a call to this routine replaces that module with the
   supplied *NLS* object.



.. _ARKODE.Usage.MRIStep.RootFinding:

Rootfinding initialization function
--------------------------------------

As described in the section :numref:`ARKODE.Mathematics.Rootfinding`, while
solving the IVP, ARKODE's time-stepping modules have the capability to
find the roots of a set of user-defined functions.  In the MRIStep module root
finding is performed between slow solution time steps only (i.e., it is not
performed within the sub-stepping a fast time scales).  To activate the
root-finding algorithm, call the following function.  This is normally
called only once, prior to the first call to
:c:func:`MRIStepEvolve()`, but if the rootfinding problem is to be
changed during the solution, :c:func:`MRIStepRootInit()` can also be
called prior to a continuation call to :c:func:`MRIStepEvolve()`.


.. c:function:: int MRIStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)

   Initializes a rootfinding problem to be solved during the
   integration of the ODE system.  It must be called after
   :c:func:`MRIStepCreate()`, and before :c:func:`MRIStepEvolve()`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *nrtfn* -- number of functions :math:`g_i`, an integer :math:`\ge` 0.
      * *g* -- name of user-supplied function, of type :c:func:`ARKRootFn()`,
        defining the functions :math:`g_i` whose roots are sought.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the MRIStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if there was a memory allocation failure
      * *ARK_ILL_INPUT* if *nrtfn* is greater than zero but *g* = ``NULL``.

   **Notes:** To disable the rootfinding feature after it has already
   been initialized, or to free memory associated with MRIStep's
   rootfinding module, call *MRIStepRootInit* with *nrtfn = 0*.

   Similarly, if a new IVP is to be solved with a call to
   :c:func:`MRIStepReInit()`, where the new IVP has no rootfinding
   problem but the prior one did, then call *MRIStepRootInit* with
   *nrtfn = 0*.

   Rootfinding is only supported for the slow (outer) integrator and should not
   be actived for the fast (inner) integrator.



.. _ARKODE.Usage.MRIStep.Integration:

MRIStep solver function
-------------------------

This is the central step in the solution process -- the call to perform
the integration of the IVP.  The input argument *itask* specifies one of two
modes as to where MRIStep is to return a solution.  These modes are modified if
the user has set a stop time (with a call to the optional input function
:c:func:`MRIStepSetStopTime()`) or has requested rootfinding.


.. c:function:: int MRIStepEvolve(void* arkode_mem, realtype tout, N_Vector yout, realtype *tret, int itask)

   Integrates the ODE over an interval in :math:`t`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *tout* -- the next time at which a computed solution is desired.
      * *yout* -- the computed solution vector.
      * *tret* -- the time corresponding to *yout* (output).
      * *itask* -- a flag indicating the job of the solver for the next
        user step.

        The *ARK_NORMAL* option causes the solver to take internal
        steps until it has just overtaken a user-specified output
        time, *tout*, in the direction of integration,
        i.e. :math:`t_{n-1} <` *tout* :math:`\le t_{n}` for forward
        integration, or :math:`t_{n} \le` *tout* :math:`< t_{n-1}` for
        backward integration.  It will then compute an approximation
        to the solution :math:`y(tout)` by interpolation (as described
        in :numref:`ARKODE.Mathematics.Interpolation`).

        The *ARK_ONE_STEP* option tells the solver to only take a
        single internal step :math:`y_{n-1} \to y_{n}` and then return
        control back to the calling program.  If this step will
        overtake *tout* then the solver will again return an
        interpolated result; otherwise it will return a copy of the
        internal solution :math:`y_{n}` in the vector *yout*.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_ROOT_RETURN* if :c:func:`MRIStepEvolve()` succeeded, and
        found one or more roots.  If the number of root functions,
        *nrtfn*, is greater than 1, call
        :c:func:`MRIStepGetRootInfo()` to see which :math:`g_i` were
        found to have a root at (*\*tret*).
      * *ARK_TSTOP_RETURN* if :c:func:`MRIStepEvolve()` succeeded and
        returned at *tstop*.
      * *ARK_MEM_NULL* if the *arkode_mem* argument was ``NULL``.
      * *ARK_NO_MALLOC* if *arkode_mem* was not allocated.
      * *ARK_ILL_INPUT* if one of the inputs to
        :c:func:`MRIStepEvolve()` is illegal, or some other input to
        the solver was either illegal or missing.  Details will be
        provided in the error message.  Typical causes of this failure:

        (a) A component of the error weight vector became zero during
            internal time-stepping.

        (b) The linear solver initialization function (called by the
            user after calling :c:func:`ARKStepCreate`) failed to set
            the linear solver-specific *lsolve* field in
            *arkode_mem*.

        (c) A root of one of the root functions was found both at a
            point :math:`t` and also very near :math:`t`.

      * *ARK_TOO_MUCH_WORK* if the solver took *mxstep* internal steps
        but could not reach *tout*.  The default value for *mxstep* is
        *MXSTEP_DEFAULT = 500*.
      * *ARK_CONV_FAILURE* if convergence test failures occurred
        too many times (*ark_maxncf*) during one internal time step.
      * *ARK_LINIT_FAIL* if the linear solver's initialization
        function failed.
      * *ARK_LSETUP_FAIL* if the linear solver's setup routine failed in
        an unrecoverable manner.
      * *ARK_LSOLVE_FAIL* if the linear solver's solve routine failed in
        an unrecoverable manner.
      * *ARK_VECTOROP_ERR* a vector operation error occurred.
      * *ARK_INNERSTEP_FAILED* if the inner stepper returned with an
        unrecoverable error. The value returned from the inner stepper can be
        obtained with :c:func:`MRIStepGetLastInnerStepFlag()`.
      * *ARK_INVALID_TABLE* if an invalid coupling table was provided.

   **Notes:**
      The input vector *yout* can use the same memory as the
      vector *y0* of initial conditions that was passed to
      :c:func:`MRIStepCreate`.

      In *ARK_ONE_STEP* mode, *tout* is used only on the first call, and
      only to get the direction and a rough scale of the independent
      variable.

      All failure return values are negative and so testing the return argument for
      negative values will trap all :c:func:`MRIStepEvolve()` failures.

      Since interpolation may reduce the accuracy in the reported
      solution, if full method accuracy is desired the user should issue
      a call to :c:func:`MRIStepSetStopTime()` before the call to
      :c:func:`MRIStepEvolve()` to specify a fixed stop time to
      end the time step and return to the user.  Upon return from
      :c:func:`MRIStepEvolve()`, a copy of the internal solution
      :math:`y_{n}` will be returned in the vector *yout*.  Once the
      integrator returns at a *tstop* time, any future testing for
      *tstop* is disabled (and can be re-enabled only though a new call
      to :c:func:`MRIStepSetStopTime()`).

      On any error return in which one or more internal steps were taken
      by :c:func:`MRIStepEvolve()`, the returned values of *tret* and
      *yout* correspond to the farthest point reached in the integration.
      On all other error returns, *tret* and *yout* are left unchanged
      from those provided to the routine.



.. _ARKODE.Usage.MRIStep.OptionalInputs:

Optional input functions
-------------------------

There are numerous optional input parameters that control the behavior
of MRIStep, each of which may be modified from its default value through
calling an appropriate input function.  The following tables list all
optional input functions, grouped by which aspect of MRIStep they control.
Detailed information on the calling syntax and arguments for each
function are then provided following each table.

The optional inputs are grouped into the following categories:

* General MRIStep options (:numref:`ARKODE.Usage.MRIStep.MRIStepInput`),

* IVP method solver options (:numref:`ARKODE.Usage.MRIStep.MRIStepMethodInput`),

* Implicit stage solver options (:numref:`ARKODE.Usage.MRIStep.MRIStepSolverInput`),

* Linear solver interface options (:numref:`ARKODE.Usage.MRIStep.ARKLsInputs`), and

* Rootfinding options (:numref:`ARKODE.Usage.MRIStep.MRIStepRootfindingInput`).

For the most casual use of MRIStep, relying on the default set of
solver parameters, the reader can skip to the section on user-supplied
functions, :numref:`ARKODE.Usage.UserSupplied`.

We note that, on an error return, all of the optional input functions send an
error message to the error handler function. All error return values are
negative, so a test on the return arguments for negative values will catch all
errors. Finally, a call to an ``MRIStepSet***`` function can generally be made
from the user's calling program at any time and, if successful, takes effect
immediately. ``MRIStepSet***`` functions that cannot be called at any time note
this in the "**Notes**:" section of the function documentation.



.. _ARKODE.Usage.MRIStep.MRIStepInput:

Optional inputs for MRIStep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.MRIStep.MRIStepInput.Table:
.. table:: Optional inputs for MRIStep

   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Optional input                                                | Function name                             | Default                |
   +===============================================================+===========================================+========================+
   | Return MRIStep solver parameters to their defaults            | :c:func:`MRIStepSetDefaults()`            | internal               |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Set dense output interpolation type                           | :c:func:`MRIStepSetInterpolantType()`     | ``ARK_INTERP_HERMITE`` |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Set dense output polynomial degree                            | :c:func:`MRIStepSetInterpolantDegree()`   | 5                      |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Supply a pointer to a diagnostics output file                 | :c:func:`MRIStepSetDiagnostics()`         | ``NULL``               |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Supply a pointer to an error output file                      | :c:func:`MRIStepSetErrFile()`             | ``stderr``             |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Supply a custom error handler function                        | :c:func:`MRIStepSetErrHandlerFn()`        | internal fn            |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Run with fixed-step sizes                                     | :c:func:`MRIStepSetFixedStep()`           | required               |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Maximum no. of warnings for :math:`t_n+h = t_n`               | :c:func:`MRIStepSetMaxHnilWarns()`        | 10                     |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Maximum no. of internal steps before *tout*                   | :c:func:`MRIStepSetMaxNumSteps()`         | 500                    |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Set a value for :math:`t_{stop}`                              | :c:func:`MRIStepSetStopTime()`            | undefined              |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Interpolate at :math:`t_{stop}`                               | :c:func:`MRIStepSetInterpolateStopTime()` | ``SUNFALSE``           |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Disable the stop time                                         | :c:func:`MRIStepClearStopTime`            | N/A                    |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Supply a pointer for user data                                | :c:func:`MRIStepSetUserData()`            | ``NULL``               |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Supply a function to be called prior to the inner integration | :c:func:`MRIStepSetPreInnerFn()`          | ``NULL``               |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+
   | Supply a function to be called after the inner integration    | :c:func:`MRIStepSetPostInnerFn()`         | ``NULL``               |
   +---------------------------------------------------------------+-------------------------------------------+------------------------+




.. c:function:: int MRIStepSetDefaults(void* arkode_mem)

   Resets all optional input parameters to MRIStep's original
   default values.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

  **Notes:** This function does not change problem-defining function pointers
  *fs* and *ff* or the *user_data* pointer. It also does not affect any data
  structures or options related to root-finding (those can be reset using
  :c:func:`MRIStepRootInit()`).



.. c:function:: int MRIStepSetInterpolantType(void* arkode_mem, int itype)

   Specifies use of the Lagrange or Hermite interpolation modules (used for
   dense output -- interpolation of solution output values and implicit
   method predictors).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *itype* -- requested interpolant type (``ARK_INTERP_HERMITE`` or ``ARK_INTERP_LAGRANGE``)

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_MEM_FAIL* if the interpolation module cannot be allocated

   * *ARK_ILL_INPUT* if the *itype* argument is not recognized or the
     interpolation module has already been initialized

   **Notes:** The Hermite interpolation module is described in
   :numref:`ARKODE.Mathematics.Interpolation.Hermite`, and the Lagrange interpolation module
   is described in :numref:`ARKODE.Mathematics.Interpolation.Lagrange`.

   This routine frees any previously-allocated interpolation module, and re-creates
   one according to the specified argument.  Thus any previous calls to
   :c:func:`MRIStepSetInterpolantDegree()` will be nullified.

   This routine must be called *after* the call to :c:func:`MRIStepCreate()`.
   After the first call to :c:func:`MRIStepEvolve()` the interpolation type may
   not be changed without first calling :c:func:`MRIStepReInit()`.

   If this routine is not called, the Hermite interpolation module will be used.



.. c:function:: int MRIStepSetInterpolantDegree(void* arkode_mem, int degree)

   Specifies the degree of the polynomial interpolant
   used for dense output (i.e. interpolation of solution output values
   and implicit method predictors).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *degree* -- requested polynomial degree.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory or interpolation module are ``NULL``

   * *ARK_INTERP_FAIL* if this is called after :c:func:`MRIStepEvolve()`

   * *ARK_ILL_INPUT* if an argument has an illegal value or the
     interpolation module has already been initialized

   **Notes:** Allowed values are between 0 and 5.

   This routine should be called *after* :c:func:`MRIStepCreate()` and *before*
   :c:func:`MRIStepEvolve()`. After the first call to :c:func:`MRIStepEvolve()`
   the interpolation degree may not be changed without first calling
   :c:func:`MRIStepReInit()`.


   If a user calls both this routine and :c:func:`MRIStepSetInterpolantType()`, then
   :c:func:`MRIStepSetInterpolantType()` must be called first.

   Since the accuracy of any polynomial interpolant is limited by the accuracy
   of the time-step solutions on which it is based, the *actual* polynomial
   degree that is used by MRIStep will be the minimum of :math:`q-1` and the
   input *degree*, for :math:`q > 1` where :math:`q` is the order of accuracy
   for the time integration method.

   .. versionchanged:: 5.5.1

      When :math:`q=1`, a linear interpolant is the default to ensure values
      obtained by the integrator are returned at the ends of the time interval.



.. c:function:: int MRIStepSetDenseOrder(void* arkode_mem, int dord)

   *This function is deprecated, and will be removed in a future release.
   Users should transition to calling* :c:func:`MRIStepSetInterpolantDegree()`
   *instead.*



.. c:function:: int MRIStepSetDiagnostics(void* arkode_mem, FILE* diagfp)

   Specifies the file pointer for a diagnostics file where
   all MRIStep step adaptivity and solver information is written.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *diagfp* -- pointer to the diagnostics output file.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** This parameter can be ``stdout`` or ``stderr``, although the
   suggested approach is to specify a pointer to a unique file opened
   by the user and returned by ``fopen``.  If not called, or if called
   with a ``NULL`` file pointer, all diagnostics output is disabled.

   When run in parallel, only one process should set a non-NULL value
   for this pointer, since statistics from all processes would be
   identical.

   .. deprecated:: 5.2.0

      Use :c:func:`SUNLogger_SetInfoFilename` instead.


.. c:function:: int MRIStepSetErrFile(void* arkode_mem, FILE* errfp)

   Specifies a pointer to the file where all MRIStep warning and error
   messages will be written if the default internal error handling
   function is used.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *errfp* -- pointer to the output file.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** The default value for *errfp* is ``stderr``.

   Passing a ``NULL`` value disables all future error message output
   (except for the case wherein the MRIStep memory pointer is
   ``NULL``).  This use of the function is strongly discouraged.

   If used, this routine should be called before any other
   optional input functions, in order to take effect for subsequent
   error messages.



.. c:function:: int MRIStepSetErrHandlerFn(void* arkode_mem, ARKErrHandlerFn ehfun, void* eh_data)

   Specifies the optional user-defined function to be used
   in handling error messages.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *ehfun* -- name of user-supplied error handler function.

   * *eh_data* -- pointer to user data passed to *ehfun* every time
     it is called.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** Error messages indicating that the MRIStep solver memory is
   ``NULL`` will always be directed to ``stderr``.




.. c:function:: int MRIStepSetFixedStep(void* arkode_mem, realtype hs)

   Set the slow step size used within MRIStep for the following internal step(s).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *hs* -- value of the outer (slow) step size.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**

   The step sizes used by the inner (fast) stepper may be controlled through calling the
   appropriate "Set" routines on the inner integrator.



..
   .. c:function:: int MRIStepSetInitStep(void* arkode_mem, realtype hin)

      Specifies the initial time step size MRIStep should use after
      initialization or re-initialization.

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *hin* -- value of the initial step to be attempted :math:`(\ne 0)`.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

      * *ARK_ILL_INPUT* if an argument has an illegal value

      **Notes:** Pass 0.0 to use the default value.

      By default, MRIStep estimates the initial step size to be the
      solution :math:`h` of the equation :math:`\left\| \frac{h^2
      \ddot{y}}{2}\right\| = 1`, where :math:`\ddot{y}` is an estimated
      value of the second derivative of the solution at *t0*.




.. c:function:: int MRIStepSetMaxHnilWarns(void* arkode_mem, int mxhnil)

   Specifies the maximum number of messages issued by the
   solver to warn that :math:`t+h=t` on the next internal step, before
   MRIStep will instead return with an error.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *mxhnil* -- maximum allowed number of warning messages :math:`(>0)`.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** The default value is 10; set *mxhnil* to zero to specify
   this default.

   A negative value indicates that no warning messages should be issued.




.. c:function:: int MRIStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)

   Specifies the maximum number of steps to be taken by the
   solver in its attempt to reach the next output time, before MRIStep
   will return with an error.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *mxsteps* -- maximum allowed number of internal steps.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** Passing *mxsteps* = 0 results in MRIStep using the
   default value (500).

   Passing *mxsteps* < 0 disables the test (not recommended).



..
   .. c:function:: int MRIStepSetMaxStep(void* arkode_mem, realtype hmax)

      Specifies the upper bound on the magnitude of the time step size.

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *hmax* -- maximum absolute value of the time step size :math:`(\ge 0)`.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

      * *ARK_ILL_INPUT* if an argument has an illegal value

      **Notes:** Pass *hmax* :math:`\le 0.0` to set the default value of :math:`\infty`.



..
   .. c:function:: int MRIStepSetMinStep(void* arkode_mem, realtype hmin)

      Specifies the lower bound on the magnitude of the time step size.

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *hmin* -- minimum absolute value of the time step size :math:`(\ge 0)`.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

      * *ARK_ILL_INPUT* if an argument has an illegal value

      **Notes:** Pass *hmin* :math:`\le 0.0` to set the default value of 0.


.. c:function:: int MRIStepSetStopTime(void* arkode_mem, realtype tstop)

   Specifies the value of the independent variable
   :math:`t` past which the solution is not to proceed.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *tstop* -- stopping time for the integrator.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**

      The default is that no stop time is imposed.

      Once the integrator returns at a stop time, any future testing for
      ``tstop`` is disabled (and can be reenabled only though a new call to
      :c:func:`MRIStepSetStopTime`).

      A stop time not reached before a call to :c:func:`MRIStepReInit` or
      :c:func:`MRIStepReset` will remain active but can be disabled by calling
      :c:func:`MRIStepClearStopTime`.


.. c:function:: int MRIStepSetInterpolateStopTime(void* arkode_mem, booleantype interp)

   Specifies that the output solution should be interpolated when the current
   :math:`t` equals the specified ``tstop`` (instead of merely copying the
   internal solution :math:`y_n`).

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *interp* -- flag indicating to use interpolation (1) or copy (0).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory is ``NULL``

   .. versionadded:: 5.6.0


.. c:function:: int MRIStepClearStopTime(void* arkode_mem)

   Disables the stop time set with :c:func:`MRIStepSetStopTime`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   **Notes:**
      The stop time can be reenabled though a new call to
      :c:func:`MRIStepSetStopTime`.

   .. versionadded:: 5.5.1


.. c:function:: int MRIStepSetUserData(void* arkode_mem, void* user_data)

   Specifies the user data block *user_data* for the outer integrator and
   attaches it to the main MRIStep memory block.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *user_data* -- pointer to the user data.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** If specified, the pointer to *user_data* is passed to all
   user-supplied functions called by the outer integrator for which it is an
   argument; otherwise ``NULL`` is passed.

   To attach a user data block to the inner integrator call the appropriate
   *SetUserData* function for the inner integrator memory structure (e.g.,
   :c:func:`ARKStepSetUserData()` if the inner stepper is ARKStep). This pointer
   may be the same as or different from the pointer attached to the outer
   integrator depending on what is required by the user code.


.. c:function:: int MRIStepSetPreInnerFn(void* arkode_mem, MRIStepPreInnerFn prefn)

   Specifies the function called *before* each inner integration.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *prefn* -- the name of the C function (of type :c:func:`MRIStepPreInnerFn()`)
     defining pre inner integration function.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``


.. c:function:: int MRIStepSetPostInnerFn(void* arkode_mem, MRIStepPostInnerFn postfn)

   Specifies the function called *after* each inner integration.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *postfn* -- the name of the C function (of type :c:func:`MRIStepPostInnerFn()`)
     defining post inner integration function.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

..
   .. c:function:: int MRIStepSetMaxErrTestFails(void* arkode_mem, int maxnef)

      Specifies the maximum number of error test failures
      permitted in attempting one step, before returning with an error.

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *maxnef* -- maximum allowed number of error test failures :math:`(>0)`.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

      * *ARK_ILL_INPUT* if an argument has an illegal value

      **Notes:** The default value is 7; set *maxnef* :math:`\le 0`
      to specify this default.



.. _ARKODE.Usage.MRIStep.MRIStepMethodInput:

Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.MRIStep.MRIStepMethodInputTable:
.. table:: Optional inputs for IVP method selection

   +--------------------------------+-------------------------------------+----------+
   | Optional input                 | Function name                       | Default  |
   +--------------------------------+-------------------------------------+----------+
   | Select the default MRI method  | :c:func:`MRIStepSetOrder()`         | 3        |
   | of a given order               |                                     |          |
   +--------------------------------+-------------------------------------+----------+
   | Set MRI coupling coefficients  | :c:func:`MRIStepSetCoupling()`      | internal |
   +--------------------------------+-------------------------------------+----------+


.. c:function:: int MRIStepSetOrder(void* arkode_mem, int ord)

   Select the default MRI method of a given order.

   The default order is 3. An order less than 3 or greater than 4 will result in
   using the default.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *ord* -- the method order.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``


.. c:function:: int MRIStepSetCoupling(void* arkode_mem, MRIStepCoupling C)

   Specifies a customized set of slow-to-fast coupling coefficients for the MRI method.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *C* -- the table of coupling coefficients for the MRI method.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**

   For a description of the :c:type:`MRIStepCoupling` type and related
   functions for creating Butcher tables see :numref:`ARKODE.Usage.MRIStep.MRIStepCoupling`.



.. _ARKODE.Usage.MRIStep.MRIStepSolverInput:

Optional inputs for implicit stage solves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation for the nonlinear solver strategies used
by MRIStep, including how each of the parameters below is used within
the code, is provided in :numref:`ARKODE.Mathematics.Nonlinear`.


.. cssclass:: table-bordered

=========================================================  =========================================  ============
Optional input                                             Function name                              Default
=========================================================  =========================================  ============
Specify linearly implicit :math:`f^I`                      :c:func:`MRIStepSetLinear()`               ``SUNFALSE``
Specify nonlinearly implicit :math:`f^I`                   :c:func:`MRIStepSetNonlinear()`            ``SUNTRUE``
Implicit predictor method                                  :c:func:`MRIStepSetPredictorMethod()`      0
Maximum number of nonlinear iterations                     :c:func:`MRIStepSetMaxNonlinIters()`       3
Coefficient in the nonlinear convergence test              :c:func:`MRIStepSetNonlinConvCoef()`       0.1
Nonlinear convergence rate constant                        :c:func:`MRIStepSetNonlinCRDown()`         0.3
Nonlinear residual divergence ratio                        :c:func:`MRIStepSetNonlinRDiv()`           2.3
User-provided implicit stage predictor                     :c:func:`MRIStepSetStagePredictFn()`       ``NULL``
RHS function for nonlinear system evaluations              :c:func:`MRIStepSetNlsRhsFn()`             ``NULL``
Specify if :math:`f^I` is deduced after a nonlinear solve  :c:func:`MRIStepSetDeduceImplicitRhs`      ``SUNFALSE``
=========================================================  =========================================  ============




.. c:function:: int MRIStepSetLinear(void* arkode_mem, int timedepend)

   Specifies that the implicit slow right-hand side function, :math:`f^I(t,y)`
   is linear in :math:`y`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *timedepend* -- flag denoting whether the Jacobian of
        :math:`f^I(t,y)` is time-dependent (1) or not (0).
        Alternately, when using a matrix-free iterative linear solver
        this flag denotes time dependence of the preconditioner.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** Tightens the linear solver tolerances and takes only a
   single Newton iteration.  Calls :c:func:`MRIStepSetDeltaGammaMax()`
   to enforce Jacobian recomputation when the step size ratio changes
   by more than 100 times the unit roundoff (since nonlinear
   convergence is not tested).  Only applicable when used in
   combination with the modified or inexact Newton iteration (not the
   fixed-point solver).

   The only SUNDIALS-provided SUNNonlinearSolver module that is compatible
   with the :c:func:`MRIStepSetLinear()` option is the Newton solver.



.. c:function:: int MRIStepSetNonlinear(void* arkode_mem)

   Specifies that the implicit slow right-hand side function, :math:`f^I(t,y)`
   is nonlinear in :math:`y`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** This is the default behavior of MRIStep, so the function
   is primarily useful to undo a previous call to
   :c:func:`MRIStepSetLinear()`.  Calls
   :c:func:`MRIStepSetDeltaGammaMax()` to reset the step size ratio
   threshold to the default value.



.. c:function:: int MRIStepSetPredictorMethod(void* arkode_mem, int method)

   Specifies the method to use for predicting implicit solutions.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *method* -- method choice (0 :math:`\le` *method* :math:`\le` 4):

        * 0 is the trivial predictor,

        * 1 is the maximum order (dense output) predictor,

        * 2 is the variable order predictor, that decreases the
          polynomial degree for more distant RK stages,

        * 3 is the cutoff order predictor, that uses the maximum order
          for early RK stages, and a first-order predictor for distant
          RK stages,

        * 4 is the bootstrap predictor, that uses a second-order
          predictor based on only information within the current step.
          **deprecated**

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** The default value is 0.  If *method* is set to an
   undefined value, this default predictor will be used.

   **The "bootstrap" predictor (option 4 above) has been deprecated, and
   will be removed from a future release.**



.. c:function:: int MRIStepSetMaxNonlinIters(void* arkode_mem, int maxcor)

   Specifies the maximum number of nonlinear solver
   iterations permitted per slow MRI stage within each time step.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *maxcor* -- maximum allowed solver iterations per stage :math:`(>0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value or if the SUNNONLINSOL module is ``NULL``
      * *ARK_NLS_OP_ERR* if the SUNNONLINSOL object returned a failure flag

   **Notes:** The default value is 3; set *maxcor* :math:`\le 0`
   to specify this default.



.. c:function:: int MRIStepSetNonlinConvCoef(void* arkode_mem, realtype nlscoef)

   Specifies the safety factor used within the nonlinear solver convergence test.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *nlscoef* -- coefficient in nonlinear solver convergence test :math:`(>0.0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** The default value is 0.1; set *nlscoef* :math:`\le 0`
   to specify this default.



.. c:function:: int MRIStepSetNonlinCRDown(void* arkode_mem, realtype crdown)

   Specifies the constant used in estimating the nonlinear solver convergence rate.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *crdown* -- nonlinear convergence rate estimation constant (default is 0.3).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** Any non-positive parameter will imply a reset to the default value.



.. c:function:: int MRIStepSetNonlinRDiv(void* arkode_mem, realtype rdiv)

   Specifies the nonlinear correction threshold beyond which the
   iteration will be declared divergent.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *rdiv* -- tolerance on nonlinear correction size ratio to
        declare divergence (default is 2.3).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** Any non-positive parameter will imply a reset to the default value.



.. c:function:: int MRIStepSetStagePredictFn(void* arkode_mem, ARKStagePredictFn PredictStage)

   Sets the user-supplied function to update the implicit stage predictor prior to
   execution of the nonlinear or linear solver algorithms that compute the implicit stage solution.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *PredictStage* -- name of user-supplied predictor function. If ``NULL``, then any
        previously-provided stage prediction function will be disabled.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   **Notes:** See :numref:`ARKODE.Usage.StagePredictFn` for more information on
   this user-supplied routine.


.. c:function:: int MRIStepSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fs)

   Specifies an alternative implicit slow right-hand side function for
   evaluating :math:`f^I(t,y)` within nonlinear system function evaluations.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *nls_fs* -- the alternative C function for computing the right-hand side
        function :math:`f^I(t,y)` in the ODE.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``.

   **Notes:** The default is to use the implicit slow right-hand side function
   provided to :c:func:`MRIStepCreate()` in nonlinear system functions. If the
   input implicit slow right-hand side function is ``NULL``, the default is
   used.

   When using a non-default nonlinear solver, this function must be called
   *after* :c:func:`MRIStepSetNonlinearSolver()`.


.. c:function:: int MRIStepSetDeduceImplicitRhs(void *arkode_mem, sunbooleantype deduce)

   Specifies if implicit stage derivatives are deduced without evaluating
   :math:`f^I`. See :numref:`ARKODE.Mathematics.Nonlinear` for more details.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *deduce* -- If ``SUNFALSE`` (default), the stage derivative is obtained
        by evaluating :math:`f^I` with the stage solution returned from the
        nonlinear solver. If ``SUNTRUE``, the stage derivative is deduced
        without an additional evaluation of :math:`f^I`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   .. versionadded:: 5.2.0


.. _ARKODE.Usage.MRIStep.ARKLsInputs:

Linear solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation of the linear solver methods
available to MRIStep is provided in :numref:`ARKODE.Mathematics.Linear`.  We
group the user-callable routines into
four categories: general routines concerning the update frequency for
matrices and/or preconditioners, optional inputs for matrix-based
linear solvers, optional inputs for matrix-free linear solvers, and
optional inputs for iterative linear solvers.  We note that the
matrix-based and matrix-free groups are mutually exclusive, whereas the
"iterative" tag can apply to either case.



.. _ARKODE.Usage.MRIStep.ARKLsInputs.General:

.. index::
   single: optional input; generic linear solver interface (MRIStep)

Optional inputs for the ARKLS linear solver interface
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

As discussed in :numref:`ARKODE.Mathematics.Linear.Setup`, ARKODE
strives to reuse matrix and preconditioner data for as many solves as
possible to amortize the high costs of matrix construction and
factorization.  To that end, MRIStep provides user-callable
routines to modify this behavior.  Recall that the
Newton system matrices that arise within an implicit stage solve are
:math:`{\mathcal A}(t,z) \approx I - \gamma J(t,z)`, where the
implicit right-hand side function has Jacobian matrix
:math:`J(t,z) = \frac{\partial f^I(t,z)}{\partial z}`.

The matrix or preconditioner for :math:`{\mathcal A}` can only be
updated within a call to the linear solver 'setup' routine.  In
general, the frequency with which the linear solver setup routine is
called may be controlled with the *msbp* argument to
:c:func:`MRIStepSetLSetupFrequency()`.  When this occurs, the
validity of :math:`{\mathcal A}` for successive time steps
intimately depends on whether the corresponding :math:`\gamma` and
:math:`J` inputs remain valid.

At each call to the linear solver setup routine the decision to update
:math:`\mathcal{A}` with a new value of :math:`\gamma`, and to reuse
or reevaluate Jacobian information, depends on several factors including:

* the success or failure of previous solve attempts,
* the success or failure of the previous time step attempts,
* the change in :math:`\gamma` from the value used when constructing :math:`\mathcal{A}`, and
* the number of steps since Jacobian information was last evaluated.

The frequency with which to update Jacobian information can be controlled
with the *msbj* argument to :c:func:`MRIStepSetJacEvalFrequency()`.
We note that this is only checked *within* calls to the linear solver setup
routine, so values *msbj* :math:`<` *msbp* do not make sense. For
linear-solvers with user-supplied preconditioning the above factors are used
to determine whether to recommend updating the Jacobian information in the
preconditioner (i.e., whether to set *jok* to ``SUNFALSE`` in calling the
user-supplied :c:type:`ARKLsPrecSetupFn()`). For matrix-based linear solvers
these factors determine whether the matrix :math:`J(t,y) = \frac{\partial f^I(t,y)}{\partial y}`
should be updated (either with an internal finite difference approximation or
a call to the user-supplied :c:type:`ARKLsJacFn`); if not then the previous
value is reused and the system matrix :math:`{\mathcal A}(t,y) \approx I - \gamma J(t,y)`
is recomputed using the current :math:`\gamma` value.



.. cssclass:: table-bordered

=============================================  =========================================  ============
Optional input                                 Function name                              Default
=============================================  =========================================  ============
Max change in step signaling new :math:`J`     :c:func:`MRIStepSetDeltaGammaMax()`        0.2
Linear solver setup frequency                  :c:func:`MRIStepSetLSetupFrequency()`      20
Jacobian / preconditioner update frequency     :c:func:`MRIStepSetJacEvalFrequency()`     51
=============================================  =========================================  ============


.. c:function:: int MRIStepSetDeltaGammaMax(void* arkode_mem, realtype dgmax)

   Specifies a scaled step size ratio tolerance, beyond which the
   linear solver setup routine will be signaled.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *dgmax* -- tolerance on step size ratio change before calling
        linear solver setup routine (default is 0.2).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**  Any non-positive parameter will imply a reset to the default value.


.. index::
   single: optional input; linear solver setup frequency (MRIStep)

.. c:function:: int MRIStepSetLSetupFrequency(void* arkode_mem, int msbp)

   Specifies the frequency of calls to the linear solver setup
   routine.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *msbp* -- the linear solver setup frequency.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   **Notes:**
   Positive values of **msbp** specify the linear solver setup frequency. For
   example, an input of 1 means the setup function will be called every time
   step while an input of 2 means it will be called called every other time
   step. If **msbp** is 0, the default value of 20 will be used. A negative
   value forces a linear solver step at each implicit stage.


.. index::
   single: optional input; Jacobian update frequency (MRIStep)
   single: optional input; preconditioner update frequency (MRIStep)

.. c:function:: int MRIStepSetJacEvalFrequency(void* arkode_mem, long int msbj)

   Specifies the frequency for recomputing the Jacobian or recommending a
   preconditioner update.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *msbj* -- the Jacobian re-computation or preconditioner update frequency.

   **Return value:**
      * *ARKLS_SUCCESS* if successful.
      * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``.
      * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``.

   **Notes:**
   The Jacobian update frequency is only checked *within* calls to the linear
   solver setup routine, as such values of *msbj* :math:`<` *msbp* will result
   in recomputing the Jacobian every *msbp* steps. See
   :c:func:`MRIStepSetLSetupFrequency()` for setting the linear solver setup
   frequency *msbp*.

   Passing a value *msbj* :math:`\le 0` indicates to use the
   default value of 50.

   This function must be called *after* the ARKLS system solver interface has
   been initialized through a call to :c:func:`MRIStepSetLinearSolver()`.







.. _ARKODE.Usage.MRIStep.ARKLsInputs.MatrixBased:

Optional inputs for matrix-based ``SUNLinearSolver`` modules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

=========================================  ===========================================  =============
Optional input                             Function name                                Default
=========================================  ===========================================  =============
Jacobian function                          :c:func:`MRIStepSetJacFn()`                  ``DQ``
Linear system function                     :c:func:`MRIStepSetLinSysFn()`               internal
Enable or disable linear solution scaling  :c:func:`MRIStepSetLinearSolutionScaling()`  on
=========================================  ===========================================  =============

When using matrix-based linear solver modules, the ARKLS solver interface needs
a function to compute an approximation to the Jacobian matrix :math:`J(t,y)` or
the linear system :math:`I - \gamma J`. The function to evaluate the Jacobian
must be of type :c:func:`ARKLsJacFn()`. The user can supply a custom Jacobian
function, or if using a dense or banded :math:`J` can use the default internal
difference quotient approximation that comes with the ARKLS interface.  At
present, we do not supply a corresponding routine to approximate Jacobian
entries in sparse matrices :math:`J`. To specify a user-supplied Jacobian
function *jac*, MRIStep provides the function :c:func:`MRIStepSetJacFn()`.
Alternatively, a function of type :c:func:`ARKLsLinSysFn()` can be provided to
evaluate the matrix :math:`I - \gamma J`. By default, ARKLS uses an
internal linear system function leveraging the SUNMATRIX API to form the matrix
:math:`I - \gamma J`. To specify a user-supplied linear system function
*linsys*, MRIStep provides the function :c:func:`MRIStepSetLinSysFn()`. In
either case the matrix information will be updated infrequently to reduce matrix
construction and, with direct solvers, factorization costs. As a result the
value of :math:`\gamma` may not be current and a scaling factor is applied to the
solution of the linear system to account for lagged value of :math:`\gamma`. See
:numref:`SUNLinSol.Lagged_matrix` for more details. The function
:c:func:`MRIStepSetLinearSolutionScaling()` can be used to disable this scaling
when necessary, e.g., when providing a custom linear solver that updates the
matrix using the current :math:`\gamma` as part of the solve.

The ARKLS interface passes the user data pointer to the Jacobian and linear
system functions. This allows the user to create an arbitrary structure with
relevant problem data and access it during the execution of the user-supplied
Jacobian or linear system functions, without using global data in the
program. The user data pointer may be specified through
:c:func:`MRIStepSetUserData()`.



.. c:function:: int MRIStepSetJacFn(void* arkode_mem, ARKLsJacFn jac)

   Specifies the Jacobian approximation routine to
   be used for the matrix-based solver with the ARKLS interface.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *jac* -- name of user-supplied Jacobian approximation function.

   **Return value:**
      * *ARKLS_SUCCESS*  if successful
      * *ARKLS_MEM_NULL*  if the MRIStep memory was ``NULL``
      * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** This routine must be called after the ARKLS linear
   solver interface has been initialized through a call to
   :c:func:`MRIStepSetLinearSolver()`.

   By default, ARKLS uses an internal difference quotient function for
   dense and band matrices.  If ``NULL`` is passed in for *jac*, this
   default is used. An error will occur if no *jac* is supplied when
   using other matrix types.

   The function type :c:func:`ARKLsJacFn()` is described in
   :numref:`ARKODE.Usage.UserSupplied`.


.. c:function:: int MRIStepSetLinSysFn(void* arkode_mem, ARKLsLinSysFn linsys)

   Specifies the linear system approximation routine to be used for the
   matrix-based solver with the ARKLS interface.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *linsys* -- name of user-supplied linear system approximation function.

   **Return value:**
      * *ARKLS_SUCCESS*  if successful
      * *ARKLS_MEM_NULL*  if the MRIStep memory was ``NULL``
      * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** This routine must be called after the ARKLS linear
   solver interface has been initialized through a call to
   :c:func:`MRIStepSetLinearSolver()`.

   By default, ARKLS uses an internal linear system function that leverages the
   SUNMATRIX API to form the system :math:`I - \gamma J`.  If ``NULL`` is passed
   in for *linsys*, this default is used.

   The function type :c:func:`ARKLsLinSysFn()` is described in
   :numref:`ARKODE.Usage.UserSupplied`.


.. c:function:: int MRIStepSetLinearSolutionScaling(void* arkode_mem, booleantype onoff)

   Enables or disables scaling the linear system solution to account for a
   change in :math:`\gamma` in the linear system. For more details see
   :numref:`SUNLinSol.Lagged_matrix`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *onoff* -- flag to enable (``SUNTRUE``) or disable (``SUNFALSE``)
        scaling

   **Return value:**
      * *ARKLS_SUCCESS* if successful
      * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``
      * *ARKLS_ILL_INPUT* if the attached linear solver is not matrix-based

   **Notes:** Linear solution scaling is enabled by default when a matrix-based
   linear solver is attached.


.. _ARKODE.Usage.MRIStep.ARKLsInputs.MatrixFree:

Optional inputs for matrix-free ``SUNLinearSolver`` modules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==================================================  =========================================  ==================
Optional input                                      Function name                              Default
==================================================  =========================================  ==================
:math:`Jv` functions (*jtimes* and *jtsetup*)       :c:func:`MRIStepSetJacTimes()`             DQ,  none
:math:`Jv` DQ rhs function (*jtimesRhsFn*)          :c:func:`MRIStepSetJacTimesRhsFn()`        fs
==================================================  =========================================  ==================


As described in :numref:`ARKODE.Mathematics.Linear`, when solving
the Newton linear systems with matrix-free methods, the ARKLS
interface requires a *jtimes* function to compute an approximation to
the product between the Jacobian matrix
:math:`J(t,y)` and a vector :math:`v`. The user can supply a custom
Jacobian-times-vector approximation function, or use the default
internal difference quotient function that comes with the ARKLS
interface.

A user-defined Jacobian-vector function must be of type
:c:type:`ARKLsJacTimesVecFn` and can be specified through a call
to :c:func:`MRIStepSetJacTimes()` (see
:numref:`ARKODE.Usage.UserSupplied` for specification details).  As with the
user-supplied preconditioner functions, the evaluation and
processing of any Jacobian-related data needed by the user's
Jacobian-times-vector function is done in the optional user-supplied
function of type :c:type:`ARKLsJacTimesSetupFn` (see
:numref:`ARKODE.Usage.UserSupplied` for specification details).  As with
the preconditioner functions, a pointer to the user-defined
data structure, *user_data*, specified through
:c:func:`MRIStepSetUserData()` (or a ``NULL`` pointer otherwise) is
passed to the Jacobian-times-vector setup and product functions each
time they are called.


.. c:function:: int MRIStepSetJacTimes(void* arkode_mem, ARKLsJacTimesSetupFn jtsetup, ARKLsJacTimesVecFn jtimes)

   Specifies the Jacobian-times-vector setup and product functions.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *jtsetup* -- user-defined Jacobian-vector setup function.
        Pass ``NULL`` if no setup is necessary.
      * *jtimes* -- user-defined Jacobian-vector product function.

   **Return value:**
      * *ARKLS_SUCCESS* if successful.
      * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``.
      * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``.
      * *ARKLS_ILL_INPUT* if an input has an illegal value.
      * *ARKLS_SUNLS_FAIL* if an error occurred when setting up
        the Jacobian-vector product in the ``SUNLinearSolver``
        object used by the ARKLS interface.

   **Notes:** The default is to use an internal finite difference
   quotient for *jtimes* and to leave out *jtsetup*.  If ``NULL`` is
   passed to *jtimes*, these defaults are used.  A user may
   specify non-``NULL`` *jtimes* and ``NULL`` *jtsetup* inputs.

   This function must be called *after* the ARKLS system solver
   interface has been initialized through a call to
   :c:func:`MRIStepSetLinearSolver()`.

   The function types :c:type:`ARKLsJacTimesSetupFn` and
   :c:type:`ARKLsJacTimesVecFn` are described in
   :numref:`ARKODE.Usage.UserSupplied`.


When using the internal difference quotient the user may optionally supply
an alternative implicit right-hand side function for use in the Jacobian-vector
product approximation by calling :c:func:`MRIStepSetJacTimesRhsFn()`. The
alternative implicit right-hand side function should compute a suitable (and
differentiable) approximation to the :math:`f^I` function provided to
:c:func:`MRIStepCreate()`. For example, as done in :cite:p:`dorr2010numerical`, the alternative
function may use lagged values when evaluating a nonlinearity in :math:`f^I` to
avoid differencing a potentially non-differentiable factor.


.. c:function:: int MRIStepSetJacTimesRhsFn(void* arkode_mem, ARKRhsFn jtimesRhsFn)

   Specifies an alternative implicit right-hand side function for use in the
   internal Jacobian-vector product difference quotient approximation.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *jtimesRhsFn* -- the name of the C function (of type
        :c:func:`ARKRhsFn()`) defining the alternative right-hand side function.

   **Return value:**
      * *ARKLS_SUCCESS* if successful.
      * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``.
      * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``.
      * *ARKLS_ILL_INPUT* if an input has an illegal value.

   **Notes:** The default is to use the implicit right-hand side function
   provided to :c:func:`MRIStepCreate()` in the internal difference quotient. If
   the input implicit right-hand side function is ``NULL``, the default is used.

   This function must be called *after* the ARKLS system solver interface has
   been initialized through a call to :c:func:`MRIStepSetLinearSolver()`.





.. _ARKODE.Usage.MRIStep.ARKLsInputs.Iterative:

Optional inputs for iterative ``SUNLinearSolver`` modules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

===============================================  =========================================  ==================
Optional input                                   Function name                              Default
===============================================  =========================================  ==================
Newton preconditioning functions                 :c:func:`MRIStepSetPreconditioner()`       ``NULL``, ``NULL``
Newton linear and nonlinear tolerance ratio      :c:func:`MRIStepSetEpsLin()`               0.05
Newton linear solve tolerance conversion factor  :c:func:`MRIStepSetLSNormFactor()`         vector length
===============================================  =========================================  ==================


As described in :numref:`ARKODE.Mathematics.Linear`, when using
an iterative linear solver the user may supply a preconditioning
operator to aid in solution of the system.  This operator consists of
two user-supplied functions, *psetup* and *psolve*, that are supplied
to MRIStep using the function :c:func:`MRIStepSetPreconditioner()`.
The *psetup* function supplied to these routines
should handle evaluation and preprocessing of any Jacobian data
needed by the user's preconditioner solve function,
*psolve*.  The user data pointer received through
:c:func:`MRIStepSetUserData()` (or a pointer to ``NULL`` if user data
was not specified) is passed to the *psetup* and *psolve* functions.
This allows the user to create an arbitrary
structure with relevant problem data and access it during the
execution of the user-supplied preconditioner functions without using
global data in the program.

Also, as described in :numref:`ARKODE.Mathematics.Error.Linear`, the
ARKLS interface requires that iterative linear solvers stop when
the norm of the preconditioned residual satisfies

.. math::
   \|r\| \le \frac{\epsilon_L \epsilon}{10}

where the default :math:`\epsilon_L = 0.05`, which may be modified by
the user through the :c:func:`MRIStepSetEpsLin()` function.


.. c:function:: int MRIStepSetPreconditioner(void* arkode_mem, ARKLsPrecSetupFn psetup, ARKLsPrecSolveFn psolve)

   Specifies the user-supplied preconditioner setup and solve functions.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *psetup* -- user defined preconditioner setup function.  Pass
        ``NULL`` if no setup is needed.
      * *psolve* -- user-defined preconditioner solve function.

   **Return value:**
      * *ARKLS_SUCCESS* if successful.
      * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``.
      * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``.
      * *ARKLS_ILL_INPUT* if an input has an illegal value.
      * *ARKLS_SUNLS_FAIL* if an error occurred when setting up
        preconditioning in the ``SUNLinearSolver`` object used
        by the ARKLS interface.

   **Notes:** The default is ``NULL`` for both arguments (i.e., no
   preconditioning).

   This function must be called *after* the ARKLS system solver
   interface has been initialized through a call to
   :c:func:`MRIStepSetLinearSolver()`.

   Both of the function types :c:func:`ARKLsPrecSetupFn()` and
   :c:func:`ARKLsPrecSolveFn()` are described in
   :numref:`ARKODE.Usage.UserSupplied`.


.. c:function:: int MRIStepSetEpsLin(void* arkode_mem, realtype eplifac)

   Specifies the factor by which the tolerance on the nonlinear
   iteration is multiplied to get a tolerance on the linear
   iteration.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *eplifac* -- linear convergence safety factor.

   **Return value:**
      * *ARKLS_SUCCESS* if successful.
      * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``.
      * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``.
      * *ARKLS_ILL_INPUT* if an input has an illegal value.

   **Notes:** Passing a value *eplifac* :math:`\le 0` indicates to use the
   default value of 0.05.

   This function must be called *after* the ARKLS system solver
   interface has been initialized through a call to
   :c:func:`MRIStepSetLinearSolver()`.


.. c:function:: int MRIStepSetLSNormFactor(void* arkode_mem, realtype nrmfac)

   Specifies the factor to use when converting from the integrator tolerance
   (WRMS norm) to the linear solver tolerance (L2 norm) for Newton linear system
   solves e.g., ``tol_L2 = fac * tol_WRMS``.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *nrmfac* -- the norm conversion factor. If *nrmfac* is:

        :math:`> 0` then the provided value is used.

        :math:`= 0` then the conversion factor is computed using the vector
        length i.e., ``nrmfac = sqrt(N_VGetLength(y))`` (*default*).

        :math:`< 0` then the conversion factor is computed using the vector dot
        product i.e., ``nrmfac = sqrt(N_VDotProd(v,v))`` where all the entries
        of ``v`` are one.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``.

   **Notes:**
   This function must be called *after* the ARKLS system solver interface has
   been initialized through a call to :c:func:`MRIStepSetLinearSolver()`.



.. _ARKODE.Usage.MRIStep.MRIStepRootfindingInput:

Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions can be called to set optional inputs to
control the rootfinding algorithm, the mathematics of which are
described in the section :numref:`ARKODE.Mathematics.Rootfinding`.


.. cssclass:: table-bordered

======================================  ========================================  ==================
Optional input                          Function name                             Default
======================================  ========================================  ==================
Direction of zero-crossings to monitor  :c:func:`MRIStepSetRootDirection()`       both
Disable inactive root warnings          :c:func:`MRIStepSetNoInactiveRootWarn()`  enabled
======================================  ========================================  ==================



.. c:function:: int MRIStepSetRootDirection(void* arkode_mem, int* rootdir)

   Specifies the direction of zero-crossings to be located and returned.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *rootdir* -- state array of length *nrtfn*, the number of root
        functions :math:`g_i`  (the value of *nrtfn* was supplied in
        the call to :c:func:`MRIStepRootInit()`).  If ``rootdir[i] ==
        0`` then crossing in either direction for :math:`g_i` should be
        reported.  A value of +1 or -1 indicates that the solver
        should report only zero-crossings where :math:`g_i` is
        increasing or decreasing, respectively.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**  The default behavior is to monitor for both zero-crossing directions.



.. c:function:: int MRIStepSetNoInactiveRootWarn(void* arkode_mem)

   Disables issuing a warning if some root function appears
   to be identically zero at the beginning of the integration.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   **Notes:** MRIStep will not report the initial conditions as a
   possible zero-crossing (assuming that one or more components
   :math:`g_i` are zero at the initial time).  However, if it appears
   that some :math:`g_i` is identically zero at the initial time
   (i.e., :math:`g_i` is zero at the initial time *and* after the
   first step), MRIStep will issue a warning which can be disabled with
   this optional input function.



.. _ARKODE.Usage.MRIStep.InterpolatedOutput:

Interpolated output function
--------------------------------

An optional function :c:func:`MRIStepGetDky()` is available to obtain
additional values of solution-related quantities.  This function
should only be called after a successful return from
:c:func:`MRIStepEvolve()`, as it provides interpolated values either of
:math:`y` or of its derivatives (up to the 3rd derivative)
interpolated to any value of :math:`t` in the last internal step taken
by :c:func:`MRIStepEvolve()`.  Internally, this "dense output" or
"continuous extension" algorithm is identical to the algorithm used for
the maximum order implicit predictors, described in
:numref:`ARKODE.Mathematics.Predictors.Max`, except that derivatives of the
polynomial model may be evaluated upon request.



.. c:function:: int MRIStepGetDky(void* arkode_mem, realtype t, int k, N_Vector dky)

   Computes the *k*-th derivative of the function
   :math:`y` at the time *t*,
   i.e. :math:`y^{(k)}(t)`, for values of the
   independent variable satisfying :math:`t_n-h_n \le t \le t_n`, with
   :math:`t_n` as current internal time reached, and :math:`h_n` is
   the last internal step size successfully used by the solver.  This
   routine uses an interpolating polynomial of degree *min(degree, 5)*,
   where *degree* is the argument provided to
   :c:func:`MRIStepSetInterpolantDegree()`.  The user may request *k* in the
   range {0,..., *min(degree, kmax)*} where *kmax* depends on the choice of
   interpolation module. For Hermite interpolants *kmax = 5* and for Lagrange
   interpolants *kmax = 3*.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *t* -- the value of the independent variable at which the
     derivative is to be evaluated.

   * *k* -- the derivative order requested.

   * *dky* -- output vector (must be allocated by the user).

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_BAD_K* if *k* is not in the range {0,..., *min(degree, kmax)*}.

   * *ARK_BAD_T* if *t* is not in the interval :math:`[t_n-h_n, t_n]`

   * *ARK_BAD_DKY* if the *dky* vector was ``NULL``

   * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``

   **Notes:** It is only legal to call this function after a successful
   return from :c:func:`MRIStepEvolve()`.

   A user may access the values :math:`t_n` and :math:`h_n` via the
   functions :c:func:`MRIStepGetCurrentTime()` and
   :c:func:`MRIStepGetLastStep()`, respectively.



.. _ARKODE.Usage.MRIStep.OptionalOutputs:

Optional output functions
------------------------------

MRIStep provides an extensive set of functions that can be used to
obtain solver performance information.  We organize these into groups:

#. General MRIStep output routines are in
   :numref:`ARKODE.Usage.MRIStep.MRIStepMainOutputs`,

#. MRIStep implicit solver output routines are in
   :numref:`ARKODE.Usage.MRIStep.MRIStepImplicitSolverOutputs`,

#. Linear solver output routines are in
   :numref:`ARKODE.Usage.MRIStep.ARKLsOutputs` and

#. General usability routines (e.g. to print the current MRIStep
   parameters, or output the current coupling table) are in
   :numref:`ARKODE.Usage.MRIStep.MRIStepExtraOutputs`.

#. Output routines regarding root-finding results are in
   :numref:`ARKODE.Usage.MRIStep.MRIStepRootOutputs`,

Following each table, we elaborate on each function.

Some of the optional outputs, especially the various counters, can be
very useful in determining the efficiency of various methods inside
MRIStep.  For example:

* The number of steps and right-hand side evaluations at both the slow and fast
  time scales provide a rough measure of the overall cost of a given run, and can
  be compared between runs with different solver options to suggest which set of
  options is the most efficient.

* The ratio *nniters/nsteps* measures the performance of the
  nonlinear iteration in solving the nonlinear systems at each implicit stage,
  providing a measure of the degree of nonlinearity in the problem.
  Typical values of this for a Newton solver on a general problem
  range from 1.1 to 1.8.

* When using a Newton nonlinear solver, the ratio *njevals/nniters*
  (when using a direct linear solver), and the ratio
  *nliters/nniters* (when using an iterative linear solver) can
  indicate the quality of the approximate Jacobian or preconditioner being
  used.  For example, if this ratio is larger for a user-supplied
  Jacobian or Jacobian-vector product routine than for the
  difference-quotient routine, it can indicate that the user-supplied
  Jacobian is inaccurate.

It is therefore recommended that users retrieve and output these
statistics following each run, and take some time to investigate
alternate solver options that will be more optimal for their
particular problem of interest.



.. _ARKODE.Usage.MRIStep.MRIStepMainOutputs:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. _ARKODE.Usage.MRIStep.MRIStepMainOutputs.Table:
.. table:: Main solver optional output functions

   +------------------------------------------------------+-------------------------------------------+
   | Optional output                                      | Function name                             |
   +------------------------------------------------------+-------------------------------------------+
   | Size of MRIStep real and integer workspaces          | :c:func:`MRIStepGetWorkSpace()`           |
   +------------------------------------------------------+-------------------------------------------+
   | Cumulative number of internal steps                  | :c:func:`MRIStepGetNumSteps()`            |
   +------------------------------------------------------+-------------------------------------------+
   | Step size used for the last successful step          | :c:func:`MRIStepGetLastStep()`            |
   +------------------------------------------------------+-------------------------------------------+
   | Current internal time reached by the solver          | :c:func:`MRIStepGetCurrentTime()`         |
   +------------------------------------------------------+-------------------------------------------+
   | Current internal solution reached by the solver      | :c:func:`MRIStepGetCurrentState()`        |
   +------------------------------------------------------+-------------------------------------------+
   | Current :math:`\gamma` value used by the solver      | :c:func:`MRIStepGetCurrentGamma()`        |
   +------------------------------------------------------+-------------------------------------------+
   | Error weight vector for state variables              | :c:func:`MRIStepGetErrWeights()`          |
   +------------------------------------------------------+-------------------------------------------+
   | Suggested factor for tolerance scaling               | :c:func:`MRIStepGetTolScaleFactor()`      |
   +------------------------------------------------------+-------------------------------------------+
   | Print all statistics                                 | :c:func:`MRIStepPrintAllStats`            |
   +------------------------------------------------------+-------------------------------------------+
   | Name of constant associated with a return flag       | :c:func:`MRIStepGetReturnFlagName()`      |
   +------------------------------------------------------+-------------------------------------------+
   | No. of calls to the :math:`f^E` and :math:`f^I`      | :c:func:`MRIStepGetNumRhsEvals()`         |
   +------------------------------------------------------+-------------------------------------------+
   | No. of failed steps due to a nonlinear solver        | :c:func:`MRIStepGetNumStepSolveFails()`   |
   | failure                                              |                                           |
   +------------------------------------------------------+-------------------------------------------+
   | Current MRI coupling tables                          | :c:func:`MRIStepGetCurrentCoupling()`     |
   +------------------------------------------------------+-------------------------------------------+
   | Last inner stepper return value                      | :c:func:`MRIStepGetLastInnerStepFlag()`   |
   +------------------------------------------------------+-------------------------------------------+
   | Retrieve a pointer for user data                     |  :c:func:`MRIStepGetUserData`             |
   +------------------------------------------------------+-------------------------------------------+

.. Functions not currently provided by MRIStep
.. No. of explicit stability-limited steps              :c:func:`MRIStepGetNumExpSteps()`
.. No. of accuracy-limited steps                        :c:func:`MRIStepGetNumAccSteps()`
.. No. of attempted steps                               :c:func:`MRIStepGetNumStepAttempts()`
.. No. of local error test failures that have occurred  :c:func:`MRIStepGetNumErrTestFails()`
.. Estimated local truncation error vector              :c:func:`MRIStepGetEstLocalErrors()`
.. Single accessor to many statistics at once           :c:func:`MRIStepGetTimestepperStats()`
.. Actual initial time step size used                   :c:func:`MRIStepGetActualInitStep()`
.. Step size to be attempted on the next step           :c:func:`MRIStepGetCurrentStep()`
.. Single accessor to many statistics at once           :c:func:`MRIStepGetStepStats()`


.. c:function:: int MRIStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)

   Returns the MRIStep real and integer workspace sizes.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *lenrw* -- the number of ``realtype`` values in the MRIStep workspace.

   * *leniw* -- the number of integer values in the MRIStep workspace.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


.. c:function:: int MRIStepGetNumSteps(void* arkode_mem, long int* nssteps, long int* nfsteps)

   Returns the cumulative number of slow and fast internal steps taken by
   the solver (so far).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *nssteps* -- number of slow steps taken in the solver.

   * *nfsteps* -- number of fast steps taken in the solver.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


..
   .. c:function:: int MRIStepGetActualInitStep(void* arkode_mem, realtype* hinused)

      Returns the value of the integration step size used on the first step.

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *hinused* -- actual value of initial step size.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


.. c:function:: int MRIStepGetLastStep(void* arkode_mem, realtype* hlast)

   Returns the integration step size taken on the last successful
   internal step.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *hlast* -- step size taken on the last internal step.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


..
   .. c:function:: int MRIStepGetCurrentStep(void* arkode_mem, realtype* hcur)

      Returns the integration step size to be attempted on the next internal step.

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *hcur* -- step size to be attempted on the next internal step.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


.. c:function:: int MRIStepGetCurrentTime(void* arkode_mem, realtype* tcur)

   Returns the current internal time reached by the solver.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *tcur* -- current internal time reached.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


.. c:function:: int MRIStepGetCurrentState(void *arkode_mem, N_Vector *ycur)

   Returns the current internal solution reached by the solver.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *ycur* -- current internal solution.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   **Notes:** Users should exercise extreme caution when using this function,
   as altering values of *ycur* may lead to undesirable behavior, depending
   on the particular use case and on when this routine is called.


.. c:function:: int MRIStepGetCurrentGamma(void *arkode_mem, realtype *gamma)

   Returns the current internal value of :math:`\gamma` used in the implicit
   solver Newton matrix (see equation :eq:`ARKODE_NewtonMatrix`).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *gamma* -- current step size scaling factor in the Newton system.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``



.. c:function:: int MRIStepGetTolScaleFactor(void* arkode_mem, realtype* tolsfac)

   Returns a suggested factor by which the user's
   tolerances should be scaled when too much accuracy has been
   requested for some internal step.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *tolsfac* -- suggested scaling factor for user-supplied tolerances.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


.. c:function:: int MRIStepGetErrWeights(void* arkode_mem, N_Vector eweight)

   Returns the current error weight vector.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *eweight* -- solution error weights at the current time.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   **Notes:** The user must allocate space for *eweight*, that will be
   filled in by this function.


..
   .. c:function:: int MRIStepGetStepStats(void* arkode_mem, long int* nssteps, long int* nfsteps, realtype* hlast, realtype* tcur)

      Returns many of the most useful optional outputs in a single call.

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *nssteps* -- number of slow steps taken in the solver.

      * *nfsteps* -- number of fast steps taken in the solver.

      * *hlast* -- step size taken on the last internal step.

      * *tcur* -- current internal time reached.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


.. c:function:: int MRIStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)

   Outputs all of the integrator, nonlinear solver, linear solver, and other
   statistics.

   **Arguments:**
     * *arkode_mem* -- pointer to the MRIStep memory block.
     * *outfile* -- pointer to output file.
     * *fmt* -- the output format:

       * :c:enumerator:`SUN_OUTPUTFORMAT_TABLE` -- prints a table of values
       * :c:enumerator:`SUN_OUTPUTFORMAT_CSV` -- prints a comma-separated list
         of key and value pairs e.g., ``key1,value1,key2,value2,...``

   **Return value:**
     * *ARK_SUCCESS* -- if the output was successfully.
     * *CV_MEM_NULL* -- if the MRIStep memory was ``NULL``.
     * *CV_ILL_INPUT* -- if an invalid formatting option was provided.

   .. note::

      The file ``scripts/sundials_csv.py`` provides python utility functions to
      read and output the data from a SUNDIALS CSV output file using the key
      and value pair format.

   .. versionadded:: 5.2.0


.. c:function:: char *MRIStepGetReturnFlagName(long int flag)

   Returns the name of the MRIStep constant corresponding to *flag*.
   See :ref:`ARKODE.Constants`.

   **Arguments:**

   * *flag* -- a return flag from an MRIStep function.

   **Return value:**
   The return value is a string containing the name of
   the corresponding constant.


..
   .. c:function:: int MRIStepGetNumExpSteps(void* arkode_mem, long int* expsteps)

      Returns the cumulative number of stability-limited steps
      taken by the solver (so far).

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *expsteps* -- number of stability-limited steps taken in the solver.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


..
   .. c:function:: int MRIStepGetNumAccSteps(void* arkode_mem, long int* accsteps)

      Returns the cumulative number of accuracy-limited steps
      taken by the solver (so far).

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *accsteps* -- number of accuracy-limited steps taken in the solver.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


..
   .. c:function:: int MRIStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts)

      Returns the cumulative number of steps attempted by the solver (so far).

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *step_attempts* -- number of steps attempted by solver.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


.. c:function:: int MRIStepGetNumRhsEvals(void* arkode_mem, long int* nfse_evals, long int* nfsi_evals)

   Returns the number of calls to the user's outer (slow) right-hand side
   functions, :math:`f^E` and :math:`f^I`, so far.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *nfse_evals* -- number of calls to the user's :math:`f^E(t,y)` function.

   * *nfsi_evals* -- number of calls to the user's :math:`f^I(t,y)` function.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


..
   .. c:function:: int MRIStepGetNumErrTestFails(void* arkode_mem, long int* netfails)

      Returns the number of local error test failures that
      have occurred (so far).

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *netfails* -- number of error test failures.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


.. c:function:: int MRIStepGetNumStepSolveFails(void* arkode_mem, long int* ncnf)

   Returns the number of failed steps due to a nonlinear solver failure (so far).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *ncnf* -- number of step failures.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``


.. c:function:: int MRIStepGetCurrentCoupling(void* arkode_mem, MRIStepCoupling *C)

   Returns the MRI coupling table currently in use by the solver.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *C* -- pointer to slow-to-fast MRI coupling structure.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   **Notes:**  The *MRIStepCoupling* data structure is defined in
   the header file ``arkode/arkode_mristep.h``.  It is defined as a
   pointer to the following C structure:

   .. code-block:: c

      struct MRIStepCouplingMem {

         int nmat;        /* number of MRI coupling matrices             */
         int stages;      /* size of coupling matrices (stages * stages) */
         int q;           /* method order of accuracy                    */
         int p;           /* embedding order of accuracy                 */
         realtype ***G;   /* coupling matrices [nmat][stages][stages]    */
         realtype *c;     /* abscissae                                   */

       };
       typedef MRIStepCouplingMem *MRIStepCoupling;

   For more details see :numref:`ARKODE.Usage.MRIStep.MRIStepCoupling`.


..
   .. c:function:: int MRIStepGetEstLocalErrors(void* arkode_mem, N_Vector ele)

      Returns the vector of estimated local truncation errors
      for the current step.

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *ele* -- vector of estimated local truncation errors.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

      **Notes:**  The user must allocate space for *ele*, that will be
      filled in by this function.

      The values returned in *ele* are valid only after a successful call
      to :c:func:`MRIStepEvolve()` (i.e., it returned a non-negative value).

      The *ele* vector, together with the *eweight* vector from
      :c:func:`MRIStepGetErrWeights()`, can be used to determine how the
      various components of the system contributed to the estimated local
      error test.  Specifically, that error test uses the WRMS norm of a
      vector whose components are the products of the components of these
      two vectors.  Thus, for example, if there were recent error test
      failures, the components causing the failures are those with largest
      values for the products, denoted loosely as ``eweight[i]*ele[i]``.


..
   .. c:function:: int MRIStepGetTimestepperStats(void* arkode_mem, long int* expsteps, long int* accsteps, long int* step_attempts, long int* nf_evals, long int* netfails)

      Returns many of the most useful time-stepper statistics in a single call.

      **Arguments:**

      * *arkode_mem* -- pointer to the MRIStep memory block.

      * *expsteps* -- number of stability-limited steps taken in the solver.

      * *accsteps* -- number of accuracy-limited steps taken in the solver.

      * *step_attempts* -- number of steps attempted by the solver.

      * *nf_evals* -- number of calls to the user's :math:`f(t,y)` function.

      * *netfails* -- number of error test failures.

      **Return value:**

      * *ARK_SUCCESS* if successful

      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

.. c:function:: int MRIStepGetLastInnerStepFlag(void* arkode_mem, int* flag)

   Returns the last return value from the inner stepper.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *flag* -- inner stepper return value.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

.. c:function:: int MRIStepGetUserData(void* arkode_mem, void** user_data)

   Returns the user data pointer previously set with
   :c:func:`MRIStepSetUserData`.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *user_data* -- memory reference to a user data pointer

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the ARKStep memory was ``NULL``

   .. versionadded:: 5.3.0


.. _ARKODE.Usage.MRIStep.MRIStepImplicitSolverOutputs:

Implicit solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.MRIStep.MRIStepImplicitSolverOutputs.Table:
.. table:: Implicit solver optional output functions

   +------------------------------------------------------+----------------------------------------------+
   | Optional output                                      | Function name                                |
   +------------------------------------------------------+----------------------------------------------+
   | No. of calls to linear solver setup function         | :c:func:`MRIStepGetNumLinSolvSetups()`       |
   +------------------------------------------------------+----------------------------------------------+
   | No. of nonlinear solver iterations                   | :c:func:`MRIStepGetNumNonlinSolvIters()`     |
   +------------------------------------------------------+----------------------------------------------+
   | No. of nonlinear solver convergence failures         | :c:func:`MRIStepGetNumNonlinSolvConvFails()` |
   +------------------------------------------------------+----------------------------------------------+
   | Single accessor to all nonlinear solver statistics   | :c:func:`MRIStepGetNonlinSolvStats()`        |
   +------------------------------------------------------+----------------------------------------------+




.. c:function:: int MRIStepGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups)

   Returns the number of calls made to the linear solver's
   setup routine (so far).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *nlinsetups* -- number of linear solver setup calls made.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   **Notes:** This is only accumulated for the "life" of the nonlinear
   solver object; the counter is reset whenever a new nonlinear solver
   module is "attached" to MRIStep, or when MRIStep is resized.


.. c:function:: int MRIStepGetNumNonlinSolvIters(void* arkode_mem, long int* nniters)

   Returns the number of nonlinear solver iterations performed (so far).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *nniters* -- number of nonlinear iterations performed.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARK_NLS_OP_ERR* if the SUNNONLINSOL object returned a failure flag

   **Notes:** This is only accumulated for the "life" of the nonlinear
   solver object; the counter is reset whenever a new nonlinear solver
   module is "attached" to MRIStep, or when MRIStep is resized.



.. c:function:: int MRIStepGetNumNonlinSolvConvFails(void* arkode_mem, long int* nncfails)

   Returns the number of nonlinear solver convergence
   failures that have occurred (so far).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *nncfails* -- number of nonlinear convergence failures.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   **Notes:** This is only accumulated for the "life" of the nonlinear
   solver object; the counter is reset whenever a new nonlinear solver
   module is "attached" to MRIStep, or when MRIStep is resized.



.. c:function:: int MRIStepGetNonlinSolvStats(void* arkode_mem, long int* nniters, long int* nncfails)

   Returns all of the nonlinear solver statistics in a single call.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *nniters* -- number of nonlinear iterations performed.

   * *nncfails* -- number of nonlinear convergence failures.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARK_NLS_OP_ERR* if the SUNNONLINSOL object returned a failure flag

   **Notes:** These are only accumulated for the "life" of the
   nonlinear solver object; the counters are reset whenever a new
   nonlinear solver module is "attached" to MRIStep, or when MRIStep is resized.



.. _ARKODE.Usage.MRIStep.MRIStepRootOutputs:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===================================================  ==========================================
Optional output                                      Function name
===================================================  ==========================================
Array showing roots found                            :c:func:`MRIStepGetRootInfo()`
No. of calls to user root function                   :c:func:`MRIStepGetNumGEvals()`
===================================================  ==========================================



.. c:function:: int MRIStepGetRootInfo(void* arkode_mem, int* rootsfound)

   Returns an array showing which functions were found to
   have a root.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *rootsfound* -- array of length *nrtfn* with the indices of the
        user functions :math:`g_i` found to have a root (the value of
        *nrtfn* was supplied in the call to
        :c:func:`MRIStepRootInit()`).  For :math:`i = 0 \ldots`
        *nrtfn*-1, ``rootsfound[i]`` is nonzero if :math:`g_i` has a
        root, and 0 if not.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   **Notes:** The user must allocate space for *rootsfound* prior to
   calling this function.

   For the components of :math:`g_i` for which a root was found, the
   sign of ``rootsfound[i]`` indicates the direction of
   zero-crossing.  A value of +1 indicates that :math:`g_i` is
   increasing, while a value of -1 indicates a decreasing :math:`g_i`.



.. c:function:: int MRIStepGetNumGEvals(void* arkode_mem, long int* ngevals)

   Returns the cumulative number of calls made to the
   user's root function :math:`g`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *ngevals* -- number of calls made to :math:`g` so far.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``



.. _ARKODE.Usage.MRIStep.ARKLsOutputs:

Linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A variety of optional outputs are available from the ARKLS interface, as
listed in the following table and elaborated below.  We note that where the
name of an output would otherwise conflict with the
name of an optional output from the main solver, a suffix LS (for
Linear Solver) has been added here (e.g. *lenrwLS*).


.. _ARKODE.Usage.MRIStep.ARKLsOutputs.Table:
.. table:: Linear solver interface optional output functions

   +--------------------------------------------------------------------+------------------------------------------+
   | Optional output                                                    | Function name                            |
   +--------------------------------------------------------------------+------------------------------------------+
   | Stored Jacobian of the ODE RHS function                            | :c:func:`MRIStepGetJac`                  |
   +--------------------------------------------------------------------+------------------------------------------+
   | Time at which the Jacobian was evaluated                           | :c:func:`MRIStepGetJacTime`              |
   +--------------------------------------------------------------------+------------------------------------------+
   | Step number at which the Jacobian was evaluated                    | :c:func:`MRIStepGetJacNumSteps`          |
   +--------------------------------------------------------------------+------------------------------------------+
   | Size of real and integer workspaces                                | :c:func:`MRIStepGetLinWorkSpace()`       |
   +--------------------------------------------------------------------+------------------------------------------+
   | No. of Jacobian evaluations                                        | :c:func:`MRIStepGetNumJacEvals()`        |
   +--------------------------------------------------------------------+------------------------------------------+
   | No. of preconditioner evaluations                                  | :c:func:`MRIStepGetNumPrecEvals()`       |
   +--------------------------------------------------------------------+------------------------------------------+
   | No. of preconditioner solves                                       | :c:func:`MRIStepGetNumPrecSolves()`      |
   +--------------------------------------------------------------------+------------------------------------------+
   | No. of linear iterations                                           | :c:func:`MRIStepGetNumLinIters()`        |
   +--------------------------------------------------------------------+------------------------------------------+
   | No. of linear convergence failures                                 | :c:func:`MRIStepGetNumLinConvFails()`    |
   +--------------------------------------------------------------------+------------------------------------------+
   | No. of Jacobian-vector setup evaluations                           | :c:func:`MRIStepGetNumJTSetupEvals()`    |
   +--------------------------------------------------------------------+------------------------------------------+
   | No. of Jacobian-vector product evaluations                         | :c:func:`MRIStepGetNumJtimesEvals()`     |
   +--------------------------------------------------------------------+------------------------------------------+
   | No. of *fs* calls for finite diff. :math:`J` or :math:`Jv` evals.  | :c:func:`MRIStepGetNumLinRhsEvals()`     |
   +--------------------------------------------------------------------+------------------------------------------+
   | Last return from a linear solver function                          | :c:func:`MRIStepGetLastLinFlag()`        |
   +--------------------------------------------------------------------+------------------------------------------+
   | Name of constant associated with a return flag                     | :c:func:`MRIStepGetLinReturnFlagName()`  |
   +--------------------------------------------------------------------+------------------------------------------+

.. c:function:: int MRIStepGetJac(void* arkode_mem, SUNMatrix* J)

   Returns the internally stored copy of the Jacobian matrix of the ODE
   implicit slow right-hand side function.

   :param arkode_mem: the MRIStep memory structure
   :param J: the Jacobian matrix

   :retval ARKLS_SUCCESS: the output value has been successfully set
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARKLS_LMEM_NULL: the linear solver interface has not been initialized

   .. warning::

      This function is provided for debugging purposes and the values in the
      returned matrix should not be altered.

.. c:function:: int MRIStepGetJacTime(void* arkode_mem, sunrealtype* t_J)

   Returns the time at which the internally stored copy of the Jacobian matrix
   of the ODE implicit slow right-hand side function was evaluated.

   :param arkode_mem: the MRIStep memory structure
   :param t_J: the time at which the Jacobian was evaluated

   :retval ARKLS_SUCCESS: the output value has been successfully set
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARKLS_LMEM_NULL: the linear solver interface has not been initialized

.. c:function:: int MRIStepGetJacNumSteps(void* arkode_mem, long int* nst_J)

   Returns the value of the internal step counter at which the internally stored copy of the
   Jacobian matrix of the ODE implicit slow right-hand side function was
   evaluated.

   :param arkode_mem: the MRIStep memory structure
   :param nst_J: the value of the internal step counter at which the Jacobian was evaluated

   :retval ARKLS_SUCCESS: the output value has been successfully set
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARKLS_LMEM_NULL: the linear solver interface has not been initialized

.. c:function:: int MRIStepGetLinWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS)

   Returns the real and integer workspace used by the ARKLS linear solver interface.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *lenrwLS* -- the number of ``realtype`` values in the ARKLS workspace.

   * *leniwLS* -- the number of integer values in the ARKLS workspace.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** The workspace requirements reported by this routine
   correspond only to memory allocated within this interface and to
   memory allocated by the ``SUNLinearSolver`` object attached
   to it.  The template Jacobian matrix allocated by the user outside
   of ARKLS is not included in this report.

   In a parallel setting, the above values are global (i.e., summed over all
   processors).


.. c:function:: int MRIStepGetNumJacEvals(void* arkode_mem, long int* njevals)

   Returns the number of Jacobian evaluations.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *njevals* -- number of Jacobian evaluations.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** This is only accumulated for the "life" of the linear
   solver object; the counter is reset whenever a new linear solver
   module is "attached" to MRIStep, or when MRIStep is resized.


.. c:function:: int MRIStepGetNumPrecEvals(void* arkode_mem, long int* npevals)

   Returns the total number of preconditioner evaluations,
   i.e., the number of calls made to *psetup* with ``jok`` = ``SUNFALSE`` and
   that returned ``*jcurPtr`` = ``SUNTRUE``.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *npevals* -- the current number of calls to *psetup*.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** This is only accumulated for the "life" of the linear
   solver object; the counter is reset whenever a new linear solver
   module is "attached" to MRIStep, or when MRIStep is resized.


.. c:function:: int MRIStepGetNumPrecSolves(void* arkode_mem, long int* npsolves)

   Returns the number of calls made to the preconditioner
   solve function, *psolve*.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *npsolves* -- the number of calls to *psolve*.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** This is only accumulated for the "life" of the linear
   solver object; the counter is reset whenever a new linear solver
   module is "attached" to MRIStep, or when MRIStep is resized.


.. c:function:: int MRIStepGetNumLinIters(void* arkode_mem, long int* nliters)

   Returns the cumulative number of linear iterations.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *nliters* -- the current number of linear iterations.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** This is only accumulated for the "life" of the linear
   solver object; the counter is reset whenever a new linear solver
   module is "attached" to MRIStep, or when MRIStep is resized.


.. c:function:: int MRIStepGetNumLinConvFails(void* arkode_mem, long int* nlcfails)

   Returns the cumulative number of linear convergence failures.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *nlcfails* -- the current number of linear convergence failures.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** This is only accumulated for the "life" of the linear
   solver object; the counter is reset whenever a new linear solver
   module is "attached" to MRIStep, or when MRIStep is resized.


.. c:function:: int MRIStepGetNumJTSetupEvals(void* arkode_mem, long int* njtsetup)

   Returns the cumulative number of calls made to the user-supplied
   Jacobian-vector setup function, *jtsetup*.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *njtsetup* -- the current number of calls to *jtsetup*.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** This is only accumulated for the "life" of the linear
   solver object; the counter is reset whenever a new linear solver
   module is "attached" to MRIStep, or when MRIStep is resized.


.. c:function:: int MRIStepGetNumJtimesEvals(void* arkode_mem, long int* njvevals)

   Returns the cumulative number of calls made to the
   Jacobian-vector product function, *jtimes*.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *njvevals* -- the current number of calls to *jtimes*.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** This is only accumulated for the "life" of the linear
   solver object; the counter is reset whenever a new linear solver
   module is "attached" to MRIStep, or when MRIStep is resized.


.. c:function:: int MRIStepGetNumLinRhsEvals(void* arkode_mem, long int* nfevalsLS)

   Returns the number of calls to the user-supplied implicit
   right-hand side function :math:`f^I` for finite difference
   Jacobian or Jacobian-vector product approximation.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *nfevalsLS* -- the number of calls to the user implicit
     right-hand side function.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** The value *nfevalsLS* is incremented only if the default
   internal difference quotient function is used.

   This is only accumulated for the "life" of the linear
   solver object; the counter is reset whenever a new linear solver
   module is "attached" to MRIStep, or when MRIStep is resized.


.. c:function:: int MRIStepGetLastLinFlag(void* arkode_mem, long int* lsflag)

   Returns the last return value from an ARKLS routine.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *lsflag* -- the value of the last return flag from an
     ARKLS function.

   **Return value:**

   * *ARKLS_SUCCESS* if successful

   * *ARKLS_MEM_NULL* if the MRIStep memory was ``NULL``

   * *ARKLS_LMEM_NULL* if the linear solver memory was ``NULL``

   **Notes:** If the ARKLS setup function failed when using the
   ``SUNLINSOL_DENSE`` or ``SUNLINSOL_BAND`` modules, then the value
   of *lsflag* is equal to the column index (numbered from one) at
   which a zero diagonal element was encountered during the LU
   factorization of the (dense or banded) Jacobian matrix.  For all
   other failures, *lsflag* is negative.

   Otherwise, if the ARKLS setup function failed
   (:c:func:`MRIStepEvolve()` returned *ARK_LSETUP_FAIL*), then
   *lsflag* will be *SUNLS_PSET_FAIL_UNREC*, *SUNLS_ASET_FAIL_UNREC*
   or *SUNLS_PACKAGE_FAIL_UNREC*.

   If the ARKLS solve function failed (:c:func:`MRIStepEvolve()`
   returned *ARK_LSOLVE_FAIL*), then *lsflag* contains the error
   return flag from the ``SUNLinearSolver`` object, which will
   be one of:
   *SUNLS_MEM_NULL*, indicating that the ``SUNLinearSolver``
   memory is ``NULL``;
   *SUNLS_ATIMES_NULL*, indicating that a matrix-free iterative solver
   was provided, but is missing a routine for the matrix-vector product
   approximation,
   *SUNLS_ATIMES_FAIL_UNREC*, indicating an unrecoverable failure in
   the :math:`Jv` function;
   *SUNLS_PSOLVE_NULL*, indicating that an iterative linear solver was
   configured to use preconditioning, but no preconditioner solve
   routine was provided,
   *SUNLS_PSOLVE_FAIL_UNREC*, indicating that the preconditioner solve
   function failed unrecoverably;
   *SUNLS_GS_FAIL*, indicating a failure in the Gram-Schmidt procedure
   (SPGMR and SPFGMR only);
   *SUNLS_QRSOL_FAIL*, indicating that the matrix :math:`R` was found
   to be singular during the QR solve phase (SPGMR and SPFGMR only); or
   *SUNLS_PACKAGE_FAIL_UNREC*, indicating an unrecoverable failure in
   an external iterative linear solver package.


.. c:function:: char *MRIStepGetLinReturnFlagName(long int lsflag)

   Returns the name of the ARKLS constant corresponding to *lsflag*.

   **Arguments:**

   * *lsflag* -- a return flag from an ARKLS function.

   **Return value:**  The return value is a string containing the name of
   the corresponding constant. If using the ``SUNLINSOL_DENSE`` or
   ``SUNLINSOL_BAND`` modules, then if  1 :math:`\le` `lsflag`
   :math:`\le n` (LU factorization failed), this routine returns "NONE".




.. _ARKODE.Usage.MRIStep.MRIStepExtraOutputs:

General usability functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional routines may be called by a user to inquire
about existing solver parameters or write the current MRI coupling table. While
neither of these would typically be called during the course of solving an
initial value problem, these may be useful for users wishing to better
understand MRIStep.


.. _ARKODE.Usage.MRIStep.MRIStepExtraOutputs.Table:
.. table:: General usability functions

   +-----------------------------------------------+-------------------------------------------+
   | Optional routine                              | Function name                             |
   +-----------------------------------------------+-------------------------------------------+
   | Output all MRIStep solver parameters          | :c:func:`MRIStepWriteParameters()`        |
   +-----------------------------------------------+-------------------------------------------+
   | Output the current MRI coupling table         | :c:func:`MRIStepWriteCoupling()`          |
   +-----------------------------------------------+-------------------------------------------+


.. c:function:: int MRIStepWriteParameters(void* arkode_mem, FILE *fp)

   Outputs all MRIStep solver parameters to the provided file pointer.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *fp* -- pointer to use for printing the solver parameters.

   **Return value:**

   * *ARKS_SUCCESS* if successful

   * *ARKS_MEM_NULL* if the MRIStep memory was ``NULL``

   **Notes:** The *fp* argument can be ``stdout`` or ``stderr``, or it
   may point to a specific file created using ``fopen``.

   When run in parallel, only one process should set a non-NULL value
   for this pointer, since parameters for all processes would be
   identical.


.. c:function:: int MRIStepWriteCoupling(void* arkode_mem, FILE *fp)

   Outputs the current MRI coupling table to the provided file pointer.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *fp* -- pointer to use for printing the Butcher tables.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   **Notes:** The *fp* argument can be ``stdout`` or ``stderr``, or it
   may point to a specific file created using ``fopen``.

   When run in parallel, only one process should set a non-NULL value
   for this pointer, since tables for all processes would be
   identical.



.. _ARKODE.Usage.MRIStep.Reinitialization:

MRIStep re-initialization function
-------------------------------------

To reinitialize the MRIStep module for the solution of a new problem,
where a prior call to :c:func:`MRIStepCreate()` has been made, the
user must call the function :c:func:`MRIStepReInit()`.  The new
problem must have the same size as the previous one.  This routine
retains the current settings for all ARKstep module options and
performs the same input checking and initializations that are done in
:c:func:`MRIStepCreate()`, but it performs no memory allocation as is
assumes that the existing internal memory is sufficient for the new
problem.  A call to this re-initialization routine deletes the
solution history that was stored internally during the previous
integration, and deletes any previously-set *tstop* value specified via a
call to :c:func:`MRIStepSetStopTime()`.  Following a successful call to
:c:func:`MRIStepReInit()`, call :c:func:`MRIStepEvolve()` again for the
solution of the new problem.

The use of :c:func:`MRIStepReInit()` requires that the number of Runge--Kutta
stages for both the slow and fast methods be no larger for the new problem than
for the previous problem.

One important use of the :c:func:`MRIStepReInit()` function is in the
treating of jump discontinuities in the RHS functions.  Except in cases
of fairly small jumps, it is usually more efficient to stop at each
point of discontinuity and restart the integrator with a readjusted
ODE model, using a call to this routine.  To stop when the location
of the discontinuity is known, simply make that location a value of
``tout``.  To stop when the location of the discontinuity is
determined by the solution, use the rootfinding feature.  In either
case, it is critical that the RHS functions *not* incorporate the
discontinuity, but rather have a smooth extension over the
discontinuity, so that the step across it (and subsequent rootfinding,
if used) can be done efficiently.  Then use a switch within the RHS
functions (communicated through ``user_data``) that can be flipped
between the stopping of the integration and the restart, so that the
restarted problem uses the new values (which have jumped).  Similar
comments apply if there is to be a jump in the dependent variable
vector.


.. c:function:: int MRIStepReInit(void* arkode_mem, ARKRhsFn fse, ARKRhsFn fsi, realtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes the
   MRIStep outer (slow) stepper.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *fse* -- the name of the function (of type :c:func:`ARKRhsFn()`)
     defining the explicit slow portion of the right-hand side function in
     :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.

   * *fsi* -- the name of the function (of type :c:func:`ARKRhsFn()`)
     defining the implicit slow portion of the right-hand side function in
     :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.

   * *t0* -- the initial value of :math:`t`.

   * *y0* -- the initial condition vector :math:`y(t_0)`.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL*  if the MRIStep memory was ``NULL``

   * *ARK_MEM_FAIL*  if a memory allocation failed

   * *ARK_ILL_INPUT* if an argument has an illegal value.

   **Notes:**
   If the inner (fast) stepper also needs to be reinitialized, its
   reinitialization function should be called before calling
   :c:func:`MRIStepReInit()` to reinitialize the outer stepper.

   All previously set options are retained but may be updated by calling
   the appropriate "Set" functions.

   If an error occurred, :c:func:`MRIStepReInit()` also
   sends an error message to the error handler function.



.. _ARKODE.Usage.MRIStep.Reset:

MRIStep reset function
----------------------

To reset the MRIStep module to a particular state :math:`(t_R,y(t_R))` for the
continued solution of a problem, where a prior
call to :c:func:`MRIStepCreate()` has been made, the user must call the function
:c:func:`MRIStepReset()`.  Like :c:func:`MRIStepReInit()` this routine retains
the current settings for all MRIStep module options and performs no memory
allocations but, unlike :c:func:`MRIStepReInit()`, this routine performs only a
*subset* of the input checking and initializations that are done in
:c:func:`MRIStepCreate()`. In particular this routine retains all internal
counter values and the step size/error history and does not reinitialize the
linear and/or nonlinear solver but it does indicate that a linear solver setup
is necessary in the next step. Like :c:func:`MRIStepReInit()`, a call to
:c:func:`MRIStepReset()` will delete any previously-set *tstop* value specified
via a call to :c:func:`MRIStepSetStopTime()`.  Following a successful call to
:c:func:`MRIStepReset()`, call :c:func:`MRIStepEvolve()` again to continue
solving the problem. By default the next call to :c:func:`MRIStepEvolve()` will
use the step size computed by MRIStep prior to calling :c:func:`MRIStepReset()`.
To set a different step size or have MRIStep estimate a new step size use
:c:func:`MRIStepSetInitStep()`.

One important use of the :c:func:`MRIStepReset()` function is in the
treating of jump discontinuities in the RHS functions.  Except in cases
of fairly small jumps, it is usually more efficient to stop at each
point of discontinuity and restart the integrator with a readjusted
ODE model, using a call to :c:func:`MRIStepReset()`.  To stop when
the location of the discontinuity is known, simply make that location
a value of ``tout``.  To stop when the location of the discontinuity
is determined by the solution, use the rootfinding feature.  In either
case, it is critical that the RHS functions *not* incorporate the
discontinuity, but rather have a smooth extension over the
discontinuity, so that the step across it (and subsequent rootfinding,
if used) can be done efficiently.  Then use a switch within the RHS
functions (communicated through ``user_data``) that can be flipped
between the stopping of the integration and the restart, so that the
restarted problem uses the new values (which have jumped).  Similar
comments apply if there is to be a jump in the dependent variable
vector.


.. c:function:: int MRIStepReset(void* arkode_mem, realtype tR, N_Vector yR)

   Resets the current MRIStep outer (slow) time-stepper module state to the
   provided independent variable value and dependent variable vector.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *tR* -- the value of the independent variable :math:`t`.

   * *yR* -- the value of the dependent variable vector :math:`y(t_R)`.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL*  if the MRIStep memory was ``NULL``

   * *ARK_MEM_FAIL*  if a memory allocation failed

   * *ARK_ILL_INPUT* if an argument has an illegal value.

   **Notes:**
   If the inner (fast) stepper also needs to be reset, its reset function should
   be called before calling :c:func:`MRIStepReset()` to reset the outer stepper.

   All previously set options are retained but may be updated by calling
   the appropriate "Set" functions.

   If an error occurred, :c:func:`MRIStepReset()` also sends an error message to
   the error handler function.

   .. versionchanged:: 5.3.0

      This now calls the corresponding :c:type:`MRIStepInnerResetFn` with the same
      (*tR*, *yR*) arguments for the :c:type:`MRIStepInnerStepper` object that is
      used to evolve the MRI "fast" time scale subproblems.



.. _ARKODE.Usage.MRIStep.Resizing:

MRIStep system resize function
-------------------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when using spatially-adaptive
PDE simulations under a method-of-lines approach), the MRIStep
integrator may be "resized" between *slow* integration steps, through calls
to the :c:func:`MRIStepResize()` function. This function modifies
MRIStep's internal memory structures to use the new problem size.

To aid in the vector resize operation, the user can supply a vector
resize function that will take as input a vector with the previous
size, and transform it in-place to return a corresponding vector of
the new size.  If this function (of type :c:func:`ARKVecResizeFn()`)
is not supplied (i.e., is set to ``NULL``), then all existing vectors
internal to MRIStep will be destroyed and re-cloned from the new input
vector.


.. c:function:: int MRIStepResize(void* arkode_mem, N_Vector yR, realtype tR, ARKVecResizeFn resize, void* resize_data)

   Re-initializes MRIStep with a different state vector.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *yR* -- the newly-sized solution vector, holding the current
     dependent variable values :math:`y(t_R)`.

   * *tR* -- the current value of the independent variable
     :math:`t_R` (this must be consistent with *yR*).

   * *resize* -- the user-supplied vector resize function (of type
     :c:func:`ARKVecResizeFn()`.

   * *resize_data* -- the user-supplied data structure to be passed
     to *resize* when modifying internal MRIStep vectors.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL*  if the MRIStep memory was ``NULL``

   * *ARK_NO_MALLOC* if *arkode_mem* was not allocated.

   * *ARK_ILL_INPUT* if an argument has an illegal value.

   **Notes:** If an error occurred, :c:func:`MRIStepResize()` also sends an error
   message to the error handler function.

   **Resizing the linear solver:**
      When using any of the SUNDIALS-provided linear solver modules, the
      linear solver memory structures must also be resized.  At present,
      none of these include a solver-specific "resize" function, so the linear
      solver memory must be destroyed and re-allocated **following** each
      call to :c:func:`MRIStepResize()`.  Moreover, the existing ARKLS
      interface should then be deleted and recreated by attaching the
      updated ``SUNLinearSolver`` (and possibly ``SUNMatrix``) object(s)
      through calls to :c:func:`MRIStepSetLinearSolver()`.

      If any user-supplied routines are provided to aid the linear solver
      (e.g. Jacobian construction, Jacobian-vector product,
      mass-matrix-vector product, preconditioning), then the corresponding
      "set" routines must be called again **following** the solver
      re-specification.

   **Resizing the absolute tolerance array:**
      If using array-valued absolute tolerances, the absolute tolerance
      vector will be invalid after the call to :c:func:`MRIStepResize()`, so
      the new absolute tolerance vector should be re-set **following** each
      call to :c:func:`MRIStepResize()` through a new call to
      :c:func:`MRIStepSVtolerances()` and possibly
      :c:func:`MRIStepResVtolerance()` if applicable.

      If scalar-valued tolerances or a tolerance function was specified
      through either :c:func:`MRIStepSStolerances()` or
      :c:func:`MRIStepWFtolerances()`, then these will remain valid and no
      further action is necessary.

      .. note::

         For an example showing usage of the similar :c:func:`ARKStepResize()`
         routine, see the supplied serial C example problem,
         ``ark_heat1D_adapt.c``.
