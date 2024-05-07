.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.MRIStep.UserCallable:

MRIStep User-callable functions
==================================

This section describes the MRIStep-specific functions that may be called
by the user to setup and then solve an IVP using the MRIStep time-stepping
module.  The large majority of these routines merely wrap :ref:`underlying
ARKODE functions <ARKODE.Usage.UserCallable>`, and are now deprecated
-- each of these are clearly marked.  However, some
of these user-callable functions are specific to ERKStep, as explained
below.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
MRIStep supports the following categories:

* implicit nonlinear and/or linear solvers



.. _ARKODE.Usage.MRIStep.Initialization:

MRIStep initialization and deallocation functions
------------------------------------------------------


.. c:function:: void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, sunrealtype t0, N_Vector y0, MRIStepInnerStepper stepper, SUNContext sunctx)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeFree` instead.



.. _ARKODE.Usage.MRIStep.Tolerances:

MRIStep tolerance specification functions
------------------------------------------------------

.. c:function:: int MRIStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSStolerances` instead.



.. c:function:: int MRIStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSVtolerances` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeWFtolerances` instead.



.. _ARKODE.Usage.MRIStep.LinearSolvers:

Linear solver interface functions
-------------------------------------------

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetLinearSolver` instead.



.. _ARKODE.Usage.MRIStep.NonlinearSolvers:

Nonlinear solver interface functions
-------------------------------------------

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetNonlinearSolver` instead.



.. _ARKODE.Usage.MRIStep.RootFinding:

Rootfinding initialization function
--------------------------------------

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeRootInit` instead.



.. _ARKODE.Usage.MRIStep.Integration:

MRIStep solver function
-------------------------

.. c:function:: int MRIStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout, sunrealtype *tret, int itask)

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
        single internal step, :math:`y_{n-1} \to y_{n}`, and return the solution
        at that point, :math:`y_{n}`, in the vector *yout*.

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

      All failure return values are negative and so testing the return argument
      for negative values will trap all :c:func:`MRIStepEvolve()` failures.

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeEvolve` instead.



.. _ARKODE.Usage.MRIStep.OptionalInputs:

Optional input functions
-------------------------


.. _ARKODE.Usage.MRIStep.MRIStepInput:

Optional inputs for MRIStep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetDefaults` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetInterpolantType` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetInterpolantDegree` instead.



.. c:function:: int MRIStepSetDenseOrder(void* arkode_mem, int dord)

   .. deprecated:: 5.2.0

      Use :c:func:`ARKodeSetInterpolantDegree` instead.


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



.. c:function:: int MRIStepSetFixedStep(void* arkode_mem, sunrealtype hs)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetFixedStep` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetMaxHnilWarns` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetMaxNumSteps` instead.



.. c:function:: int MRIStepSetStopTime(void* arkode_mem, sunrealtype tstop)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetStopTime` instead.



.. c:function:: int MRIStepSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetInterpolateStopTime` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeClearStopTime` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetUserData` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetOrder` instead.



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

   **Warning:**

   This should not be used with :c:func:`ARKodeSetOrder`.



.. _ARKODE.Usage.MRIStep.MRIStepSolverInput:

Optional inputs for implicit stage solves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetLinear` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetNonlinear` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetPredictorMethod` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetMaxNonlinIters` instead.



.. c:function:: int MRIStepSetNonlinConvCoef(void* arkode_mem, sunrealtype nlscoef)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetNonlinConvCoef` instead.



.. c:function:: int MRIStepSetNonlinCRDown(void* arkode_mem, sunrealtype crdown)

   Specifies the constant used in estimating the nonlinear solver convergence rate.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *crdown* -- nonlinear convergence rate estimation constant (default is 0.3).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:** Any non-positive parameter will imply a reset to the default value.

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetNonlinCRDown` instead.



.. c:function:: int MRIStepSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetNonlinRDiv` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetStagePredictFn` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetNlsRhsFn` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetDeduceImplicitRhs` instead.



.. _ARKODE.Usage.MRIStep.ARKLsInputs:

Linear solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. _ARKODE.Usage.MRIStep.ARKLsInputs.General:

Optional inputs for the ARKLS linear solver interface
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. c:function:: int MRIStepSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetDeltaGammaMax` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetLSetupFrequency` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetJacEvalFrequency` instead.




.. _ARKODE.Usage.MRIStep.ARKLsInputs.MatrixBased:

Optional inputs for matrix-based ``SUNLinearSolver`` modules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetJacFn` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetLinSysFn` instead.



.. c:function:: int MRIStepSetLinearSolutionScaling(void* arkode_mem, sunbooleantype onoff)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetLinearSolutionScaling` instead.



.. _ARKODE.Usage.MRIStep.ARKLsInputs.MatrixFree:

Optional inputs for matrix-free ``SUNLinearSolver`` modules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetJacTimes` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetJacTimesRhsFn` instead.





.. _ARKODE.Usage.MRIStep.ARKLsInputs.Iterative:

Optional inputs for iterative ``SUNLinearSolver`` modules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetPreconditioner` instead.



.. c:function:: int MRIStepSetEpsLin(void* arkode_mem, sunrealtype eplifac)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetEpsLin` instead.



.. c:function:: int MRIStepSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetLSNormFactor` instead.



.. _ARKODE.Usage.MRIStep.MRIStepRootfindingInput:

Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetRootDirection` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeSetNoInactiveRootWarn` instead.



.. _ARKODE.Usage.MRIStep.InterpolatedOutput:

Interpolated output function
--------------------------------

.. c:function:: int MRIStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetDky` instead.



.. _ARKODE.Usage.MRIStep.OptionalOutputs:

Optional output functions
------------------------------


.. _ARKODE.Usage.MRIStep.MRIStepMainOutputs:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: int MRIStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)

   Returns the MRIStep real and integer workspace sizes.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *lenrw* -- the number of ``sunrealtype`` values in the MRIStep workspace.

   * *leniw* -- the number of integer values in the MRIStep workspace.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetWorkSpace` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumSteps` instead.



.. c:function:: int MRIStepGetLastStep(void* arkode_mem, sunrealtype* hlast)

   Returns the integration step size taken on the last successful
   internal step.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *hlast* -- step size taken on the last internal step.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetLastStep` instead.



.. c:function:: int MRIStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur)

   Returns the current internal time reached by the solver.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *tcur* -- current internal time reached.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetCurrentTime` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetCurrentState` instead.


.. c:function:: int MRIStepGetCurrentGamma(void *arkode_mem, sunrealtype *gamma)

   Returns the current internal value of :math:`\gamma` used in the implicit
   solver Newton matrix (see equation :eq:`ARKODE_NewtonMatrix`).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *gamma* -- current step size scaling factor in the Newton system.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetCurrentGamma` instead.


.. c:function:: int MRIStepGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfac)

   Returns a suggested factor by which the user's
   tolerances should be scaled when too much accuracy has been
   requested for some internal step.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *tolsfac* -- suggested scaling factor for user-supplied tolerances.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetTolScaleFactor` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetErrWeights` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodePrintAllStats` instead.


.. c:function:: char* MRIStepGetReturnFlagName(long int flag)

   Returns the name of the MRIStep constant corresponding to *flag*.
   See :ref:`ARKODE.Constants`.

   **Arguments:**

   * *flag* -- a return flag from an MRIStep function.

   **Return value:**
   The return value is a string containing the name of
   the corresponding constant.

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetReturnFlagName` instead.



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



.. c:function:: int MRIStepGetNumStepSolveFails(void* arkode_mem, long int* ncnf)

   Returns the number of failed steps due to a nonlinear solver failure (so far).

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *ncnf* -- number of step failures.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumStepSolveFails` instead.


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
         sunrealtype ***G;   /* coupling matrices [nmat][stages][stages]    */
         sunrealtype *c;     /* abscissae                                   */

       };
       typedef MRIStepCouplingMem *MRIStepCoupling;

   For more details see :numref:`ARKODE.Usage.MRIStep.MRIStepCoupling`.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetUserData` instead.



.. _ARKODE.Usage.MRIStep.MRIStepImplicitSolverOutputs:

Implicit solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumLinSolvSetups` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumNonlinSolvIters` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumNonlinSolvConvFails` instead.



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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNonlinSolvStats` instead.



.. _ARKODE.Usage.MRIStep.MRIStepRootOutputs:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetRootInfo` instead.



.. c:function:: int MRIStepGetNumGEvals(void* arkode_mem, long int* ngevals)

   Returns the cumulative number of calls made to the
   user's root function :math:`g`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *ngevals* -- number of calls made to :math:`g` so far.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumGEvals` instead.



.. _ARKODE.Usage.MRIStep.ARKLsOutputs:

Linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetJac` instead.


.. c:function:: int MRIStepGetJacTime(void* arkode_mem, sunrealtype* t_J)

   Returns the time at which the internally stored copy of the Jacobian matrix
   of the ODE implicit slow right-hand side function was evaluated.

   :param arkode_mem: the MRIStep memory structure
   :param t_J: the time at which the Jacobian was evaluated

   :retval ARKLS_SUCCESS: the output value has been successfully set
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARKLS_LMEM_NULL: the linear solver interface has not been initialized

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetJacTime` instead.


.. c:function:: int MRIStepGetJacNumSteps(void* arkode_mem, long int* nst_J)

   Returns the value of the internal step counter at which the internally stored copy of the
   Jacobian matrix of the ODE implicit slow right-hand side function was
   evaluated.

   :param arkode_mem: the MRIStep memory structure
   :param nst_J: the value of the internal step counter at which the Jacobian was evaluated

   :retval ARKLS_SUCCESS: the output value has been successfully set
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARKLS_LMEM_NULL: the linear solver interface has not been initialized

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetJacNumSteps` instead.


.. c:function:: int MRIStepGetLinWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS)

   Returns the real and integer workspace used by the ARKLS linear solver interface.

   **Arguments:**

   * *arkode_mem* -- pointer to the MRIStep memory block.

   * *lenrwLS* -- the number of ``sunrealtype`` values in the ARKLS workspace.

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetLinWorkSpace` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumJacEvals` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumPrecEvals` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumPrecSolves` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumLinIters` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumLinConvFails` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumJTSetupEvals` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumJtimesEvals` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetNumLinRhsEvals` instead.


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
   or *SUN_ERR_EXT_FAIL*.

   If the ARKLS solve function failed (:c:func:`MRIStepEvolve()`
   returned *ARK_LSOLVE_FAIL*), then *lsflag* contains the error
   return flag from the ``SUNLinearSolver`` object, which will
   be one of:
   *SUN_ERR_ARG_CORRUPTRRUPT*, indicating that the ``SUNLinearSolver``
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
   *SUN_ERR_EXT_FAIL*, indicating an unrecoverable failure in
   an external iterative linear solver package.

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetLastLinFlag` instead.


.. c:function:: char* MRIStepGetLinReturnFlagName(long int lsflag)

   Returns the name of the ARKLS constant corresponding to *lsflag*.

   **Arguments:**

   * *lsflag* -- a return flag from an ARKLS function.

   **Return value:**  The return value is a string containing the name of
   the corresponding constant. If using the ``SUNLINSOL_DENSE`` or
   ``SUNLINSOL_BAND`` modules, then if  1 :math:`\le` `lsflag`
   :math:`\le n` (LU factorization failed), this routine returns "NONE".

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeGetLinReturnFlagName` instead.




.. _ARKODE.Usage.MRIStep.MRIStepExtraOutputs:

General usability functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeWriteParameters` instead.


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

   .. deprecated:: x.y.z

      Use :c:func:`MRIStepGetCurrentCoupling` and :c:func:`MRIStepCoupling_Write`
      instead.


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


.. c:function:: int MRIStepReInit(void* arkode_mem, ARKRhsFn fse, ARKRhsFn fsi, sunrealtype t0, N_Vector y0)

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

.. c:function:: int MRIStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)

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

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeReset` instead.



.. _ARKODE.Usage.MRIStep.Resizing:

MRIStep system resize function
-------------------------------------

.. c:function:: int MRIStepResize(void* arkode_mem, N_Vector yR, sunrealtype tR, ARKVecResizeFn resize, void* resize_data)

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
      :c:func:`MRIStepSVtolerances()`.

      If scalar-valued tolerances or a tolerance function was specified
      through either :c:func:`MRIStepSStolerances()` or
      :c:func:`MRIStepWFtolerances()`, then these will remain valid and no
      further action is necessary.

      .. note::

         For an example showing usage of the similar :c:func:`ARKStepResize()`
         routine, see the supplied serial C example problem,
         ``ark_heat1D_adapt.c``.

   .. deprecated:: x.y.z

      Use :c:func:`ARKodeResize` instead.
