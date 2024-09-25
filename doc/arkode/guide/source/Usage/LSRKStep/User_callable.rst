.. ----------------------------------------------------------------
   Programmer(s): Mustafa Aggul @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.LSRKStep.UserCallable:

LSRKStep User-callable functions
==================================

This section describes the LSRKStep-specific functions that may be called
by the user to setup and then solve an IVP using the LSRKStep time-stepping
module.  As mentioned in Section :numref:`ARKODE.Usage.UserCallable`,
shared ARKODE-level routines may be used for the large majority of LSRKStep
configuration and use.  In this section, we describe only those routines
that are specific to LSRKStep.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
LSRKStep supports the following categories:

* temporal adaptivity



.. _ARKODE.Usage.LSRKStep.Initialization:

LSRKStep initialization functions
---------------------------------


.. c:function:: void* LSRKStepCreate(ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0, N_Vector y0, SUNContext sunctx);

   This function allocates and initializes memory for a problem to
   be solved using the LSRKStep time-stepping module in ARKODE.

   **Arguments:**
      * *fe* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the explicit portion of the right-hand side function in
        :math:`y'(t) = f^E(t,y) + f^I(t,y)`.
      * *fi* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the implicit portion of the right-hand side function in
        :math:`y'(t) = f^E(t,y) + f^I(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a pointer to initialized problem memory
      of type ``void*``, to be passed to all user-facing LSRKStep routines
      listed below.  If unsuccessful, a ``NULL`` pointer will be
      returned, and an error message will be printed to ``stderr``.

   .. warning: LSRKStep does not currently support implicit or ImEx methods. The *fi* argument is included here as a placeholder for upcoming development.

.. _ARKODE.Usage.LSRKStep.OptionalInputs:

Optional input functions
-------------------------


.. c:function:: int LSRKStepSetMethod(void* arkode_mem, ARKODE_LSRKMethodType method);

   This function selects the LSRK method that should be used.  The list of allowable values for this input is below.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *method* -- Type of the method.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. typo in the method type).

   .. note:: If this routine is not called, then LSRKStep will use the Runge--Kutta--Chebyshev method by default.


Allowable Method Families

.. c:enum:: ARKODE_LSRKMethodType

   * ``ARKODE_LSRK_RKC_2`` -- Runge--Kutta--Chebyshev
   * ``ARKODE_LSRK_RKL_2`` -- Runge--Kutta--Legendre
   * ``ARKODE_LSRK_SSP_S_2`` -- SSP(s,2) -- 2nd order, s-stage
   * ``ARKODE_LSRK_SSP_S_3`` -- SSP(s,3) -- 3rd order, s-stage
   * ``ARKODE_LSRK_SSP_10_4`` -- SSP(10,4) -- 4th order, 10-stage

.. c:function:: int LSRKStepSetDomEigFn(void* arkode_mem, ARKDomEigFn dom_eig);

   Specifies the dominant eigenvalue approximation routine to
   be used for determining the number of stages that will be used by either the
   RKC or RKL methods.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *dom_eig* -- name of user-supplied dominant eigenvalue approximation function (of type :c:func:`ARKDomEigFn()`).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARKLS_MEM_NULL* if ``arkode_mem`` was ``NULL``.
      * *ARK_ILL_INPUT* ``dom_eig = NULL`` and LSRKStep does not currently estimate this internally.

   .. note:: This function is currently required when either the RKC or RKL methods are used.


.. c:function:: int LSRKStepSetDomEigFrequency(void* arkode_mem, int nsteps);

   Specifies the number of steps after which the dominant eigenvalue information is
   considered out-of-date, and should be recomputed. This only applies to RKL and RKC methods.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *nsteps* -- the dominant eigenvalue re-computation update frequency.  A value ``nsteps = 0`` indicates that the dominant eigenvalue will not change throughout the simulation.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARKLS_MEM_NULL* if ``arkode_mem`` was ``NULL``.
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. ``nsteps < 0``)


.. c:function:: int LSRKStepSetMaxNumStages(void* arkode_mem, int stage_max_limit);

   Specifies the maximum number of stages allowed within each time step.  This bound only applies to
   RKL and RKC methods.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *stage_max_limit* -- maximum allowed number of stages :math:`(>1)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARKLS_MEM_NULL* if ``arkode_mem`` was ``NULL``.
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. ``stage_max_limit < 2``)


.. c:function:: int LSRKStepSetDomEigSafetyFactor(void* arkode_mem, sunrealtype dom_eig_sfty);

   Specifies a safety factor to use for the result of the dominant eigenvalue estimation function.  This value is used to scale the magnitude of the dominant eigenvalue, in the hope of ensuring a sufficient number of stages for the method to be stable.  This input is only used for RKC and RKL methods.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *dom_eig_sfty* -- safety factor :math:`(\ge 1)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARKLS_MEM_NULL* if ``arkode_mem`` was ``NULL``.
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. ``dom_eig_sfty < 1``)


.. c:function:: int LSRKStepSetSSPStageNum(void* arkode_mem, int num_of_stages);

   Sets the number of stages, ``s`` in ``SSP(s, p)`` methods. This input is only utilized by SSPRK methods.

      * ``ARKODE_LSRK_SSP_S_2``  -- ``num_of_stages`` must be greater than or equal to 2
      * ``ARKODE_LSRK_SSP_S_3``  -- ``num_of_stages`` must be a perfect-square greater than or equal to 9
      * ``ARKODE_LSRK_SSP_10_4`` -- ``num_of_stages`` cannot be modified from 10, so this function should not be called.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *num_of_stages* -- number of stages :math:`(>1)` for ``SSP(s,2)`` and :math:`(n^2 = s \geq 9)` for ``SSP(s,3)``.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARKLS_MEM_NULL* if ``arkode_mem`` was ``NULL``.
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. ``num_of_stages < 2``)


.. _ARKODE.Usage.LSRKStep.OptionalOutputs:

Optional output functions
------------------------------


.. c:function:: int LSRKStepGetNumRhsEvals(void* arkode_mem, long int* fe_evals, long int* fi_evals);

   Returns the number of calls to the user's right-hand
   side functions, :math:`f^E` and :math:`f^I` (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *fe_evals* -- number of calls to the user's :math:`f^E(t,y)` function.
      * *fi_evals* -- number of calls to the user's :math:`f^I(t,y)` function.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``


.. c:function:: int LSRKStepGetNumDomEigUpdates(void* arkode_mem, long int* num_dom_eig_updates);

   Returns the number of dominant eigenvalue evaluations (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *num_dom_eig_updates* -- number of calls to the user's ``dom_eig`` function.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``


.. c:function:: int LSRKStepGetMaxNumStages(void* arkode_mem, int* stage_max);

   Returns the max number of stages used in any single step (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *stage_max* -- max number of stages used.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``


.. c:function:: int LSRKStepGetAverageStageNum(void* arkode_mem, sunrealtype* avg_stage);

   Returns the average number of stages per step (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *avg_stage* -- average number of stages.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``


.. _ARKODE.Usage.LSRKStep.Reinitialization:

LSRKStep re-initialization function
-------------------------------------

To reinitialize the LSRKStep module for the solution of a new problem,
where a prior call to :c:func:`LSRKStepCreate` has been made, the
user must call the function :c:func:`LSRKStepReInit()`.  The new
problem must have the same size as the previous one.  This routine
retains the current settings for all LSRKstep module options and
performs the same input checking and initializations that are done in
:c:func:`LSRKStepCreate`, but it performs no memory allocation as it
assumes that the existing internal memory is sufficient for the new
problem.  A call to this re-initialization routine deletes the
solution history that was stored internally during the previous
integration, and deletes any previously-set *tstop* value specified via a
call to :c:func:`ARKodeSetStopTime()`.  Following a successful call to
:c:func:`LSRKStepReInit()`, call :c:func:`ARKodeEvolve()` again for the
solution of the new problem.

One important use of the :c:func:`LSRKStepReInit()` function is in the
treating of jump discontinuities in the RHS function.  Except in cases
of fairly small jumps, it is usually more efficient to stop at each
point of discontinuity and restart the integrator with a readjusted
ODE model, using a call to this routine.  To stop when the location
of the discontinuity is known, simply make that location a value of
``tout``.  To stop when the location of the discontinuity is
determined by the solution, use the rootfinding feature.  In either
case, it is critical that the RHS function *not* incorporate the
discontinuity, but rather have a smooth extension over the
discontinuity, so that the step across it (and subsequent rootfinding,
if used) can be done efficiently.  Then use a switch within the RHS
function (communicated through ``user_data``) that can be flipped
between the stopping of the integration and the restart, so that the
restarted problem uses the new values (which have jumped).  Similar
comments apply if there is to be a jump in the dependent variable
vector.


.. c:function:: int LSRKStepReInit(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0, N_Vector y0);

   Provides required problem specifications and re-initializes the
   LSRKStep time-stepper module.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *fe* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the explicit right-hand side function in :math:`\dot{y} = f^E(t,y)`.
      * *fi* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the implicit right-hand side function in :math:`\dot{y} = f^I(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the LSRKStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if memory allocation failed
      * *ARK_NO_MALLOC*  if memory allocation failed
      * *ARK_CONTROLLER_ERR*  if unable to reset error controller object
      * *ARK_ILL_INPUT* if an argument had an illegal value.

   **Notes:**
      All previously set options are retained but may be updated by calling
      the appropriate "Set" functions.

      If an error occurred, :c:func:`LSRKStepReInit()` also
      sends an error message to the error handler function.
