..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdjoint.Stepper:

The SUNAdjointStepper API
=========================

.. versionadded:: x.y.z

The :c:type:`SUNAdjointStepper` API provides a package-agnostic interface to SUNDIALS ASA
capabilities. It currently only supports the discrete ASA capabilities in the ARKODE
package, but in the future this support may be expanded.

A :c:type:`SUNAdjointStepper` is a pointer to the
:c:struct:`SUNAdjointStepper_` structure:

.. c:type:: struct SUNAdjointStepper_ *SUNAdjointStepper

.. c:struct:: SUNAdjointStepper_

   .. c:member:: SUNStepper adj_sunstepper

      The :c:type:`SUNStepper` object used for backwards time stepping of the adjoint ODE system.

   .. c:member:: SUNStepper fwd_sunstepper

      The :c:type:`SUNStepper` object used for forward time stepping of the original ODE system if any recomputation of missing
      state data is required during the backwards integration.

   .. c:member:: sunbooleantype own_adj_sunstepper

      If true, then the :c:type:`SUNAdjointStepper` will be responsible for destroying the ``adj_sunstepper``.

   .. c:member:: sunbooleantype own_fwd_sunstepper

      If true, then the :c:type:`SUNAdjointStepper` will be responsible for destroying the ``fwd_sunstepper``.

   .. c:member:: sunrealtype tf

      The terminal time of the backwards adjoint ODE.

   .. c:member:: suncountertype step_idx

      The index of the current backward integration step with respect to the forward integration.

   .. c:member:: suncountertype final_step_idx

      The index of the final step in the forward integration (corresponds to ``tf``).

   .. c:member:: SUNMatrix Jac

      Matrix data for the Jacobian :math:`\partial f / \partial y`.

   .. c:member:: SUNMatrix JacP

      Matrix data for the Jacobian :math:`\partial f / \partial p`.

   .. c:member:: SUNRhsJacFn JacFn

      Jacobian function pointer to evaluate :math:`\partial f / \partial y`.

   .. c:member:: SUNRhsJacFn JacPFn

      Jacobian function pointer to evaluate :math:`\partial f / \partial p`.

   .. c:member:: SUNRhsJacTimesFn JvpFn

      Jacobian-times-vector function pointer to evaluate :math:`(\partial f/\partial y)^* v`
      or :math:`v^*(\partial f/\partial y)`.

   .. c:member:: SUNRhsJacTimesFn JPvpFn

      Jacobian-times-vector function pointer to evaluate :math:`(\partial f/\partial p)^* v`,
      or :math:`v^*(\partial f/\partial p)`.

   .. c:member:: suncountertype nst

      Holds the count of the number of backwards steps taken.

   .. c:member:: suncountertype njeval

      Holds the count of the number of :math:`\partial f / \partial y` evaluations.

   .. c:member:: suncountertype njpeval

      Holds the count of the number of :math:`\partial f / \partial p` evaluations.

   .. c:member:: suncountertype njtimesv

      Holds the count of the number of :math:`(\partial f/\partial y)^* v`, or
      :math:`v^*(\partial f/\partial y)` evaluations.

   .. c:member:: suncountertype njptimesv

      Holds the count of the number of :math:`(\partial f/\partial p)^* v`, or
      :math:`v^*(\partial f/\partial p)`evaluations.

   .. c:member:: suncountertype nrecompute

      Holds the count of the number of partial recomputations of the forward problem.

   .. c:member:: void* user_data

      A pointer that is passed back to user-supplied functions

   .. c:member:: void* content

      Pointer to derived class specific member data

   .. c:member:: SUNContext sunctx

      The SUNDIALS simulation context


The :c:type:`SUNAdjointStepper` class has the following functions:

.. c:function:: SUNErrCode SUNAdjointStepper_Create(SUNStepper fwd_sunstepper, sunbooleantype own_fwd, \
   SUNStepper adj_sunstepper, sunbooleantype own_adj, suncountertype final_step_idx, N_Vector sf, \
   sunrealtype tf, SUNAdjointCheckpointScheme checkpoint_scheme, SUNContext sunctx, SUNAdjointStepper* adj_stepper)

   Creates the ``SUNAdjointStepper`` object needed to solve the adjoint problem.

   :param fwd_sunstepper: The :c:type:`SUNStepper` to be used for forward computations of the original ODE.
   :param own_fwd: Should `fwd_sunstepper` be owned (and destroyed) by the `SUNAdjointStepper` or not.
   :param adj_sunstepper: The :c:type:`SUNStepper` to be used for the backward integration of the adjoint ODE.
   :param own_adj: Should `adj_sunstepper` be owned (and destroyed) by the `SUNAdjointStepper` or not.
   :param final_step_idx: The index (step number) of the step corresponding to ``t_f`` for the forward ODE.
   :param sf: The terminal condition for the adjoint ODE.
   :param tf: The terminal time for the forward ODE (the initial time for the adjoint ODE).
   :param checkpoint_scheme: The :c:type:`SUNAdjointCheckpointScheme` object that determines the checkpointing strategy to use. This should be the same object provided to the forward integrator/stepper.
   :param sunctx: The :c:type:`SUNContext` for the simulation.
   :param adj_stepper: The :c:type:`SUNAdjointStepper` to construct (will be ``NULL`` on failure).

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper adj, N_Vector sf, sunrealtype tf)

   Reinitializes the adjoint stepper to solve a new problem of the same size.

   :param adj_stepper: The adjoint solver object.
   :param sf: The terminal condition vector of sensitivity solutions :math:`\partial g/\partial y_0` and :math:`\partial g/\partial p`.
   :param tf: The time to start integrating the adjoint system from.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper adj_stepper, sunrealtype tout,\
   N_Vector sens, sunrealtype* tret)

   Integrates the adjoint system.

   :param adj_stepper: The adjoint solver object.
   :param tout: The time at which the adjoint solution is desired.
   :param sens: The vector of sensitivity solutions :math:`\partial g/\partial y_0` and :math:`\partial g/\partial p`.
   :param tret: On return, the time reached by the adjoint solver.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper adj_stepper, sunrealtype tout,\
   N_Vector sens, sunrealtype* tret)

   Evolves the adjoint system backwards one step.

   :param adj_stepper: The adjoint solver object.
   :param tout: The time at which the adjoint solution is desired.
   :param sens: The vector of sensitivity solutions :math:`\partial g/\partial y_0` and :math:`\partial g/\partial p`.
   :param tret: On return, the time reached by the adjoint solver.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper adj_stepper, suncountertype start_idx,\
                                                          sunrealtype t0, sunrealtype tf, N_Vector y0)

   Evolves the forward system in time from (``start_idx``, ``t0``) to (``stop_idx``, ``tf``) with dense checkpointing.

   :param adj_stepper: The SUNAdjointStepper object.
   :param start_idx: the index of the step, w.r.t. the original forward integration, to begin forward integration from.
   :param t0: the initial time, w.r.t. the original forward integration, to start forward integration from.
   :param tf: the final time, w.r.t. the original forward integration, to stop forward integration.
   :param y0: the initial state, w.r.t. the original forward integration, to start forward integration.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_SetJacFn(SUNAdjointStepper adj_stepper, SUNRhsJacFn JacFn, \
      SUNMatrix Jac, SUNRhsJacFn JacPFn, SUNMatrix JacP)

   Sets the function pointers and matrices needed to evaluate and store :math:`\partial f / \partial y` and
   :math:`\partial f / \partial p`. ``Jac`` should have dimensions ``neq x neq`` where ``neq`` is the number of states
   in the forward problem. ``JacP`` should have dimensions ``nparams x neq`` where ``nparams`` is the
   number of parameters in the model to get sensitivities for.

   :param adj_stepper: The SUNAdjointStepper object.
   :param JacFn: the function that evaluates :math:`\partial f / \partial y`.
   :param Jac: a :c:type:`SUNMatrix` that will hold :math:`\partial f / \partial y`.
   :param JacPFn: the function that evaluates :math:`\partial f / \partial p`.
   :param JacP: a :c:type:`SUNMatrix` that will hold :math:`\partial f / \partial p`.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_SetJacHermitianTransposeVecFn(SUNAdjointStepper adj_stepper, SUNRhsJacTimesFn Jvp, SUNRhsJacTimesFn JPvp)


   Sets the function pointers to evaluate :math:`(\partial f/\partial y)^* v`  and :math:`(\partial f/\partial p)^* v`

   :param adj_stepper: The SUNAdjointStepper object.
   :param Jvp: function that evaluates :math:`(\partial f/\partial y)^* v`.
   :param JPvp: function that evaluates :math:`(\partial f/\partial p)^* v`.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_SetUserData(SUNAdjointStepper adj_stepper, void* user_data)

   Sets the user data pointer.

   :param adj_stepper: The SUNAdjointStepper object.
   :param user_data: the user data pointer that will be passed back to user-supplied callback functions.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumSteps(SUNAdjointStepper adj_stepper, suncountertype* num_steps)

   Retrieves the number of steps taken by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_steps: Pointer to store the number of steps.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumJacEvals(SUNAdjointStepper adj_stepper, suncountertype* num_jac_evals)

   Retrieves the number of Jacobian evaluations performed by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_jac_evals: Pointer to store the number of Jacobian evaluations.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumJacPEvals(SUNAdjointStepper adj_stepper, suncountertype* num_jac_p_evals)

   Retrieves the number of Jacobian parameter evaluations performed by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_jac_p_evals: Pointer to store the number of Jacobian parameter evaluations.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumJacTimesVecEvals(SUNAdjointStepper adj_stepper, suncountertype* num_jac_times_vec_evals)

   Retrieves the number of Jacobian-times-vector evaluations performed by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_jac_times_vec_evals: Pointer to store the number of Jacobian-times-vector evaluations.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumJacPTimesVecEvals(SUNAdjointStepper adj_stepper, suncountertype* num_jac_p_times_vec_evals)

   Retrieves the number of parameter Jacobian-times-vector evaluations performed by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_jac_p_times_vec_evals: Pointer to store the number of parameter Jacobian-times-vector evaluations.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumVecTimesJacEvals(SUNAdjointStepper adj_stepper, suncountertype* num_vec_times_jac_evals)

   Retrieves the number of vector-times-Jacobian evaluations performed by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_vec_times_jac_evals: Pointer to store the number of vector-times-Jacobian evaluations.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumVecTimesJacPEvals(SUNAdjointStepper adj_stepper, suncountertype* num_vec_times_jac_p_evals)

   Retrieves the number of vector-times-Jacobian parameter evaluations performed by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_vec_times_jac_p_evals: Pointer to store the number of vector-times-Jacobian parameter evaluations.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumRecompute(SUNAdjointStepper adj_stepper, suncountertype* num_recompute)

   Retrieves the number of recomputations performed by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_recompute: Pointer to store the number of recomputations.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_PrintAllStats(SUNAdjointStepper adj_stepper, \
                                                           FILE* outfile, SUNOutputFormat fmt)

   Prints the adjoint stepper statistics/counters in a human-readable table format or CSV format.

   :param adj_stepper: The SUNAdjointStepper object.
   :param outfile: A file to write the output to.
   :param fmt: the format to write in (:c:type:`SUN_OUTPUTFORMAT_TABLE` or :c:type:`SUN_OUTPUTFORMAT_CSV`).

   :return: A :c:type:`SUNErrCode` indicating failure or success.
