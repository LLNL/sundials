..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdjointStepper:

The SUNAdjointStepper API and Module
====================================

The :c:type:`SUNAdjointStepper` API and module provides a package-agnostic
interface to SUNDIALS ASA capabilities. It currently only supports the discrete
ASA capabilities in the ARKODE package, but in the future this support may be expanded.

A :c:type:`SUNAdjointStepper` is a pointer to the
:c:struct:`SUNAdjointStepper_` structure:

.. c:type:: struct SUNAdjointStepper_ *SUNAdjointStepper

.. c:struct:: SUNAdjointStepper_

   .. c:member:: SUNStepper adj_sunstepper

      The :c:type:`SUNStepper` object used for backwards time stepping of the adjoint ODE system.

   .. c:member:: SUNStepper fwd_sunstepper

      The :c:type:`SUNStepper` object used for forward time stepping of the original ODE system if any recomputation of missing
      state data is required during the backwards integration.

   .. c:member:: sunrealtype tf

      The terminal time of the backwards adjoint ODE.

   .. c:member:: int64_t step_idx

      The index of the current backward integration step with respect to the forward integration.

   .. c:member:: int64_t final_step_idx

      The index of the final step in the forward integration (corresponds to ``tf``).

   .. c:member:: SUNMatrix Jac

      Matrix data for the Jacobian :math:`df/dy`.

   .. c:member:: SUNMatrix JacP

      Matrix data for the Jacobian :math:`df/dp`.

   .. c:member:: SUNRhsJacFn JacFn

      Jacobian function pointer to evaluate :math:`df/dy`.

   .. c:member:: SUNRhsJacFn JacPFn

      Jacobian function pointer to evaluate :math:`df/dp`.

   .. c:member:: SUNRhsJacTimesFn JvpFn

      Jacobian-times-vector function pointer to evaluate :math:`(df/dy)^T v`.

   .. c:member:: SUNRhsJacTimesFn JPvpFn

      Jacobian-times-vector function pointer to evaluate :math:`(df/dp)^T v`.

   .. c:member:: SUNRhsJacTimesFn vJpFn

      Jacobian-times-vector function pointer to evaluate :math:`v^T(df/dy)`.

   .. c:member:: SUNRhsJacTimesFn vJPpFn

      Jacobian-times-vector function pointer to evaluate :math:`v^T(df/dp)`.

   .. c:member:: int64_t nst

      Holds the count of the number of backwards steps taken.

   .. c:member:: int64_t njeval

      Holds the count of the number of :math:`df/dy` evaluations.

   .. c:member:: int64_t njpeval

      Holds the count of the number of :math:`df/dp` evaluations.

   .. c:member:: int64_t njtimesv

      Holds the count of the number of :math:`(df/dy)^T v` evaluations.

   .. c:member:: int64_t njptimesv

      Holds the count of the number of :math:`(df/dp)^T v` evaluations.

   .. c:member:: int64_t nvtimesj

      Holds the count of the number of :math:`v^T(df/dy)` evaluations.

   .. c:member:: int64_t nvtimesjp

      Holds the count of the number of :math:`v^T(df/dp)` evaluations.

   .. c:member:: int64_t nrecompute

   .. c:member:: void* user_data

      A pointer that is passed back to user-supplied functions

   .. c:member:: void* content

      Pointer to derived class specific member data

   .. c:member:: SUNContext sunctx

      The SUNDIALS simulation context


The :c:type:`SUNAdjointStepper` class has the following functions:

.. c:function:: SUNErrCode SUNAdjointStepper_Create(SUNStepper fwd_sunstepper, SUNStepper adj_sunstepper, \
   int64_t final_step_idx, N_Vector sf, sunrealtype tf, SUNAdjointCheckpointScheme checkpoint_scheme, \
   SUNContext sunctx, SUNAdjointStepper* adj_stepper)

   Creates the ``SUNAdjointStepper`` object needed to solve the adjoint problem.

   :param fwd_sunstepper: The :c:type:`SUNStepper` to be used for forward computations of the original ODE.
   :param adj_sunstepper: The :c:type:`SUNStepper` to be usef for the backward integration of the adjoint ODE.
   :param final_step_idx: The index (step number) of the step corresponding to ``t_f`` for the forward ODE.
   :param sf: The terminal condition for the adjoint ODE.
   :param tf: The terminal time for the forward ODE and (which is the initial time for the adjoint ODE).
   :param checkpoint_scheme: The :c:type:`SUNAdjointCheckpointScheme` object that determines the checkpointing strategy to use. This should be the same scheme provided to the forward integrator/stepper.
   :param sunctx: The :c:type:`SUNContext` for the simulation context.


.. c:function:: SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper adj, N_Vector sf, sunrealtype tf)

   Reinitializes the adjoint stepper to solve a new problem of the same size.

   :param adj_stepper: The adjoint solver object.
   :param sf: The terminal condition vector of sensitivity solutions :math:`dg/dy_0`` and :math:`dg/dp`.
   :param tf: The time to start integrating the adjoint system from.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper adj_stepper, sunrealtype tout,\
   N_Vector sens, sunrealtype* tret)

   Integrates the adjoint system.

   :param adj_stepper: The adjoint solver object.
   :param tout: The time at which the adjoint solution is desired.
   :param sens: The vector of sensitivity solutions :math:`dg/dy_0`` and :math:`dg/dp`.
   :param tret: On return, the time reached by the adjoint solver.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper adj_stepper, sunrealtype tout,\
   N_Vector sens, sunrealtype* tret)

   Evolves the adjoint system backwards one step.

   :param adj_stepper: The adjoint solver object.
   :param tout: The time at which the adjoint solution is desired.
   :param sens: The vector of sensitivity solutions :math:`dg/dy_0` and :math:`dg/dp`.
   :param tret: On return, the time reached by the adjoint solver.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper adj_stepper, int64_t start_idx,\
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

   Sets the function pointers and matrices needed to evaluate and store :math:`df/dy` and
   :math:`df/dp`. ``Jac`` should have dimensions ``neq x neq`` where ``neq`` is the number of states
   in the forward problem. ``JacP`` should have dimensions ``nparams x neq`` where ``nparams`` is the
   number of parameters in the model to get sensitivities for.

   :param adj_stepper: The SUNAdjointStepper object.
   :param JacFn: the function that evaluates :math:`df/dy`.
   :param Jac: a :c:type:`SUNMatrix` that will hold :math:`df/dy`.
   :param JacPFn: the function that evaluates :math:`df/dp`.
   :param JacP: a :c:type:`SUNMatrix` that will hold :math:`df/dp`.

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointStepper_SetVecTimesJacFn(SUNAdjointStepper adj_stepper, SUNRhsJacTimesFn Jvp, SUNRhsJacTimesFn JPvp)


   Sets the function pointers to evaluate :math:`(df/dy)^T v`  and :math:`(df/dp)^T v`

   :param adj_stepper: The SUNAdjointStepper object.
   :param Jvp: function that evaluates :math:`(df/dy)^T v`.
   :param JPvp: function that evaluates :math:`(df/dp)^T v`.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_SetJacTimesVecFn(SUNAdjointStepper adj_stepper, SUNRhsJacTimesFn Jvp, SUNRhsJacTimesFn JPvp)


   Sets the function pointers to evaluate :math:`v^T (df/dy)`  and :math:`v^T (df/dp)`

   :param adj_stepper: The SUNAdjointStepper object.
   :param Jvp: function that evaluates :math:`v^T (df/dy)`.
   :param JPvp: function that evaluates :math:`v^T (df/dp)`.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_SetUserData(SUNAdjointStepper adj_stepper, void* user_data)

   Sets the user data pointer.

   :param adj_stepper: The SUNAdjointStepper object.
   :param user_data: the user data pointer that will be passed back to user-supplied callback functions.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_PrintAllStats(SUNAdjointStepper adj_stepper, \
                                                           FILE* outfile, SUNOutputFormat fmt)

   Prints the adjoint stepper statistics/counters in a human-readable table format or CSV format.

   :param adj_stepper: The SUNAdjointStepper object.
   :param outfile: A file to write the output to.
   :param fmt: the format to write in (:c:type:`SUN_OUTPUTFORMAT_TABLE` or :c:type:`SUN_OUTPUTFORMAT_CSV`).

   :return: A :c:type:`SUNErrCode` indicating failure or success.

