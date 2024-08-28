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

The SUNAdjointStepper class
===========================

The :c:type:`SUNAdjointStepper` class provides a package-agnostic
interface to SUNDIALS ASA capabilties. It currently only supports the discrete
ASA capabilties in the ARKODE package, but in the future this support may be expanded.

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

   .. c:member:: SUNJacFn JacFn

      Jacobian function pointer to evaluate :math:`df/dy`.

   .. c:member:: SUNJacFn JacPFn

      Jacobian function pointer to evaluate :math:`df/dp`.

   .. c:member:: SUNJacTimesFn JvpFn

      Jacobian-times-vector function pointer to evaluate :math:`(df/dy)^T v`.

   .. c:member:: SUNJacTimesFn JPvpFn

      Jacobian-times-vector function pointer to evaluate :math:`(df/dp)^T v`.

   .. c:member:: SUNJacTimesFn vJpFn

      Jacobian-times-vector function pointer to evaluate :math:`v^T(df/dy)`.

   .. c:member:: SUNJacTimesFn vJPpFn

      Jacobian-times-vector function pointer to evaluate :math:`v^T(df/dp)`.

   .. c:member:: uint64_t nst

      Holds the count of the number of backwards steps taken.

   .. c:member:: uint64_t njeval

      Holds the count of the number of :math:`df/dy` evaluations.

   .. c:member:: uint64_t njpeval

      Holds the count of the number of :math:`df/dp` evaluations.

   .. c:member:: uint64_t njtimesv

      Holds the count of the number of :math:`(df/dy)^T v` evaluations.

   .. c:member:: uint64_t njptimesv

      Holds the count of the number of :math:`(df/dp)^T v` evaluations.

   .. c:member:: uint64_t nvtimesj

      Holds the count of the number of :math:`v^T(df/dy)` evaluations.

   .. c:member:: uint64_t nvtimesjp

      Holds the count of the number of :math:`v^T(df/dp)` evaluations.

   .. c:member:: uint64_t nrecompute

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



.. c:function:: SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper adj, N_Vector sf, sunrealtype tf)

   Reinitializes the adjoint stepper to solve a new problem of the same size.

   :param adj_stepper: The adjoint solver object.
   :param tf: The time to start integrating the adjoint system from.
   :param sf: The terminal condition vector of sensitivity solutions dg/dy0 and dg/dp.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper adj_stepper, sunrealtype tout,\
   N_Vector sens, sunrealtype* tret, int* stop_reason)

   Integrates the adjoint system.

   :param adj_stepper: The adjoint solver object.
   :param tout: The time at which the adjoint solution is desired.
   :param sens: The vector of sensitivity solutions dg/dy0 and dg/dp.
   :param tret: On return, the time reached by the adjoint solver.
   :param stop_reason: On return, an integer code that indicates why the adjoint solver stopped.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper adj_stepper, sunrealtype tout,\
   N_Vector sens, sunrealtype* tret, int* stop_reason)

   Evolves the adjoint system backwards one step.

   :param adj_stepper: The adjoint solver object.
   :param tout: The time at which the adjoint solution is desired.
   :param sens: The vector of sensitivity solutions dg/dy0 and dg/dp.
   :param tret: On return, the time reached by the adjoint solver.
   :param stop_reason: On return, an integer code that indicates why the adjoint solver stopped.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper adj_stepper, \
                                                          int64_t start_idx, int64_t stop_idx, \
                                                          sunrealtype t0, sunrealtype tf, N_Vector y0)

   Evolves the forward system in time from start_idx/t0 to stop_idx/tf with dense checkpointing.

   :param adj_stepper: The SUNAdjointStepper object.
   :param start_idx: the index of the step, w.r.t. the original forward integration, to begin forward integration from.
   :param stop_idx: the index of the step, w.r.t. the original forward integration, to end forward integration.
   :param t0: the initial time, w.r.t. the original forward integration, to start forward integration from.
   :param tf: the final time, w.r.t. the original forward integration, to stop forward integration.
   :param y0: the initial state, w.r.t. the original forward integration, to start forward integration.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_SetJacFn(SUNAdjointStepper adj_stepper, SUNJacFn JacFn, \
      SUNMatrix Jac, SUNJacFn JacPFn, SUNMatrix JacP)

   Sets the function pointers and matrices needed to evluate and store :math:`df/dy` and
   :math:`df/dp`. ``Jac`` should have dimensions ``neq x neq`` where ``neq`` is the number of states
   in the forward problem. ``JacP`` should have dimensions ``nparams x neq`` where ``nparams`` is the
   number of parameters in the model to get sensitivities for.

   :param adj_stepper: The SUNAdjointStepper object.
   :param JacFn: the function that evaluates :math:`df/dy`.
   :param Jac: a :c:type:`SUNMatrix` that will hold :math:`df/dy`.
   :param JacPFn: the function that evaluates :math:`df/dp`.
   :param JacP: a :c:type:`SUNMatrix` that will hold :math:`df/dp`.

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointStepper_SetJacTimesVecFn(SUNAdjointStepper adj_stepper, SUNJacTimesFn Jvp, SUNJacTimesFn JPvp)


   Sets the function pointers to evaluate :math:`(df/dy)^T v`  and :math:`(df/dp)^T v`

   :param adj_stepper: The SUNAdjointStepper object.
   :param Jvp: function that evaluates :math:`(df/dy)^T v`.
   :param JPvp: function that evaluates :math:`(df/dp)^T v`.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_SetJacTimesVecFn(SUNAdjointStepper adj_stepper, SUNJacTimesFn Jvp, SUNJacTimesFn JPvp)


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


.. _SUNAdjointCheckpointScheme:

The SUNAdjointCheckpointScheme Class
====================================

The :c:type:`SUNAdjointCheckpointScheme` base class provides an inteface for checkpointing
states during forward integration and accessing them as needed during the backwards integration
of the adjoint model.

A :c:type:`SUNAdjointCheckpointScheme` is a pointer to the
:c:struct:`SUNAdjointCheckpointScheme_` structure:

.. c:type:: struct SUNAdjointCheckpointScheme_ *SUNAdjointCheckpointScheme

.. c:struct:: SUNAdjointCheckpointScheme_

   .. c:member:: SUNAdjointCheckpointScheme_Ops ops

      The ops structure holds the vtable of function pointers for the base class.

   .. c:member:: void* content

      Pointer to derived class specific member data.

   .. c:member:: SUNContext sunctx

      The SUNDIALS simulation context.


.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(SUNContext sunctx, \
   SUNAdjointCheckpointScheme* cs_ptr)

   Allocates a new object but without any content.

   :param sunctx: The SUNDIALS simulation context
   :param cs_ptr: on output, the pointer to the new :c:type:`SUNAdjointCheckpointScheme` object

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave(SUNAdjointCheckpointScheme cs, \
   sunindextype step_num, sunindextype stage_num, sunrealtype t, sunbooleantype* yes_or_no)

   Determines if the (step_num, stage_num) should be checkpointed or not.

   :param cs: The :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param yes_or_no: boolean indicating if the checkpoint should be saved or not

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete(SUNAdjointCheckpointScheme cs, \
   sunindextype step_num, sunindextype stage_num, sunbooleantype* yes_or_no)

   Determines if the (step_num, stage_num) checkpoint should be deleted or not.

   :param cs: The :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param yes_or_no: boolean indicating if the checkpoint should be deleted or not

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_InsertVector(SUNAdjointCheckpointScheme cs, \
   sunindextype step_num, sunindextype stage_num, sunrealtype t, N_Vector state)

   Inserts the vector as the checkpoint for (step_num, stage_num).

   :param cs: The :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param state: the state vector to checkpoint

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_LoadVector(SUNAdjointCheckpointScheme cs, \
   sunindextype step_num, sunindextype stage_num, sunbooleantype peek, N_Vector* out, sunrealtype* tout)

   Loads the checkpointed vector for (step_num, stage_num).

   :param cs: The :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param peek: if true, then the checkpoint will be loaded but not deleted regardless
      of other implementation-specific settings. If false, then the checkpoint may be
      deleted depending on the implementation.
   :param out: the loaded state vector
   :param tout: on output, the time of the checkpoint

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_RemoveVector(SUNAdjointCheckpointScheme cs, \
   sunindextype step_num, sunindextype stage_num, N_Vector* out)

   Removes the checkpointed vector for (step_num, stage_num).

   :param cs: The :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param out: the loaded state vector

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_EnableDense(SUNAdjointCheckpointScheme cs, \
   sunbooleantype on_or_off)

   Enables or disables dense checkpointing (checkpointing every step/stage).

   :param cs: The :c:type:`SUNAdjointCheckpointScheme` object
   :param on_or_off: if true, dense checkpointing will be turned on, ifalse it will be turned off.

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_Destroy(SUNAdjointCheckpointScheme* cs_ptr)

   Destroys (deallocates) the SUNAdjointCheckpointScheme object.

   :param cs_ptr: pointer to a :c:type:`SUNAdjointCheckpointScheme` object

   :return: A :c:type:`SUNErrCode` indicating failure or success.