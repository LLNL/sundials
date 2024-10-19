.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdjointCheckpointScheme:

The SUNAdjointCheckpointScheme API
==================================

The :c:type:`SUNAdjointCheckpointScheme` base class provides an interface for checkpointing
states during forward integration and accessing them as needed during the backwards integration
of the adjoint model.

.. c:enum:: SUNDataIOMode

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


.. c:type:: struct SUNAdjointCheckpointScheme_Ops_ *SUNAdjointCheckpointScheme_Ops


.. c:struct:: SUNAdjointCheckpointScheme_Ops_

   .. c:member:: SUNErrCode (*shouldWeSave)(SUNAdjointCheckpointScheme cs, int64_t step_num, int64_t stage_num, sunrealtype t, sunbooleantype* yes_or_no)

      Function pointer to determine if a checkpoint should be saved at the current timestep.

   .. c:member:: SUNErrCode (*shouldWeDelete)(SUNAdjointCheckpointScheme cs, int64_t step_num, int64_t stage_num, sunbooleantype* yes_or_no)

      Function pointer to determine if a checkpoint should be deleted at the current timestep.

   .. c:member:: SUNErrCode (*insertVector)(SUNAdjointCheckpointScheme cs, int64_t step_num, int64_t stage_num, sunrealtype t, N_Vector state)

      Function pointer to insert a checkpoint state represented as a :c:type:`N_Vector`.

   .. c:member:: SUNErrCode (*loadVector)(SUNAdjointCheckpointScheme cs, int64_t step_num, int64_t stage_num, sunbooleantype peek, N_Vector* out, sunrealtype* tout)

      Function pointer to load a checkpoint state represented as a :c:type:`N_Vector`.

   .. c:member:: SUNErrCode (*removeVector)(SUNAdjointCheckpointScheme cs, int64_t step_num, int64_t stage_num, N_Vector* out)

      Function pointer to remove a checkpoint state represented as a :c:type:`N_Vector`.

   .. c:member:: SUNErrCode (*destroy)(SUNAdjointCheckpointScheme*)

      Function pointer to destroy and free the memory for the :c:type:`SUNAdjointCheckpointScheme` object.

   .. c:member:: SUNErrCode (*enableDense)(SUNAdjointCheckpointScheme cs, sunbooleantype on_or_off)

      Function pointer to enable or disable dense checkpointing, saving all steps.


.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(SUNContext sunctx, \
   SUNAdjointCheckpointScheme* cs_ptr)

   Allocates a new object but without any content.

   :param sunctx: The SUNDIALS simulation context
   :param cs_ptr: on output, the pointer to the new :c:type:`SUNAdjointCheckpointScheme` object

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave(SUNAdjointCheckpointScheme cs, \
   int64_t step_num, int64_t stage_num, sunrealtype t, sunbooleantype* yes_or_no)

   Determines if the (step_num, stage_num) should be checkpointed or not.

   :param cs: The :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param yes_or_no: boolean indicating if the checkpoint should be saved or not

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete(SUNAdjointCheckpointScheme cs, \
   int64_t step_num, int64_t stage_num, sunbooleantype* yes_or_no)

   Determines if the (step_num, stage_num) checkpoint should be deleted or not.

   :param cs: The :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param yes_or_no: boolean indicating if the checkpoint should be deleted or not

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_InsertVector(SUNAdjointCheckpointScheme cs, \
   int64_t step_num, int64_t stage_num, sunrealtype t, N_Vector state)

   Inserts the vector as the checkpoint for (step_num, stage_num).

   :param cs: The :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param state: the state vector to checkpoint

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_LoadVector(SUNAdjointCheckpointScheme cs, \
   int64_t step_num, int64_t stage_num, sunbooleantype peek, N_Vector* out, sunrealtype* tout)

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
   int64_t step_num, int64_t stage_num, N_Vector* out)

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
   :param on_or_off: if true, dense checkpointing will be turned on, if false it will be turned off.

   :return: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_Destroy(SUNAdjointCheckpointScheme* cs_ptr)

   Destroys (deallocates) the SUNAdjointCheckpointScheme object.

   :param cs_ptr: pointer to a :c:type:`SUNAdjointCheckpointScheme` object

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. _SUNAdjointCheckpointScheme.Basic:

The SUNAdjointCheckpointScheme_Basic Module
===========================================

The SUNAdjointCheckpointScheme_Basic module implements a scheme where a checkpoint is saved at some
fixed interval (in timesteps). The module supports checkpointing of time step states only, or time
step stages with intermediate stage states as well (for multistage methods). When used with a
fixed timestep size then the number of checkpoints that will be saved is fixed. However, with
adaptive timesteps the number of checkpoints stored with this scheme is unbounded.

The diagram below illustrates how checkpoints are stored with this scheme:



The SUNAdjointCheckpointScheme_Basic module has the following user-callable functions:

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_Create_Basic(SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, int64_t interval, int64_t estimate, sunbooleantype save_stages, sunbooleantype keep, SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr)

   Creates a new :c:type:`SUNAdjointCheckpointScheme` object that checkpoints at a fixed interval.

   :param io_mode: The IO mode used for storing the checkpoints.
   :param mem_helper: Memory helper for managing memory.
   :param interval: The interval (in steps) between checkpoints.
   :param estimate: An estimate of the total number of checkpoints needed.
   :param save_stages: If using a multistage method, should stages be saved with the step.
   :param keep: Keep data stored even after it is not needed anymore.
   :param sunctx: The :c:type:`SUNContext` for the simulation.
   :param check_scheme_ptr: Pointer to the newly constructed object.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave_Basic(SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num, sunrealtype t, sunbooleantype* yes_or_no)

   Queries the checkpointing scheme to determine if a checkpoint should be saved at this timestep.

   :param check_scheme: The `SUNAdjointCheckpointScheme` object.
   :param step_num: The current time step number.
   :param stage_num: The current stage number (only nonzero for multistage methods).
   :param t: The current time.
   :param yes_or_no: On output, will be 1 if you should save, 0 otherwise.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_InsertVector_Basic(SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num, sunrealtype t, N_Vector state)

   Inserts a checkpoint state represented as a `N_Vector`.

   :param check_scheme: The `SUNAdjointCheckpointScheme` object.
   :param step_num: The current time step number.
   :param stage_num: The current stage number (only nonzero for multistage methods).
   :param t: The current time.
   :param state: A `N_Vector` object that holds the current state to be inserted.
   :return: A `SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete_Basic(SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num, sunrealtype t, sunbooleantype* yes_or_no)

   Queries the checkpointing scheme to determine if a checkpoint should be deleted at this timestep.

   :param check_scheme: The `SUNAdjointCheckpointScheme` object.
   :param step_num: The current time step number.
   :param stage_num: The current stage number (only nonzero for multistage methods).
   :param t: The current time.
   :param yes_or_no: On output, will be 1 if you should delete, 0 otherwise.
   :return: A `SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_RemoveVector_Basic(SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num, N_Vector* out)

   Removes a checkpoint state represented as a `N_Vector`.

   :param check_scheme: The `SUNAdjointCheckpointScheme` object.
   :param step_num: The current time step number.
   :param stage_num: The current stage number (only nonzero for multistage methods).
   :param out: Pointer to the `N_Vector` object that holds the current state to be removed.
   :return: A `SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Basic(SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num, sunbooleantype peek, N_Vector* out, sunrealtype* tout)

   Loads a checkpoint state represented as a `N_Vector`.

   :param check_scheme: The `SUNAdjointCheckpointScheme` object.
   :param step_num: The current time step number.
   :param stage_num: The current stage number (only nonzero for multistage methods).
   :param peek: Load the checkpointed vector without removing it regardless of the "keep" setting.
   :param out: Pointer to the `N_Vector` object that holds loaded state.
   :param tout: Pointer to the time associated with the loaded state.
   :return: A `SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_Destroy_Basic(SUNAdjointCheckpointScheme* check_scheme_ptr)

   Destroys and frees the memory for the `SUNAdjointCheckpointScheme` object.

   :param check_scheme_ptr: Pointer to the `SUNAdjointCheckpointScheme` object.
   :return: A `SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_EnableDense_Basic(SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off)

   Enables dense checkpointing, saving all steps.

   :param check_scheme: The `SUNAdjointCheckpointScheme` object.
   :param on_or_off: Turn dense checkpoints on or off.
   :return: A `SUNErrCode` indicating success or failure.