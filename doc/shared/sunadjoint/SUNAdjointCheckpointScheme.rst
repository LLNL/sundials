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

.. c:enumerator:: SUNDATAIOMODE_INMEM

   The IO mode for data that is stored in addressible random access memory.
   The location of the memory (e.g., CPU or GPU) is not specified by this mode.


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

The SUNAdjointCheckpointScheme_Fixed Module
===========================================

The SUNAdjointCheckpointScheme_Fixed module implements a scheme where a checkpoint is saved at some
fixed interval (in timesteps). The module supports checkpointing of time step states only, or timestep
states with intermediate stage states as well (for multistage methods). When used with a
fixed timestep size then the number of checkpoints that will be saved is fixed. However, with
adaptive timesteps the number of checkpoints stored with this scheme is unbounded.

The diagram below illustrates how checkpoints are stored with this scheme:

.. figure:: /figs/sunadjoint_ckpt_fixed.png
   :width: 75 %


Base-class Method Overrides
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``SUNAdjointCheckpointScheme_Fixed`` module implements the following :c:type:`SUNAdjointCheckpointScheme` functions:

* :c:func:`SUNAdjointCheckpointScheme_ShouldWeSave`
* :c:func:`SUNAdjointCheckpointScheme_InsertVector`
* :c:func:`SUNAdjointCheckpointScheme_ShouldWeDelete`
* :c:func:`SUNAdjointCheckpointScheme_RemoveVector`
* :c:func:`SUNAdjointCheckpointScheme_LoadVector`
* :c:func:`SUNAdjointCheckpointScheme_Destroy`
* :c:func:`SUNAdjointCheckpointScheme_EnableDense`


Implementation Specific Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``SUNAdjointCheckpointScheme_Fixed`` module also implements the following module-specific functions:

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_Create_Fixed(SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, int64_t interval, int64_t estimate, sunbooleantype save_stages, sunbooleantype keep, SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr)

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
