.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdjoint.CheckpointScheme:

The SUNAdjointCheckpointScheme API
==================================

.. versionadded:: x.y.z

The :c:type:`SUNAdjointCheckpointScheme` base class provides an interface for checkpointing
states during forward integration and accessing them as needed during the backwards integration
of the adjoint model.

.. c:enum:: SUNDataIOMode

   .. c:enumerator:: SUNDATAIOMODE_INMEM

      The IO mode for data that is stored in addressable random access memory.
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


The virtual table structure is defined as

.. c:type:: struct SUNAdjointCheckpointScheme_Ops_ *SUNAdjointCheckpointScheme_Ops


.. c:struct:: SUNAdjointCheckpointScheme_Ops_

   .. c:member:: SUNErrCode (*needssaving)(SUNAdjointCheckpointScheme cs, suncountertype step_num, suncountertype stage_num, sunrealtype t, sunbooleantype* yes_or_no)

      Function pointer to determine if a checkpoint should be saved at the current time step.

   .. c:member:: SUNErrCode (*needsdeleting)(SUNAdjointCheckpointScheme cs, suncountertype step_num, suncountertype stage_num, sunrealtype t, sunbooleantype* yes_or_no)

      Function pointer to determine if a checkpoint should be deleted at the current time step.

   .. c:member:: SUNErrCode (*insertvector)(SUNAdjointCheckpointScheme cs, suncountertype step_num, suncountertype stage_num, sunrealtype t, N_Vector y)

      Function pointer to insert a checkpoint state represented as a :c:type:`N_Vector`.

   .. c:member:: SUNErrCode (*loadvector)(SUNAdjointCheckpointScheme cs, suncountertype step_num, suncountertype stage_num, sunrealtype t, sunbooleantype peek, N_Vector* yout, sunrealtype* tout)

      Function pointer to load a checkpoint state represented as a :c:type:`N_Vector`.

   .. c:member:: SUNErrCode (*deleteVector)(SUNAdjointCheckpointScheme cs, suncountertype step_num, suncountertype stage_num, N_Vector* out)

      Function pointer to remove a checkpoint state represented as a :c:type:`N_Vector`.

   .. c:member:: SUNErrCode (*destroy)(SUNAdjointCheckpointScheme*)

      Function pointer to destroy and free the memory for the :c:type:`SUNAdjointCheckpointScheme` object.

   .. c:member:: SUNErrCode (*enableDense)(SUNAdjointCheckpointScheme cs, sunbooleantype on_or_off)

      Function pointer to enable or disable dense checkpointing, saving all steps.


.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(SUNContext sunctx, \
   SUNAdjointCheckpointScheme* cs_ptr)

   :param sunctx: The SUNDIALS simulation context
   :param cs_ptr: on output, the pointer to the new :c:type:`SUNAdjointCheckpointScheme` object

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_NeedsSaving(SUNAdjointCheckpointScheme cs, \
   suncountertype step_num, suncountertype stage_num, sunrealtype t, sunbooleantype* yes_or_no)

   Determines if the (step_num, stage_num) should be checkpointed or not.

   :param cs: the :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param yes_or_no: boolean indicating if the checkpoint should be saved or not

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_NeedsDeleting(SUNAdjointCheckpointScheme cs, \
   suncountertype step_num, suncountertype stage_num, sunbooleantype* yes_or_no)

   Determines if the (step_num, stage_num) checkpoint should be deleted or not.

   :param cs: the :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param yes_or_no: boolean indicating if the checkpoint should be deleted or not

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_InsertVector(SUNAdjointCheckpointScheme cs, \
   suncountertype step_num, suncountertype stage_num, sunrealtype t, N_Vector y)

   Inserts the vector as the checkpoint for (step_num, stage_num).

   :param cs: the :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param y: the state vector to checkpoint

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_LoadVector(SUNAdjointCheckpointScheme cs, \
   suncountertype step_num, suncountertype stage_num, sunrealtype t, sunbooleantype peek, N_Vector* yout, sunrealtype* tout)

   Loads the checkpointed vector for (step_num, stage_num).

   :param cs: the :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the desired time of the checkpoint
   :param peek: if true, then the checkpoint will be loaded but not deleted regardless
      of other implementation-specific settings. If false, then the checkpoint may be
      deleted depending on the implementation.
   :param yout: the loaded state vector
   :param tout: on output, the time of the checkpoint

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_RemoveVector(SUNAdjointCheckpointScheme cs, \
   suncountertype step_num, suncountertype stage_num, N_Vector* out)

   Removes the checkpointed vector for (step_num, stage_num).

   :param cs: the :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param out: the loaded state vector

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_EnableDense(SUNAdjointCheckpointScheme cs, \
   sunbooleantype on_or_off)

   Enables or disables dense checkpointing (checkpointing every step/stage). When dense checkpointing
   is disabled, the checkpointing interval that was set when the object was created is restored.

   :param cs: the :c:type:`SUNAdjointCheckpointScheme` object
   :param on_or_off: if true, dense checkpointing will be turned on, if false it will be turned off.

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_Destroy(SUNAdjointCheckpointScheme* cs_ptr)

   Destroys (deallocates) the SUNAdjointCheckpointScheme object.

   :param cs_ptr: pointer to a :c:type:`SUNAdjointCheckpointScheme` object

   :returns: A :c:type:`SUNErrCode` indicating failure or success.


.. _SUNAdjoint.CheckpointScheme.Fixed:

The SUNAdjointCheckpointScheme_Fixed Module
===========================================

The ``SUNAdjointCheckpointScheme_Fixed`` module implements a scheme where a checkpoint is saved at some
fixed interval (in time steps). The module supports checkpointing of time step states only, or time step
states with intermediate stage states as well (for multistage methods). When used with a
fixed time step size then the number of checkpoints that will be saved is fixed. However, with
adaptive time steps the number of checkpoints stored with this scheme is unbounded.

The diagram below illustrates how checkpoints are stored with this scheme:

.. figure:: /figs/sunadjoint_ckpt_fixed.png
   :width: 75 %
   :align: center


Base-class Method Overrides
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``SUNAdjointCheckpointScheme_Fixed`` module implements the following :c:type:`SUNAdjointCheckpointScheme` functions:

* :c:func:`SUNAdjointCheckpointScheme_NeedsSaving`
* :c:func:`SUNAdjointCheckpointScheme_InsertVector`
* :c:func:`SUNAdjointCheckpointScheme_NeedsDeleting`
* :c:func:`SUNAdjointCheckpointScheme_RemoveVector`
* :c:func:`SUNAdjointCheckpointScheme_LoadVector`
* :c:func:`SUNAdjointCheckpointScheme_Destroy`
* :c:func:`SUNAdjointCheckpointScheme_EnableDense`


Implementation Specific Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``SUNAdjointCheckpointScheme_Fixed`` module also implements the following module-specific functions:

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_Create_Fixed(SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, suncountertype interval, suncountertype estimate, sunbooleantype save_stages, sunbooleantype keep, SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr)

   Creates a new :c:type:`SUNAdjointCheckpointScheme` object that checkpoints at a fixed interval.

   :param io_mode: The IO mode used for storing the checkpoints.
   :param mem_helper: Memory helper for managing memory.
   :param interval: The interval (in steps) between checkpoints.
   :param estimate: An estimate of the total number of checkpoints needed.
   :param save_stages: If using a multistage method, should stages be saved with the step.
   :param keep: Keep data stored even after it is not needed anymore.
   :param sunctx: The :c:type:`SUNContext` for the simulation.
   :param check_scheme_ptr: Pointer to the newly constructed object.
   :returns: A :c:type:`SUNErrCode` indicating success or failure.
