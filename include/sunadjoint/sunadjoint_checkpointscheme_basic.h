/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNAdjointCheckpointScheme_Basic class declaration.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_CHECKPOINTSCHEME_BASIC_H
#define _SUNADJOINT_CHECKPOINTSCHEME_BASIC_H

#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_datanode.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_export.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/**

  This function creates a new SUNAdjointCheckPointScheme object that checkpoints at a fixed interval.

  :param io_mode: the IO mode that will be used for storing the checkpoints
  :param interval: the interval (in steps) between checkpoints
  :param estimate: an estimate of the total number of checkpoints needed - underestimating may result in slower performance but overestimating may result in unnecessarily high memory usage
  :param save_stages: if using a multistage method, should stages be saved with the step
  :param keep: keep data stored even after it is not needed anymore
  :param sunctx: the SUNContext for the simulation
  :param check_scheme_ptr: the newly constructed object

  :returns: a :c:type:`SUNErrCode` indicating success or failure

 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Create_Basic(
  SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, uint64_t interval,
  uint64_t estimate, sunbooleantype save_stages, sunbooleantype keep,
  SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr);

/**

  This function queries the checkpointing scheme to determine if you should save a checkpoint at this timestep.

  :param check_scheme: the SUNAdjointCheckpointScheme object
  :param step_num: the current time step number
  :param stage_num: the current stage number (only nonzero for multistage methods)
  :param t: the current time
  :param yes_or_no: on output, will be 1 if you should save, 0 otherwise.

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, sunrealtype t, sunbooleantype* yes_or_no);

/**
  This function inserts a checkpoint state (represented as a SUNDataNode).

  :param check_scheme: the SUNAdjointCheckpointScheme object
  :param step_num: the current time step number
  :param stage_num: the current stage number (only nonzero for multistage methods)
  :param t: the current time
  :param state: a :c:type:`SUNDataNode` object that holds the current state to be inserted

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Insert_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, sunrealtype t, SUNDataNode state);

/**
  This function inserts a checkpoint state (represented as a N_Vector).

  :param check_scheme: the SUNAdjointCheckpointScheme object
  :param step_num: the current time step number
  :param stage_num: the current stage number (only nonzero for multistage methods)
  :param t: the current time
  :param state: a :c:type:`N_Vector` object that holds the current state to be inserted

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_InsertVector_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, sunrealtype t, N_Vector state);

/**
  This function queries the checkpointing scheme to determine if you should delete a checkpoint at this timestep.

  :param check_scheme: the SUNAdjointCheckpointScheme object
  :param step_num: the current time step number
  :param stage_num: the current stage number (only nonzero for multistage methods)
  :param t: the current time
  :param yes_or_no: on output, will be 1 if you should delete, 0 otherwise.

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, sunbooleantype* yes_or_no);

/**
  This function removes a checkpoint state (represented as a SUNDataNode).

  :param check_scheme: the SUNAdjointCheckpointScheme object
  :param step_num: the current time step number
  :param stage_num: the current stage number (only nonzero for multistage methods)
  :param t: the current time
  :param state: a :c:type:`SUNDataNode` object that holds the current state to be removed

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Remove_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, SUNDataNode* out);

/**
  This function removes a checkpoint state (represented as a N_Vector).

  :param check_scheme: the SUNAdjointCheckpointScheme object
  :param step_num: the current time step number
  :param stage_num: the current stage number (only nonzero for multistage methods)
  :param t: the current time
  :param state: a :c:type:`N_Vector` object that holds the current state to be removed

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_RemoveVector_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, N_Vector* out);

/**
  This function removes multiple checkpoints.

  :param check_scheme: the SUNAdjointCheckpointScheme object

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_RemoveRange_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num_start,
  sunindextype step_num_end, sunindextype stage_num_start,
  sunindextype stage_num_end);

/**
  This function loads a checkpoint state (represented as a N_Vector).

  :param check_scheme: the SUNAdjointCheckpointScheme object
  :param step_num: the current time step number
  :param stage_num: the current stage number (only nonzero for multistage methods)
  :param t: the current time
  :param state: a :c:type:`N_Vector` object that holds loaded state

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Load_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, SUNDataNode* out);

/**
  This function loads a checkpoint state (represented as a N_Vector).

  :param check_scheme: the SUNAdjointCheckpointScheme object
  :param step_num: the current time step number
  :param stage_num: the current stage number (only nonzero for multistage methods)
  :param t: the current time
  :param state: a :c:type:`N_Vector` object that holds loaded state

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, N_Vector* out);

/**
  This function destroys and frees the memory for the SUNAdjointCheckpointScheme object.

  :param check_scheme_ptr: the SUNAdjointCheckpointScheme object

  :returns: a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Destroy_Basic(
  SUNAdjointCheckpointScheme* check_scheme_ptr);

/**
  This function enables dense checkpointing, that is it saves
  all steps. 

  :param check_scheme: the SUNAdjointCheckpointScheme object
  :param on_or_off: turn dense checkpoints on or off

  :returns a :c:type:`SUNErrCode` indicating success or failure
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_EnableDense_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off);

#ifdef __cplusplus
}
#endif

#endif /* _SUNADJOINT_CHECKPOINTSCHEME_BASIC_H */
