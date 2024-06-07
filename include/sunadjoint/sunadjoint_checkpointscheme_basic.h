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
  :param SUNAdjointCheckpointScheme: the newly constructed object

  :returns: a :c:type:`SUNErrCode` indicating success or failure

 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Create_Basic(
  SUNDataIOMode io_mod, uint64_t interval, uint64_t estimate,
  sunbooleantype save_stages, sunbooleantype keep, SUNContext sunctx,
  SUNAdjointCheckpointScheme*);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave_Basic(
  SUNAdjointCheckpointScheme, sunindextype step_num, sunindextype stage_num,
  sunrealtype t, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Insert_Basic(SUNAdjointCheckpointScheme,
                                                   sunindextype step_num,
                                                   sunindextype stage_num,
                                                   sunrealtype t,
                                                   SUNDataNode state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_InsertVector_Basic(
  SUNAdjointCheckpointScheme, sunindextype step_num, sunindextype stage_num,
  sunrealtype t, N_Vector state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete_Basic(
  SUNAdjointCheckpointScheme, sunindextype step_num, sunindextype stage_num,
  sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Remove_Basic(SUNAdjointCheckpointScheme,
                                                   sunindextype step_num,
                                                   sunindextype stage_num,
                                                   SUNDataNode* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_RemoveVector_Basic(SUNAdjointCheckpointScheme,
                                                         sunindextype step_num,
                                                         sunindextype stage_num,
                                                         N_Vector* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_RemoveRange_Basic(
  SUNAdjointCheckpointScheme, sunindextype step_num_start,
  sunindextype step_num_end, sunindextype stage_num_start,
  sunindextype stage_num_end);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Load_Basic(SUNAdjointCheckpointScheme,
                                                 sunindextype step_num,
                                                 sunindextype stage_num,
                                                 SUNDataNode* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Basic(SUNAdjointCheckpointScheme,
                                                       sunindextype step_num,
                                                       sunindextype stage_num,
                                                       N_Vector* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Destroy_Basic(SUNAdjointCheckpointScheme*);

#ifdef __cpluplus
}
#endif

#endif /* _SUNADJOINT_CHECKPOINTSCHEME_BASIC_H */
