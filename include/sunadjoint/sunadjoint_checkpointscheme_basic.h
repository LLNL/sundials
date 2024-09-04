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
#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Create_Basic(
  SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, uint64_t interval,
  uint64_t estimate, sunbooleantype save_stages, sunbooleantype keep,
  SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, sunrealtype t, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Insert_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, sunrealtype t, SUNDataNode state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_InsertVector_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, sunrealtype t, N_Vector state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Remove_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, SUNDataNode* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_RemoveVector_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, N_Vector* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Load_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, SUNDataNode* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunindextype step_num,
  sunindextype stage_num, sunbooleantype peek, N_Vector* out, sunrealtype* tout);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Destroy_Basic(
  SUNAdjointCheckpointScheme* check_scheme_ptr);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_EnableDense_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off);

#ifdef __cplusplus
}
#endif

#endif /* _SUNADJOINT_CHECKPOINTSCHEME_BASIC_H */
