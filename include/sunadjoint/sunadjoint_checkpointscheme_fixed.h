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
 * SUNAdjointCheckpointScheme_Fixed class declaration.
 * ----------------------------------------------------------------*/

#ifndef _sunadjoint_checkpointscheme_fixed_H
#define _sunadjoint_checkpointscheme_fixed_H

#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_export.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Create_Fixed(
  SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, int64_t interval,
  int64_t estimate, sunbooleantype save_stages, sunbooleantype keep,
  SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave_Fixed(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  sunrealtype t, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_InsertVector_Fixed(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  sunrealtype t, N_Vector state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete_Fixed(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_RemoveVector_Fixed(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  N_Vector* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Fixed(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  sunbooleantype peek, N_Vector* out, sunrealtype* tout);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Destroy_Fixed(
  SUNAdjointCheckpointScheme* check_scheme_ptr);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_EnableDense_Fixed(
  SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off);

#ifdef __cplusplus
}
#endif

#endif /* _sunadjoint_checkpointscheme_fixed_H */
