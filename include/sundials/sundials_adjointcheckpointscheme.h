/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNAdjointCheckpointScheme class declaration.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_CHECKPOINTSCHEME_H
#define _SUNADJOINT_CHECKPOINTSCHEME_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef _SUNDIALS_STRUCT_ SUNAdjointCheckpointScheme_Ops_* SUNAdjointCheckpointScheme_Ops;
typedef _SUNDIALS_STRUCT_ SUNAdjointCheckpointScheme_* SUNAdjointCheckpointScheme;

struct SUNAdjointCheckpointScheme_Ops_
{
  SUNErrCode (*needsSaving)(SUNAdjointCheckpointScheme, int64_t step_num,
                            int64_t stage_num, sunrealtype t,
                            sunbooleantype* yes_or_no);

  SUNErrCode (*needsDeleting)(SUNAdjointCheckpointScheme, int64_t step_num,
                              int64_t stage_num, sunbooleantype* yes_or_no);

  SUNErrCode (*insertVector)(SUNAdjointCheckpointScheme, int64_t step_num,
                             int64_t stage_num, sunrealtype t, N_Vector state);

  SUNErrCode (*loadVector)(SUNAdjointCheckpointScheme, int64_t step_num,
                           int64_t stage_num, sunbooleantype peek,
                           N_Vector* out, sunrealtype* tout);

  SUNErrCode (*removeVector)(SUNAdjointCheckpointScheme, int64_t step_num,
                             int64_t stage_num, N_Vector* out);

  SUNErrCode (*destroy)(SUNAdjointCheckpointScheme*);

  SUNErrCode (*enableDense)(SUNAdjointCheckpointScheme, sunbooleantype on_or_off);
};

struct SUNAdjointCheckpointScheme_
{
  SUNAdjointCheckpointScheme_Ops ops;
  void* content;
  SUNContext sunctx;
};

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(SUNContext sunctx,
                                               SUNAdjointCheckpointScheme*);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave(SUNAdjointCheckpointScheme,
                                                   int64_t step_num,
                                                   int64_t stage_num,
                                                   sunrealtype t,
                                                   sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete(SUNAdjointCheckpointScheme,
                                                     int64_t step_num,
                                                     int64_t stage_num,
                                                     sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_InsertVector(SUNAdjointCheckpointScheme,
                                                   int64_t step_num,
                                                   int64_t stage_num,
                                                   sunrealtype t, N_Vector state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_LoadVector(
  SUNAdjointCheckpointScheme, int64_t step_num, int64_t stage_num,
  sunbooleantype peek, N_Vector* out, sunrealtype* tout);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_RemoveVector(SUNAdjointCheckpointScheme,
                                                   int64_t step_num,
                                                   int64_t stage_num,
                                                   N_Vector* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Destroy(SUNAdjointCheckpointScheme*);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_EnableDense(SUNAdjointCheckpointScheme,
                                                  sunbooleantype on_or_off);

#ifdef __cplusplus
}
#endif

#endif /*_SUNADJOINT_CHECKPOINTSCHEME_H*/
