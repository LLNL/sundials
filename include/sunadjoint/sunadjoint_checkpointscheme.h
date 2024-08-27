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
 * SUNAdjointCheckpointScheme class declaration.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_CHECKPOINTSCHEME_H
#define _SUNADJOINT_CHECKPOINTSCHEME_H

#include <stdint.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_datanode.h>
#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef struct SUNAdjointCheckpointScheme_Ops_* SUNAdjointCheckpointScheme_Ops;
typedef struct SUNAdjointCheckpointScheme_* SUNAdjointCheckpointScheme;

struct SUNAdjointCheckpointScheme_Ops_
{
  SUNErrCode (*shouldWeSave)(SUNAdjointCheckpointScheme, sunindextype step_num,
                             sunindextype stage_num, sunrealtype t,
                             sunbooleantype* yes_or_no);

  SUNErrCode (*shouldWeDelete)(SUNAdjointCheckpointScheme, sunindextype step_num,
                               sunindextype stage_num, sunbooleantype* yes_or_no);

  SUNErrCode (*insertVector)(SUNAdjointCheckpointScheme, sunindextype step_num,
                             sunindextype stage_num, sunrealtype t,
                             N_Vector state);

  SUNErrCode (*loadVector)(SUNAdjointCheckpointScheme, sunindextype step_num,
                           sunindextype stage_num, sunbooleantype peek,
                           N_Vector* out, sunrealtype* tout);

  SUNErrCode (*removeVector)(SUNAdjointCheckpointScheme, sunindextype step_num,
                             sunindextype stage_num, N_Vector* out);

  SUNErrCode (*removeRange)(SUNAdjointCheckpointScheme,
                            sunindextype step_num_start,
                            sunindextype step_num_end,
                            sunindextype stage_num_start,
                            sunindextype stage_num_end);

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
                                                   sunindextype step_num,
                                                   sunindextype stage_num,
                                                   sunrealtype t,
                                                   sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete(SUNAdjointCheckpointScheme,
                                                     sunindextype step_num,
                                                     sunindextype stage_num,
                                                     sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_InsertVector(SUNAdjointCheckpointScheme,
                                                   sunindextype step_num,
                                                   sunindextype stage_num,
                                                   sunrealtype t, N_Vector state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_LoadVector(
  SUNAdjointCheckpointScheme, sunindextype step_num, sunindextype stage_num,
  sunbooleantype peek, N_Vector* out, sunrealtype* tout);

SUNErrCode SUNAdjointCheckpointScheme_RemoveVector(SUNAdjointCheckpointScheme,
                                                   sunindextype step_num,
                                                   sunindextype stage_num,
                                                   N_Vector* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_RemoveRange(SUNAdjointCheckpointScheme,
                                                  sunindextype step_num_start,
                                                  sunindextype step_num_end,
                                                  sunindextype stage_num_start,
                                                  sunindextype stage_num_end);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Destroy(SUNAdjointCheckpointScheme*);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_EnableDense(SUNAdjointCheckpointScheme,
                                                  sunbooleantype on_or_off);

#ifdef __cplusplus
}
#endif

#endif /*_SUNADJOINT_CHECKPOINTSCHEME_H*/
