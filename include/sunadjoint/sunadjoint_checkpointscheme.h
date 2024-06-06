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

#include <sundials/sundials_core.h>
#include <sundials/sundials_datanode.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef struct SUNAdjointCheckpointScheme_Ops_* SUNAdjointCheckpointScheme_Ops;
typedef struct SUNAdjointCheckpointScheme_* SUNAdjointCheckpointScheme;

struct SUNAdjointCheckpointScheme_Ops_
{
  // This function lets the caller know if they should checkpoint the state
  // at (step_num, stage_num).
  SUNErrCode (*shouldWeSave)(SUNAdjointCheckpointScheme, sunindextype step_num,
                             sunindextype stage_num, sunrealtype t,
                             sunbooleantype* yes_or_no);

  // This function inserts the state (step_num, stage_num) in the checkpoint list.
  SUNErrCode (*insert)(SUNAdjointCheckpointScheme, sunindextype step_num,
                       sunindextype stage_num, sunrealtype t, SUNDataNode state);

  SUNErrCode (*insertVector)(SUNAdjointCheckpointScheme, sunindextype step_num,
                             sunindextype stage_num, sunrealtype t,
                             N_Vector state);

  // This function lets the caller know if they should remove the checkpoint state.
  SUNErrCode (*shouldWeDelete)(SUNAdjointCheckpointScheme, sunindextype step_num,
                               sunindextype stage_num, sunbooleantype* yes_or_no);

  // This function removes the state (step_num, stage_num) from the list.
  // Optionally, the removed state will be in the SUNDataNode 'out'.
  SUNErrCode (*remove)(SUNAdjointCheckpointScheme, sunindextype step_num,
                       sunindextype stage_num, SUNDataNode* out);

  SUNErrCode (*removeVector)(SUNAdjointCheckpointScheme, sunindextype step_num,
                             sunindextype stage_num, N_Vector* out);

  // This function removes the states in the given range from the list.
  // This is primarily useful for removing stage states when a step is unsuccessful.
  SUNErrCode (*removeRange)(SUNAdjointCheckpointScheme,
                            sunindextype step_num_start,
                            sunindextype step_num_end,
                            sunindextype stage_num_start,
                            sunindextype stage_num_end);

  // This function loads the state (step_num, stage_num) into 'out'.
  // Passing in 0 for stage_num will return the step vector,
  // but data for the stages could be fetched too depending on the implementation.
  SUNErrCode (*load)(SUNAdjointCheckpointScheme, sunindextype step_num,
                     sunindextype stage_num, SUNDataNode* out);

  SUNErrCode (*loadVector)(SUNAdjointCheckpointScheme, sunindextype step_num,
                           sunindextype stage_num, N_Vector* out);

  SUNErrCode (*destroy)(SUNAdjointCheckpointScheme*);
};

struct SUNAdjointCheckpointScheme_
{
  SUNAdjointCheckpointScheme_Ops ops;
  SUNDataNode root_data_node;
  void* content;
  SUNContext sunctx;
};

SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(SUNContext sunctx,
                                               SUNAdjointCheckpointScheme*);

SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave(SUNAdjointCheckpointScheme,
                                                   sunindextype step_num,
                                                   sunindextype stage_num,
                                                   sunrealtype t,
                                                   sunbooleantype* yes_or_no);

SUNErrCode SUNAdjointCheckpointScheme_Insert(SUNAdjointCheckpointScheme,
                                             sunindextype step_num,
                                             sunindextype stage_num,
                                             sunrealtype t, SUNDataNode state);

SUNErrCode SUNAdjointCheckpointScheme_InsertVector(SUNAdjointCheckpointScheme,
                                                   sunindextype step_num,
                                                   sunindextype stage_num,
                                                   sunrealtype t, N_Vector state);

SUNErrCode SUNAdjointCheckpointScheme_ShouldWeDelete(SUNAdjointCheckpointScheme,
                                                     sunindextype step_num,
                                                     sunindextype stage_num,
                                                     sunbooleantype* yes_or_no);

SUNErrCode SUNAdjointCheckpointScheme_Remove(SUNAdjointCheckpointScheme,
                                             sunindextype step_num,
                                             sunindextype stage_num,
                                             SUNDataNode* out);

SUNErrCode SUNAdjointCheckpointScheme_Load(SUNAdjointCheckpointScheme,
                                           sunindextype step_num,
                                           sunindextype stage_num,
                                           SUNDataNode* out);

SUNErrCode SUNAdjointCheckpointScheme_LoadVector(SUNAdjointCheckpointScheme,
                                                 sunindextype step_num,
                                                 sunindextype stage_num,
                                                 N_Vector* out);

SUNErrCode SUNAdjointCheckpointScheme_Destroy(SUNAdjointCheckpointScheme*);

#ifdef __cplusplus
}
#endif

#endif /*_SUNADJOINT_CHECKPOINTSCHEME_H*/
