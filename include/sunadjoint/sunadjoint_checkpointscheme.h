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
 * SUNAdjointCheckpointScheme class definition.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_CHECKPOINTSCHEME_H
#define _SUNADJOINT_CHECKPOINTSCHEME_H

#include <sundials/sundials_core.h>
#include <sundials/sundials_datanode.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

struct SUNAdjointCheckpointScheme_s
{
  // This function lets the caller know if they should checkpoint the state
  // at (step_num, stage_num).
  SUNErrCode (*shouldWeSave)(SUNAdjointCheckpointScheme, sunindextype step_num,
                             sunindextype stage_num, sunrealtype t,
                             sunbooleantype* yes_or_no);

  // This function inserts the state (step_num, stage_num) in the checkpoint list.
  SUNErrCode (*insert)(SUNAdjointCheckpointScheme, sunindextype step_num,
                       sunindextype stage_num, sunrealtype t, N_Vector state);

  // This function lets the caller know if they should remove the checkpoint state.
  SUNErrCode (*shouldWeDelete)(SUNAdjointCheckpointScheme, sunindextype step_num,
                               sunindextype stage_num, sunbooleantype* yes_or_no);

  // This function removes the state (step_num, stage_num) from the list.
  // Optionally, the removed state will be in the vector 'out'.
  SUNErrCode (*remove)(SUNAdjointCheckpointScheme, sunindextype step_num,
                       sunindextype stage_num, SUNDataNode* out);

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
                     sunindextype stage_num, N_vector* out);

  SUNErrCode (*destroy)(SUNAdjointCheckpointScheme*);

  SUNDataNode root_data_node;

  void* impl;

  SUNContext sunctx;
};

typedef SUNAdjointCheckpointScheme_s* SUNAdjointCheckpointScheme;

SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(SUNContext);

#ifdef __cplusplus
}
#endif

#endif /*_SUNADJOINT_CHECKPOINTSCHEME_H*/
