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
 * SUNAdjointCheckpointScheme_Basic class definition.
 * ----------------------------------------------------------------*/

#include <sunadjoint/sunadjoint_checkpointscheme_basic.h>
#include <sundials/sundatanode_inmem.h>
#include <sundials/sundials_core.h>

#include "sunadjoint/sunadjoint_checkpointscheme.h"
#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_datanode.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_logger.h"
#include "sundials/sundials_types.h"
#include "sundials_utils.h"

struct SUNAdjointCheckpointScheme_Basic_Content_
{
  uint64_t interval;
  sunbooleantype save_stages;
  sunbooleantype keep;
  SUNDataIOMode io_mode;
  SUNDataNode root_node;
  uint64_t stepnum_of_current_insert;
  SUNDataNode current_insert_step_node;
  uint64_t stepnum_of_current_load;
  SUNDataNode current_load_step_node;
};

typedef struct SUNAdjointCheckpointScheme_Basic_Content_*
  SUNAdjointCheckpointScheme_Basic_Content;

#define GET_CONTENT(S)     ((SUNAdjointCheckpointScheme_Basic_Content)S->content)
#define IMPL_PROP(S, prop) (GET_CONTENT(S)->prop)

SUNErrCode SUNAdjointCheckpointScheme_Create_Basic(
  SUNDataIOMode io_mode, uint64_t interval, uint64_t estimate,
  sunbooleantype save_stages, sunbooleantype keep, SUNContext sunctx,
  SUNAdjointCheckpointScheme* check_scheme_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointCheckpointScheme check_scheme = NULL;
  SUNCheckCall(SUNAdjointCheckpointScheme_NewEmpty(sunctx, &check_scheme));

  check_scheme->ops->shouldWeSave = SUNAdjointCheckpointScheme_ShouldWeSave_Basic;
  // check_scheme->ops->shouldWeDelete = SUNAdjointCheckpointScheme_ShouldWeDelete;
  // check_scheme->ops->insert = SUNAdjointCheckpointScheme_Insert_Basic;
  check_scheme->ops->insertVector = SUNAdjointCheckpointScheme_InsertVector_Basic;
  // check_scheme->ops->load       = SUNAdjointCheckpointScheme_Load_Basic;
  check_scheme->ops->loadVector = SUNAdjointCheckpointScheme_LoadVector_Basic;
  // check_scheme->ops->remove     = SUNAdjointCheckpointScheme_Remove_Basic;
  // check_scheme->ops->removeVector = SUNAdjointCheckpointScheme_RemoveVector_Basic;
  // check_scheme->ops->removeRange = SUNAdjointCheckpointScheme_RemoveRange_Basic;
  // check_scheme->ops->destroy     = SUNAdjointCheckpointScheme_Destroy_Basic;

  SUNAdjointCheckpointScheme_Basic_Content content = NULL;

  content = malloc(sizeof(*content));
  SUNAssert(content, SUN_ERR_MALLOC_FAIL);

  content->interval                  = interval;
  content->save_stages               = save_stages;
  content->keep                      = keep;
  content->root_node                 = NULL;
  content->current_insert_step_node  = NULL;
  content->stepnum_of_current_insert = 0;
  content->current_load_step_node    = NULL;
  content->stepnum_of_current_load   = 0;
  content->io_mode                   = io_mode;

  SUNCheckCall(
    SUNDataNode_CreateObject(io_mode, estimate, sunctx, &content->root_node));

  check_scheme->content = content;
  *check_scheme_ptr     = check_scheme;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave_Basic(
  SUNAdjointCheckpointScheme cs, sunindextype step_num, sunindextype stage_num,
  sunrealtype t, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(cs->sunctx);

  if (!(step_num % IMPL_PROP(cs, interval)))
  {
    if (stage_num == 0) { *yes_or_no = SUNTRUE; }
    else if (IMPL_PROP(cs, save_stages)) { *yes_or_no = SUNTRUE; }
    else { *yes_or_no = SUNFALSE; }
  }
  else { *yes_or_no = SUNFALSE; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_InsertVector_Basic(
  SUNAdjointCheckpointScheme cs, sunindextype step_num, sunindextype stage_num,
  sunrealtype t, N_Vector state)
{
  SUNFunctionBegin(cs->sunctx);

  /* If this is the step solution, then we need to create a list node first to
     store the step and stage solutions in. We keep a pointer to the list node
     until this step is over for fast access when inserting stages. */
  SUNDataNode step_data_node = NULL;
  if (step_num != IMPL_PROP(cs, stepnum_of_current_insert) || stage_num == -1)
  {
    SUNCheckCall(SUNDataNode_CreateList(IMPL_PROP(cs, io_mode), 0, SUNCTX_,
                                        &step_data_node));
    IMPL_PROP(cs, current_insert_step_node)  = step_data_node;
    IMPL_PROP(cs, stepnum_of_current_insert) = step_num;

    /* Store the step node in the root node object. */
    char* key = sunUnsignedToString(step_num);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                       "insert-new-step", "step_num = %d, key = %s", step_num,
                       key);
#endif
    SUNCheckCall(
      SUNDataNode_AddNamedChild(IMPL_PROP(cs, root_node), key, step_data_node));
    free(key);
  }
  else { step_data_node = IMPL_PROP(cs, current_insert_step_node); }

  /* Add the solution data as a node in the step list. */
  SUNDataNode solution_node = NULL;
  SUNCheckCall(SUNDataNode_CreateLeaf(IMPL_PROP(cs, io_mode), state, 1,
                                      sizeof(N_Vector), SUNCTX_, &solution_node));

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                     "insert-stage", "step_num = %d, stage_num = %d", step_num,
                     stage_num);
#endif
  SUNDataNode_AddChild(step_data_node, solution_node);

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Basic(
  SUNAdjointCheckpointScheme cs, sunindextype step_num, sunindextype stage_num,
  N_Vector* out)
{
  SUNFunctionBegin(cs->sunctx);

  /* If we are trying to load the step solution, we need to load the list which holds
     the step and stage solutions. We keep a pointer to the list node until
     this step is over for fast access when loading stages. */
  SUNDataNode step_data_node = NULL;
  if (step_num != IMPL_PROP(cs, stepnum_of_current_load) || stage_num == -1)
  {
    char* key = sunUnsignedToString(step_num);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                       "load-new-step", "step_num = %d, key = %s", step_num, key);
#endif
    SUNCheckCall(SUNDataNode_GetNamedChild(IMPL_PROP(cs, root_node), key,
                                           &step_data_node));
    free(key);
    IMPL_PROP(cs, current_load_step_node)  = step_data_node;
    IMPL_PROP(cs, stepnum_of_current_load) = step_num;
  }
  else { step_data_node = IMPL_PROP(cs, current_load_step_node); }

  SUNAssert(step_data_node, SUN_ERR_CORRUPT);

  SUNDataNode solution_node = NULL;
  SUNDataNode_GetChild(step_data_node, stage_num + 1, &solution_node);

  void* data = NULL;
  SUNCheckCall(SUNDataNode_GetData(solution_node, &data));
  SUNAssert(data, SUN_ERR_CORRUPT);

  *out = (N_Vector)(data);

  return SUN_SUCCESS;
}