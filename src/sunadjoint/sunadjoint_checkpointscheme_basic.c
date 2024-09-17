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
#include <sundials/sundials_core.h>

#include "sunadjoint/sunadjoint_checkpointscheme.h"
#include "sundatanode/sundatanode_inmem.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_logger.h"
#include "sundials/sundials_memory.h"
#include "sundials/sundials_types.h"
#include "sundials_datanode.h"
#include "sundials_macros.h"
#include "sundials_utils.h"

struct SUNAdjointCheckpointScheme_Basic_Content_
{
  SUNMemoryHelper mem_helper;
  uint64_t backup_interval;
  uint64_t interval;
  sunbooleantype save_stages;
  sunbooleantype keep;
  SUNDataIOMode io_mode;
  SUNDataNode root_node;
  int64_t stepnum_of_current_insert;
  SUNDataNode current_insert_step_node;
  int64_t stepnum_of_current_load;
  SUNDataNode current_load_step_node;
};

typedef struct SUNAdjointCheckpointScheme_Basic_Content_*
  SUNAdjointCheckpointScheme_Basic_Content;

#define GET_CONTENT(S)    ((SUNAdjointCheckpointScheme_Basic_Content)S->content)
#define PROPERTY(S, prop) (GET_CONTENT(S)->prop)

SUNErrCode SUNAdjointCheckpointScheme_Create_Basic(
  SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, uint64_t interval,
  uint64_t estimate, sunbooleantype save_stages, sunbooleantype keep,
  SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointCheckpointScheme check_scheme = NULL;
  SUNCheckCall(SUNAdjointCheckpointScheme_NewEmpty(sunctx, &check_scheme));

  check_scheme->ops->shouldWeSave = SUNAdjointCheckpointScheme_ShouldWeSave_Basic;
  check_scheme->ops->insertVector = SUNAdjointCheckpointScheme_InsertVector_Basic;
  check_scheme->ops->loadVector  = SUNAdjointCheckpointScheme_LoadVector_Basic;
  check_scheme->ops->enableDense = SUNAdjointCheckpointScheme_EnableDense_Basic;
  check_scheme->ops->destroy     = SUNAdjointCheckpointScheme_Destroy_Basic;

  SUNAdjointCheckpointScheme_Basic_Content content = NULL;

  content = malloc(sizeof(*content));
  SUNAssert(content, SUN_ERR_MALLOC_FAIL);

  content->mem_helper                = mem_helper;
  content->interval                  = interval;
  content->save_stages               = save_stages;
  content->keep                      = keep;
  content->root_node                 = NULL;
  content->current_insert_step_node  = NULL;
  content->stepnum_of_current_insert = -2;
  content->current_load_step_node    = NULL;
  content->stepnum_of_current_load   = -2;
  content->io_mode                   = io_mode;

  SUNCheckCall(
    SUNDataNode_CreateObject(io_mode, estimate, sunctx, &content->root_node));

  check_scheme->content = content;
  *check_scheme_ptr     = check_scheme;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave_Basic(
  SUNAdjointCheckpointScheme self, sunindextype step_num, sunindextype stage_num,
  SUNDIALS_MAYBE_UNUSED sunrealtype t, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  if (!(step_num % PROPERTY(self, interval)))
  {
    if (stage_num == 0) { *yes_or_no = SUNTRUE; }
    else if (PROPERTY(self, save_stages)) { *yes_or_no = SUNTRUE; }
    else { *yes_or_no = SUNFALSE; }
  }
  else { *yes_or_no = SUNFALSE; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_InsertVector_Basic(
  SUNAdjointCheckpointScheme self, sunindextype step_num,
  sunindextype stage_num, sunrealtype t, N_Vector state)
{
  SUNFunctionBegin(self->sunctx);

  /* If this is the first state for a step, then we need to create a
     list node first to store the step and all stage solutions in.
     We keep a pointer to the list node until this step is over for
     fast access when inserting stages. */
  SUNDataNode step_data_node = NULL;
  if (step_num != PROPERTY(self, stepnum_of_current_insert))
  {
    SUNCheckCall(SUNDataNode_CreateList(PROPERTY(self, io_mode), 0, SUNCTX_,
                                        &step_data_node));
    PROPERTY(self, current_insert_step_node)  = step_data_node;
    PROPERTY(self, stepnum_of_current_insert) = step_num;

    /* Store the step node in the root node object. */
    char* key = sunSignedToString(step_num);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                       "insert-new-step", "step_num = %d", step_num);

#endif
    SUNCheckCall(SUNDataNode_AddNamedChild(PROPERTY(self, root_node), key,
                                           step_data_node));
    free(key);
  }
  else { step_data_node = PROPERTY(self, current_insert_step_node); }

  /* Add the state data as a leaf node in the step node's list of children. */
  SUNDataNode solution_node = NULL;
  SUNCheckCall(SUNDataNode_CreateLeaf(PROPERTY(self, io_mode),
                                      PROPERTY(self, mem_helper), SUNCTX_,
                                      &solution_node));
  SUNCheckCall(SUNDataNode_SetDataNvector(solution_node, state, t));

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                     "insert-stage", "step_num = %d, stage_num = %d, t = %g",
                     step_num, stage_num, t);
#endif
  SUNCheckCall(SUNDataNode_AddChild(step_data_node, solution_node));

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Basic(
  SUNAdjointCheckpointScheme self, sunindextype step_num, sunindextype stage_num,
  sunbooleantype peek, N_Vector* loaded_state, sunrealtype* t)
{
  SUNFunctionBegin(self->sunctx);

  SUNErrCode errcode = SUN_SUCCESS;

  /* If we are trying to load the step solution, we need to load the list which holds
     the step and stage solutions. We keep a pointer to the list node until
     this step is over for fast access when loading stages. */
  SUNDataNode step_data_node = NULL;
  if (step_num != PROPERTY(self, stepnum_of_current_load))
  {
    char* key = sunSignedToString(step_num);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                       "try-load-new-step", "step_num = %d, stage_num = %d",
                       step_num, stage_num);
#endif
    errcode = SUNDataNode_GetNamedChild(PROPERTY(self, root_node), key,
                                        &step_data_node);
    if (errcode == SUN_SUCCESS)
    {
      PROPERTY(self, current_load_step_node)  = step_data_node;
      PROPERTY(self, stepnum_of_current_load) = step_num;
    }
    else if (errcode == SUN_ERR_DATANODE_NODENOTFOUND)
    {
      step_data_node = NULL;
    }
    else { SUNCheckCall(errcode); }
    free(key);
  }
  else { step_data_node = PROPERTY(self, current_load_step_node); }

  if (!step_data_node)
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                       "step-not-found", "step_num = %d, stage_num = %d",
                       step_num, stage_num);
#endif
    return SUN_ERR_CHECKPOINT_NOT_FOUND;
  }

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__, "step-loaded",
                     "step_num = %d, stage_num = %d", step_num, stage_num);
#endif

  SUNDataNode solution_node = NULL;
  if (PROPERTY(self, keep) || peek)
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                       "try-load-stage",
                       "keep = 1, step_num = %d, stage_num = %d", step_num,
                       stage_num);
#endif
    errcode = SUNDataNode_GetChild(step_data_node, stage_num, &solution_node);
    if (errcode == SUN_ERR_DATANODE_NODENOTFOUND) { solution_node = NULL; }
    else { SUNCheckCall(errcode); }
  }
  else
  {
    sunbooleantype has_children = SUNFALSE;
    SUNCheckCall(SUNDataNode_HasChildren(step_data_node, &has_children));

    if (has_children)
    {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
      SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                         "try-load-stage",
                         "keep = 0, step_num = %d, stage_num = %d", step_num,
                         stage_num);
#endif
      errcode = SUNDataNode_RemoveChild(step_data_node, stage_num,
                                        &solution_node);
      if (errcode == SUN_ERR_DATANODE_NODENOTFOUND) { solution_node = NULL; }
      else { SUNCheckCall(errcode); }
    }

    SUNCheckCall(SUNDataNode_HasChildren(step_data_node, &has_children));
    if (!has_children)
    {
      char* key = sunSignedToString(step_num);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
      SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                         "remove-step", "step_num = %d", step_num);
#endif
      SUNCheckCall(SUNDataNode_RemoveNamedChild(PROPERTY(self, root_node), key,
                                                &step_data_node));
      free(key);
      SUNCheckCall(SUNDataNode_Destroy(&step_data_node));
    }
  }

  if (!solution_node)
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                       "stage-not-found", "step_num = %d, stage_num = %d",
                       step_num, stage_num);
#endif
    return SUN_ERR_CHECKPOINT_NOT_FOUND;
  }

  SUNCheckCall(SUNDataNode_GetDataNvector(solution_node, *loaded_state, t));
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                     "stage-loaded", "step_num = %d, stage_num = %d, t = %g",
                     step_num, stage_num, *t);
#endif

  /* Cleanup the checkpoint memory if need be */
  if (!(PROPERTY(self, keep) || peek))
  {
    SUNCheckCall(SUNDataNode_Destroy(&solution_node));
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_Destroy_Basic(
  SUNAdjointCheckpointScheme* check_scheme_ptr)
{
  SUNFunctionBegin((*check_scheme_ptr)->sunctx);

  SUNAdjointCheckpointScheme self = *check_scheme_ptr;
  SUNAdjointCheckpointScheme_Basic_Content content =
    (SUNAdjointCheckpointScheme_Basic_Content)self->content;

  SUNCheckCall(SUNDataNode_Destroy(&content->root_node));
  free(content);
  free(self);

  *check_scheme_ptr = NULL;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_EnableDense_Basic(
  SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off)
{
  SUNFunctionBegin(check_scheme->sunctx);

  if (on_or_off)
  {
    PROPERTY(check_scheme, backup_interval) = PROPERTY(check_scheme, interval);
    PROPERTY(check_scheme, interval)        = 1;
  }
  else
  {
    PROPERTY(check_scheme, interval) = PROPERTY(check_scheme, backup_interval);
  }

  return SUN_SUCCESS;
}
