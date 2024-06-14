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
#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_datanode.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_logger.h"
#include "sundials/sundials_memory.h"
#include "sundials/sundials_types.h"
#include "sundials_utils.h"

struct SUNAdjointCheckpointScheme_Basic_Content_
{
  SUNMemoryHelper mem_helper;
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
  check_scheme->ops->loadVector = SUNAdjointCheckpointScheme_LoadVector_Basic;
  check_scheme->ops->destroy    = SUNAdjointCheckpointScheme_Destroy_Basic;

  SUNAdjointCheckpointScheme_Basic_Content content = NULL;

  content = malloc(sizeof(*content));
  SUNAssert(content, SUN_ERR_MALLOC_FAIL);

  content->mem_helper                = mem_helper;
  content->interval                  = interval;
  content->save_stages               = save_stages;
  content->keep                      = keep;
  content->root_node                 = NULL;
  content->current_insert_step_node  = NULL;
  content->stepnum_of_current_insert = -1;
  content->current_load_step_node    = NULL;
  content->stepnum_of_current_load   = -1;
  content->io_mode                   = io_mode;

  SUNCheckCall(
    SUNDataNode_CreateObject(io_mode, estimate, sunctx, &content->root_node));

  check_scheme->content = content;
  *check_scheme_ptr     = check_scheme;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_ShouldWeSave_Basic(
  SUNAdjointCheckpointScheme self, sunindextype step_num,
  sunindextype stage_num, sunrealtype t, sunbooleantype* yes_or_no)
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

  /* If this is the step solution, then we need to create a list node first to
     store the step and stage solutions in. We keep a pointer to the list node
     until this step is over for fast access when inserting stages. */
  SUNDataNode step_data_node = NULL;
  if (step_num != PROPERTY(self, stepnum_of_current_insert))
  {
    SUNCheckCall(SUNDataNode_CreateList(PROPERTY(self, io_mode), 0, SUNCTX_,
                                        &step_data_node));
    PROPERTY(self, current_insert_step_node)  = step_data_node;
    PROPERTY(self, stepnum_of_current_insert) = step_num;

    /* Store the step node in the root node object. */
    char* key = sunUnsignedToString(step_num);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                       "insert-new-step", "step_num = %d, key = %s", step_num,
                       key);
#endif
    SUNCheckCall(SUNDataNode_AddNamedChild(PROPERTY(self, root_node), key,
                                           step_data_node));
    free(key);
  }
  else { step_data_node = PROPERTY(self, current_insert_step_node); }

  /* Add the solution data as a node in the step list. */
  SUNDataNode solution_node = NULL;
  SUNCheckCall(SUNDataNode_CreateLeaf(PROPERTY(self, io_mode),
                                      PROPERTY(self, mem_helper), SUNCTX_,
                                      &solution_node));
  SUNCheckCall(SUNDataNode_SetDataNvector(solution_node, state));

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                     "insert-stage", "step_num = %d, stage_num = %d", step_num,
                     stage_num);
#endif
  SUNCheckCall(SUNDataNode_AddChild(step_data_node, solution_node));

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Basic(
  SUNAdjointCheckpointScheme self, sunindextype step_num,
  sunindextype stage_num, N_Vector* loaded_state)
{
  SUNFunctionBegin(self->sunctx);

  /* If we are trying to load the step solution, we need to load the list which holds
     the step and stage solutions. We keep a pointer to the list node until
     this step is over for fast access when loading stages. */
  SUNDataNode step_data_node = NULL;
  if (step_num != PROPERTY(self, stepnum_of_current_load))
  {
    char* key = sunUnsignedToString(step_num);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(SUNCTX_->logger, SUN_LOGLEVEL_DEBUG, __func__,
                       "load-new-step", "step_num = %d, key = %s", step_num, key);
#endif
    SUNCheckCall(SUNDataNode_GetNamedChild(PROPERTY(self, root_node), key,
                                           &step_data_node));
    free(key);
    PROPERTY(self, current_load_step_node)  = step_data_node;
    PROPERTY(self, stepnum_of_current_load) = step_num;
  }
  else { step_data_node = PROPERTY(self, current_load_step_node); }

  if (!step_data_node) { return SUN_ERR_CHECKPOINT_NOT_FOUND; }

  SUNDataNode solution_node = NULL;
  if (PROPERTY(self, keep))
  {
    SUNCheckCall(
      SUNDataNode_GetChild(step_data_node, stage_num + 1, &solution_node));
  }
  else
  {
    SUNCheckCall(
      SUNDataNode_RemoveChild(step_data_node, stage_num + 1, &solution_node));
  }

  if (!solution_node) { return SUN_ERR_CHECKPOINT_NOT_FOUND; }

  SUNCheckCall(SUNDataNode_GetDataNvector(solution_node, *loaded_state));

  /* Cleanup the checkpoint memory if need be */
  if (!PROPERTY(self, keep))
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