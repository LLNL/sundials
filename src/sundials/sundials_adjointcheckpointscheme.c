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
 * SUNAdjointCheckpointScheme class definition.
 * ----------------------------------------------------------------*/

#include <stdint.h>
#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_core.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_types.h"

SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(
  SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointCheckpointScheme check_scheme = NULL;
  check_scheme                            = malloc(sizeof(*check_scheme));
  SUNAssert(check_scheme, SUN_ERR_MALLOC_FAIL);

  check_scheme->sunctx  = sunctx;
  check_scheme->content = NULL;
  check_scheme->ops     = NULL;

  SUNAdjointCheckpointScheme_Ops ops = NULL;
  ops                                = malloc(sizeof(*ops));
  SUNAssert(ops, SUN_ERR_MALLOC_FAIL);

  ops->needssaving   = NULL;
  ops->needsdeleting = NULL;
  ops->insertvector  = NULL;
  ops->loadvector    = NULL;
  ops->deleteVector  = NULL;
  ops->enableDense   = NULL;
  ops->destroy       = NULL;

  check_scheme->ops = ops;
  *check_scheme_ptr = check_scheme;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_NeedsSaving(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  sunrealtype t, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(check_scheme->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (check_scheme->ops->needssaving)
  {
    SUNErrCode err = check_scheme->ops->needssaving(check_scheme, step_num,
                                                    stage_num, t, yes_or_no);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNAdjointCheckpointScheme_NeedsDeleting(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(check_scheme->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (check_scheme->ops->needsdeleting)
  {
    SUNErrCode err = check_scheme->ops->needsdeleting(check_scheme, step_num,
                                                      stage_num, yes_or_no);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNAdjointCheckpointScheme_InsertVector(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  sunrealtype t, N_Vector state)
{
  SUNFunctionBegin(check_scheme->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);
  if (check_scheme->ops->insertvector)
  {
    SUNErrCode err = check_scheme->ops->insertvector(check_scheme, step_num,
                                                     stage_num, t, state);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNAdjointCheckpointScheme_LoadVector(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  sunbooleantype peek, N_Vector* out, sunrealtype* tout)
{
  SUNFunctionBegin(check_scheme->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);
  if (check_scheme->ops->loadvector)
  {
    SUNErrCode err = check_scheme->ops->loadvector(check_scheme, step_num,
                                                   stage_num, peek, out, tout);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNAdjointCheckpointScheme_RemoveVector(
  SUNAdjointCheckpointScheme check_scheme, int64_t step_num, int64_t stage_num,
  N_Vector* out)
{
  SUNFunctionBegin(check_scheme->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);
  if (check_scheme->ops->deleteVector)
  {
    SUNErrCode err = check_scheme->ops->deleteVector(check_scheme, step_num,
                                                     stage_num, out);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNAdjointCheckpointScheme_Destroy(
  SUNAdjointCheckpointScheme* check_scheme_ptr)
{
  SUNFunctionBegin((*check_scheme_ptr)->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);
  if ((*check_scheme_ptr)->ops->destroy)
  {
    SUNErrCode err = (*check_scheme_ptr)->ops->destroy(check_scheme_ptr);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  else if (*check_scheme_ptr)
  {
    free((*check_scheme_ptr)->ops);
    free(*check_scheme_ptr);
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_EnableDense(
  SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off)
{
  SUNFunctionBegin(check_scheme->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);
  if (check_scheme->ops->enableDense)
  {
    SUNErrCode err = check_scheme->ops->enableDense(check_scheme, on_or_off);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}
