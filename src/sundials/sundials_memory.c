/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS memory helper.
 * ----------------------------------------------------------------*/

#include <string.h>
#include <sundials/impl/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sundials_debug.h"

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(SUNMemoryHelper H)
{
  return (H->sunctx->profiler);
}
#endif

SUNMemory SUNMemoryNewEmpty()
{
  SUNMemory mem = NULL;

  mem = (SUNMemory)malloc(sizeof(struct _SUNMemory));

  mem->bytes = 0;

  return (mem);
}

SUNMemoryHelper SUNMemoryHelper_NewEmpty(SUNContext sunctx)
{
  SUNAssignSUNCTX(sunctx);
  SUNMemoryHelper helper = NULL;

  helper = (SUNMemoryHelper)malloc(sizeof(struct _SUNMemoryHelper));
  SUNAssert(helper, SUN_ERR_MALLOC_FAIL);

  helper->ops = (SUNMemoryHelper_Ops)malloc(sizeof(struct _SUNMemoryHelper_Ops));
  SUNAssert(helper->ops, SUN_ERR_MALLOC_FAIL);

  /* Set all ops to NULL */
  memset(helper->ops, 0, sizeof(struct _SUNMemoryHelper_Ops));
  helper->content = NULL;
  helper->sunctx  = sunctx;

  return helper;
}

SUNErrCode SUNMemoryHelper_CopyOps(SUNMemoryHelper src, SUNMemoryHelper dst)
{
  memcpy(dst->ops, src->ops, sizeof(struct _SUNMemoryHelper_Ops));
  return SUN_SUCCESS;
}

booleantype SUNMemoryHelper_ImplementsRequiredOps(SUNMemoryHelper helper)
{
  if (helper->ops->alloc == NULL || helper->ops->dealloc == NULL ||
      helper->ops->copy == NULL)
  {
    return SUNFALSE;
  }
  return SUNTRUE;
}

SUNMemory SUNMemoryHelper_Alias(SUNMemory mem)
{
  SUNMemory alias = SUNMemoryNewEmpty();

  alias->ptr  = mem->ptr;
  alias->type = mem->type;
  alias->own  = SUNFALSE;

  return alias;
}

SUNMemory SUNMemoryHelper_Wrap(void* ptr, SUNMemoryType mem_type)
{
  SUNMemory mem = SUNMemoryNewEmpty();

  mem->ptr = ptr;
  mem->own = SUNFALSE;

  switch (mem_type)
  {
  case SUNMEMTYPE_HOST: mem->type = SUNMEMTYPE_HOST; break;
  case SUNMEMTYPE_PINNED: mem->type = SUNMEMTYPE_PINNED; break;
  case SUNMEMTYPE_DEVICE: mem->type = SUNMEMTYPE_DEVICE; break;
  case SUNMEMTYPE_UVM: mem->type = SUNMEMTYPE_UVM; break;
  default:
    free(mem);
    /* TODO(CJB): We dont have access to SUNContext to handle this error */
    SUNDIALS_DEBUG_PRINT(
      "ERROR in SUNMemoryHelper_Wrap: unknown memory type\n");
    return (NULL);
  }

  return mem;
}

SUNErrCode SUNMemoryHelper_GetAllocStats(SUNMemoryHelper helper,
                                         SUNMemoryType mem_type,
                                         unsigned long* num_allocations,
                                         unsigned long* num_deallocations,
                                         size_t* bytes_allocated,
                                         size_t* bytes_high_watermark)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (helper->ops->getallocstats)
  {
    ier = helper->ops->getallocstats(helper, mem_type, num_allocations,
                                     num_deallocations, bytes_allocated,
                                     bytes_high_watermark);
  }
  else { ier = SUN_ERR_NOT_IMPLEMENTED; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return (ier);
}

SUNErrCode SUNMemoryHelper_Alloc(SUNMemoryHelper helper, SUNMemory* memptr,
                                 size_t mem_size, SUNMemoryType mem_type,
                                 void* queue)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (helper->ops->alloc)
  {
    ier = helper->ops->alloc(helper, memptr, mem_size, mem_type, queue);
  }
  else { ier = SUN_ERR_NOT_IMPLEMENTED; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return ier;
}

SUNErrCode SUNMemoryHelper_Dealloc(SUNMemoryHelper helper, SUNMemory mem,
                                   void* queue)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (helper->ops->dealloc == NULL) { ier = -1; }
  if (!mem) { ier = SUN_SUCCESS; }
  else { ier = helper->ops->dealloc(helper, mem, queue); }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return ier;
}

SUNErrCode SUNMemoryHelper_Copy(SUNMemoryHelper helper, SUNMemory dst,
                                SUNMemory src, size_t memory_size, void* queue)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (!helper->ops->copy) { ier = SUN_ERR_NOT_IMPLEMENTED; }
  else { ier = helper->ops->copy(helper, dst, src, memory_size, queue); }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return ier;
}

SUNErrCode SUNMemoryHelper_CopyAsync(SUNMemoryHelper helper, SUNMemory dst,
                                     SUNMemory src, size_t memory_size,
                                     void* queue)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (!helper->ops->copyasync)
  {
    ier = SUNMemoryHelper_Copy(helper, dst, src, memory_size, queue);
  }
  else { ier = helper->ops->copyasync(helper, dst, src, memory_size, queue); }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return ier;
}

SUNErrCode SUNMemoryHelper_Destroy(SUNMemoryHelper helper)
{
  if (!helper) return SUN_SUCCESS;

  if (helper->ops->destroy)
  {
    /* user helper defined destroy */
    return helper->ops->destroy(helper);
  }
  else if (helper->content)
  {
    /* helper should have defined destroy */
    return SUN_ERR_NOT_IMPLEMENTED;
  }
  else
  {
    /* default destroy */
    free(helper->ops);
    free(helper);
    return SUN_SUCCESS;
  }

  return SUN_SUCCESS;
}

SUNMemoryHelper SUNMemoryHelper_Clone(SUNMemoryHelper helper)
{
  SUNAssignSUNCTX(helper->sunctx);
  if (!helper->ops->clone)
  {
    if (helper->content)
    {
      SUNCheck(!helper->content, SUN_ERR_NOT_IMPLEMENTED);
      return (NULL);
    }
    else
    {
      SUNMemoryHelper hclone = SUNMemoryHelper_NewEmpty(helper->sunctx);
      if (hclone) SUNMemoryHelper_CopyOps(helper, hclone);
      return (hclone);
    }
  }
  else { return (helper->ops->clone(helper)); }
}
