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

#include <sundials/sundials_math.h>
#include <sundials/sundials_memory.h>
#include "sundials_debug.h"
#include "sundials_context_impl.h"

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(SUNMemoryHelper H)
{
  return(H->sunctx->profiler);
}
#endif

SUNMemory SUNMemoryNewEmpty()
{
  SUNMemory mem = NULL;

  mem = (SUNMemory) malloc(sizeof(struct _SUNMemory));
  if (mem == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryNewEmpty: malloc failed\n");
    return(NULL);
  }

  mem->bytes = 0;

  return(mem);
}


SUNMemoryHelper SUNMemoryHelper_NewEmpty(SUNContext sunctx)
{
  SUNMemoryHelper helper = NULL;

  if (sunctx == NULL) return(NULL);

  helper = (SUNMemoryHelper) malloc(sizeof(struct _SUNMemoryHelper));
  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_NewEmpty: malloc failed\n");
    return(NULL);
  }

  helper->ops = (SUNMemoryHelper_Ops) malloc(sizeof(struct _SUNMemoryHelper_Ops));
  if (helper->ops == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_NewEmpty: malloc failed\n");
    free(helper);
    return(NULL);
  }

  /* Set all ops to NULL */
  memset(helper->ops, 0, sizeof(struct _SUNMemoryHelper_Ops));
  helper->content = NULL;
  helper->sunctx = sunctx;

  return(helper);
}


int SUNMemoryHelper_CopyOps(SUNMemoryHelper src, SUNMemoryHelper dst)
{
  /* Check that ops structures exist */
  if (src == NULL || dst == NULL || src->ops == NULL || dst->ops == NULL)
    return(-1);
  memcpy(dst->ops, src->ops, sizeof(struct _SUNMemoryHelper_Ops));
  return(0);
}


booleantype SUNMemoryHelper_ImplementsRequiredOps(SUNMemoryHelper helper)
{
  if (helper->ops->alloc == NULL || helper->ops->dealloc == NULL ||
      helper->ops->copy == NULL)
  {
    return(SUNFALSE);
  }
  return(SUNTRUE);
}


SUNMemory SUNMemoryHelper_Alias(SUNMemory mem)
{
  SUNMemory alias = SUNMemoryNewEmpty();

  alias->ptr  = mem->ptr;
  alias->type = mem->type;
  alias->own  = SUNFALSE;

  return(alias);
}


SUNMemory SUNMemoryHelper_Wrap(void* ptr, SUNMemoryType mem_type)
{
  SUNMemory mem = SUNMemoryNewEmpty();

  mem->ptr = ptr;
  mem->own = SUNFALSE;

  switch(mem_type)
  {
    case SUNMEMTYPE_HOST:
      mem->type = SUNMEMTYPE_HOST;
      break;
    case SUNMEMTYPE_PINNED:
      mem->type = SUNMEMTYPE_PINNED;
      break;
    case SUNMEMTYPE_DEVICE:
      mem->type = SUNMEMTYPE_DEVICE;
      break;
    case SUNMEMTYPE_UVM:
      mem->type = SUNMEMTYPE_UVM;
      break;
    default:
      free(mem);
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Wrap: unknown memory type\n");
      return(NULL);
  }

  return(mem);
}

int SUNMemoryHelper_GetAllocStats(SUNMemoryHelper helper, SUNMemoryType mem_type, unsigned long* num_allocations,
                                  unsigned long* num_deallocations, size_t* bytes_allocated,
                                  size_t* bytes_high_watermark)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (helper->ops->getallocstats) {
    return helper->ops->getallocstats(helper, mem_type, num_allocations, num_deallocations, bytes_allocated, bytes_high_watermark);
  } else {
    ier = helper->ops->getallocstats(helper, mem_type, num_allocations, num_deallocations, bytes_allocated, bytes_high_watermark);
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return(ier);
}


int SUNMemoryHelper_Alloc(SUNMemoryHelper helper, SUNMemory* memptr,
                          size_t mem_size, SUNMemoryType mem_type, void* queue)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (helper->ops->alloc == NULL) {
    ier = -1;
  } else {
    ier = helper->ops->alloc(helper, memptr, mem_size, mem_type, queue);
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return(ier);
}


int SUNMemoryHelper_Dealloc(SUNMemoryHelper helper, SUNMemory mem, void* queue)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (helper->ops->dealloc == NULL) { ier = -1; }
  if (!mem) {
    ier = 0;
  } else {
    ier = helper->ops->dealloc(helper, mem, queue);
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return(ier);
}


int SUNMemoryHelper_Copy(SUNMemoryHelper helper, SUNMemory dst,
                         SUNMemory src, size_t memory_size, void* queue)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (helper->ops->copy == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Copy: function pointer is NULL\n");
    ier = -1;
  }
  else
  {
    ier = helper->ops->copy(helper, dst, src, memory_size, queue);
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return(ier);
}


int SUNMemoryHelper_CopyAsync(SUNMemoryHelper helper, SUNMemory dst,
                              SUNMemory src, size_t memory_size,
                              void* queue)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (helper->ops->copyasync == NULL)
    ier = SUNMemoryHelper_Copy(helper, dst, src, memory_size, queue);
  else
    ier = helper->ops->copyasync(helper, dst, src, memory_size, queue);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return(ier);
}


int SUNMemoryHelper_Destroy(SUNMemoryHelper helper)
{
  if (!helper) return 0;

  if (helper->ops->destroy)
  {
    /* user helper defined destroy */
    return helper->ops->destroy(helper);
  }
  else if (helper->content)
  {
    /* helper should have defined destroy */
    return -1;
  }
  else
  {
    /* default destroy */
    free(helper->ops);
    free(helper);
    return 0;
  }

  return 0;
}


SUNMemoryHelper SUNMemoryHelper_Clone(SUNMemoryHelper helper)
{
  if (helper->ops->clone == NULL)
  {
    if (helper->content != NULL)
    {
      return(NULL);
    }
    else
    {
      SUNMemoryHelper hclone = SUNMemoryHelper_NewEmpty(helper->sunctx);
      if (hclone) SUNMemoryHelper_CopyOps(helper, hclone);
      return(hclone);
    }
  }
  else
  {
    return(helper->ops->clone(helper));
  }
}
