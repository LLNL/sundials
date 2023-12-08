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
 * SUNDIALS memory helper implementation that uses the standard
 * system memory allocators.
 * ----------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sunmemory/sunmemory_system.h>

#include "sundials_debug.h"

struct SUNMemoryHelper_Content_Sys_
{
  unsigned long num_allocations;
  unsigned long num_deallocations;
  size_t bytes_allocated;
  size_t bytes_high_watermark;
};

typedef struct SUNMemoryHelper_Content_Sys_ SUNMemoryHelper_Content_Sys;

#define SUNHELPER_CONTENT(h) ((SUNMemoryHelper_Content_Sys*)h->content)

SUNMemoryHelper SUNMemoryHelper_Sys(SUNContext sunctx)
{
  SUNMemoryHelper helper;

  /* Allocate the helper */
  helper = SUNMemoryHelper_NewEmpty(sunctx);

  /* Set the ops */
  helper->ops->alloc         = SUNMemoryHelper_Alloc_Sys;
  helper->ops->dealloc       = SUNMemoryHelper_Dealloc_Sys;
  helper->ops->copy          = SUNMemoryHelper_Copy_Sys;
  helper->ops->getallocstats = SUNMemoryHelper_GetAllocStats_Sys;
  helper->ops->clone         = SUNMemoryHelper_Clone_Sys;
  helper->ops->destroy       = SUNMemoryHelper_Destroy_Sys;

  /* Attach content and ops */
  helper->content =
    (SUNMemoryHelper_Content_Sys*)malloc(sizeof(SUNMemoryHelper_Content_Sys));
  SUNHELPER_CONTENT(helper)->num_allocations      = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations    = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated      = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark = 0;

  return helper;
}

int SUNMemoryHelper_Alloc_Sys(SUNMemoryHelper helper, SUNMemory* memptr,
                              size_t mem_size, SUNMemoryType mem_type,
                              void* queue)
{
  SUNMemory mem = SUNMemoryNewEmpty(helper->sunctx);

  mem->ptr   = NULL;
  mem->own   = SUNTRUE;
  mem->type  = mem_type;
  mem->bytes = mem_size;

  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr = malloc(mem_size);
    if (mem->ptr == NULL)
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Sys: malloc returned NULL\n");
      free(mem);
      return (-1);
    }
    SUNHELPER_CONTENT(helper)->bytes_allocated += mem_size;
    SUNHELPER_CONTENT(helper)->num_allocations++;
    SUNHELPER_CONTENT(helper)->bytes_high_watermark =
      SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated,
             SUNHELPER_CONTENT(helper)->bytes_high_watermark);
  }
  else
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in SUNMemoryHelper_Alloc_Sys: unsupported memory type\n");
    free(mem);
    return (-1);
  }

  *memptr = mem;
  return (0);
}

int SUNMemoryHelper_Dealloc_Sys(SUNMemoryHelper helper, SUNMemory mem, void* queue)
{
  if (mem == NULL) { return (0); }

  if (mem->ptr != NULL && mem->own)
  {
    if (mem->type == SUNMEMTYPE_HOST)
    {
      SUNHELPER_CONTENT(helper)->num_deallocations++;
      SUNHELPER_CONTENT(helper)->bytes_allocated -= mem->bytes;
      free(mem->ptr);
      mem->ptr = NULL;
    }
    else
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Dealloc_Sys: unsupported memory type\n");
      return (-1);
    }
  }

  free(mem);
  return (0);
}

int SUNMemoryHelper_Copy_Sys(SUNMemoryHelper helper, SUNMemory dst,
                             SUNMemory src, size_t memory_size, void* queue)
{
  int retval = 0;
  memcpy(dst->ptr, src->ptr, memory_size);
  return (retval);
}

int SUNMemoryHelper_GetAllocStats_Sys(SUNMemoryHelper helper,
                                      SUNMemoryType mem_type,
                                      unsigned long* num_allocations,
                                      unsigned long* num_deallocations,
                                      size_t* bytes_allocated,
                                      size_t* bytes_high_watermark)
{
  if (mem_type == SUNMEMTYPE_HOST)
  {
    *num_allocations      = SUNHELPER_CONTENT(helper)->num_allocations;
    *num_deallocations    = SUNHELPER_CONTENT(helper)->num_deallocations;
    *bytes_allocated      = SUNHELPER_CONTENT(helper)->bytes_allocated;
    *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark;
  }
  else { return -1; }
  return 0;
}

SUNMemoryHelper SUNMemoryHelper_Clone_Sys(SUNMemoryHelper helper)
{
  SUNMemoryHelper hclone = SUNMemoryHelper_Sys(helper->sunctx);
  return hclone;
}

int SUNMemoryHelper_Destroy_Sys(SUNMemoryHelper helper)
{
  if (helper)
  {
    if (helper->content) { free(helper->content); }
    if (helper->ops) { free(helper->ops); }
    free(helper);
  }
  return 0;
}
