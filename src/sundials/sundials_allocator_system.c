/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Implementation of the SUNDIALS allocator class wrapping malloc() and free()
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "sundials/sundials_errors.h"
#include "sundials/sundials_math.h"
#include "sundials_allocator_system.h"

/* -----------------------------------------------------------------------------
 * Definition of class member data
 * ---------------------------------------------------------------------------*/

struct SUNAllocator_Content_Sys_
{
  unsigned long num_allocations;
  unsigned long num_deallocations;
  size_t bytes_allocated;
  size_t bytes_max;
};

typedef struct SUNAllocator_Content_Sys_* SUNAllocator_Content_Sys;

#define GET_CONTENT(x) ((SUNAllocator_Content_Sys)(x->content))

/* -----------------------------------------------------------------------------
 * Constructor and destructor
 * ---------------------------------------------------------------------------*/

SUNErrCode SUNAllocator_Create_System(SUNAllocator* sun_alloc)
{
  if (sun_alloc == NULL) { return SUN_ERR_ARG_CORRUPT; }

  if (SUNAllocator_NewEmpty(sun_alloc)) { return SUN_ERR_MEM_FAIL; }

  /* Attach method implementations */
  (*sun_alloc)->ops->allocate   = SUNAllocator_Allocate_System;
  (*sun_alloc)->ops->deallocate = SUNAllocator_Deallocate_System;
  (*sun_alloc)->ops->printstats = SUNAllocator_PrintStats_System;
  (*sun_alloc)->ops->destroy    = SUNAllocator_Destroy_System;

  /* Allocate, initialize, and attach member data */
  SUNAllocator_Content_Sys content =
    (SUNAllocator_Content_Sys)malloc(sizeof(*content));
  if (!content)
  {
    (void)SUNAllocator_DestroyEmpty(sun_alloc);
    return SUN_ERR_MALLOC_FAIL;
  }

  content->num_allocations   = 0;
  content->num_deallocations = 0;
  content->bytes_allocated   = 0;
  content->bytes_max         = 0;

  (*sun_alloc)->content = content;

  return SUN_SUCCESS;
}

SUNErrCode SUNAllocator_Destroy_System(SUNAllocator* sun_alloc)
{
  if (sun_alloc == NULL) { return SUN_SUCCESS; }
  if (*sun_alloc == NULL) { return SUN_SUCCESS; }
  if ((*sun_alloc)->content) { free((*sun_alloc)->content); }
  return SUNAllocator_DestroyEmpty(sun_alloc);
}

/* -----------------------------------------------------------------------------
 * Implementations of class methods
 * ---------------------------------------------------------------------------*/

void* SUNAllocator_Allocate_System(SUNAllocator sun_alloc, size_t mem_size,
                                   SUNMemoryType mem_type)
{
  if (mem_type != SUNMEMTYPE_HOST) { return NULL; }

  void* ptr_out = malloc(mem_size);
  if (!ptr_out) { return NULL; }

  SUNAllocator_Content_Sys content = GET_CONTENT(sun_alloc);
  content->num_allocations++;
  content->bytes_allocated += mem_size;
  content->bytes_max = SUNMAX(content->bytes_allocated, content->bytes_max);

  return ptr_out;
}

void SUNAllocator_Deallocate_System(SUNAllocator sun_alloc, void* ptr,
                                    size_t mem_size, SUNMemoryType mem_type)
{
  if (mem_type != SUNMEMTYPE_HOST) { return; }

  free(ptr);
  SUNAllocator_Content_Sys content = GET_CONTENT(sun_alloc);
  content->num_deallocations++;
  content->bytes_allocated -= mem_size;
}

SUNErrCode SUNAllocator_PrintStats_System(SUNAllocator sun_alloc, FILE* outfile,
                                          SUNOutputFormat fmt)
{
  SUNAllocator_Content_Sys content = GET_CONTENT(sun_alloc);

  switch (fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "Number of allocations   = %lu\n", content->num_allocations);
    fprintf(outfile, "Number of deallocations = %lu\n",
            content->num_deallocations);
    fprintf(outfile, "Current bytes allocated = %zu\n", content->bytes_allocated);
    fprintf(outfile, "Max bytes allocated     = %zu\n", content->bytes_max);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, "Number of allocations,%lu", content->num_allocations);
    fprintf(outfile, ",Number of deallocations,%lu", content->num_deallocations);
    fprintf(outfile, ",Current bytes allocated,%zu", content->bytes_allocated);
    fprintf(outfile, ",Max bytes allocated,%zu", content->bytes_max);
    break;
  default: return SUN_ERR_ARG_OUTOFRANGE;
  }

  return SUN_SUCCESS;
}
