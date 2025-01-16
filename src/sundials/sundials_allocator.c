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
 * SUNDIALS allocator class
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "sundials/sundials_allocator.h"
#include "sundials/sundials_errors.h"

/* -----------------------------------------------------------------------------
 * Exported SUNAllocator functions
 * ---------------------------------------------------------------------------*/

void* SUNAllocator_Allocate(SUNAllocator sun_alloc, size_t mem_size,
                            SUNMemoryType mem_type)
{
  return sun_alloc->ops->allocate(sun_alloc, mem_size, mem_type);
}

void SUNAllocator_Deallocate(SUNAllocator sun_alloc, void* mem_ptr,
                             size_t mem_size, SUNMemoryType mem_type)
{
  sun_alloc->ops->deallocate(sun_alloc, mem_ptr, mem_size, mem_type);
}

SUNErrCode SUNAllocator_PrintStats(SUNAllocator sun_alloc, FILE* outfile,
                                   SUNOutputFormat fmt)
{
  if (sun_alloc->ops->printstats)
  {
    return sun_alloc->ops->printstats(sun_alloc, outfile, fmt);
  }
  else { return SUN_SUCCESS; }
}

SUNErrCode SUNAllocator_Destroy(SUNAllocator* sun_alloc)
{
  if (sun_alloc == NULL) { return SUN_SUCCESS; }
  if (*sun_alloc == NULL) { return SUN_SUCCESS; }

  if ((*sun_alloc)->ops->destroy)
  {
    return (*sun_alloc)->ops->destroy(sun_alloc);
  }
  else { return SUNAllocator_DestroyEmpty(sun_alloc); }

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Utilities for creating a SUNAllocator implementation
 * ---------------------------------------------------------------------------*/

SUNErrCode SUNAllocator_NewEmpty(SUNAllocator* sun_alloc)
{
  if (sun_alloc == NULL) { return SUN_ERR_ARG_CORRUPT; }

  (*sun_alloc) = (SUNAllocator)malloc(sizeof(struct SUNAllocator_));
  if (*sun_alloc == NULL) { return SUN_ERR_MALLOC_FAIL; }
  memset((*sun_alloc), 0, sizeof(struct SUNAllocator_));

  (*sun_alloc)->ops = (SUNAllocatorOps)malloc(sizeof(struct SUNAllocatorOps_));
  if ((*sun_alloc)->ops == NULL)
  {
    (void)SUNAllocator_DestroyEmpty(sun_alloc);
    return SUN_ERR_MALLOC_FAIL;
  }
  memset((*sun_alloc)->ops, 0, sizeof(struct SUNAllocatorOps_));

  return SUN_SUCCESS;
}

SUNErrCode SUNAllocator_DestroyEmpty(SUNAllocator* sun_alloc)
{
  if (sun_alloc == NULL) { return SUN_SUCCESS; }
  if (*sun_alloc == NULL) { return SUN_SUCCESS; }

  /* not empty, implementations should free their own content */
  if ((*sun_alloc)->content) { return SUN_ERR_NOT_IMPLEMENTED; }

  if ((*sun_alloc)->ops) { free((*sun_alloc)->ops); }
  (*sun_alloc)->ops     = NULL;
  (*sun_alloc)->content = NULL;
  free(*sun_alloc);
  *sun_alloc = NULL;

  return SUN_SUCCESS;
}
