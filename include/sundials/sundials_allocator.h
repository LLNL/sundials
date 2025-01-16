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
 * SUNDIALS Allocator Class
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ALLOCATOR_H
#define _SUNDIALS_ALLOCATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_export.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------------------
 * SUNAllocator definition
 * ---------------------------------------------------------------------------*/

typedef _SUNDIALS_STRUCT_ SUNAllocatorOps_* SUNAllocatorOps;

struct SUNAllocator_
{
  void* content;       /* class member data */
  SUNAllocatorOps ops; /* vtable of class methods */
};

struct SUNAllocatorOps_
{
  void* (*allocate)(SUNAllocator sun_alloc, size_t mem_size,
                    SUNMemoryType mem_type);
  void (*deallocate)(SUNAllocator sun_alloc, void* mem_ptr, size_t mem_size,
                     SUNMemoryType mem_type);
  SUNErrCode (*printstats)(SUNAllocator sun_alloc, FILE* outfile,
                           SUNOutputFormat fmt);
  SUNErrCode (*destroy)(SUNAllocator* sun_alloc);
};

/* -----------------------------------------------------------------------------
 * Exported SUNAllocator functions
 * ---------------------------------------------------------------------------*/

SUNDIALS_EXPORT
void* SUNAllocator_Allocate(SUNAllocator sun_alloc, size_t mem_size,
                            SUNMemoryType mem_type);

SUNDIALS_EXPORT
void SUNAllocator_Deallocate(SUNAllocator sun_alloc, void* mem_ptr,
                             size_t mem_size, SUNMemoryType mem_type);

SUNDIALS_EXPORT
SUNErrCode SUNAllocator_PrintStats(SUNAllocator sun_alloc, FILE* outfile,
                                   SUNOutputFormat fmt);

SUNDIALS_EXPORT
SUNErrCode SUNAllocator_Destroy(SUNAllocator* sun_alloc);

/* -----------------------------------------------------------------------------
 * Utilities for creating a SUNAllocator implementation
 * ---------------------------------------------------------------------------*/

SUNDIALS_EXPORT
SUNErrCode SUNAllocator_NewEmpty(SUNAllocator* sun_alloc);

SUNDIALS_EXPORT
SUNErrCode SUNAllocator_DestroyEmpty(SUNAllocator* sun_alloc);

#ifdef __cplusplus
}
#endif

#endif
