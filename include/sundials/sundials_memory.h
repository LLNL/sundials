/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * SUNDIALS memory helpers and types.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_MEMORY_H
#define _SUNDIALS_MEMORY_H

#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef enum
{
  SUNMEMTYPE_HOST,   /* pageable memory accessible on the host     */
  SUNMEMTYPE_PINNED, /* page-locked memory accesible on the host   */
  SUNMEMTYPE_DEVICE, /* memory accessible from the device          */
  SUNMEMTYPE_UVM     /* memory accessible from the host or device  */
} SUNMemoryType;

/*
 * SUNMemory is a simple abstraction of a pointer to some
 * contiguos memory, so that we can keep track of its type
 * and its ownership.
 */

typedef struct SUNMemory_* SUNMemory;

struct SUNMemory_
{
  void* ptr;
  SUNMemoryType type;
  sunbooleantype own;
  size_t bytes;
};

/* Creates a new SUNMemory object with a NULL ptr */
SUNDIALS_EXPORT SUNMemory SUNMemoryNewEmpty(SUNContext sunctx);

/*
 * SUNMemoryHelper holds ops which can allocate, deallocate,
 * and copy SUNMemory.
 */

typedef struct SUNMemoryHelper_Ops_* SUNMemoryHelper_Ops;
typedef struct SUNMemoryHelper_* SUNMemoryHelper;

struct SUNMemoryHelper_
{
  void* content;
  SUNMemoryHelper_Ops ops;
  SUNContext sunctx;
};

struct SUNMemoryHelper_Ops_
{
  /* operations that implementations are required to provide */
  SUNErrCode (*alloc)(SUNMemoryHelper, SUNMemory* memptr, size_t mem_size,
                      SUNMemoryType mem_type, void* queue);
  SUNErrCode (*dealloc)(SUNMemoryHelper, SUNMemory mem, void* queue);
  SUNErrCode (*copy)(SUNMemoryHelper, SUNMemory dst, SUNMemory src,
                     size_t mem_size, void* queue);

  /* operations that provide default implementations */
  SUNErrCode (*copyasync)(SUNMemoryHelper, SUNMemory dst, SUNMemory src,
                          size_t mem_size, void* queue);
  SUNErrCode (*getallocstats)(SUNMemoryHelper, SUNMemoryType mem_type,
                              unsigned long* num_allocations,
                              unsigned long* num_deallocations,
                              size_t* bytes_allocated,
                              size_t* bytes_high_watermark);
  SUNMemoryHelper (*clone)(SUNMemoryHelper);
  SUNErrCode (*destroy)(SUNMemoryHelper);
};

/*
 * Generic SUNMemoryHelper functions that work without a SUNMemoryHelper object.
 */

/* Creates a new SUNMemory object which points to the same data as another
 * SUNMemory object.
 * The SUNMemory returned will not own the ptr, therefore, it will not free
 * the ptr in Dealloc. */
SUNDIALS_EXPORT
SUNMemory SUNMemoryHelper_Alias(SUNMemoryHelper, SUNMemory mem);

/* Creates a new SUNMemory object with ptr set to the user provided pointer
 * The SUNMemory returned will not own the ptr, therefore, it will not free
 * the ptr in Dealloc. */
SUNDIALS_EXPORT
SUNMemory SUNMemoryHelper_Wrap(SUNMemoryHelper, void* ptr,
                               SUNMemoryType mem_type);

/*
 * Required SUNMemoryHelper operations.
 */

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Alloc(SUNMemoryHelper, SUNMemory* memptr,
                                 size_t mem_size, SUNMemoryType mem_type,
                                 void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Dealloc(SUNMemoryHelper, SUNMemory mem, void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Copy(SUNMemoryHelper, SUNMemory dst, SUNMemory src,
                                size_t mem_size, void* queue);

/*
 * Optional SUNMemoryHelper operations.
 */

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_CopyAsync(SUNMemoryHelper, SUNMemory dst,
                                     SUNMemory src, size_t mem_size, void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_GetAllocStats(SUNMemoryHelper, SUNMemoryType mem_type,
                                         unsigned long* num_allocations,
                                         unsigned long* num_deallocations,
                                         size_t* bytes_allocated,
                                         size_t* bytes_high_watermark);

/* Clones the SUNMemoryHelper */
SUNDIALS_EXPORT
SUNMemoryHelper SUNMemoryHelper_Clone(SUNMemoryHelper);

/* Frees the SUNMemoryHelper */
SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Destroy(SUNMemoryHelper);

/*
 * Utility SUNMemoryHelper functions.
 */

/* Creates an empty SUNMemoryHelper object */
SUNDIALS_EXPORT
SUNMemoryHelper SUNMemoryHelper_NewEmpty(SUNContext sunctx);

/* Copyies the SUNMemoryHelper ops structure from src->ops to dst->ops. */
SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_CopyOps(SUNMemoryHelper src, SUNMemoryHelper dst);

/* Checks that all required SUNMemoryHelper ops are provided */
SUNDIALS_EXPORT
sunbooleantype SUNMemoryHelper_ImplementsRequiredOps(SUNMemoryHelper);

#ifdef __cplusplus
}
#endif

#endif
