/* ----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * SUNDIALS SYCL memory helper implementation.
 * ----------------------------------------------------------------*/

#include <cstdlib>
#include <sundials/sundials_math.h>
#include <sunmemory/sunmemory_sycl.h>

#include "sundials_debug.h"

struct SUNMemoryHelper_Content_Sycl_
{
  unsigned long num_allocations_host;
  unsigned long num_deallocations_host;
  unsigned long num_allocations_device;
  unsigned long num_deallocations_device;
  unsigned long num_allocations_pinned;
  unsigned long num_deallocations_pinned;
  unsigned long num_allocations_uvm;
  unsigned long num_deallocations_uvm;
  size_t bytes_allocated_host;
  size_t bytes_high_watermark_host;
  size_t bytes_allocated_device;
  size_t bytes_high_watermark_device;
  size_t bytes_allocated_pinned;
  size_t bytes_high_watermark_pinned;
  size_t bytes_allocated_uvm;
  size_t bytes_high_watermark_uvm;
};

typedef struct SUNMemoryHelper_Content_Sycl_ SUNMemoryHelper_Content_Sycl;

#define SUNHELPER_CONTENT(h) ((SUNMemoryHelper_Content_Sycl*)h->content)

SUNMemoryHelper SUNMemoryHelper_Sycl(SUNContext sunctx)
{
  // Allocate the helper
  SUNMemoryHelper helper = SUNMemoryHelper_NewEmpty(sunctx);
  if (!helper)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: "
                         "SUNMemoryHelper_NewEmpty returned NULL\n");
    return NULL;
  }

  // Set the ops
  helper->ops->alloc         = SUNMemoryHelper_Alloc_Sycl;
  helper->ops->dealloc       = SUNMemoryHelper_Dealloc_Sycl;
  helper->ops->clone         = SUNMemoryHelper_Clone_Sycl;
  helper->ops->getallocstats = SUNMemoryHelper_GetAllocStats_Sycl;
  helper->ops->copy          = SUNMemoryHelper_Copy_Sycl;
  helper->ops->copyasync     = SUNMemoryHelper_CopyAsync_Sycl;
  helper->ops->destroy       = SUNMemoryHelper_Destroy_Sycl;

  // Attach content
  helper->content =
    (SUNMemoryHelper_Content_Sycl*)malloc(sizeof(SUNMemoryHelper_Content_Sycl));
  SUNHELPER_CONTENT(helper)->num_allocations_host        = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations_host      = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated_host        = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark_host   = 0;
  SUNHELPER_CONTENT(helper)->num_allocations_device      = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations_device    = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated_device      = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark_device = 0;
  SUNHELPER_CONTENT(helper)->num_allocations_pinned      = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations_pinned    = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated_pinned      = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark_pinned = 0;
  SUNHELPER_CONTENT(helper)->num_allocations_uvm         = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations_uvm       = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated_uvm         = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark_uvm    = 0;

  return helper;
}

SUNMemoryHelper SUNMemoryHelper_Clone_Sycl(SUNMemoryHelper helper)
{
  SUNMemoryHelper hclone = SUNMemoryHelper_Sycl(helper->sunctx);
  return hclone;
}

int SUNMemoryHelper_Alloc_Sycl(SUNMemoryHelper helper, SUNMemory* memptr,
                               size_t mem_size, SUNMemoryType mem_type,
                               void* queue)
{
  // Check inputs
  if (!queue)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in SUNMemoryHelper_Alloc_Sycl: queue is NULL\n");
    return -1;
  }
  ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

  // Allocate the memory struct
  SUNMemory mem = SUNMemoryNewEmpty(helper->sunctx);
  if (!mem)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in SUNMemoryHelper_Sycl: SUNMemoryNewEmpty returned NULL\n");
    return -1;
  }

  // Initialize the memory content
  mem->ptr   = nullptr;
  mem->own   = SUNTRUE;
  mem->type  = mem_type;
  mem->bytes = mem_size;

  // Allocate the data pointer
  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr = malloc(mem_size);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Sycl: malloc returned NULL\n");
      free(mem);
      return -1;
    }
    else
    {
      SUNHELPER_CONTENT(helper)->bytes_allocated_host += mem_size;
      SUNHELPER_CONTENT(helper)->num_allocations_host++;
      SUNHELPER_CONTENT(helper)->bytes_high_watermark_host =
        SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated_host,
               SUNHELPER_CONTENT(helper)->bytes_high_watermark_host);
    }
  }
  else if (mem_type == SUNMEMTYPE_PINNED)
  {
    mem->ptr = ::sycl::malloc_host(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_host returned NULL\n");
      free(mem);
      return -1;
    }
    else
    {
      SUNHELPER_CONTENT(helper)->bytes_allocated_pinned += mem_size;
      SUNHELPER_CONTENT(helper)->num_allocations_pinned++;
      SUNHELPER_CONTENT(helper)->bytes_high_watermark_pinned =
        SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated_pinned,
               SUNHELPER_CONTENT(helper)->bytes_high_watermark_pinned);
    }
  }
  else if (mem_type == SUNMEMTYPE_DEVICE)
  {
    mem->ptr = ::sycl::malloc_device(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_device returned NULL\n");
      free(mem);
      return -1;
    }
    else
    {
      SUNHELPER_CONTENT(helper)->bytes_allocated_device += mem_size;
      SUNHELPER_CONTENT(helper)->num_allocations_device++;
      SUNHELPER_CONTENT(helper)->bytes_high_watermark_device =
        SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated_device,
               SUNHELPER_CONTENT(helper)->bytes_high_watermark_device);
    }
  }
  else if (mem_type == SUNMEMTYPE_UVM)
  {
    mem->ptr = ::sycl::malloc_shared(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_shared returned NULL\n");
      free(mem);
      return -1;
    }
    else
    {
      SUNHELPER_CONTENT(helper)->bytes_allocated_uvm += mem_size;
      SUNHELPER_CONTENT(helper)->num_allocations_uvm++;
      SUNHELPER_CONTENT(helper)->bytes_high_watermark_uvm =
        SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated_uvm,
               SUNHELPER_CONTENT(helper)->bytes_high_watermark_uvm);
    }
  }
  else
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in SUNMemoryHelper_Alloc_Sycl: unknown memory type\n");
    free(mem);
    return -1;
  }

  *memptr = mem;
  return 0;
}

int SUNMemoryHelper_Dealloc_Sycl(SUNMemoryHelper helper, SUNMemory mem,
                                 void* queue)
{
  if (!mem) { return 0; }

  if (mem->ptr && mem->own)
  {
    if (!queue)
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Dealloc_Sycl: queue is NULL\n");
      return -1;
    }
    ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

    if (mem->type == SUNMEMTYPE_HOST)
    {
      SUNHELPER_CONTENT(helper)->num_deallocations_host++;
      SUNHELPER_CONTENT(helper)->bytes_allocated_host -= mem->bytes;
      free(mem->ptr);
      mem->ptr = nullptr;
    }
    else if (mem->type == SUNMEMTYPE_PINNED || mem->type == SUNMEMTYPE_DEVICE ||
             mem->type == SUNMEMTYPE_UVM)
    {
      if (mem->type == SUNMEMTYPE_PINNED)
      {
        SUNHELPER_CONTENT(helper)->num_deallocations_pinned++;
        SUNHELPER_CONTENT(helper)->bytes_allocated_pinned -= mem->bytes;
      }
      else if (mem->type == SUNMEMTYPE_DEVICE)
      {
        SUNHELPER_CONTENT(helper)->num_deallocations_device++;
        SUNHELPER_CONTENT(helper)->bytes_allocated_device -= mem->bytes;
      }
      else if (mem->type == SUNMEMTYPE_UVM)
      {
        SUNHELPER_CONTENT(helper)->num_deallocations_uvm++;
        SUNHELPER_CONTENT(helper)->bytes_allocated_uvm -= mem->bytes;
      }
      ::sycl::free(mem->ptr, *sycl_queue);
      mem->ptr = nullptr;
    }
    else
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Dealloc_Sycl: unknown memory type\n");
      return -1;
    }
  }

  free(mem);
  return 0;
}

int SUNMemoryHelper_Copy_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                              SUNMemory src, size_t memory_size, void* queue)
{
  if (!queue)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Copy_Sycl: queue is NULL\n");
    return -1;
  }
  ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

  if (SUNMemoryHelper_CopyAsync_Sycl(helper, dst, src, memory_size, queue))
  {
    return -1;
  }
  sycl_queue->wait_and_throw();
  return 0;
}

int SUNMemoryHelper_CopyAsync_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                                   SUNMemory src, size_t memory_size, void* queue)
{
  if (!queue)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in SUNMemoryHelper_CopyAsync_Sycl: queue is NULL\n");
    return -1;
  }
  ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

  if (src->type == SUNMEMTYPE_HOST && dst->type == SUNMEMTYPE_HOST)
  {
    memcpy(dst->ptr, src->ptr, memory_size);
  }
  else { sycl_queue->memcpy(dst->ptr, src->ptr, memory_size); }
  return 0;
}

int SUNMemoryHelper_Destroy_Sycl(SUNMemoryHelper helper)
{
  if (helper)
  {
    if (helper->content) { free(helper->content); }
    if (helper->ops) { free(helper->ops); }
    free(helper);
  }
  return 0;
}

int SUNMemoryHelper_GetAllocStats_Sycl(SUNMemoryHelper helper,
                                       SUNMemoryType mem_type,
                                       unsigned long* num_allocations,
                                       unsigned long* num_deallocations,
                                       size_t* bytes_allocated,
                                       size_t* bytes_high_watermark)
{
  if (mem_type == SUNMEMTYPE_HOST)
  {
    *num_allocations   = SUNHELPER_CONTENT(helper)->num_allocations_host;
    *num_deallocations = SUNHELPER_CONTENT(helper)->num_deallocations_host;
    *bytes_allocated   = SUNHELPER_CONTENT(helper)->bytes_allocated_host;
    *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark_host;
  }
  else if (mem_type == SUNMEMTYPE_PINNED)
  {
    *num_allocations   = SUNHELPER_CONTENT(helper)->num_allocations_pinned;
    *num_deallocations = SUNHELPER_CONTENT(helper)->num_deallocations_pinned;
    *bytes_allocated   = SUNHELPER_CONTENT(helper)->bytes_allocated_pinned;
    *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark_pinned;
  }
  else if (mem_type == SUNMEMTYPE_DEVICE)
  {
    *num_allocations   = SUNHELPER_CONTENT(helper)->num_allocations_device;
    *num_deallocations = SUNHELPER_CONTENT(helper)->num_deallocations_device;
    *bytes_allocated   = SUNHELPER_CONTENT(helper)->bytes_allocated_device;
    *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark_device;
  }
  else if (mem_type == SUNMEMTYPE_UVM)
  {
    *num_allocations      = SUNHELPER_CONTENT(helper)->num_allocations_uvm;
    *num_deallocations    = SUNHELPER_CONTENT(helper)->num_deallocations_uvm;
    *bytes_allocated      = SUNHELPER_CONTENT(helper)->bytes_allocated_uvm;
    *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark_uvm;
  }
  else { return -1; }
  return 0;
}
