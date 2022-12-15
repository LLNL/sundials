/* ----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
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

#include <sunmemory/sunmemory_sycl.h>
#include "sundials_debug.h"

struct SUNMemoryHelper_Content_Sycl_ {
  unsigned long long  num_allocations_host;
  unsigned long long  num_deallocations_host;
  unsigned long long  num_allocations_device;
  unsigned long long  num_deallocations_device;
  unsigned long long  num_allocations_pinned;
  unsigned long long  num_deallocations_pinned;
  unsigned long long  num_allocations_uvm;
  unsigned long long  num_deallocations_uvm;
  size_t              bytes_allocated_host;
  size_t              bytes_high_watermark_host;
  size_t              bytes_allocated_device;
  size_t              bytes_high_watermark_device;
  size_t              bytes_allocated_pinned;
  size_t              bytes_high_watermark_pinned;
  size_t              bytes_allocated_uvm;
  size_t              bytes_high_watermark_uvm;
};

typedef struct SUNMemoryHelper_Content_Sycl_ SUNMemoryHelper_Content_Sycl;

SUNMemoryHelper SUNMemoryHelper_Sycl(SUNContext sunctx)
{
  // Allocate the helper
  SUNMemoryHelper helper = SUNMemoryHelper_NewEmpty(sunctx);
  if (!helper)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: SUNMemoryHelper_NewEmpty returned NULL\n");
    return NULL;
  }

  // Set the ops
  helper->ops->alloc     = SUNMemoryHelper_Alloc_Sycl;
  helper->ops->dealloc   = SUNMemoryHelper_Dealloc_Sycl;
  helper->ops->copy      = SUNMemoryHelper_Copy_Sycl;
  helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Sycl;

  // Attach content
  helper->content = (SUNMemoryHelper_Content_Sycl*) malloc(sizeof(SUNMemoryHelper_Content_Sycl));
  helper->content->num_allocations_host = 0;
  helper->content->num_deallocations_host = 0;
  helper->content->bytes_allocated_host = 0;
  helper->content->bytes_high_watermark_host = 0;
  helper->content->num_allocations_device = 0;
  helper->content->num_deallocations_device = 0;
  helper->content->bytes_allocated_device = 0;
  helper->content->bytes_high_watermark_device = 0;
  helper->content->num_allocations_pinned = 0;
  helper->content->num_deallocations_pinned = 0;
  helper->content->bytes_allocated_pinned = 0;
  helper->content->bytes_high_watermark_pinned = 0;
  helper->content->num_allocations_uvm = 0;
  helper->content->num_deallocations_uvm = 0;
  helper->content->bytes_allocated_uvm = 0;
  helper->content->bytes_high_watermark_uvm = 0;

  return helper;
}

int SUNMemoryHelper_Alloc_Sycl(SUNMemoryHelper helper, SUNMemory* memptr,
                               size_t mem_size, SUNMemoryType mem_type,
                               void* queue)
{
  // Check inputs
  if (!queue) {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: queue is NULL\n");
    return -1;
  }
  ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

  // Allocate the memory struct
  SUNMemory mem = SUNMemoryNewEmpty();
  if (!mem)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: SUNMemoryNewEmpty returned NULL\n");
    return -1;
  }

  // Initialize the memory content
  mem->ptr  = nullptr;
  mem->own  = SUNTRUE;
  mem->type = mem_type;

  // Allocate the data pointer
  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr = malloc(mem_size);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc returned NULL\n");
      free(mem);
      return -1;
    }
    else
    {
      helper->content->bytes_allocated_host += memsize;
      helper->content->num_allocations_host++;
      helper->content->bytes_high_watermark_host = SUNMAX(helper->bytes_allocated_host, helper->bytes_high_watermark_host);
    }
  }
  else if (mem_type == SUNMEMTYPE_PINNED)
  {
    mem->ptr = ::sycl::malloc_host(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_host returned NULL\n");
      free(mem);
      return -1;
    }
    else
    {
      helper->content->bytes_allocated_pinned += memsize;
      helper->content->num_allocations_pinned++;
      helper->content->bytes_high_watermark_pinned = SUNMAX(helper->bytes_allocated_pinned, helper->bytes_high_watermark_pinned);
    }
  }
  else if (mem_type == SUNMEMTYPE_DEVICE)
  {
    mem->ptr = ::sycl::malloc_device(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_device returned NULL\n");
      free(mem);
      return -1;
    }
    else
    {
      helper->content->bytes_allocated_device += memsize;
      helper->content->num_allocations_device++;
      helper->content->bytes_high_watermark_device = SUNMAX(helper->bytes_allocated_device, helper->bytes_high_watermark_device);
    }
  }
  else if (mem_type == SUNMEMTYPE_UVM)
  {
    mem->ptr = ::sycl::malloc_shared(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_shared returned NULL\n");
      free(mem);
      return -1;
    }
    else
    {
      helper->content->bytes_allocated_uvm += memsize;
      helper->content->num_allocations_uvm++;
      helper->content->bytes_high_watermark_uvm = SUNMAX(helper->bytes_allocated_uvm, helper->bytes_high_watermark_uvm);
    }
  }
  else
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: unknown memory type\n");
    free(mem);
    return -1;
  }

  *memptr = mem;
  return 0;
}

int SUNMemoryHelper_Dealloc_Sycl(SUNMemoryHelper helper, SUNMemory mem,
                                 void* queue)
{
  if (!mem) return 0;

  if (mem->ptr && mem->own)
  {
    if (!queue) {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Sycl: queue is NULL\n");
      return -1;
    }
    ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

    if (mem->type == SUNMEMTYPE_HOST)
    {
      free(mem->ptr);
      mem->ptr = nullptr;
      helper->num_deallocations_host++;
      helper->bytes_allocated_host -= mem->bytes;
    }
    else if (mem->type == SUNMEMTYPE_PINNED ||
             mem->type == SUNMEMTYPE_DEVICE ||
             mem->type == SUNMEMTYPE_UVM)
    {
      ::sycl::free(mem->ptr, *sycl_queue);
      mem->ptr = nullptr;
      if (mem->type == SUNMEMTYPE_PINNED) 
      {
        helper->num_deallocations_pinned++;
        helper->bytes_allocated_pinned -= mem->bytes;
      }
      else if (mem->type == SUNMEMTYPE_DEVICE)
      {
        helper->num_deallocations_device++;
        helper->bytes_allocated_device -= mem->bytes;
      }
      else if (mem->type == SUNMEMTYPE_UVM)
      {
        helper->num_deallocations_uvm++;
        helper->bytes_allocated_uvm -= mem->bytes;
      }
    }
    else
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Sycl: unknown memory type\n");
      return -1;
    }
  }

  free(mem);
  return 0;
}


int SUNMemoryHelper_Copy_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                              SUNMemory src, size_t memory_size, void* queue)
{
  if (!queue) {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Copy_Sycl: queue is NULL\n");
    return -1;
  }
  ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

  if (SUNMemoryHelper_CopyAsync_Sycl(helper, dst, src, memory_size, queue))
    return -1;
  sycl_queue->wait_and_throw();
  return 0;
}


int SUNMemoryHelper_CopyAsync_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                                   SUNMemory src, size_t memory_size,
                                   void* queue)
{
  if (!queue) {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_CopyAsync_Sycl: queue is NULL\n");
    return -1;
  }
  ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

  if (src->type == SUNMEMTYPE_HOST && dst->type == SUNMEMTYPE_HOST)
  {
    memcpy(dst->ptr, src->ptr, memory_size);
  }
  else
  {
    sycl_queue->memcpy(dst->ptr, src->ptr, memory_size);
  }
  return 0;
}

int SUNMemoryHelper_Destroy_Sycl(SUNMemoryHelper helper)
{
  if (helper)
  {
    free(helper->content);
    free(helper);
  }
  return 0;
}

int SUNMemoryHelper_GetHostAllocStats_Sycl(SUNMemoryHelper helper, unsigned long long* num_allocations_host,
                                           unsigned long long* num_deallocations_host, size_t* bytes_allocated_host,
                                           size_t* bytes_high_watermark_host)
{
  *num_allocations_host = helper->num_allocations_host;
  *num_deallocations_host = helper->num_deallocations_host;
  *bytes_allocated_host = helper->bytes_allocated_host;
  *bytes_high_watermark_host = helper->bytes_high_watermark_host;
  return 0;
}

int SUNMemoryHelper_GetPinnedAllocStats_Sycl(SUNMemoryHelper helper, unsigned long long* num_allocations_pinned,
                                             unsigned long long* num_deallocations_pinned, size_t* bytes_allocated_pinned,
                                             size_t* bytes_high_watermark_pinned)
{
  *num_allocations_pinned = helper->num_allocations_pinned;
  *num_deallocations_pinned = helper->num_deallocations_pinned;
  *bytes_allocated_pinned = helper->bytes_allocated_pinned;
  *bytes_high_watermark_pinned = helper->bytes_high_watermark_pinned;
  return 0;
}

int SUNMemoryHelper_GetDeviceAllocStats_Sycl(SUNMemoryHelper helper, unsigned long long* num_allocations_device,
                                             unsigned long long* num_deallocations_device, size_t* bytes_allocated_device,
                                             size_t* bytes_high_watermark_device)
{
  *num_allocations_device = helper->num_allocations_device;
  *num_deallocations_device = helper->num_deallocations_device;
  *bytes_allocated_device = helper->bytes_allocated_device;
  *bytes_high_watermark_device = helper->bytes_high_watermark_device;
  return 0;
}

int SUNMemoryHelper_GetUVMAllocStats_Sycl(SUNMemoryHelper helper, unsigned long long* num_allocations_uvm,
                                          unsigned long long* num_deallocations_uvm, size_t* bytes_allocated_uvm,
                                          size_t* bytes_high_watermark_uvm)
{
  *num_allocations_uvm = helper->num_allocations_uvm;
  *num_deallocations_uvm = helper->num_deallocations_uvm;
  *bytes_allocated_uvm = helper->bytes_allocated_uvm;
  *bytes_high_watermark_uvm = helper->bytes_high_watermark_uvm;
  return 0;
}
