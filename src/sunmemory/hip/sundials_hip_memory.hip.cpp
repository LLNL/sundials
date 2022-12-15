/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Sehiprity
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS HIP memory helper implementation.
 * ----------------------------------------------------------------*/

#include <cstdlib>

#include <sunmemory/sunmemory_hip.h>
#include "sundials_debug.h"
#include "sundials_hip.h"

struct SUNMemoryHelper_Content_Hip_ {
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

typedef struct SUNMemoryHelper_Content_Hip_ SUNMemoryHelper_Content_Hip;

SUNMemoryHelper SUNMemoryHelper_Hip(SUNContext sunctx)
{
  SUNMemoryHelper helper;

  /* Allocate the helper */
  helper = SUNMemoryHelper_NewEmpty(sunctx);

  /* Set the ops */
  helper->ops->alloc     = SUNMemoryHelper_Alloc_Hip;
  helper->ops->dealloc   = SUNMemoryHelper_Dealloc_Hip;
  helper->ops->copy      = SUNMemoryHelper_Copy_Hip;
  helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Hip;
  helper->ops->destroy   = SUNMemoryHelper_Destroy_Hip;

  /* Attach content and ops */
  helper->content = (SUNMemoryHelper_Content_Hip*) malloc(sizeof(SUNMemoryHelper_Content_Hip));
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

int SUNMemoryHelper_Alloc_Hip(SUNMemoryHelper helper, SUNMemory* memptr,
                              size_t mem_size, SUNMemoryType mem_type,
                              void* queue)
{
  SUNMemory mem = SUNMemoryNewEmpty();

  mem->ptr  = NULL;
  mem->own  = SUNTRUE;
  mem->type = mem_type;

  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr = malloc(mem_size);
    if (mem->ptr == NULL)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: malloc returned NULL\n");
      free(mem);
      return(-1);
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
    if (!SUNDIALS_HIP_VERIFY(hipMallocHost(&(mem->ptr), mem_size)))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: hipMallocHost failed\n");
      free(mem);
      return(-1);
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
    if (!SUNDIALS_HIP_VERIFY(hipMalloc(&(mem->ptr), mem_size)))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: hipMalloc failed\n");
      free(mem);
      return(-1);
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
    if (!SUNDIALS_HIP_VERIFY(hipMallocManaged(&(mem->ptr), mem_size)))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: hipMallocManaged failed\n");
      free(mem);
      return(-1);
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
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: unknown memory type\n");
    free(mem);
    return(-1);
  }

  *memptr = mem;
  return(0);
}

int SUNMemoryHelper_Dealloc_Hip(SUNMemoryHelper helper, SUNMemory mem,
                                void *queue)
{
  if (mem == NULL) return(0);

  if (mem->ptr != NULL && mem->own)
  {
    if (mem->type == SUNMEMTYPE_HOST)
    {
      free(mem->ptr);
      mem->ptr = NULL;
      helper->num_deallocations_host++;
      helper->bytes_allocated_host -= mem->bytes;
    }
    else if (mem->type == SUNMEMTYPE_PINNED)
    {
      if (!SUNDIALS_HIP_VERIFY(hipFreeHost(mem->ptr)))
      {
        SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Hip: hipFreeHost failed\n");
        return(-1);
      }
      else 
      {
        helper->num_deallocations_pinned++;
        helper->bytes_allocated_pinned -= mem->bytes;
      }
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_DEVICE)
    {
      if (!SUNDIALS_HIP_VERIFY(hipFree(mem->ptr)))
      {
        SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Hip: hipFree failed\n");
        return(-1);
      }
      else 
      {
        helper->num_deallocations_device++;
        helper->bytes_allocated_device -= mem->bytes;
      }
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_UVM)
    {
      if (!SUNDIALS_HIP_VERIFY(hipFree(mem->ptr)))
      {
        SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Hip: hipFree failed\n");
        return(-1);
      }
      else 
      {
        helper->num_deallocations_uvm++;
        helper->bytes_allocated_uvm -= mem->bytes;
      }
      mem->ptr = NULL;
    }
    else
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Hip: unknown memory type\n");
      return(-1);
    }
  }

  free(mem);
  return(0);
}

int SUNMemoryHelper_Copy_Hip(SUNMemoryHelper helper, SUNMemory dst,
                             SUNMemory src, size_t memory_size, void* queue)
{
  int retval = 0;
  hipError_t hiperr = hipSuccess;

  switch(src->type)
  {
    case SUNMEMTYPE_HOST:
    case SUNMEMTYPE_PINNED:
      if (dst->type == SUNMEMTYPE_HOST ||
          dst->type == SUNMEMTYPE_PINNED)
      {
        memcpy(dst->ptr, src->ptr, memory_size);
      }
      else if (dst->type == SUNMEMTYPE_DEVICE ||
               dst->type == SUNMEMTYPE_UVM)
      {
        hiperr = hipMemcpy(dst->ptr, src->ptr,
                           memory_size,
                           hipMemcpyHostToDevice);
      }
      if (!SUNDIALS_HIP_VERIFY(hiperr)) retval = -1;
      break;
    case SUNMEMTYPE_UVM:
    case SUNMEMTYPE_DEVICE:
      if (dst->type == SUNMEMTYPE_HOST ||
          dst->type == SUNMEMTYPE_PINNED)
      {
        hiperr = hipMemcpy(dst->ptr, src->ptr,
                           memory_size,
                           hipMemcpyDeviceToHost);
      }
      else if (dst->type == SUNMEMTYPE_DEVICE ||
               dst->type == SUNMEMTYPE_UVM)
      {
        hiperr = hipMemcpy(dst->ptr, src->ptr,
                           memory_size,
                           hipMemcpyDeviceToDevice);
      }
      if (!SUNDIALS_HIP_VERIFY(hiperr)) retval = -1;
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_CopyAsync_Hip: unknown memory type\n");
      retval = -1;
  }

  return(retval);
}

int SUNMemoryHelper_CopyAsync_Hip(SUNMemoryHelper helper, SUNMemory dst,
                                  SUNMemory src, size_t memory_size,
                                  void* queue)
{
  int retval = 0;
  hipError_t hiperr = hipSuccess;
  hipStream_t stream = 0;

  if (queue != NULL)
  {
    stream = *((hipStream_t*) queue);
  }

  switch(src->type)
  {
    case SUNMEMTYPE_HOST:
    case SUNMEMTYPE_PINNED:
      if (dst->type == SUNMEMTYPE_HOST ||
          dst->type == SUNMEMTYPE_PINNED)
      {
        memcpy(dst->ptr, src->ptr, memory_size);
      }
      else if (dst->type == SUNMEMTYPE_DEVICE ||
               dst->type == SUNMEMTYPE_UVM)
      {
        hiperr = hipMemcpyAsync(dst->ptr, src->ptr,
                                memory_size,
                                hipMemcpyHostToDevice,
                                stream);
      }
      if (!SUNDIALS_HIP_VERIFY(hiperr)) retval = -1;
      break;
    case SUNMEMTYPE_UVM:
    case SUNMEMTYPE_DEVICE:
      if (dst->type == SUNMEMTYPE_HOST ||
          dst->type == SUNMEMTYPE_PINNED)
      {
        hiperr = hipMemcpyAsync(dst->ptr, src->ptr,
                                memory_size,
                                hipMemcpyDeviceToHost,
                                stream);
      }
      else if (dst->type == SUNMEMTYPE_DEVICE ||
              dst->type == SUNMEMTYPE_UVM)
      {
        hiperr = hipMemcpyAsync(dst->ptr, src->ptr,
                                memory_size,
                                hipMemcpyDeviceToDevice,
                                stream);
      }
      if (!SUNDIALS_HIP_VERIFY(hiperr)) retval = -1;
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_CopyAsync_Hip: unknown memory type\n");
      retval = -1;
  }

  return(retval);
}

int SUNMemoryHelper_Destroy_Hip(SUNMemoryHelper helper)
{
  if (helper)
  {
    free(helper->content);
    free(helper);
  }
  return 0;
}

int SUNMemoryHelper_GetHostAllocStatsHip(SUNMemoryHelper helper, unsigned long long* num_allocations_host,
                                         unsigned long long* num_deallocations_host, size_t* bytes_allocated_host,
                                         size_t* bytes_high_watermark_host)
{
  *num_allocations_host = helper->num_allocations_host;
  *num_deallocations_host = helper->num_deallocations_host;
  *bytes_allocated_host = helper->bytes_allocated_host;
  *bytes_high_watermark_host = helper->bytes_high_watermark_host;
  return 0;
}

int SUNMemoryHelper_GetPinnedAllocStatsHip(SUNMemoryHelper helper, unsigned long long* num_allocations_pinned,
                                           unsigned long long* num_deallocations_pinned, size_t* bytes_allocated_pinned,
                                           size_t* bytes_high_watermark_pinned)
{
  *num_allocations_pinned = helper->num_allocations_pinned;
  *num_deallocations_pinned = helper->num_deallocations_pinned;
  *bytes_allocated_pinned = helper->bytes_allocated_pinned;
  *bytes_high_watermark_pinned = helper->bytes_high_watermark_pinned;
  return 0;
}

int SUNMemoryHelper_GetDeviceAllocStatsHip(SUNMemoryHelper helper, unsigned long long* num_allocations_device,
                                           unsigned long long* num_deallocations_device, size_t* bytes_allocated_device,
                                           size_t* bytes_high_watermark_device)
{
  *num_allocations_device = helper->num_allocations_device;
  *num_deallocations_device = helper->num_deallocations_device;
  *bytes_allocated_device = helper->bytes_allocated_device;
  *bytes_high_watermark_device = helper->bytes_high_watermark_device;
  return 0;
}

int SUNMemoryHelper_GetUVMAllocStatsHip(SUNMemoryHelper helper, unsigned long long* num_allocations_uvm,
                                        unsigned long long* num_deallocations_uvm, size_t* bytes_allocated_uvm,
                                        size_t* bytes_high_watermark_uvm)
{
  *num_allocations_uvm = helper->num_allocations_uvm;
  *num_deallocations_uvm = helper->num_deallocations_uvm;
  *bytes_allocated_uvm = helper->bytes_allocated_uvm;
  *bytes_high_watermark_uvm = helper->bytes_high_watermark_uvm;
  return 0;
}
