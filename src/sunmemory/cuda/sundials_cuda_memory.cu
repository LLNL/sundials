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
 * SUNDIALS CUDA memory helper implementation.
 * ----------------------------------------------------------------*/

#include <cstdlib>
#include <sundials/sundials_math.h>
#include <sunmemory/sunmemory_cuda.h>

#include "sundials_cuda.h"
#include "sundials_debug.h"

struct SUNMemoryHelper_Content_Cuda_
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

typedef struct SUNMemoryHelper_Content_Cuda_ SUNMemoryHelper_Content_Cuda;

#define SUNHELPER_CONTENT(h) ((SUNMemoryHelper_Content_Cuda*)h->content)

SUNMemoryHelper SUNMemoryHelper_Cuda(SUNContext sunctx)
{
  SUNMemoryHelper helper;

  /* Allocate the helper */
  helper = SUNMemoryHelper_NewEmpty(sunctx);

  /* Set the ops */
  helper->ops->alloc         = SUNMemoryHelper_Alloc_Cuda;
  helper->ops->dealloc       = SUNMemoryHelper_Dealloc_Cuda;
  helper->ops->copy          = SUNMemoryHelper_Copy_Cuda;
  helper->ops->copyasync     = SUNMemoryHelper_CopyAsync_Cuda;
  helper->ops->getallocstats = SUNMemoryHelper_GetAllocStats_Cuda;
  helper->ops->clone         = SUNMemoryHelper_Clone_Cuda;
  helper->ops->destroy       = SUNMemoryHelper_Destroy_Cuda;

  /* Attach content */
  helper->content =
    (SUNMemoryHelper_Content_Cuda*)malloc(sizeof(SUNMemoryHelper_Content_Cuda));
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

SUNMemoryHelper SUNMemoryHelper_Clone_Cuda(SUNMemoryHelper helper)
{
  SUNMemoryHelper hclone = SUNMemoryHelper_Cuda(helper->sunctx);
  return hclone;
}

int SUNMemoryHelper_Alloc_Cuda(SUNMemoryHelper helper, SUNMemory* memptr,
                               size_t mem_size, SUNMemoryType mem_type,
                               void* queue)
{
  SUNMemory mem = SUNMemoryNewEmpty();

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
        "ERROR in SUNMemoryHelper_Alloc_Cuda: malloc returned NULL\n");
      free(mem);
      return (-1);
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
    if (!SUNDIALS_CUDA_VERIFY(cudaMallocHost(&(mem->ptr), mem_size)))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Cuda: cudaMallocHost failed\n");
      free(mem);
      return (-1);
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
    if (!SUNDIALS_CUDA_VERIFY(cudaMalloc(&(mem->ptr), mem_size)))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Cuda: cudaMalloc failed\n");
      free(mem);
      return (-1);
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
    if (!SUNDIALS_CUDA_VERIFY(cudaMallocManaged(&(mem->ptr), mem_size)))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Cuda: cudaMallocManaged failed\n");
      free(mem);
      return (-1);
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
      "ERROR in SUNMemoryHelper_Alloc_Cuda: unknown memory type\n");
    free(mem);
    return (-1);
  }

  *memptr = mem;
  return (0);
}

int SUNMemoryHelper_Dealloc_Cuda(SUNMemoryHelper helper, SUNMemory mem,
                                 void* queue)
{
  if (mem == NULL) { return (0); }

  if (mem->ptr != NULL && mem->own)
  {
    if (mem->type == SUNMEMTYPE_HOST)
    {
      SUNHELPER_CONTENT(helper)->num_deallocations_host++;
      SUNHELPER_CONTENT(helper)->bytes_allocated_host -= mem->bytes;
      free(mem->ptr);
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_PINNED)
    {
      SUNHELPER_CONTENT(helper)->num_deallocations_pinned++;
      SUNHELPER_CONTENT(helper)->bytes_allocated_pinned -= mem->bytes;
      if (!SUNDIALS_CUDA_VERIFY(cudaFreeHost(mem->ptr)))
      {
        SUNDIALS_DEBUG_PRINT(
          "ERROR in SUNMemoryHelper_Dealloc_Cuda: cudaFreeHost failed\n");
        return (-1);
      }
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_DEVICE)
    {
      SUNHELPER_CONTENT(helper)->num_deallocations_device++;
      SUNHELPER_CONTENT(helper)->bytes_allocated_device -= mem->bytes;
      if (!SUNDIALS_CUDA_VERIFY(cudaFree(mem->ptr)))
      {
        SUNDIALS_DEBUG_PRINT(
          "ERROR in SUNMemoryHelper_Dealloc_Cuda: cudaFree failed\n");
        return (-1);
      }
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_UVM)
    {
      SUNHELPER_CONTENT(helper)->num_deallocations_uvm++;
      SUNHELPER_CONTENT(helper)->bytes_allocated_uvm -= mem->bytes;
      if (!SUNDIALS_CUDA_VERIFY(cudaFree(mem->ptr)))
      {
        SUNDIALS_DEBUG_PRINT(
          "ERROR in SUNMemoryHelper_Dealloc_Cuda: cudaFree failed\n");
        return (-1);
      }
      mem->ptr = NULL;
    }
    else
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Dealloc_Cuda: unknown memory type\n");
      return (-1);
    }
  }

  free(mem);
  return (0);
}

int SUNMemoryHelper_Copy_Cuda(SUNMemoryHelper helper, SUNMemory dst,
                              SUNMemory src, size_t memory_size, void* queue)
{
  int retval        = 0;
  cudaError_t cuerr = cudaSuccess;

  switch (src->type)
  {
  case SUNMEMTYPE_HOST:
  case SUNMEMTYPE_PINNED:
    if (dst->type == SUNMEMTYPE_HOST || dst->type == SUNMEMTYPE_PINNED)
    {
      memcpy(dst->ptr, src->ptr, memory_size);
    }
    else if (dst->type == SUNMEMTYPE_DEVICE || dst->type == SUNMEMTYPE_UVM)
    {
      cuerr = cudaMemcpy(dst->ptr, src->ptr, memory_size, cudaMemcpyHostToDevice);
    }
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) { retval = -1; }
    break;
  case SUNMEMTYPE_UVM:
  case SUNMEMTYPE_DEVICE:
    if (dst->type == SUNMEMTYPE_HOST || dst->type == SUNMEMTYPE_PINNED)
    {
      cuerr = cudaMemcpy(dst->ptr, src->ptr, memory_size, cudaMemcpyDeviceToHost);
    }
    else if (dst->type == SUNMEMTYPE_DEVICE || dst->type == SUNMEMTYPE_UVM)
    {
      cuerr = cudaMemcpy(dst->ptr, src->ptr, memory_size,
                         cudaMemcpyDeviceToDevice);
    }
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) { retval = -1; }
    break;
  default:
    SUNDIALS_DEBUG_PRINT(
      "ERROR in SUNMemoryHelper_CopyAsync_Cuda: unknown memory type\n");
    retval = -1;
  }

  return (retval);
}

int SUNMemoryHelper_CopyAsync_Cuda(SUNMemoryHelper helper, SUNMemory dst,
                                   SUNMemory src, size_t memory_size, void* queue)
{
  int retval          = 0;
  cudaError_t cuerr   = cudaSuccess;
  cudaStream_t stream = 0;

  if (queue != NULL) { stream = *((cudaStream_t*)queue); }

  switch (src->type)
  {
  case SUNMEMTYPE_HOST:
  case SUNMEMTYPE_PINNED:
    if (dst->type == SUNMEMTYPE_HOST || dst->type == SUNMEMTYPE_PINNED)
    {
      memcpy(dst->ptr, src->ptr, memory_size);
    }
    else if (dst->type == SUNMEMTYPE_DEVICE || dst->type == SUNMEMTYPE_UVM)
    {
      cuerr = cudaMemcpyAsync(dst->ptr, src->ptr, memory_size,
                              cudaMemcpyHostToDevice, stream);
    }
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) { retval = -1; }
    break;
  case SUNMEMTYPE_UVM:
  case SUNMEMTYPE_DEVICE:
    if (dst->type == SUNMEMTYPE_HOST || dst->type == SUNMEMTYPE_PINNED)
    {
      cuerr = cudaMemcpyAsync(dst->ptr, src->ptr, memory_size,
                              cudaMemcpyDeviceToHost, stream);
    }
    else if (dst->type == SUNMEMTYPE_DEVICE || dst->type == SUNMEMTYPE_UVM)
    {
      cuerr = cudaMemcpyAsync(dst->ptr, src->ptr, memory_size,
                              cudaMemcpyDeviceToDevice, stream);
    }
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) { retval = -1; }
    break;
  default:
    SUNDIALS_DEBUG_PRINT(
      "ERROR in SUNMemoryHelper_CopyAsync_Cuda: unknown memory type\n");
    retval = -1;
  }

  return (retval);
}

int SUNMemoryHelper_Destroy_Cuda(SUNMemoryHelper helper)
{
  if (helper)
  {
    if (helper->content) { free(helper->content); }
    if (helper->ops) { free(helper->ops); }
    free(helper);
  }
  return 0;
}

int SUNMemoryHelper_GetAllocStats_Cuda(SUNMemoryHelper helper,
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
