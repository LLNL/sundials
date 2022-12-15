/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS CUDA memory helper header file.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_CUDAMEMORY_H
#define _SUNDIALS_CUDAMEMORY_H

#include <cuda_runtime.h>
#include <sundials/sundials_memory.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* Implementation specific functions */

SUNDIALS_EXPORT
SUNMemoryHelper SUNMemoryHelper_Cuda(SUNContext sunctx);

/* SUNMemoryHelper functions */

SUNDIALS_EXPORT
int SUNMemoryHelper_Alloc_Cuda(SUNMemoryHelper helper, SUNMemory* memptr,
                               size_t mem_size, SUNMemoryType mem_type,
                               void* queue);

SUNDIALS_EXPORT
SUNMemoryHelper SUNMemoryHelper_Clone_Cuda(SUNMemoryHelper helper);

SUNDIALS_EXPORT
int SUNMemoryHelper_Dealloc_Cuda(SUNMemoryHelper helper, SUNMemory mem,
                                 void* queue);

SUNDIALS_EXPORT
int SUNMemoryHelper_Copy_Cuda(SUNMemoryHelper helper, SUNMemory dst,
                              SUNMemory src, size_t memory_size, void* queue);

SUNDIALS_EXPORT
int SUNMemoryHelper_CopyAsync_Cuda(SUNMemoryHelper helper, SUNMemory dst,
                                   SUNMemory src, size_t memory_size,
                                   void* queue);

SUNDIALS_EXPORT
int SUNMemoryHelper_Destroy_Cuda(SUNMemoryHelper helper);

SUNDIALS_EXPORT
int SUNMemoryHelper_GetHostAllocStats_Cuda(SUNMemoryHelper helper, unsigned long long* num_allocations_host,
                                           unsigned long long* num_deallocations_host, size_t* bytes_allocated_host,
                                           size_t* bytes_high_watermark_host);

SUNDIALS_EXPORT
int SUNMemoryHelper_GetPinnedAllocStats_Cuda(SUNMemoryHelper helper, unsigned long long* num_allocations_pinned,
                                             unsigned long long* num_deallocations_pinned, size_t* bytes_allocated_pinned,
                                             size_t* bytes_high_watermark_pinned);

SUNDIALS_EXPORT
int SUNMemoryHelper_GetDeviceAllocStats_Cuda(SUNMemoryHelper helper, unsigned long long* num_allocations_device,
                                             unsigned long long* num_deallocations_device, size_t* bytes_allocated_device,
                                             size_t* bytes_high_watermark_device);

SUNDIALS_EXPORT                                            
int SUNMemoryHelper_GetUVMAllocStats_Cuda(SUNMemoryHelper helper, unsigned long long* num_allocations_uvm,
                                          unsigned long long* num_deallocations_uvm, size_t* bytes_allocated_uvm,
                                          size_t* bytes_high_watermark_uvm);                                            

#ifdef __cplusplus
}
#endif

#endif
