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

#ifndef _SUNDIALS_ALLOCATOR_SYSTEM_H
#define _SUNDIALS_ALLOCATOR_SYSTEM_H

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_allocator.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNErrCode SUNAllocator_Create_System(SUNAllocator* sun_alloc);

void* SUNAllocator_Allocate_System(SUNAllocator sun_alloc, size_t mem_size,
                                   SUNMemoryType mem_type);

void SUNAllocator_Deallocate_System(SUNAllocator sun_alloc, void* mem_ptr,
                                    size_t mem_size, SUNMemoryType mem_type);

SUNErrCode SUNAllocator_PrintStats_System(SUNAllocator sun_alloc, FILE* outfile,
                                          SUNOutputFormat fmt);

SUNErrCode SUNAllocator_Destroy_System(SUNAllocator* sun_alloc);

#ifdef __cplusplus
}
#endif

#endif
