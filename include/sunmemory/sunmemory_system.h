/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS system memory helper header file.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_SYSMEMORY_H
#define _SUNDIALS_SYSMEMORY_H

#include <sundials/sundials_memory.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* Implementation specific functions */

SUNMemoryHelper SUNMemoryHelper_Sys(SUNContext sunctx);

/* SUNMemoryHelper functions */

SUNDIALS_EXPORT int SUNMemoryHelper_Alloc_Sys(SUNMemoryHelper helper, SUNMemory* memptr,
                                              size_t memsize, SUNMemoryType mem_type);

SUNDIALS_EXPORT int SUNMemoryHelper_Dealloc_Sys(SUNMemoryHelper helper, SUNMemory mem);

SUNDIALS_EXPORT int SUNMemoryHelper_Copy_Sys(SUNMemoryHelper helper, SUNMemory dst,
                                             SUNMemory src, size_t memory_size);

#ifdef __cplusplus
}
#endif

#endif
