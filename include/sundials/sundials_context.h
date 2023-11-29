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
 * SUNDIALS context class. A context object holds data that all
 * SUNDIALS objects in a simulation share. It is thread-safe provided
 * that each thread has its own context object.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_CONTEXT_H
#define _SUNDIALS_CONTEXT_H

#include "sundials/sundials_logger.h"
#include "sundials/sundials_profiler.h"
#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef struct _SUNContext* SUNContext;

SUNDIALS_EXPORT int SUNContext_Create(SUNComm comm, SUNContext* ctx);
SUNDIALS_EXPORT int SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler* profiler);
SUNDIALS_EXPORT int SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler);
SUNDIALS_EXPORT int SUNContext_GetLogger(SUNContext sunctx, SUNLogger* logger);
SUNDIALS_EXPORT int SUNContext_SetLogger(SUNContext sunctx, SUNLogger logger);
SUNDIALS_EXPORT int SUNContext_Free(SUNContext* ctx);


#ifdef __cplusplus
}

/* We include this here for backwards compatibility
   (the contents used to be defined here directly) */
#include <sundials/sundials_context.hpp>

#endif
#endif
