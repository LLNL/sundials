/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
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

#include <stdio.h>
#include <sundials/priv/sundials_context_impl.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT
SUNErrCode SUNContext_Create(SUNComm comm, SUNContext* sunctx_out);

SUNDIALS_EXPORT
SUNErrCode SUNContext_GetLastError(SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNContext_PeekLastError(SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNContext_PushErrHandler(SUNContext sunctx, SUNErrHandlerFn err_fn,
                                     void* err_user_data);

SUNDIALS_EXPORT
SUNErrCode SUNContext_PopErrHandler(SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNContext_ClearErrHandlers(SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler* profiler);

SUNDIALS_EXPORT
SUNErrCode SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler);

SUNDIALS_EXPORT
SUNErrCode SUNContext_GetLogger(SUNContext sunctx, SUNLogger* logger);

SUNDIALS_EXPORT
SUNErrCode SUNContext_SetLogger(SUNContext sunctx, SUNLogger logger);

SUNDIALS_EXPORT
SUNErrCode SUNContext_CreateTempVecStack(SUNContext sunctx, N_Vector tmpl);

SUNDIALS_EXPORT
SUNErrCode SUNContext_PopTempVector(SUNContext sunctx, N_Vector* vec);

SUNDIALS_EXPORT
SUNErrCode SUNContext_PushTempVector(SUNContext sunctx, N_Vector* vec);

SUNDIALS_EXPORT
SUNErrCode SUNContext_PrintAllocStats(SUNContext sunctx, FILE* outfile,
                                      SUNOutputFormat fmt);

SUNDIALS_EXPORT
SUNErrCode SUNContext_Free(SUNContext* sunctx);

#ifdef __cplusplus
}

#endif
#endif
