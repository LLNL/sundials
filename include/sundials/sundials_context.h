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

struct SUNErrHandler_;

struct SUNContext_ {
  SUNProfiler profiler;
  sunbooleantype own_profiler;
  SUNLogger logger;
  sunbooleantype own_logger;
  SUNErrCode last_err;
  struct SUNErrHandler_ *err_handler;
  void *comm;
};

typedef struct SUNContext_ *SUNContext;

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT
SUNErrCode SUNContext_Create(void *comm, SUNContext *ctx);

SUNDIALS_EXPORT
SUNErrCode SUNContext_GetLastError(SUNContext sunctx, SUNErrCode *last_err);

SUNDIALS_EXPORT
struct SUNErrHandler_ *SUNContext_PushErrHandler(
    SUNContext sunctx,
    int (*err_fn)(int line, const char *func, const char *file, const char *msg,
                  SUNErrCode err_code, void *err_user_data,
                  struct SUNContext_ *sunctx),
    void *err_user_data);

SUNDIALS_EXPORT
SUNErrCode SUNContext_PopErrHandler(SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler *profiler);

SUNDIALS_EXPORT
SUNErrCode SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler);

SUNDIALS_EXPORT
SUNErrCode SUNContext_GetLogger(SUNContext sunctx, SUNLogger *logger);

SUNDIALS_EXPORT
SUNErrCode SUNContext_SetLogger(SUNContext sunctx, SUNLogger logger);

SUNDIALS_EXPORT
SUNErrCode SUNContext_Free(SUNContext *ctx);

#ifdef __cplusplus
}

/* We include this here for backwards compatibility
   (the contents used to be defined here directly) */
#include <sundials/sundials_context.hpp>

#endif
#endif
