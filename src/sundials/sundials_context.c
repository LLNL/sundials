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
 * SUNDIALS objects in a simulation share.
 * ----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_profiler.h>

#include "sundials_context_impl.h"
#include "sundials_debug.h"

SUNErrCode SUNContext_Create(void* comm, SUNContext* sunctx_ptr)
{
  SUNProfiler profiler = NULL;
  SUNLogger logger     = NULL;

  SUNContext sunctx;

  sunctx = NULL;
  sunctx = (SUNContext)malloc(sizeof(struct SUNContext_));

  /* SUNContext_Create cannot assert or log since the SUNContext is not yet
   * created */
  if (!sunctx) {
    return SUN_ERR_MALLOC_FAIL;
  }

#if SUNDIALS_LOGGING_LEVEL > 0
#if defined(SUNDIALS_LOGGING_ENABLE_MPI)
  if (SUNLogger_CreateFromEnv(comm, &logger)) {
    return SUN_ERR_LOGGER_CORRUPT;
  }
#else
  if (SUNLogger_CreateFromEnv(NULL, &logger)) {
    return SUN_ERR_LOGGER_CORRUPT;
  }
#endif
#else
  if (SUNLogger_Create(NULL, 0, &logger)) {
    return SUN_ERR_LOGGER_CORRUPT;
  }
  SUNCheckCallReturn(SUNLogger_SetErrorFilename(logger, ""), sunctx);
  SUNCheckCallReturn(SUNLogger_SetWarningFilename(logger, ""), sunctx);
  SUNCheckCallReturn(SUNLogger_SetInfoFilename(logger, ""), sunctx);
  SUNCheckCallReturn(SUNLogger_SetDebugFilename(logger, ""), sunctx);
#endif

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  SUNCheckCallReturn(SUNProfiler_Create(comm, "SUNContext Default", &profiler), sunctx);
#endif

  sunctx->logger        = logger;
  sunctx->own_logger    = logger != NULL;
  sunctx->profiler      = profiler;
  sunctx->own_profiler  = profiler != NULL;
  sunctx->last_err      = 0;
  sunctx->err_handler   = SUNErrHandler_Create(SUNLogErrHandlerFn, NULL);
  sunctx->comm          = comm;

  *sunctx_ptr = sunctx;

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_GetLastError(SUNContext sunctx, SUNErrCode* last_err)
{
  *last_err = sunctx->last_err;
  return SUN_SUCCESS;
}

SUNErrHandler SUNContext_PushErrHandler(SUNContext sunctx, SUNErrHandlerFn err_fn, void* err_user_data)
{
  SUNErrHandler new_err_handler = SUNErrHandler_Create(err_fn, err_user_data);
  new_err_handler->previous = sunctx->err_handler;
  sunctx->err_handler = new_err_handler;
  return new_err_handler;
}

SUNErrCode SUNContext_PopErrHandler(SUNContext sunctx)
{
  if (sunctx->err_handler) {
    sunctx->err_handler = sunctx->err_handler->previous;
  }
  free(sunctx->err_handler);
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler* profiler)
{
#ifdef SUNDIALS_BUILD_WITH_PROFILING
  /* get profiler */
  *profiler = sunctx->profiler;
#else
  *profiler = NULL;
#endif

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler)
{
#ifdef SUNDIALS_BUILD_WITH_PROFILING
  /* free any existing profiler */
  if (sunctx->profiler && sunctx->own_profiler) {
    if (SUNProfiler_Free(&(sunctx->profiler))) return -1;
    sunctx->profiler = NULL;
  }

  /* set profiler */
  sunctx->profiler     = profiler;
  sunctx->own_profiler = SUNFALSE;
#endif

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_GetLogger(SUNContext sunctx, SUNLogger* logger)
{
  /* get logger */
  *logger = sunctx->logger;
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_SetLogger(SUNContext sunctx, SUNLogger logger)
{
  /* free any existing logger */
  if (sunctx->logger && sunctx->own_logger) {
    if (SUNLogger_Destroy(&(sunctx->logger))) {
      return -1;
    }
    sunctx->logger = NULL;
  }

  /* set logger */
  sunctx->logger     = logger;
  sunctx->own_logger = SUNFALSE;

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_Free(SUNContext* sunctx)
{
#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  FILE* fp;
  char* sunprofiler_print_env;
#endif

  if (!sunctx || !(*sunctx)) {
    return SUN_SUCCESS;
  }

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  /* Find out where we are printing to */
  sunprofiler_print_env = getenv("SUNPROFILER_PRINT");
  fp                    = NULL;
  if (sunprofiler_print_env) {
    if (!strcmp(sunprofiler_print_env, "0")) fp = NULL;
    else if (!strcmp(sunprofiler_print_env, "1") || !strcmp(sunprofiler_print_env, "TRUE") ||
             !strcmp(sunprofiler_print_env, "stdout"))
      fp = stdout;
    else
      fp = fopen(sunprofiler_print_env, "a");
  }

  /* Enforce that the profiler is freed before finalizing,
     if it is not owned by the sunctx. */
  if ((*sunctx)->profiler) {
    if (fp) SUNProfiler_Print((*sunctx)->profiler, fp);
    if (fp) fclose(fp);
    if ((*sunctx)->own_profiler) SUNProfiler_Free(&(*sunctx)->profiler);
  }
#endif

  if ((*sunctx)->logger && (*sunctx)->own_logger) {
    SUNLogger_Destroy(&(*sunctx)->logger);
  }

  SUNErrHandler_Destroy((*sunctx)->err_handler);

  free(*sunctx);
  *sunctx = NULL;

  return SUN_SUCCESS;
}
