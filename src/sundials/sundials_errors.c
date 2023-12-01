/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include "sundials/sundials_errors.h"

#include <stdlib.h>
#include <string.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

static inline char* combineFileAndLine(int line, const char* file)
{
  size_t total_str_len = strlen(file) + 6;
  char* file_and_line  = malloc(total_str_len * sizeof(char));
  snprintf(file_and_line, total_str_len, "%s:%d", file, line);
  return file_and_line;
}

SUNErrCode SUNErrHandler_Create(SUNErrHandlerFn eh_fn, void* eh_data,
                                SUNErrHandler* eh_out)
{
  SUNErrHandler eh = NULL;

  eh = (SUNErrHandler)malloc(sizeof(struct SUNErrHandler_));
  if (!eh) { return SUN_ERR_MALLOC_FAIL; }

  eh->previous = NULL;
  eh->call     = eh_fn;
  eh->data     = eh_data;

  *eh_out = eh;
  return SUN_SUCCESS;
}

void SUNErrHandler_Destroy(SUNErrHandler eh) { free(eh); }

const char* SUNGetErrMsg(SUNErrCode code, SUNContext sunctx)
{
#define SUN_EXPAND_TO_CASES(name, description) \
  case name: return description; break;

  switch (code)
  {
    SUN_ERR_CODE_LIST(SUN_EXPAND_TO_CASES)
  default: return "unknown error";
  }

  return NULL;
}

void SUNLogErrHandlerFn(int line, const char* func, const char* file,
                        const char* msg, SUNErrCode err_code,
                        void* err_user_data, SUNContext sunctx)
{
  char* file_and_line = combineFileAndLine(line, file);
  if (msg == NULL) { msg = SUNGetErrMsg(err_code, sunctx); }
  SUNLogger_QueueMsg(sunctx->logger, SUN_LOGLEVEL_ERROR, file_and_line, func,
                     msg);
  free(file_and_line);
}

void SUNAbortErrHandlerFn(int line, const char* func, const char* file,
                          const char* msg, SUNErrCode err_code,
                          void* err_user_data, SUNContext sunctx)
{
  char* file_and_line = combineFileAndLine(line, file);
  SUNLogger_QueueMsg(sunctx->logger, SUN_LOGLEVEL_ERROR, file_and_line, func,
                     "SUNAbortErrHandler: Calling abort now, use a different "
                     "error handler to avoid program termination.\n");
  free(file_and_line);
  abort();
}

void SUNAssertErrHandlerFn(int line, const char* func, const char* file,
                           const char* stmt, SUNErrCode err_code,
                           void* err_user_data, SUNContext sunctx)
{
  char* file_and_line = combineFileAndLine(line, file);
  SUNLogger_QueueMsg(sunctx->logger, SUN_LOGLEVEL_ERROR, file_and_line, func,
                     "SUNAssertErrHandler: assert(%s) failed... terminating\n",
                     stmt);
  free(file_and_line);
  abort();
}
