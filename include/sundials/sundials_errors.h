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
 * -----------------------------------------------------------------
 * Contains all error handling interfaces for SUNDIALS that do
 * not depend on MPI.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_ERRORS_H
#define _SUNDIALS_ERRORS_H

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>

/* ----------------------------------------------------------------------------
 * Error code definitions
 * ---------------------------------------------------------------------------*/

/* SUN_ERR_CODE_LIST is an X macro that can be expanded in various ways */
#define SUN_ERR_CODE_LIST(ENTRY)                                               \
  ENTRY(SUN_ERR_ARG_CORRUPT, "argument provided is NULL or corrupted")         \
  ENTRY(SUN_ERR_ARG_OUTOFRANGE, "argument is out of the valid range")          \
  ENTRY(SUN_ERR_PROFILER_MAPFULL,                                              \
        "too many SUNProfiler entries, try setting SUNPROFILER_MAX_ENTRIES "   \
        "in environment to a bigger number")                                   \
  ENTRY(SUN_ERR_PROFILER_MAPGET, "unknown error getting SUNProfiler timer")    \
  ENTRY(SUN_ERR_PROFILER_MAPINSERT,                                            \
        "unknown error inserting SUNProfiler timer")                           \
  ENTRY(SUN_ERR_PROFILER_MAPKEYNOTFOUND, "timer was not found in SUNProfiler") \
  ENTRY(SUN_ERR_PROFILER_MAPSORT, "error sorting SUNProfiler map")             \
  ENTRY(SUN_ERR_LOGGER_CORRUPT, "SUNLogger is NULL or corrupt")                \
  ENTRY(SUN_ERR_LOGGER_CANNOTOPENFILE,                                         \
        "File provided to SUNLogger could not be opened")                      \
  ENTRY(SUN_ERR_MALLOC_FAIL, "malloc returned NULL")                           \
  ENTRY(SUN_ERR_MANYVECTOR_COMMNOTSAME,                                        \
        "not all subvectors have the same MPI_Comm")                           \
  ENTRY(SUN_ERR_MANYVECTOR_COMMNULL, "MPI_Comm is NULL")                       \
  ENTRY(SUN_ERR_MPI_FAIL,                                                      \
        "the MPI call returned something other than MPI_SUCCESS")              \
  ENTRY(SUN_ERR_NOT_IMPLEMENTED,                                               \
        "operation is not implemented: function pointer is NULL")              \
  ENTRY(SUN_ERR_SUNCTX_CORRUPT, "SUNContext is NULL or corrupt")               \
  ENTRY(SUN_ERR_GENERIC, "")                                                   \
  ENTRY(SUN_ERR_UNKNOWN, "Unknown error occured: open an issue at "            \
                         "https://github.com/LLNL/sundials")

/* Expand SUN_ERR_CODE_LIST to enum */
#define SUN_EXPAND_TO_ENUM(name, description) name,

/* SUNErrorCode range is [-1000, -2000] to avoid conflicts wiht package error
   codes, and old/deprecated codes for matrix and (non)linear solvers. */
enum
{
  SUN_ERR_MINIMUM                                       = -2000,
  SUN_ERR_CODE_LIST(SUN_EXPAND_TO_ENUM) SUN_ERR_MAXIMUM = -1000,
  SUN_SUCCESS                                           = 0
};

/* ----------------------------------------------------------------------------
 * Error handler definitions
 * ---------------------------------------------------------------------------*/

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef int (*SUNErrHandlerFn)(int line, const char* func, const char* file,
                               const char* msg, SUNErrCode err_code,
                               void* err_user_data, SUNContext sunctx);

struct SUNErrHandler_
{
  struct SUNErrHandler_* previous; /* next error handler to call (singly
                                      linked-list) */
  SUNErrHandlerFn call;
  void* data;
};

typedef struct SUNErrHandler_* SUNErrHandler;

SUNErrHandler SUNErrHandler_Create(SUNErrHandlerFn eh_fn, void* eh_data);

void SUNErrHandler_Destroy(SUNErrHandler eh);

int SUNLogErrHandlerFn(int line, const char* func, const char* file,
                       const char* msg, SUNErrCode err_code,
                       void* err_user_data, SUNContext sunctx);

int SUNAbortErrHandlerFn(int line, const char* func, const char* file,
                         const char* msg, SUNErrCode err_code,
                         void* err_user_data, SUNContext sunctx);

int SUNAssertErrHandlerFn(int line, const char* func, const char* file,
                          const char* msg, SUNErrCode err_code,
                          void* err_user_data, SUNContext sunctx);

int SUNMPIAbortErrHandlerFn(int line, const char* func, const char* file,
                            const char* msg, SUNErrCode err_code,
                            void* err_user_data, SUNContext sunctx);

int SUNMPIAssertErrHandlerFn(int line, const char* func, const char* file,
                             const char* msg, SUNErrCode err_code,
                             void* err_user_data, SUNContext sunctx);

/* ----------------------------------------------------------------------------
 * Error functions
 * ---------------------------------------------------------------------------*/

/* Turn error code into error message */
const char* SUNGetErrMsg(SUNErrCode code, SUNContext sunctx);

/* Alternative function to SUNContext_GetLastError that is more concise. */
static inline SUNErrCode SUNLastErr(SUNContext sunctx)
{
  SUNErrCode code;
  SUNContext_GetLastError(sunctx, &code);
  return code;
}

static inline int SUNHandleErr(int line, const char* func, const char* file,
                               SUNErrCode code, SUNContext sunctx)
{
  sunctx->last_err = code;
  SUNErrHandler eh = sunctx->err_handler;
  while (eh != NULL)
  {
    eh->call(line, func, file, NULL, code, eh->data, sunctx);
    eh = eh->previous;
  }
  return 0;
}

static inline int SUNHandleErrWithMsg(int line, const char* func,
                                      const char* file, const char* msg,
                                      SUNErrCode code, SUNContext sunctx)
{
  sunctx->last_err = code;
  SUNErrHandler eh = sunctx->err_handler;
  while (eh != NULL)
  {
    eh->call(line, func, file, msg, code, eh->data, sunctx);
    eh = eh->previous;
  }
  return 0;
}

static inline int SUNHandleErrWithFmtMsg(int line, const char* func,
                                         const char* file, const char* msgfmt,
                                         SUNErrCode code, SUNContext sunctx, ...)
{
  size_t msglen;
  char* msg;
  va_list values;
  va_start(values, sunctx);
  msglen = vsnprintf(NULL, (size_t)0, msgfmt, values); /* determine size of
                                                          buffer needed */
  msg = (char*)malloc(msglen);
  vsnprintf(msg, msglen, msgfmt, values);
  SUNHandleErrWithMsg(line, func, file, msg, code, sunctx);
  va_end(values);
  free(msg);
  return 0;
}

#define SUNError(code, sunctx) \
  SUNHandleErr(__LINE__, __func__, __FILE__, code, sunctx)

#define SUNErrorWithMsg(code, msg, sunctx) \
  SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, code, sunctx)

#define SUNErrorWithFmtMsg(code, msg, sunctx, ...)                        \
  SUNHandleErrWithFmtMsg(__LINE__, __func__, __FILE__, msg, code, sunctx, \
                         __VA_ARGS__)

#ifdef __cplusplus /* wrapper to enable C++ usage */
} /* extern "C" */
#endif

/* ----------------------------------------------------------------------------
 * Error checking macros
 * ---------------------------------------------------------------------------*/

/* SUNCheckCall performs the SUNDIALS function call, and checks the returned
   error code. If an error occured, then it will log the error, set the last_err
   value, and call the error handler. */
#define SUNCheckCall(call, sunctx)                                       \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_ = call;                            \
    if (sun_chk_call_err_code_ < 0)                                      \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   sunctx);                                              \
    }                                                                    \
  }                                                                      \
  while (0);

/* Same as SUNCheckCall, but returns with the error code. */
#define SUNCheckCallReturn(call, sunctx)                                 \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_ = call;                            \
    if (sun_chk_call_err_code_ < 0)                                      \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   sunctx);                                              \
      return sun_chk_call_err_code_;                                     \
    }                                                                    \
  }                                                                      \
  while (0);

/* Same as SUNCheckCall, but returns with NULL. */
#define SUNCheckCallReturnNull(call, sunctx)                             \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_ = call;                            \
    if (sun_chk_call_err_code_ < 0)                                      \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   sunctx);                                              \
      return NULL;                                                       \
    }                                                                    \
  }                                                                      \
  while (0);

/* SUNCheckLastErr checks the last_err value in the SUNContext.
   If an error occured, then it will log the error, set the last_err
   value, and calls the error handler. */
#define SUNCheckLastErr(sunctx)                                          \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_;                                   \
    SUNContext_GetLastError(sunctx, &sun_chk_call_err_code_);            \
    if (sun_chk_call_err_code_ < 0)                                      \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   sunctx);                                              \
    }                                                                    \
  }                                                                      \
  while (0);

/* Same as SUNCheckLastErr, but returns with the error code. */
#define SUNCheckLastErrReturn(sunctx)                                    \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_;                                   \
    SUNContext_GetLastError(sunctx, &sun_chk_call_err_code_);            \
    if (sun_chk_call_err_code_ < 0)                                      \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   sunctx);                                              \
    }                                                                    \
  }                                                                      \
  while (0);

/* Same as SUNCheckLastErr, but returns void. */
#define SUNCheckLastErrReturnVoid(sunctx)                                \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_;                                   \
    SUNContext_GetLastError(sunctx, &sun_chk_call_err_code_);            \
    if (sun_chk_call_err_code_ < 0)                                      \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   sunctx);                                              \
      return;                                                            \
    }                                                                    \
  }                                                                      \
  while (0);

/* Same as SUNCheckLastErr, but returns with the error code. */
#define SUNCheckLastErrReturnNull(sunctx)                                \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_;                                   \
    SUNContext_GetLastError(sunctx, &sun_chk_call_err_code_);            \
    if (sun_chk_call_err_code_ < 0)                                      \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   sunctx);                                              \
      return NULL;                                                       \
    }                                                                    \
  }                                                                      \
  while (0);

/* TODO(CJB): should asserts really need a code even though expr is printed? */
/* SUNAssert checks if an expression is true.
   If the expression is false, it calls the SUNAssertErrHandler. */
#if !defined(NDEBUG)
#define SUNAssert(expr, code, sunctx)                                  \
  do {                                                                 \
    if (!(expr))                                                       \
    {                                                                  \
      SUNAssertErrHandlerFn(__LINE__, __func__, __FILE__, #expr, code, \
                            sunctx->err_handler->data, sunctx);        \
    }                                                                  \
  }                                                                    \
  while (0);
#endif

/* TODO(CJB): should asserts really need a code even though expr is printed? */
/* SUNMPIAssert checks if an expression is true.
   If the expression is false, it calls the SUNMPIAbortErrHandler. */
#if !defined(NDEBUG)
#define SUNMPIAssert(expr, code, sunctx)                                  \
  do {                                                                    \
    if (!(expr))                                                          \
    {                                                                     \
      SUNMPIAssertErrHandlerFn(__LINE__, __func__, __FILE__, #expr, code, \
                               sunctx->err_handler->data, sunctx);        \
    }                                                                     \
  }                                                                       \
  while (0);
#else
#define SUNMPIAssert(expr, code, sunctx)
#endif

/* TODO(CJB): the error code scheme used isnt yet properly handled by
  SUNHandleErr Also, should these be asserts? */
#define SUNCheckMPICall(call, sunctx)                                      \
  do {                                                                     \
    int sun_chk_mpi_call_err_code_ = call;                                 \
    if (sun_chk_mpi_call_err_code_ != MPI_SUCCESS)                         \
    {                                                                      \
      SUNHandleErr(__LINE__, __func__, __FILE__,                           \
                   SUN_ERR_MPI_FAIL + sun_chk_mpi_call_err_code_, sunctx); \
    }                                                                      \
  }                                                                        \
  while (0);

/* Ensure SUNContext is not NULL by trying to access the last_err member.
   This purposely causes a segmentation fault if the context is NULL so
   that the error occurs on the correct line in the correct file.  */
#if !defined(NDEBUG)
#define SUNAssertContext(sunctx)                                               \
  do {                                                                         \
    sunctx->last_err = sunctx->last_err; /* trigger segfault if sunctx is NULL \
                                          */                                   \
  }                                                                            \
  while (0);
#else
#define SUNAssertContext(sunctx)
#endif

#endif /* _SUNDIALS_ERRORS_H */
