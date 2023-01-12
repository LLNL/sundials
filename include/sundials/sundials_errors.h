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
#include <stdarg.h>
#include <stdio.h>
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
  ENTRY(SUN_ERR_ARG_INCOMPATIBLE, "argument provided is not comptaible")       \
  ENTRY(SUN_ERR_ARG_OUTOFRANGE, "argument is out of the valid range")          \
  ENTRY(SUN_ERR_ARG_WRONGTYPE, "argument provided is not the right type")      \
  ENTRY(SUN_ERR_ARG_DIMSMISMATCH, "argument dimensions do not agree")          \
  ENTRY(SUN_ERR_HYPRE_GENERIC, "a generic error ocurred in hypre")             \
  ENTRY(SUN_ERR_HYPRE_MEMORY, "a memory error ocurred in hypre")               \
  ENTRY(SUN_ERR_HYPRE_ARG, "an argument error ocurred in hypre")               \
  ENTRY(SUN_ERR_HYPRE_CONV, "a convergence error ocurred in hypre")            \
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
  ENTRY(SUN_ERR_PROFILER_MAPFULL,                                              \
        "too many SUNProfiler entries, try setting SUNPROFILER_MAX_ENTRIES "   \
        "in environment to a bigger number")                                   \
  ENTRY(SUN_ERR_PROFILER_MAPGET, "unknown error getting SUNProfiler timer")    \
  ENTRY(SUN_ERR_PROFILER_MAPINSERT,                                            \
        "unknown error inserting SUNProfiler timer")                           \
  ENTRY(SUN_ERR_PROFILER_MAPKEYNOTFOUND, "timer was not found in SUNProfiler") \
  ENTRY(SUN_ERR_PROFILER_MAPSORT, "error sorting SUNProfiler map")             \
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

SUNDIALS_EXPORT
SUNErrHandler SUNErrHandler_Create(SUNErrHandlerFn eh_fn, void* eh_data);

SUNDIALS_EXPORT
void SUNErrHandler_Destroy(SUNErrHandler eh);

SUNDIALS_EXPORT
int SUNLogErrHandlerFn(int line, const char* func, const char* file,
                       const char* msg, SUNErrCode err_code,
                       void* err_user_data, SUNContext sunctx);

SUNDIALS_EXPORT
int SUNAbortErrHandlerFn(int line, const char* func, const char* file,
                         const char* msg, SUNErrCode err_code,
                         void* err_user_data, SUNContext sunctx);

SUNDIALS_EXPORT
int SUNAssertErrHandlerFn(int line, const char* func, const char* file,
                          const char* msg, SUNErrCode err_code,
                          void* err_user_data, SUNContext sunctx);

/* ----------------------------------------------------------------------------
 * Error functions
 * ---------------------------------------------------------------------------*/

/* Turn error code into error message */
SUNDIALS_EXPORT
const char* SUNGetErrMsg(SUNErrCode code, SUNContext sunctx);

/* Alternative function to SUNContext_GetLastError that is more concise. */
static inline SUNErrCode SUNGetLastErr(SUNContext sunctx)
{
  SUNErrCode code = SUN_SUCCESS;
  (void)SUNContext_GetLastError(sunctx, &code);
  return code;
}

/* Alternative function to SUNContext_SetLastError that is more concise. */
static inline SUNErrCode SUNSetLastErr(SUNErrCode code, SUNContext sunctx)
{
  sunctx->last_err = code;
  return SUN_SUCCESS;
}

/* Alternative function to SUNContext_PeekLastError that is more concise. */
static inline SUNErrCode SUNPeekLastErr(SUNContext sunctx)
{
  SUNErrCode code = SUN_SUCCESS;
  (void)SUNContext_PeekLastError(sunctx, &code);
  return code;
}

static inline void SUNHandleErr(int line, const char* func, const char* file,
                                SUNErrCode code, SUNContext sunctx)
{
  sunctx->last_err = code;
  SUNErrHandler eh = sunctx->err_handler;
  while (eh != NULL)
  {
    eh->call(line, func, file, NULL, code, eh->data, sunctx);
    eh = eh->previous;
  }
}

static inline void SUNHandleErrWithMsg(int line, const char* func,
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
}

static inline void SUNHandleErrWithFmtMsg(int line, const char* func,
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
}

#ifdef __cplusplus /* wrapper to enable C++ usage */
} /* extern "C" */
#endif

/* ----------------------------------------------------------------------------
 * Error checking macros
 * ---------------------------------------------------------------------------*/

/* SUNCheck is like SUNAssert except it is compiled even in release mode
   unless SUNDIALS_DISABLE_ERROR_CHECKS is defined  */
#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheck(expr, code)                                                  \
  do {                                                                        \
    if (SUNHintFalse(!(expr)))                                                \
    {                                                                         \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, #expr, code, SUNCTX); \
    }                                                                         \
  }                                                                           \
  while (0)
#else
#define SUNCheck(expr, code)
#endif

/* SUNCheckCallNoRet performs the SUNDIALS function call, and checks the
   returned error code. If an error occured, then it will log the error, set the
   last_err value, and call the error handler. */
#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheckCallNoRet(call)                                          \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_ = call;                            \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))                        \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   SUNCTX);                                              \
      (void)SUNGetLastErr(SUNCTX);                                       \
    }                                                                    \
  }                                                                      \
  while (0)
#else
#define SUNCheckCallNoRet(call) (void)call
#endif

/* Same as SUNCheckCallNoRet, but returns with the error code if an error
 * occured. */
#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheckCall(call)                                               \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_ = call;                            \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))                        \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   SUNCTX);                                              \
      return sun_chk_call_err_code_;                                     \
    }                                                                    \
  }                                                                      \
  while (0)
#else
#define SUNCheckCall(call) (void)call
#endif

/* Same as SUNCheckCall, but returns with NULL. */
#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheckCallNull(call)                                           \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_ = call;                            \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))                        \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   SUNCTX);                                              \
      return NULL;                                                       \
    }                                                                    \
  }                                                                      \
  while (0)
#else
#define SUNCheckCallNull(call) (void)call;
#endif

/* SUNCheckLastErr checks the last_err value in the SUNContext.
   If an error occured, then it will log the error, set the last_err
   value, and calls the error handler. */

#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheckCallLastErrNoRet(call) \
  call;                                \
  SUNCheckCallNoRet(SUNGetLastErr(SUNCTX))

/* Same as SUNCheckCallLastErrNoRet, but returns with the error code. */
#define SUNCheckCallLastErr(call) \
  call;                           \
  SUNCheckCall(SUNGetLastErr(SUNCTX))

/* Same as SUNCheckCallLastErrNoRet, but returns void. */
#define SUNCheckCallLastErrVoid(call) \
  call;                               \
  SUNCheckCallVoid(SUNGetLastErr(SUNCTX))

/* Same as SUNCheckCallLastErrNoRet, but returns with the error code. */
#define SUNCheckCallLastErrNull(call) \
  call;                               \
  SUNCheckCallNull(SUNGetLastErr(SUNCTX))
#else
#define SUNCheckCallLastErrNoRet(call) call
#define SUNCheckCallLastErr(call)      call
#define SUNCheckCallLastErrVoid(call)  call
#define SUNCheckCallLastErrNull(call)  call
#endif

/* SUNAssert checks if an expression is true.
   If the expression is false, it calls the SUNAssertErrHandler. */
#if !defined(NDEBUG)
#define SUNAssert(expr, code)                                          \
  do {                                                                 \
    if (!(expr))                                                       \
    {                                                                  \
      SUNAssertErrHandlerFn(__LINE__, __func__, __FILE__, #expr, code, \
                            SUNCTX->err_handler->data, SUNCTX);        \
    }                                                                  \
  }                                                                    \
  while (0)
#else
#define SUNAssert(expr, code)
#endif

/* TODO(CJB): probably should create a global handler that does not need
   a SUNContext for cases where SUNContext is not available.
   It won't be thread safe, but that should be OK since this is a catastropic
   error. */

#endif /* _SUNDIALS_ERRORS_H */
