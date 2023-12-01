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
 * !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * This is a 'private' header file and should not be used in user
 * code. It is subject to change without warning.
 * !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * -----------------------------------------------------------------
 * Contains all error checking macros and private error handling API.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_ERRORS_IMPL_H
#define _SUNDIALS_ERRORS_IMPL_H

#include <sundials/sundials_errors.h>

/* ----------------------------------------------------------------------------
 * SUNErrHandler_ definition
 * ---------------------------------------------------------------------------*/

struct SUNErrHandler_
{
  SUNErrHandler previous; /* next error handler to call (singly linked-list) */
  SUNErrHandlerFn call;
  void* data;
};

SUNErrCode SUNErrHandler_Create(SUNErrHandlerFn eh_fn, void* eh_data,
                                SUNErrHandler* eh_out);

void SUNErrHandler_Destroy(SUNErrHandler eh);

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
  msglen = (size_t)vsnprintf(NULL, (size_t)0, msgfmt, values); /* determine size
                                                                  of buffer
                                                                  needed */
  msg = (char*)malloc(msglen + 1);
  vsnprintf(msg, msglen + 1, msgfmt, values);
  SUNHandleErrWithMsg(line, func, file, msg, code, sunctx);
  va_end(values);
  free(msg);
}

/* ----------------------------------------------------------------------------
 * Error checking macros
 * ---------------------------------------------------------------------------*/

#define SUNCTX_ sunctx_local_scope_

/* The SUNFunctionBegin macro is used to declare the SUNContext
   object to be used a function. */
#define SUNFunctionBegin(sunctx) \
  SUNContext SUNCTX_ = sunctx;   \
  (void)SUNCTX_

/* SUNCheck evaluates the given expression and calls the error handler if it is
 * false. */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheck(expr, code)                                                   \
  do {                                                                         \
    if (SUNHintFalse(!(expr)))                                                 \
    {                                                                          \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, #expr, code, SUNCTX_); \
    }                                                                          \
  }                                                                            \
  while (0)
#else
#define SUNCheck(expr, code)
#endif

/* SUNCheckCallNoRet performs the SUNDIALS function call, and checks the
   returned error code. If an error occured, then it will log the error, set the
   last_err value, and call the error handler. */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckCallNoRetMsg(call, msg)                      \
  do {                                                       \
    SUNErrCode sun_chk_call_err_code_ = call;                \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))            \
    {                                                        \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, \
                          sun_chk_call_err_code_, SUNCTX_);  \
      (void)SUNContext_GetLastError(SUNCTX_);                \
    }                                                        \
  }                                                          \
  while (0)
#else
#define SUNCheckCallNoRetMsg(call, msg) (void)call
#endif

/* Same as SUNCheckCallNoRet, but returns with the error code if an error
 * occured. */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckCallMsg(call, msg)                           \
  do {                                                       \
    SUNErrCode sun_chk_call_err_code_ = call;                \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))            \
    {                                                        \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, \
                          sun_chk_call_err_code_, SUNCTX_);  \
      return sun_chk_call_err_code_;                         \
    }                                                        \
  }                                                          \
  while (0)
#else
#define SUNCheckCallMsg(call, msg) (void)call
#endif

/* Same as SUNCheckCall, but returns with NULL. */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckCallNullMsg(call, msg)                       \
  do {                                                       \
    SUNErrCode sun_chk_call_err_code_ = call;                \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))            \
    {                                                        \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, \
                          sun_chk_call_err_code_, SUNCTX_);  \
      return NULL;                                           \
    }                                                        \
  }                                                          \
  while (0)
#else
#define SUNCheckCallNullMsg(call, msg) (void)call
#endif

/* Same as SUNCheckCall, but returns void. */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckCallVoidMsg(call, msg)                       \
  do {                                                       \
    SUNErrCode sun_chk_call_err_code_ = call;                \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))            \
    {                                                        \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, \
                          sun_chk_call_err_code_, SUNCTX_);  \
      return;                                                \
    }                                                        \
  }                                                          \
  while (0)
#else
#define SUNCheckCallVoidMsg(call, msg) (void)call
#endif

/* These versions of SUNCheckCall do not take a custom message so the
   default message associated with the error code will be used. */
#define SUNCheckCallNoRet(call) SUNCheckCallNoRetMsg(call, NULL)
#define SUNCheckCall(call)      SUNCheckCallMsg(call, NULL)
#define SUNCheckCallNull(call)  SUNCheckCallNullMsg(call, NULL)
#define SUNCheckCallVoid(call)  SUNCheckCallVoidMsg(call, NULL)

/* SUNCheckLastErrMoRetMsg checks the last_err value in the SUNContext.
   If an error occured, then it will log the error, set the last_err
   value, and calls the error handler. */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckLastErrNoRetMsg(msg) \
  SUNCheckCallNoRetMsg(SUNContext_GetLastError(SUNCTX_), msg)
#define SUNCheckLastErrMsg(msg) \
  SUNCheckCallMsg(SUNContext_GetLastError(SUNCTX_), msg)
#define SUNCheckLastErrNullMsg(msg) \
  SUNCheckCallNullMsg(SUNContext_GetLastError(SUNCTX_), msg)
#define SUNCheckLastErrVoidMsg(msg) \
  SUNCheckCallVoidMsg(SUNContext_GetLastError(SUNCTX_), msg)
#else
#define SUNCheckLastErrNoRetMsg(msg)
#define SUNCheckLastErrMsg(msg)
#define SUNCheckLastErrVoidMsg(msg)
#define SUNCheckLastErrNullMsg(msg)
#endif

/* These versions of SUNCheckLastErr do not take a custom message so the
   default message associated with the error code will be used. */
#define SUNCheckLastErr()      SUNCheckLastErrMsg(NULL)
#define SUNCheckLastErrNoRet() SUNCheckLastErrNoRetMsg(NULL)
#define SUNCheckLastErrVoid()  SUNCheckLastErrVoidMsg(NULL)
#define SUNCheckLastErrNull()  SUNCheckLastErrNullMsg(NULL)

/* SUNAssert checks if an expression is true.
   It expands to SUNCheck when error checks are enabled.
   If error checks are disabled, then we try to expand it to an assumption,
   if the compiler supoprts, so that the compiler can make optimizations based
   on the assumption.
*/

#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNAssert(expr, code) SUNCheck(expr, code)
#else
#define SUNAssert(expr, code) SUNAssume(expr)
#endif

#endif /* _SUNDIALS_ERRORS_IMPL_H */
