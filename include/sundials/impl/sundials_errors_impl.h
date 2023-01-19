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

/* The SUNAssignSUNCTX macro is used to declare the SUNContext
   object to be used a function. */
#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNAssignSUNCTX(sunctx) \
  SUNContext sunctx_ = sunctx;  \
  (void)sunctx_
#else
#define SUNAssignSUNCTX(sunctx) \
  SUNContext sunctx_ = sunctx;  \
  (void)sunctx_
#endif

/* ----------------------------------------------------------------------------
 * Error checking macros
 * ---------------------------------------------------------------------------*/

/* SUNCheck is like SUNAssert except it is compiled even in release mode
   unless SUNDIALS_DISABLE_ERROR_CHECKS is defined  */
#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheck(expr, code)                                                   \
  do {                                                                         \
    if (SUNHintFalse(!(expr)))                                                 \
    {                                                                          \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, #expr, code, sunctx_); \
    }                                                                          \
  }                                                                            \
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
                   sunctx_);                                             \
      (void)SUNGetLastErr(sunctx_);                                      \
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
                   sunctx_);                                             \
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
                   sunctx_);                                             \
      return NULL;                                                       \
    }                                                                    \
  }                                                                      \
  while (0)
#else
#define SUNCheckCallNull(call) (void)call;
#endif

/* Same as SUNCheckCall, but returns void. */
#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheckCallVoid(call)                                           \
  do {                                                                   \
    SUNErrCode sun_chk_call_err_code_ = call;                            \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))                        \
    {                                                                    \
      SUNHandleErr(__LINE__, __func__, __FILE__, sun_chk_call_err_code_, \
                   sunctx_);                                             \
      return;                                                            \
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
  SUNCheckCallNoRet(SUNGetLastErr(sunctx_))

/* Same as SUNCheckCallLastErrNoRet, but returns with the error code. */
#define SUNCheckCallLastErr(call) \
  call;                           \
  SUNCheckCall(SUNGetLastErr(sunctx_))

/* Same as SUNCheckCallLastErrNoRet, but returns void. */
#define SUNCheckCallLastErrVoid(call) \
  call;                               \
  SUNCheckCallVoid(SUNGetLastErr(sunctx_))

/* Same as SUNCheckCallLastErrNoRet, but returns NULL. */
#define SUNCheckCallLastErrNull(call) \
  call;                               \
  SUNCheckCallNull(SUNGetLastErr(sunctx_))
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
                            sunctx_->err_handler->data, sunctx_);      \
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

#endif /* _SUNDIALS_ERRORS_IMPL_H */
