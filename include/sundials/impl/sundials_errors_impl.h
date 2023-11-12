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

#define SUNCTX sunctx_local_scope_

/* The SUNFunctionBegin macro is used to declare the SUNContext
   object to be used a function. */
#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNFunctionBegin(sunctx) \
  SUNContext SUNCTX = sunctx;  \
  (void)SUNCTX
#else
#define SUNFunctionBegin(sunctx) \
  SUNContext SUNCTX = sunctx;  \
  (void)SUNCTX
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
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, #expr, code, SUNCTX); \
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
                   SUNCTX);                                             \
      (void)SUNGetLastErr(SUNCTX);                                      \
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
                   SUNCTX);                                             \
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
                   SUNCTX);                                             \
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
                   SUNCTX);                                             \
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
#define SUNCheckLastErrNoRet() \
  SUNCheckCallNoRet(SUNGetLastErr(SUNCTX))

/* Same as SUNCheckLastErrNoRet, but returns with the error code. */
#define SUNCheckLastErr() \
  SUNCheckCall(SUNGetLastErr(SUNCTX))

/* Same as SUNCheckLastErrNoRet, but returns void. */
#define SUNCheckLastErrVoid() \
  SUNCheckCallVoid(SUNGetLastErr(SUNCTX))

/* Same as SUNCheckLastErrNoRet, but returns NULL. */
#define SUNCheckLastErrNull() \
  SUNCheckCallNull(SUNGetLastErr(SUNCTX))
#else
#define SUNCheckLastErrNoRet() 
#define SUNCheckLastErr()      
#define SUNCheckLastErrVoid()  
#define SUNCheckLastErrNull()  
#endif

/* SUNAssert checks if an expression is true.
   If the expression is false, it calls the SUNAssertErrHandler. */
#if !defined(NDEBUG)
#define SUNAssert(expr, code)                                          \
  do {                                                                 \
    if (!(expr))                                                       \
    {                                                                  \
      SUNAssertErrHandlerFn(__LINE__, __func__, __FILE__, #expr, code, \
                            SUNCTX->err_handler->data, SUNCTX);      \
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
