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

#ifndef _SUNDIALS_MPI_ERRORS_H
#define _SUNDIALS_MPI_ERRORS_H

#include <sundials/sundials_errors.h>

int SUNMPIAbortErrHandlerFn(int line, const char* func, const char* file,
                            const char* msg, SUNErrCode err_code,
                            void* err_user_data, SUNContext sunctx);

int SUNMPIAssertErrHandlerFn(int line, const char* func, const char* file,
                             const char* msg, SUNErrCode err_code,
                             void* err_user_data, SUNContext sunctx);


/* TODO(CJB): the error code scheme used isnt yet properly handled by
    SUNHandleErr Also, should these be asserts? */
#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheckMPICallNoRet(call, sunctx)                                 \
  do {                                                                     \
    int sun_chk_mpi_call_err_code_ = call;                                 \
    if (sun_chk_mpi_call_err_code_ != MPI_SUCCESS)                         \
    {                                                                      \
      SUNHandleErr(__LINE__, __func__, __FILE__,                           \
                   SUN_ERR_MPI_FAIL + sun_chk_mpi_call_err_code_, sunctx); \
    }                                                                      \
  }                                                                        \
  while (0)
#else
#define SUNCheckMPICallNoRet(call, sunctx) \
  call;                                    \
  (void)sunctx
#endif

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
  while (0)
#else
#define SUNMPIAssert(expr, code, sunctx)
#endif

#endif /* _SUNDIALS_MPI_ERRORS_H */
