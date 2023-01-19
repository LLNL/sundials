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

#ifndef _SUNDIALS_MPI_ERRORS_IMPL_H
#define _SUNDIALS_MPI_ERRORS_IMPL_H

#include <mpi.h>
#include <sundials/sundials_mpi_errors.h>

#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheckMPICall(call)                                                \
  do {                                                                       \
    int sun_chk_mpi_call_err_code_ = call;                                   \
    if (sun_chk_mpi_call_err_code_ != MPI_SUCCESS)                           \
    {                                                                        \
      SUNHandleErr(__LINE__, __func__, __FILE__, SUN_ERR_MPI_FAIL, sunctx_); \
      return SUN_ERR_MPI_FAIL;                                               \
    }                                                                        \
  }                                                                          \
  while (0)
#else
#define SUNCheckMPICall(call) call;
#endif

#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheckMPICallNull(call)                                            \
  do {                                                                       \
    int sun_chk_mpi_call_err_code_ = call;                                   \
    if (sun_chk_mpi_call_err_code_ != MPI_SUCCESS)                           \
    {                                                                        \
      SUNHandleErr(__LINE__, __func__, __FILE__, SUN_ERR_MPI_FAIL, sunctx_); \
      return NULL;                                                           \
    }                                                                        \
  }                                                                          \
  while (0)
#else
#define SUNCheckMPICallNull(call) call;
#endif

#if !defined(SUNDIALS_DISABLE_ERROR_CHECKS)
#define SUNCheckMPICallNoRet(call)                                           \
  do {                                                                       \
    int sun_chk_mpi_call_err_code_ = call;                                   \
    if (sun_chk_mpi_call_err_code_ != MPI_SUCCESS)                           \
    {                                                                        \
      SUNHandleErr(__LINE__, __func__, __FILE__, SUN_ERR_MPI_FAIL, sunctx_); \
    }                                                                        \
  }                                                                          \
  while (0)
#else
#define SUNCheckMPICallNoRet(call) call;
#endif

/* SUNMPIAssert checks if an expression is true.
   If the expression is false, it calls the SUNMPIAssertErrHandler. */
#if !defined(NDEBUG)
#define SUNMPIAssert(expr, code)                                          \
  do {                                                                    \
    if (!(expr))                                                          \
    {                                                                     \
      SUNMPIAssertErrHandlerFn(__LINE__, __func__, __FILE__, #expr, code, \
                               sunctx_->err_handler->data, sunctx_);      \
    }                                                                     \
  }                                                                       \
  while (0)
#else
#define SUNMPIAssert(expr, code)
#endif

#endif /* _SUNDIALS_MPI_ERRORS_IMPL_H */
