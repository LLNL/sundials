/* -----------------------------------------------------------------
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
 * !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * This is a 'private' header file and should not be used in user
 * code. It is subject to change without warning.
 * !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * -----------------------------------------------------------------
 * Contains error checking macros and prototypes for MPI calls.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_MPI_ERRORS_IMPL_H
#define _SUNDIALS_MPI_ERRORS_IMPL_H

#include <mpi.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_mpi_errors.h>

/*
   SUNCheckMPICallMsg performs the MPI function call, and checks the
   returned error code. If an error occured, then it will log the error, set the
   last_err value, call the error handler, **and then returns SUN_ERR_MPI_FAIL**.

   :param call: the MPI function call
   :param msg: an error message
 */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckMPICallMsg(call, msg)                                          \
  do {                                                                         \
    int sun_chk_mpi_call_err_code_ = call;                                     \
    if (sun_chk_mpi_call_err_code_ != MPI_SUCCESS)                             \
    {                                                                          \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, SUN_ERR_MPI_FAIL, \
                          SUNCTX_);                                            \
      return SUN_ERR_MPI_FAIL;                                                 \
    }                                                                          \
  }                                                                            \
  while (0)
#else
#define SUNCheckMPICallMsg(call, msg) call
#endif

/*
   SUNCheckMPICallNullMsg performs the MPI function call, and checks the
   returned error code. If an error occured, then it will log the error, set the
   last_err value, call the error handler, **and then returns NULL**.

   :param call: the MPI function call
   :param msg: an error message
 */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckMPICallNullMsg(call, msg)                                      \
  do {                                                                         \
    int sun_chk_mpi_call_err_code_ = call;                                     \
    if (sun_chk_mpi_call_err_code_ != MPI_SUCCESS)                             \
    {                                                                          \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, SUN_ERR_MPI_FAIL, \
                          SUNCTX_);                                            \
      return NULL;                                                             \
    }                                                                          \
  }                                                                            \
  while (0)
#else
#define SUNCheckMPICallNullMsg(call, msg) call
#endif

/*
   SUNCheckMPICallVoidMsg performs the MPI function call, and checks the
   returned error code. If an error occured, then it will log the error, set the
   last_err value, call the error handler, **and then returns void**.

   :param call: the MPI function call
   :param msg: an error message
 */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckMPICallVoidMsg(call, msg)                                      \
  do {                                                                         \
    int sun_chk_mpi_call_err_code_ = call;                                     \
    if (sun_chk_mpi_call_err_code_ != MPI_SUCCESS)                             \
    {                                                                          \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, SUN_ERR_MPI_FAIL, \
                          SUNCTX_);                                            \
      return;                                                                  \
    }                                                                          \
  }                                                                            \
  while (0)
#else
#define SUNCheckMPICallVoidMsg(call, msg) call
#endif

/*
   SUNCheckMPICallNoRetMsg performs the MPI function call, and checks the
   returned error code. If an error occured, then it will log the error, set the
   last_err value, call the error handler. **It does not return**.

   :param call: the MPI function call
   :param msg: an error message
 */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckMPICallNoRetMsg(call)                                          \
  do {                                                                         \
    int sun_chk_mpi_call_err_code_ = call;                                     \
    if (sun_chk_mpi_call_err_code_ != MPI_SUCCESS)                             \
    {                                                                          \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, SUN_ERR_MPI_FAIL, \
                          SUNCTX_);                                            \
    }                                                                          \
  }                                                                            \
  while (0)
#else
#define SUNCheckMPICallNoRetMsg(call) call
#endif

/* These versions of SUNCheckMPICall do not take a custom message so a
   default message associated with the error code will be used. */
#define SUNCheckMPICall(call)      SUNCheckMPICallMsg(call, NULL)
#define SUNCheckMPICallNoRet(call) SUNCheckMPICallNoRetMsg(call, NULL)
#define SUNCheckMPICallNull(call)  SUNCheckMPICallNullMsg(call, NULL)
#define SUNCheckMPICallVoid(call)  SUNCheckMPICallVoidMsg(call, NULL)

#endif /* _SUNDIALS_MPI_ERRORS_IMPL_H */
