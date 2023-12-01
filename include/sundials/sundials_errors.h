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
#include <sundials/sundials_types.h>

/* ----------------------------------------------------------------------------
 * Error code definitions
 * ---------------------------------------------------------------------------*/

/* SUN_ERR_CODE_LIST is an X macro that can be expanded in various ways */
#define SUN_ERR_CODE_LIST(ENTRY)                                               \
  ENTRY(SUN_ERR_ARG_CORRUPT, "argument provided is NULL or corrupted")         \
  ENTRY(SUN_ERR_ARG_INCOMPATIBLE, "argument provided is not compatible")       \
  ENTRY(SUN_ERR_ARG_OUTOFRANGE, "argument is out of the valid range")          \
  ENTRY(SUN_ERR_ARG_WRONGTYPE, "argument provided is not the right type")      \
  ENTRY(SUN_ERR_ARG_DIMSMISMATCH, "argument dimensions do not agree")          \
                                                                               \
  ENTRY(SUN_ERR_CORRUPT, "Object is NULL or corrupt")                          \
  ENTRY(SUN_ERR_FILE_OPEN, "Unable to open file")                              \
  ENTRY(SUN_ERR_MALLOC_FAIL, "malloc returned NULL")                           \
  ENTRY(SUN_ERR_DESTROY_FAIL, "a destroy function returned an error")          \
  ENTRY(SUN_ERR_NOT_IMPLEMENTED,                                               \
        "operation is not implemented: function pointer is NULL")              \
                                                                               \
  ENTRY(SUN_ERR_PROFILER_MAPFULL,                                              \
        "the number of profiler entries exceeded SUNPROFILER_MAX_ENTRIES")     \
  ENTRY(SUN_ERR_PROFILER_MAPGET, "unknown error getting SUNProfiler timer")    \
  ENTRY(SUN_ERR_PROFILER_MAPINSERT,                                            \
        "unknown error inserting SUNProfiler timer")                           \
  ENTRY(SUN_ERR_PROFILER_MAPKEYNOTFOUND, "timer was not found in SUNProfiler") \
  ENTRY(SUN_ERR_PROFILER_MAPSORT, "error sorting SUNProfiler map")             \
                                                                               \
  ENTRY(SUN_ERR_SUNCTX_CORRUPT, "SUNContext is NULL or corrupt")               \
                                                                               \
  ENTRY(SUN_ERR_MPI_FAIL,                                                      \
        "an MPI call returned something other than MPI_SUCCESS")               \
                                                                               \
  ENTRY(SUN_ERR_UNREACHABLE,                                                   \
        "Reached code that should be unreachable: open an issue at: "          \
        "https://github.com/LLNL/sundials")                                    \
  ENTRY(SUN_ERR_UNKNOWN, "Unknown error occured: open an issue at "            \
                         "https://github.com/LLNL/sundials")

/* Expand SUN_ERR_CODE_LIST to enum */
#define SUN_EXPAND_TO_ENUM(name, description) name,

/* SUNErrorCode range is [-10000, -1000] to avoid conflicts with package error
   codes, and old/deprecated codes for matrix and (non)linear solvers. */

/* clang-format off */
enum
{
  SUN_ERR_MINIMUM                                       = -10000,
  SUN_ERR_CODE_LIST(SUN_EXPAND_TO_ENUM) 
  SUN_ERR_MAXIMUM                                       = -1000,
  SUN_SUCCESS                                           = 0
};

/* clang-format on */

/* ----------------------------------------------------------------------------
 * Error handler definitions
 * ---------------------------------------------------------------------------*/

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT
void SUNLogErrHandlerFn(int line, const char* func, const char* file,
                        const char* msg, SUNErrCode err_code,
                        void* err_user_data, SUNContext sunctx);

SUNDIALS_EXPORT
void SUNAbortErrHandlerFn(int line, const char* func, const char* file,
                          const char* msg, SUNErrCode err_code,
                          void* err_user_data, SUNContext sunctx);

SUNDIALS_EXPORT
void SUNAssertErrHandlerFn(int line, const char* func, const char* file,
                           const char* stmt, SUNErrCode err_code,
                           void* err_user_data, SUNContext sunctx);

/* ----------------------------------------------------------------------------
 * Error functions
 * ---------------------------------------------------------------------------*/

/* Turn error code into error message */
SUNDIALS_EXPORT
const char* SUNGetErrMsg(SUNErrCode code, SUNContext sunctx);

#ifdef __cplusplus /* wrapper to enable C++ usage */
} /* extern "C" */
#endif
#endif /* _SUNDIALS_ERRORS_H */
