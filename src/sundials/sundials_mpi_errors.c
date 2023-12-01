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

#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/priv/sundials_mpi_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_mpi_types.h>
#include <unistd.h>

static inline char* combineFileAndLine(int line, const char* file)
{
  size_t total_str_len = strlen(file) + 6;
  char* file_and_line  = malloc(total_str_len * sizeof(char));
  snprintf(file_and_line, total_str_len, "%s:%d", file, line);
  return file_and_line;
}

void SUNMPIAbortErrHandlerFn(int line, const char* func, const char* file,
                             const char* msg, SUNErrCode err_code,
                             void* err_user_data, SUNContext sunctx)
{
  char* file_and_line = combineFileAndLine(line, file);
  SUNLogger_QueueMsg(sunctx->logger, SUN_LOGLEVEL_ERROR, file_and_line, func,
                     "SUNMPIAbortErrHandler: Calling MPI_Abort now, use a "
                     "different "
                     "error handler to avoid program termination.\n");
  free(file_and_line);
  sleep(1);
  MPI_Abort(sunctx->comm, err_code);
}

void SUNMPIAssertErrHandlerFn(int line, const char* func, const char* file,
                              const char* stmt, SUNErrCode err_code,
                              void* err_user_data, SUNContext sunctx)
{
  char* file_and_line = combineFileAndLine(line, file);
  SUNLogger_QueueMsg(sunctx->logger, SUN_LOGLEVEL_ERROR, file_and_line,
                     func, "SUNMPIAssertErrHandler: assert(%s) failed... terminating.\n",
                     stmt);
  free(file_and_line);
  sleep(1);
  MPI_Abort(sunctx->comm, err_code);
}
