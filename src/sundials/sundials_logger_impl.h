/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * SUNDIALS logging class implementation.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_LOGGER_IMPL_H
#define _SUNDIALS_LOGGER_IMPL_H

#include <stdarg.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>

#include "sundials_hashmap.h"
#include "sundials_utils.h"

#define SUNDIALS_LOGGING_ERROR   1
#define SUNDIALS_LOGGING_WARNING 2
#define SUNDIALS_LOGGING_INFO    3
#define SUNDIALS_LOGGING_DEBUG   4
#if SUNDIALS_LOGGING_LEVEL > SUNDIALS_LOGGING_DEBUG
#define SUNDIALS_LOGGING_EXTRA_DEBUG
#endif

struct SUNLogger_
{
  /* MPI information */
  SUNComm comm;
  int output_rank;

  /* Ouput files */
  FILE* debug_fp;
  FILE* warning_fp;
  FILE* info_fp;
  FILE* error_fp;

  /* Hashmap used to store filename, FILE* pairs */
  SUNHashMap filenames;

  /* Slic-style format string */
  const char* format;

  /* Content for custom implementations */
  void* content;

  /* Overridable operations */
  SUNErrCode (*queuemsg)(SUNLogger logger, SUNLogLevel lvl, const char* scope,
                         const char* label, const char* msg_txt, va_list args);
  SUNErrCode (*flush)(SUNLogger logger, SUNLogLevel lvl);
  SUNErrCode (*destroy)(SUNLogger* logger);
};

static void sunCreateLogMessage(SUNLogLevel lvl, int rank, const char* scope,
                                const char* label, const char* txt,
                                va_list args, char** log_msg)
{
  char* prefix;
  char* formatted_txt;
  int msg_length;

  prefix        = NULL;
  formatted_txt = NULL;
  msg_length    = 0;
  *log_msg      = NULL;

  msg_length = sunvasnprintf(&formatted_txt, txt, args);
  if (msg_length < 0)
  {
    fprintf(stderr, "[FATAL LOGGER ERROR] %s\n",
            "SUNDIALS_MAX_SPRINTF_SIZE is too small");
  }

  if (lvl == SUN_LOGLEVEL_DEBUG) { prefix = (char*)"DEBUG"; }
  else if (lvl == SUN_LOGLEVEL_WARNING) { prefix = (char*)"WARNING"; }
  else if (lvl == SUN_LOGLEVEL_INFO) { prefix = (char*)"INFO"; }
  else if (lvl == SUN_LOGLEVEL_ERROR) { prefix = (char*)"ERROR"; }

  msg_length = sunsnprintf(NULL, 0, "[%s][rank %d][%s][%s] %s\n", prefix, rank,
                           scope, label, formatted_txt);
  *log_msg   = (char*)malloc(msg_length + 1);
  sunsnprintf(*log_msg, msg_length + 1, "[%s][rank %d][%s][%s] %s\n", prefix,
              rank, scope, label, formatted_txt);
  free(formatted_txt);
}

#endif /* _SUNDIALS_LOGGER_IMPL_H */
