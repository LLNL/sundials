/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * SUNDIALS logging class implementation.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_LOGGER_IMPL_H
#define _SUNDIALS_LOGGER_IMPL_H

#include <stdarg.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>

#include "sundials_hashmap.h"

#define SUNDIALS_LOGGING_ERROR 1
#define SUNDIALS_LOGGING_WARNING 2
#define SUNDIALS_LOGGING_INFO 3
#define SUNDIALS_LOGGING_DEBUG 4
#if SUNDIALS_LOGGING_LEVEL > SUNDIALS_LOGGING_DEBUG
#define SUNDIALS_LOGGING_EXTRA_DEBUG
#endif

struct SUNLogger_ {
  /* MPI information */
  void* commptr;
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
  int (*queuemsg)(SUNLogger logger, SUNLogLevel lvl, const char* scope,
                  const char* label, const char* msg_txt, va_list args);
  int (*flush)(SUNLogger logger, SUNLogLevel lvl);
  int (*destroy)(SUNLogger* logger);
};

#endif /* _SUNDIALS_LOGGER_IMPL_H */
