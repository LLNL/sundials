/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
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
 * -----------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_logger.h>

#ifdef SUNDIALS_LOGGING_ENABLE_MPI
#include <mpi.h>
#endif

#include "sundials_logger_impl.h"
#include "sundials_utils.h"

/* max number of files that can be opened */
#define SUN_MAX_LOGFILE_HANDLES_ 8

/* shortcut */
#define SUNLOGGER_MPICOMM(logger) (*((MPI_Comm*)logger->commptr))

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

  if (lvl == SUN_LOGLEVEL_DEBUG)
  {
    prefix = "DEBUG";
  }
  else if (lvl == SUN_LOGLEVEL_WARNING)
  {
    prefix = "WARNING";
  }
  else if (lvl == SUN_LOGLEVEL_INFO)
  {
    prefix = "INFO";
  }
  else if (lvl == SUN_LOGLEVEL_ERROR)
  {
    prefix = "ERROR";
  }

  msg_length = sunsnprintf(NULL, 0, "[%s][rank::%d][%s][%s] %s\n", prefix,
                           rank, scope, label, formatted_txt);
  *log_msg   = (char*)malloc(msg_length + 1);
  sunsnprintf(*log_msg, msg_length + 1, "[%s][rank::%d][%s][%s] %s\n", prefix,
              rank, scope, label, formatted_txt);
  free(formatted_txt);
}

static FILE* sunOpenLogFile(const char* fname, const char* mode)
{
  FILE* fp = NULL;

  if (fname)
  {
    if (!strcmp(fname, "stdout"))
    {
      fp = stdout;
    }
    else if (!strcmp(fname, "stderr"))
    {
      fp = stderr;
    }
    else
    {
      fp = fopen(fname, mode);
    }
  }

  return fp;
}

static void sunCloseLogFile(void* fp)
{
  if (fp && fp != stdout && fp != stderr)
  {
    fclose((FILE*)fp);
  }
}

static sunbooleantype sunLoggerIsOutputRank(SUNLogger logger, int* rank_ref)
{
  sunbooleantype retval;

#ifdef SUNDIALS_LOGGING_ENABLE_MPI
  int rank = 0;

  if (logger->commptr)
  {
    MPI_Comm_rank(SUNLOGGER_MPICOMM(logger), &rank);

    if (logger->output_rank < 0)
    {
      if (rank_ref)
      {
        *rank_ref = rank;
      }
      retval = SUNTRUE; /* output all ranks */
    }
    else
    {
      if (rank_ref)
      {
        *rank_ref = rank;
      }
      retval = logger->output_rank == rank;
    }
  }
  else
  {
    retval = SUNTRUE; /* output all ranks */
  }
#else
  if (rank_ref)
  {
    *rank_ref = -1;
  }
  retval = SUNTRUE;
#endif

  return retval;
}

int SUNLogger_Create(void* comm, int output_rank, SUNLogger* logger_ptr)
{
  SUNLogger logger = NULL;

  *logger_ptr = logger = (SUNLogger)malloc(sizeof(struct SUNLogger_));
  if (logger == NULL)
  {
    return -1;
  }

  /* Attach the comm, duplicating it if MPI is used. */
#ifdef SUNDIALS_LOGGING_ENABLE_MPI
  logger->commptr = NULL;
  if (comm != NULL)
  {
    logger->commptr = malloc(sizeof(MPI_Comm));
    MPI_Comm_dup(*((MPI_Comm*) comm), (MPI_Comm*) logger->commptr);
  }
#else
  if (comm != NULL)
  {
    return -1;
  }
  logger->commptr = NULL;
#endif
  logger->output_rank = output_rank;
  logger->content     = NULL;

  /* use default routines */
  logger->queuemsg = NULL;
  logger->flush    = NULL;
  logger->destroy  = NULL;

  /* set the output file handles */
  logger->filenames  = NULL;
  logger->error_fp   = NULL;
  logger->warning_fp = NULL;
  logger->debug_fp   = NULL;
  logger->info_fp    = NULL;
  if (sunLoggerIsOutputRank(logger, NULL))
  {
    /* We store the FILE* in a hash map so that we can ensure
       that we do not open a file twice if the same file is used
       for multiple output levels */
    SUNHashMap_New(SUN_MAX_LOGFILE_HANDLES_, &logger->filenames);
  }

  return 0;
}

int SUNLogger_CreateFromEnv(void* comm, SUNLogger* logger)
{
  int retval = 0;

  const char* output_rank_env   = getenv("SUNLOGGER_OUTPUT_RANK");
  int output_rank               = (output_rank_env) ? atoi(output_rank_env) : 0;
  const char* error_fname_env   = getenv("SUNLOGGER_ERROR_FILENAME");
  const char* warning_fname_env = getenv("SUNLOGGER_WARNING_FILENAME");
  const char* info_fname_env    = getenv("SUNLOGGER_INFO_FILENAME");
  const char* debug_fname_env   = getenv("SUNLOGGER_DEBUG_FILENAME");

  retval += SUNLogger_Create(comm, output_rank, logger);
  retval += SUNLogger_SetErrorFilename(*logger, error_fname_env);
  retval += SUNLogger_SetWarningFilename(*logger, warning_fname_env);
  retval += SUNLogger_SetDebugFilename(*logger, debug_fname_env);
  retval += SUNLogger_SetInfoFilename(*logger, info_fname_env);

  return (retval < 0) ? -1 : 0;
}

int SUNLogger_SetErrorFilename(SUNLogger logger, const char* error_filename)
{
  if (logger == NULL)
  {
    return -1;
  }

  if (!sunLoggerIsOutputRank(logger, NULL))
  {
    return 0;
  }

  if (error_filename && strcmp(error_filename, ""))
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_ERROR
    FILE* fp = NULL;
    if (!SUNHashMap_GetValue(logger->filenames, error_filename, (void*)&fp))
    {
      logger->error_fp = fp;
    }
    else
    {
      logger->error_fp = sunOpenLogFile(error_filename, "w+");
      if (logger->error_fp)
      {
        SUNHashMap_Insert(logger->filenames, error_filename,
                          (void*)logger->error_fp);
      }
      else
      {
        return -1;
      }
    }
#else
    fprintf(stderr,
            "[LOGGER WARNING] "
            "SUNDIALS_LOGGING_LEVEL=%d (build time option) "
            "is set too low for ERROR, but a ERROR file was provided. "
            "Set the logging level to >= %d and recompile if ERROR output level "
            "is desired.\n", SUN_LOGLEVEL_ERROR, SUNDIALS_LOGGING_LEVEL);
#endif
  }

  return 0;
}

int SUNLogger_SetWarningFilename(SUNLogger logger, const char* warning_filename)
{
  if (logger == NULL)
  {
    return -1;
  }

  if (!sunLoggerIsOutputRank(logger, NULL))
  {
    return 0;
  }

  if (warning_filename && strcmp(warning_filename, ""))
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_WARNING
    FILE* fp = NULL;
    if (!SUNHashMap_GetValue(logger->filenames, warning_filename, (void*)&fp))
    {
      logger->warning_fp = fp;
    }
    else
    {
      logger->warning_fp = sunOpenLogFile(warning_filename, "w+");
      if (logger->warning_fp)
      {
        SUNHashMap_Insert(logger->filenames, warning_filename,
                          (void*)logger->warning_fp);
      }
      else
      {
        return -1;
      }
    }
#else
    fprintf(stderr,
            "[LOGGER WARNING] "
            "SUNDIALS_LOGGING_LEVEL=%d (build time option) "
            "is set too low for WARNING, but a WARNING file was provided. "
            "Set the logging level to >= %d and recompile if WARNING output "
            "level is desired.\n", SUN_LOGLEVEL_WARNING, SUNDIALS_LOGGING_LEVEL);
#endif
  }

  return 0;
}

int SUNLogger_SetInfoFilename(SUNLogger logger, const char* info_filename)
{
  if (logger == NULL)
  {
    return -1;
  }

  if (!sunLoggerIsOutputRank(logger, NULL))
  {
    return 0;
  }

  if (info_filename && strcmp(info_filename, ""))
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    FILE* fp = NULL;
    if (!SUNHashMap_GetValue(logger->filenames, info_filename, (void*)&fp))
    {
      logger->info_fp = fp;
    }
    else
    {
      logger->info_fp = sunOpenLogFile(info_filename, "w+");
      if (logger->info_fp)
      {
        SUNHashMap_Insert(logger->filenames, info_filename,
                          (void*)logger->info_fp);
      }
      else
      {
        return -1;
      }
    }
#else
    fprintf(stderr,
            "[LOGGER WARNING] "
            "SUNDIALS_LOGGING_LEVEL=%d (build time option) "
            "is set too low for INFO, but a INFO file was provided. Set the "
            "logging level to >= %d and recompile if INFO output level is "
            "desired.\n", SUN_LOGLEVEL_INFO, SUNDIALS_LOGGING_LEVEL);
#endif
  }

  return 0;
}

int SUNLogger_SetDebugFilename(SUNLogger logger, const char* debug_filename)
{
  if (logger == NULL)
  {
    return -1;
  }

  if (!sunLoggerIsOutputRank(logger, NULL))
  {
    return 0;
  }

  if (debug_filename && strcmp(debug_filename, ""))
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
    FILE* fp = NULL;
    if (!SUNHashMap_GetValue(logger->filenames, debug_filename, (void*)&fp))
    {
      logger->debug_fp = fp;
    }
    else
    {
      logger->debug_fp = sunOpenLogFile(debug_filename, "w+");
      if (logger->debug_fp)
      {
        SUNHashMap_Insert(logger->filenames, debug_filename,
                          (void*)logger->debug_fp);
      }
      else
      {
        return -1;
      }
    }
#else
    fprintf(stderr,
            "[LOGGER WARNING] "
            "SUNDIALS_LOGGING_LEVEL=%d (build time option) "
            "is set too low for DEBUG output, but a DEBUG file was provided. "
            "Set the logging level to >= %d and recompile if DEBUG output level "
            "is desired.\n", SUN_LOGLEVEL_DEBUG, SUNDIALS_LOGGING_LEVEL);
#endif
  }

  return 0;
}


int SUNLogger_QueueMsg(SUNLogger logger, SUNLogLevel lvl, const char* scope,
                       const char* label, const char* msg_txt, ...)
{
  int retval = 0;

  if (logger == NULL)
  {
    return -1;
  }

#if SUNDIALS_LOGGING_LEVEL > 0
  va_list args;
  va_start(args, msg_txt);

  if (logger->queuemsg)
  {
    retval = logger->queuemsg(logger, lvl, scope, label, msg_txt, args);
  }
  else
  {
    /* Default implementation */
    int rank = 0;
    if (sunLoggerIsOutputRank(logger, &rank))
    {
      char* log_msg = NULL;
      sunCreateLogMessage(lvl, rank, scope, label, msg_txt, args, &log_msg);

      switch (lvl)
      {
      case (SUN_LOGLEVEL_DEBUG):
        if (logger->debug_fp)
        {
          fprintf(logger->debug_fp, "%s", log_msg);
        }
        break;
      case (SUN_LOGLEVEL_WARNING):
        if (logger->warning_fp)
        {
          fprintf(logger->warning_fp, "%s", log_msg);
        }
        break;
      case (SUN_LOGLEVEL_INFO):
        if (logger->info_fp)
        {
          fprintf(logger->info_fp, "%s", log_msg);
        }
        break;
      case (SUN_LOGLEVEL_ERROR):
        if (logger->error_fp)
        {
          fprintf(logger->error_fp, "%s", log_msg);
        }
        break;
      default:
        retval = -1;
      }

      free(log_msg);
    }
  }

  va_end(args);
#endif

  return retval;
}

int SUNLogger_Flush(SUNLogger logger, SUNLogLevel lvl)
{
  int retval = 0;

  if (logger == NULL)
  {
    return -1;
  }

#if SUNDIALS_LOGGING_LEVEL > 0
  if (logger->flush)
  {
    retval = logger->flush(logger, lvl);
  }
  else
  {
    /* Default implementation */
    if (sunLoggerIsOutputRank(logger, NULL))
    {
      switch (lvl)
      {
      case (SUN_LOGLEVEL_DEBUG):
        if (logger->debug_fp)
        {
          fflush(logger->debug_fp);
        }
        break;
      case (SUN_LOGLEVEL_WARNING):
        if (logger->warning_fp)
        {
          fflush(logger->warning_fp);
        }
        break;
      case (SUN_LOGLEVEL_INFO):
        if (logger->info_fp)
        {
          fflush(logger->info_fp);
        }
        break;
      case (SUN_LOGLEVEL_ERROR):
        if (logger->error_fp)
        {
          fflush(logger->error_fp);
        }
        break;
      case (SUN_LOGLEVEL_ALL):
        if (logger->debug_fp)
        {
          fflush(logger->debug_fp);
        }
        if (logger->warning_fp)
        {
          fflush(logger->warning_fp);
        }
        if (logger->info_fp)
        {
          fflush(logger->info_fp);
        }
        if (logger->error_fp)
        {
          fflush(logger->error_fp);
        }
        break;
      default:
        retval = -1;
      }
    }
  }
#endif

  return retval;
}

int SUNLogger_GetOutputRank(SUNLogger logger, int* output_rank)
{
  int retval = 0;
  if (logger == NULL)
  {
    retval = -1;
  }
  else
  {
    *output_rank = logger->output_rank;
    retval       = 0;
  }
  return retval;
}

int SUNLogger_Destroy(SUNLogger* logger)
{
  int retval = 0;

  if ((*logger)->destroy)
  {
    retval = (*logger)->destroy(logger);
  }
  else
  {
    if (logger && (*logger))
    {
      /* Default implementation */
      if (sunLoggerIsOutputRank(*logger, NULL))
      {
        SUNHashMap_Destroy(&(*logger)->filenames, sunCloseLogFile);
      }

      free(*logger);
      *logger = NULL;
    }
  }

  return retval;
}
