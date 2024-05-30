/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This header file contains common utility functions.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_UTILS_H
#define _SUNDIALS_UTILS_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

static inline char* sunCombineFileAndLine(int line, const char* file)
{
  size_t total_str_len = strlen(file) + 6;
  char* file_and_line  = (char*)malloc(total_str_len * sizeof(char));
  snprintf(file_and_line, total_str_len, "%s:%d", file, line);
  return file_and_line;
}

static inline int sunvsnprintf(char* buffer, size_t bufsz, const char* format,
                               va_list vlist)
{
  int size = 0;
  va_list tmp;
  va_copy(tmp, vlist);
  size = vsnprintf(buffer, bufsz, format, tmp);
  va_end(tmp);
  return size;
}

static inline int sunsnprintf(char* buffer, size_t bufsz, const char* format, ...)
{
  int size = 0;
  va_list args;
  va_start(args, format);
  size = sunvsnprintf(buffer, bufsz, format, args);
  va_end(args);
  return size;
}

/*
 * Implementation of the GNU extension function vasprintf which
 * is itself an analog for vsprintf, except it allocates a string
 * large enough to hold the output byte ('\0').
 */
static inline int sunvasnprintf(char** str, const char* fmt, va_list args)
{
  int size = 0;

  /* compute string length */
  size = sunvsnprintf(NULL, 0, fmt, args);

  if (size < 0) { return -1; }

  /* add one to size for the null terminator*/
  *str = (char*)malloc(size + 1);
  if (NULL == *str) { return -1; }

  size = vsnprintf(*str, size + 1, fmt, args);

  return size;
}

static inline void sunCompensatedSum(sunrealtype base, sunrealtype inc,
                                     sunrealtype* sum, sunrealtype* error)
{
  sunrealtype err           = *error;
  volatile sunrealtype tmp1 = inc - err;
  volatile sunrealtype tmp2 = base + tmp1;
  *error                    = (tmp2 - base) - tmp1;
  *sum                      = tmp2;
}

#endif /* _SUNDIALS_UTILS_H */
