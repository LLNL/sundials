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

/*
 * Implementation of the GNU extension function vasprintf which
 * is itself an analog for vsprintf, except it allocates a string
 * large enough to hold the output byte ('\0').
 */
static inline int sunvasnprintf(char** str, const char* fmt, va_list args)
{
  int size = 0;

  /* compute string length */
  va_list tmp1;
  va_copy(tmp1, args);
  size = vsnprintf(NULL, 0, fmt, tmp1);
  va_end(tmp1);

  if (size < 0) { return -1; }

  /* add one to size for the null terminator*/
  *str = (char*)malloc(size + 1);
  if (NULL == *str) { return -1; }

  va_list tmp2;
  va_copy(tmp2, args);
  size = vsnprintf(*str, size + 1, fmt, tmp2);
  va_end(tmp2);

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
