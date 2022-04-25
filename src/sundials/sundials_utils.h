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
 * This header file contains common utility functions.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_UTILS_H
#define _SUNDIALS_UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

static int sunvsnprintf(char* buffer, size_t bufsz, const char* format, va_list vlist)
{
  int size = 0;
#ifdef SUNDIALS_COMPILER_HAS_SNPRINTF_AND_VA_COPY
  va_list tmp;
  va_copy(tmp, vlist);
  size = vsnprintf(buffer, bufsz, format, tmp);
  va_end(tmp);
#else
  size = SUNDIALS_MAX_SPRINTF_SIZE;
  if ((int) strlen(format) > size)
  {
    /* buffer is definitely not big enough */
    size = -1;
  }
  else if (buffer != NULL)
  {
    vsprintf(buffer, format, vlist);
  }
#endif
return size;
}


static int sunsnprintf(char* buffer, size_t bufsz, const char* format, ...)
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
static int sunvasnprintf(char** str, const char* fmt, va_list args)
{
  int size = 0;

  /* compute string length */
  size = sunvsnprintf(NULL, 0, fmt, args);

  if (size < 0)
  {
    return -1;
  }

  /* add one to size for the null terminator*/
  *str = (char*) malloc(size + 1);
  if (NULL == *str)
  {
    return -1;
  }

  size = vsprintf(*str, fmt, args);

  return size;
}


#endif /* _SUNDIALS_UTILS_H */
