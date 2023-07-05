/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Utility functions for C examples
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include <sundials/sundials_config.h>

/* Check return flags */
int check_flag(int flag, const char* funcname)
{
  if (flag < 0)
  {
    fprintf(stderr, "ERROR: %s() returned %d\n", funcname, flag);
    return 1;
  }
  return 0;
}

/* Check return pointers */
int check_ptr(void* ptr, const char* funcname)
{
  if (!ptr)
  {
    fprintf(stderr, "ERROR: %s() returned NULL\n", funcname);
    return 1;
  }
  return 0;
}

/* Precision specific math function macros */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define EXP(x)  (exp((x)))
#define SIN(x)  (sin((x)))
#define COS(x)  (cos((x)))
#define SQRT(x) (sqrt((x)))
#define ABS(x)  (fabs((x)))
#define LOG(x)  (log((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define EXP(x)  (expf((x)))
#define SIN(x)  (sinf((x)))
#define COS(x)  (cosf((x)))
#define SQRT(x) (sqrtf((x)))
#define ABS(x)  (fabsf((x)))
#define LOG(x)  (logf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define EXP(x)  (expl((x)))
#define SIN(x)  (sinl((x)))
#define COS(x)  (cosl((x)))
#define SQRT(x) (sqrtl((x)))
#define ABS(x)  (fabsl((x)))
#define LOG(x)  (logl((x)))
#endif

/* Precision specific output macros */
#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif
