/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Aaron Collier @ LLNL
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
 * This is the implementation file for a simple C-language math
 * library.
 * -----------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

static long double sunNextafterl(long double from, long double to);

static booleantype sunIsInf(sunrealtype a)
{
#if defined(__cplusplus) || defined(SUNDIALS_C_COMPILER_HAS_ISINF_ISNAN)
  return (isinf(a));
#else
  return (a < -BIG_REAL || a > BIG_REAL);
#endif
}

static booleantype sunIsNaN(sunrealtype a)
{
#if defined(__cplusplus) || defined(SUNDIALS_C_COMPILER_HAS_ISINF_ISNAN)
  return (isnan(a));
#else
  /* Most compilers/platforms follow NaN != a,
   * but since C89 does not require this, it is
   * possible some platforms might not follow it.
   */
  return (a != a);
#endif
}

sunrealtype SUNRpowerI(sunrealtype base, int exponent)
{
  int i, expt;
  sunrealtype prod;

  prod = RCONST(1.0);
  expt = abs(exponent);
  for (i = 1; i <= expt; i++) { prod *= base; }
  if (exponent < 0) { prod = RCONST(1.0) / prod; }
  return (prod);
}

sunrealtype SUNRpowerR(sunrealtype base, sunrealtype exponent)
{
  if (base <= RCONST(0.0)) { return (RCONST(0.0)); }

#if defined(__cplusplus) || defined(SUNDIALS_C_COMPILER_HAS_MATH_PRECISIONS)
#if defined(SUNDIALS_DOUBLE_PRECISION)
  return (pow(base, exponent));
#elif defined(SUNDIALS_SINGLE_PRECISION)
  return (powf(base, exponent));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  return (powl(base, exponent));
#endif
#else
  return ((sunrealtype)pow((double)base, (double)exponent));
#endif
}

booleantype SUNRCompare(sunrealtype a, sunrealtype b)
{
  return (SUNRCompareTol(a, b, 10 * UNIT_ROUNDOFF));
}

booleantype SUNRCompareTol(sunrealtype a, sunrealtype b, sunrealtype tol)
{
  sunrealtype diff;
  sunrealtype norm;

  /* If a and b are exactly equal.
   * This also covers the case where a and b are both inf under IEEE 754.
   */
  if (a == b) { return (SUNFALSE); }

  /* If a or b are NaN */
  if (sunIsNaN(a) || sunIsNaN(b)) { return (SUNTRUE); }

  /* If one of a or b are Inf (since we handled both being inf above) */
  if (sunIsInf(a) || sunIsInf(b)) { return (SUNTRUE); }

  diff = SUNRabs(a - b);
  norm = SUNMIN(SUNRabs(a + b), BIG_REAL);

  /* When |a + b| is very small (less than 10*UNIT_ROUNDOFF) or zero, we use an
   * absolute difference:
   *    |a - b| >= 10*UNIT_ROUNDOFF
   * Otherwise we use a relative difference:
   *    |a - b| < tol * |a + b|
   * The choice to use |a + b| over max(a, b)
   * is arbitrary, as is the choice to use
   * 10*UNIT_ROUNDOFF.
   */
  return (diff >= SUNMAX(10 * UNIT_ROUNDOFF, tol * norm));
}

long double sunNextafterl(long double from, long double to)
{
#if defined(__cplusplus) || defined(SUNDIALS_C_COMPILER_HAS_MATH_PRECISIONS)
  return nextafterl(from, to);
#else
  union
  {
    long double f;
    int i;
  } u;

  u.i = 0;
  u.f = from;

  /* if either are NaN, then return NaN via the sum */
  if (sunIsNaN((sunrealtype)from) || sunIsNaN((sunrealtype)to))
  {
    return from + to;
  }

  if (from == to) { return to; }

  /* ordering is -0.0, +0.0 so nextafter(-0.0, 0.0) should give +0.0
     and nextafter(0.0, -0.0) should give -0.0 */
  if (from == 0)
  {
    u.i = 1;
    return to > 0 ? u.f : -u.f;
  }

  if ((from > 0) == (to > from)) { u.i++; }
  else { u.i--; }

  return u.f;
#endif
}

sunrealtype SUNStrToReal(const char* str)
{
  char* end;
#if defined(__cplusplus) || defined(SUNDIALS_C_COMPILER_HAS_MATH_PRECISIONS)
#if defined(SUNDIALS_EXTENDED_PRECISION)
  return strtold(str, &end);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  return strtod(str, &end);
#elif defined(SUNDIALS_SINGLE_PRECISION)
  return strtof(str, &end);
#else
#error \
  "Should not be here, no SUNDIALS precision defined, report to github.com/LLNL/sundials/issues"
#endif
#else
#if defined(SUNDIALS_EXTENDED_PRECISION)
  /* Use strtod, but then round down to the closest double value
     since strtod will effectively round up to the closest long double. */
  double val = strtod(str, &end);
  return (sunrealtype)sunNextafterl(val, -0.0);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  return strtod(str, &end);
#elif defined(SUNDIALS_SINGLE_PRECISION)
  return strtod(str, &end);
#endif
#endif
}
