/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Aaron Collier @ LLNL
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
 * This is the implementation file for a simple C-language math
 * library.
 * -----------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

int SUNIpowerI(int base, int exponent)
{
  int i;
  int prod = 1;

  for (i = 1; i <= exponent; i++) { prod *= base; }
  return (prod);
}

sunrealtype SUNRpowerI(sunrealtype base, int exponent)
{
  int i, expt;
  sunrealtype prod;

  prod = SUN_RCONST(1.0);
  expt = abs(exponent);
  for (i = 1; i <= expt; i++) { prod *= base; }
  if (exponent < 0) { prod = SUN_RCONST(1.0) / prod; }
  return (prod);
}

sunrealtype SUNRpowerR(sunrealtype base, sunrealtype exponent)
{
  if (base <= SUN_RCONST(0.0)) { return (SUN_RCONST(0.0)); }

#if defined(SUNDIALS_DOUBLE_PRECISION)
  return (pow(base, exponent));
#elif defined(SUNDIALS_SINGLE_PRECISION)
  return (powf(base, exponent));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  return (powl(base, exponent));
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
}

sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b)
{
  return (SUNRCompareTol(a, b, 10 * SUN_UNIT_ROUNDOFF));
}

sunbooleantype SUNRCompareTol(sunrealtype a, sunrealtype b, sunrealtype tol)
{
  sunrealtype diff;
  sunrealtype norm;

  /* If a and b are exactly equal.
   * This also covers the case where a and b are both inf under IEEE 754.
   */
  if (a == b) { return (SUNFALSE); }

  /* If a or b are NaN */
  if (isnan(a) || isnan(b)) { return (SUNTRUE); }

  /* If one of a or b are Inf (since we handled both being inf above) */
  if (isinf(a) || isinf(b)) { return (SUNTRUE); }

  diff = SUNRabs(a - b);
  norm = SUNMIN(SUNRabs(a + b), SUN_BIG_REAL);

  /* When |a + b| is very small (less than 10*SUN_UNIT_ROUNDOFF) or zero, we use an
   * absolute difference:
   *    |a - b| >= 10*SUN_UNIT_ROUNDOFF
   * Otherwise we use a relative difference:
   *    |a - b| < tol * |a + b|
   * The choice to use |a + b| over max(a, b)
   * is arbitrary, as is the choice to use
   * 10*SUN_UNIT_ROUNDOFF.
   */
  return (diff >= SUNMAX(10 * SUN_UNIT_ROUNDOFF, tol * norm));
}

sunrealtype SUNStrToReal(const char* str)
{
  char* end;
#if defined(SUNDIALS_EXTENDED_PRECISION)
  return strtold(str, &end);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  return strtod(str, &end);
#elif defined(SUNDIALS_SINGLE_PRECISION)
  return strtof(str, &end);
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
}
