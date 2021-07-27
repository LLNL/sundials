/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sundials/sundials_math.h>

realtype SUNRpowerI(realtype base, int exponent)
{
  int i, expt;
  realtype prod;

  prod = RCONST(1.0);
  expt = abs(exponent);
  for(i = 1; i <= expt; i++) prod *= base;
  if (exponent < 0) prod = RCONST(1.0)/prod;
  return(prod);
}

realtype SUNRpowerR(realtype base, realtype exponent)
{
  if (base <= RCONST(0.0)) return(RCONST(0.0));

#if defined(SUNDIALS_USE_GENERIC_MATH)
  return((realtype) pow((double) base, (double) exponent));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  return(pow(base, exponent));
#elif defined(SUNDIALS_SINGLE_PRECISION)
  return(powf(base, exponent));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  return(powl(base, exponent));
#endif
}

booleantype SUNRCompare(realtype a, realtype b)
{
  return(SUNRCompareTol(a, b, 10*UNIT_ROUNDOFF));
}

booleantype SUNRCompareTol(realtype a, realtype b, realtype tol)
{
  realtype diff;
  realtype norm;

  /* If a and b are exactly equal */
  if (a == b) return(SUNFALSE);

#if ( __STDC_VERSION__ >= 199901L)
  /* If a or b are NaN */
  if (isnan(a) || isnan(b)) return(SUNTRUE);

  /* If a or b are Inf */
  if (isinf(a) || isinf(b)) return(SUNTRUE);
#endif
  
  diff = SUNRabs(a - b);
  norm = SUNMIN(SUNRabs(a + b), BIG_REAL);
  
  /* When a + b is very small or zero, we use an absolute difference:
   *    |a - b| >= SMALL_REAL
   * Otherwise we use a relative difference:
   *    |a - b| < tol * |a + b|
   * The choice to use |a + b| over max(a, b)
   * is arbitrary, as is the choice to use
   * SMALL_REAL. 
   */
  return(diff >= SUNMAX(SMALL_REAL, tol*norm));
}
