/*
 * -----------------------------------------------------------------
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
 * This is the header file for a simple C-language math library. The
 * routines listed here work with the type sunrealtype as defined in
 * the header file sundials_types.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALSMATH_H
#define _SUNDIALSMATH_H

#include <math.h>

#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Macros
 * -----------------------------------------------------------------
 * SUNMIN(A,B) returns the minimum of A and B
 *
 * SUNMAX(A,B) returns the maximum of A and B
 *
 * SUNSQR(A) returns A^2
 *
 * SUNRsqrt calls the appropriate version of sqrt
 *
 * SUNRabs calls the appropriate version of abs
 *
 * SUNRexp calls the appropriate version of exp
 *
 * SUNRceil calls the appropriate version of ceil
 * -----------------------------------------------------------------
 */

#ifndef SUNMIN
#define SUNMIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef SUNMAX
#define SUNMAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef SUNSQR
#define SUNSQR(A) ((A)*(A))
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRsqrt
 * -----------------------------------------------------------------
 * Usage : sunrealtype sqrt_x;
 *         sqrt_x = SUNRsqrt(x);
 * -----------------------------------------------------------------
 * SUNRsqrt(x) returns the square root of x. If x < ZERO, then
 * SUNRsqrt returns ZERO.
 * -----------------------------------------------------------------
 */

#ifndef SUNRsqrt
#if defined(__cplusplus) || defined(SUNDIALS_C_COMPILER_HAS_MATH_PRECISIONS)
#  if defined(SUNDIALS_DOUBLE_PRECISION)
#    define SUNRsqrt(x) ((x) <= RCONST(0.0) ? (RCONST(0.0)) : (sqrt((x))))
#  elif defined(SUNDIALS_SINGLE_PRECISION)
#    define SUNRsqrt(x) ((x) <= RCONST(0.0) ? (RCONST(0.0)) : (sqrtf((x))))
#  elif defined(SUNDIALS_EXTENDED_PRECISION)
#    define SUNRsqrt(x) ((x) <= RCONST(0.0) ? (RCONST(0.0)) : (sqrtl((x))))
#  endif
#else
#  define SUNRsqrt(x) ((x) <= RCONST(0.0) ? (RCONST(0.0)) : ((sunrealtype) sqrt((double) (x))))
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRabs
 * -----------------------------------------------------------------
 * Usage : sunrealtype abs_x;
 *         abs_x = SUNRabs(x);
 * -----------------------------------------------------------------
 * SUNRabs(x) returns the absolute value of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRabs
#if defined(__cplusplus) || defined(SUNDIALS_C_COMPILER_HAS_MATH_PRECISIONS)
#  if defined(SUNDIALS_DOUBLE_PRECISION)
#    define SUNRabs(x) (fabs((x)))
#  elif defined(SUNDIALS_SINGLE_PRECISION)
#    define SUNRabs(x) (fabsf((x)))
#  elif defined(SUNDIALS_EXTENDED_PRECISION)
#    define SUNRabs(x) (fabsl((x)))
#  endif
#else
#  define SUNRabs(x) ((sunrealtype) fabs((double) (x)))
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRexp
 * -----------------------------------------------------------------
 * Usage : sunrealtype exp_x;
 *         exp_x = SUNRexp(x);
 * -----------------------------------------------------------------
 * SUNRexp(x) returns e^x (base-e exponential function).
 * -----------------------------------------------------------------
 */

#ifndef SUNRexp
#if defined(__cplusplus) || defined(SUNDIALS_C_COMPILER_HAS_MATH_PRECISIONS)
#  if defined(SUNDIALS_DOUBLE_PRECISION)
#    define SUNRexp(x) (exp((x)))
#  elif defined(SUNDIALS_SINGLE_PRECISION)
#    define SUNRexp(x) (expf((x)))
#  elif defined(SUNDIALS_EXTENDED_PRECISION)
#    define SUNRexp(x) (expl((x)))
#  endif
#else
#  define SUNRexp(x) ((sunrealtype) exp((double) (x)))
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRceil
 * -----------------------------------------------------------------
 * Usage : sunrealtype ceil_x;
 *         ceil_x = SUNRceil(x);
 * -----------------------------------------------------------------
 * SUNRceil(x) returns the smallest integer value not less than x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRceil
#if defined(__cplusplus) || defined(SUNDIALS_C_COMPILER_HAS_MATH_PRECISIONS)
#  if defined(SUNDIALS_DOUBLE_PRECISION)
#    define SUNRceil(x) (ceil((x)))
#  elif defined(SUNDIALS_SINGLE_PRECISION)
#    define SUNRceil(x) (ceilf((x)))
#  elif defined(SUNDIALS_EXTENDED_PRECISION)
#    define SUNRceil(x) (ceill((x)))
#  endif
#else
#  define SUNRceil(x) ((sunrealtype) ceil((double) (x)))
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRpowerI
 * -----------------------------------------------------------------
 * Usage : int exponent;
 *         sunrealtype base, ans;
 *         ans = SUNRpowerI(base,exponent);
 * -----------------------------------------------------------------
 * SUNRpowerI returns the value of base^exponent, where base is of type
 * sunrealtype and exponent is of type int.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunrealtype SUNRpowerI(sunrealtype base, int exponent);

/*
 * -----------------------------------------------------------------
 * Function : SUNRpowerR
 * -----------------------------------------------------------------
 * Usage : sunrealtype base, exponent, ans;
 *         ans = SUNRpowerR(base,exponent);
 * -----------------------------------------------------------------
 * SUNRpowerR returns the value of base^exponent, where both base and
 * exponent are of type sunrealtype. If base < ZERO, then SUNRpowerR
 * returns ZERO.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunrealtype SUNRpowerR(sunrealtype base, sunrealtype exponent);

/*
 * -----------------------------------------------------------------
 * Function : SUNRCompare
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         sunrealtype a, b;
 *         isNotEqual = SUNRCompare(a, b);
 * -----------------------------------------------------------------
 * SUNRCompareTol returns 0 if the relative difference of a and b is
 * less than or equal to 10*machine epsilon. If the relative
 * difference is greater than 10*machine epsilon, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT booleantype SUNRCompare(sunrealtype a, sunrealtype b);

/*
 * -----------------------------------------------------------------
 * Function : SUNRCompareTol
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         sunrealtype a, b, tol;
 *         isNotEqual = SUNRCompareTol(a, b, tol);
 * -----------------------------------------------------------------
 * SUNRCompareTol returns 0 if the relative difference of a and b is
 * less than or equal to the provided tolerance. If the relative
 * difference is greater than the tolerance, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT booleantype SUNRCompareTol(sunrealtype a, sunrealtype b, sunrealtype tol);

/*
 * -----------------------------------------------------------------
 * Function : SUNStrToReal
 * -----------------------------------------------------------------
 * Usage : realtype a = SUNStrToReal(const char* str)
 * -----------------------------------------------------------------
 * SUNStrToReal parses str into the realtype variable. Uses standard
 * strtod variants when they are available.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunrealtype SUNStrToReal(const char* str);

#ifdef __cplusplus
}
#endif

#endif
