/*
 * -----------------------------------------------------------------
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
 * This is the header file for a simple C-language math library. The
 * routines listed here work with the type sunrealtype as defined in
 * the header file sundials_types.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALSMATH_H
#define _SUNDIALSMATH_H

#include <math.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
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
 * SUNRlog calls the appropriate version of log
 * 
 * SUNRsinh calls the appropriate version of sinh
 * 
 * SUNRcosh calls the appropriate version of cosh
 *
 * SUNRceil calls the appropriate version of ceil
 * 
 * SUNRfloor calls the appropriate version of floor
 * 
 * SUNRround calls the appropriate version of round
 * 
 * SUNIceil calls the appropriate version of ceil and returns sunindextype
 * 
 * SUNIfloor calls the appropriate version of floor and returns sunindextype
 * 
 * SUNIround calls the appropriate version of round and returns sunindextype
 * -----------------------------------------------------------------
 */

#ifndef SUNMIN
#define SUNMIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef SUNMAX
#define SUNMAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef SUNSQR
#define SUNSQR(A) ((A) * (A))
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
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRsqrt(x) ((x) <= SUN_RCONST(0.0) ? (SUN_RCONST(0.0)) : (sqrt((x))))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRsqrt(x) ((x) <= SUN_RCONST(0.0) ? (SUN_RCONST(0.0)) : (sqrtf((x))))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRsqrt(x) ((x) <= SUN_RCONST(0.0) ? (SUN_RCONST(0.0)) : (sqrtl((x))))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
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
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRabs(x) (fabs((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRabs(x) (fabsf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRabs(x) (fabsl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
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
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRexp(x) (exp((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRexp(x) (expf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRexp(x) (expl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRlog
 * -----------------------------------------------------------------
 * Usage : sunrealtype log_x;
 *         log_x = SUNRlog(x);
 * -----------------------------------------------------------------
 * SUNRlog(x) returns log(x) (base-e logarithmic function).
 * -----------------------------------------------------------------
 */

#ifndef SUNRlog
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRlog(x) (log((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRlog(x) (logf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRlog(x) (logl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRsinh
 * -----------------------------------------------------------------
 * Usage : sunrealtype sinh_x;
 *         sinh_x = SUNRsinh(x);
 * -----------------------------------------------------------------
 * SUNRsinh(x) returns sinh(x) (the hyperbolic sine of x).
 * -----------------------------------------------------------------
 */

#ifndef SUNRsinh
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRsinh(x) (sinh((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRsinh(x) (sinhf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRsinh(x) (sinhl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRcosh
 * -----------------------------------------------------------------
 * Usage : sunrealtype cosh_x;
 *         cosh_x = SUNRcosh(x);
 * -----------------------------------------------------------------
 * SUNRcosh(x) returns cosh(x) (the hyperbolic cosine of x).
 * -----------------------------------------------------------------
 */

#ifndef SUNRcosh
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRcosh(x) (cosh((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRcosh(x) (coshf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRcosh(x) (coshl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
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
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRceil(x) (ceil((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRceil(x) (ceilf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRceil(x) (ceill((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRfloor
 * -----------------------------------------------------------------
 * Usage : sunrealtype floor_x;
 *         floor_x = SUNRfloor(x);
 * -----------------------------------------------------------------
 * SUNRfloor(x) returns the largest integer value not greater than x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRfloor
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRfloor(x) (floor((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRfloor(x) (floorf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRfloor(x) (floorl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRround
 * -----------------------------------------------------------------
 * Usage : sunrealtype round_x;
 *         round_x = SUNRround(x);
 * -----------------------------------------------------------------
 * SUNRround(x) returns the nearest integer value to x (in floating-point format), 
 * rounding halfway cases away from zero, regardless of the current rounding mode.
 * -----------------------------------------------------------------
 */

#ifndef SUNRround
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRround(x) (round((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRround(x) (roundf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRround(x) (roundl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNIceil
 * -----------------------------------------------------------------
 * Usage : sunindextype ceil_x;
 *         ceil_x = SUNIceil(x);
 * -----------------------------------------------------------------
 * SUNIceil(x) returns the smallest sunindextype value not less than x.
 * -----------------------------------------------------------------
 */

#ifndef SUNIceil
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNIceil(x) ((sunindextype)(ceil((x)) + SUN_RCONST(0.5)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNIceil(x) ((sunindextype)(ceilf((x)) + SUN_RCONST(0.5)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNIceil(x) ((sunindextype)(ceill((x)) + SUN_RCONST(0.5)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNIfloor
 * -----------------------------------------------------------------
 * Usage : sunindextype floor_x;
 *         floor_x = SUNIfloor(x);
 * -----------------------------------------------------------------
 * SUNIfloor(x) returns the largest sunindextype value not greater than x.
 * -----------------------------------------------------------------
 */

#ifndef SUNIfloor
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNIfloor(x) ((sunindextype)(floor((x)) + SUN_RCONST(0.5)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNIfloor(x) ((sunindextype)(floorf((x)) + SUN_RCONST(0.5)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNIfloor(x) ((sunindextype)(floorl((x)) + SUN_RCONST(0.5)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNIround
 * -----------------------------------------------------------------
 * Usage : sunindextype round_x;
 *         round_x = SUNIround(x);
 * -----------------------------------------------------------------
 * SUNIround(x) returns the nearest sunindextype value to x, 
 * rounding halfway cases away from zero, regardless of the current rounding mode.
 * -----------------------------------------------------------------
 */

#ifndef SUNIround
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNIround(x) ((sunindextype)(round((x)) + SUN_RCONST(0.5)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNIround(x) ((sunindextype)(roundf((x)) + SUN_RCONST(0.5)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNIround(x) ((sunindextype)(roundl((x)) + SUN_RCONST(0.5)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
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

SUNDIALS_EXPORT sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b);

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

SUNDIALS_EXPORT sunbooleantype SUNRCompareTol(sunrealtype a, sunrealtype b,
                                              sunrealtype tol);

/*
 * -----------------------------------------------------------------
 * Function : SUNStrToReal
 * -----------------------------------------------------------------
 * Usage : sunrealtype a = SUNStrToReal(const char* str)
 * -----------------------------------------------------------------
 * SUNStrToReal parses str into the sunrealtype variable. Uses standard
 * strtod variants when they are available.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunrealtype SUNStrToReal(const char* str);

#ifdef __cplusplus
}
#endif

#endif
