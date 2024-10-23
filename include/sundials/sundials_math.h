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


  /*** NOTE(DRR): For now, I've added complex-specific definitions
       of all functions in this file (denoted as "SUNCxxx", e.g., SUNCsqrt
       and SUNCMIN).  Should these instead REPLACE the real-valued
       functions?  ***/


/*
 * -----------------------------------------------------------------
 * Function : SUNCCONJ
 * -----------------------------------------------------------------
 * Usage : suncomplextype sqrt_x;
 *         sqrt_x = SUNCCONJ(x);
 * -----------------------------------------------------------------
 * SUNCCONJ(x) returns the complex conjugate of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNCCONJ
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCCONJ(x) (conj((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCCONJ(x) (conjf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCCONJ(x) (conjl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Macros
 * -----------------------------------------------------------------
 * SUNMIN(A,B) returns the minimum of A and B
 * SUNCMIN(A,B) returns the whichever of A and B has minimum real part
 *
 * SUNMAX(A,B) returns the maximum of A and B
 * SUNCMAX(A,B) returns whichever of A and B has maximum real part
 *
 * SUNSQR(A) returns A^2
 * SUNCSQR(A) returns A*conj(A)
 *
 * SUNRsqrt calls the appropriate version of sqrt
 * SUNCsqrt calls the appropriate version of csqrt
 *
 * SUNRabs calls the appropriate version of abs
 * SUNCabs calls the appropriate version of cabs
 *
 * SUNRexp calls the appropriate version of exp
 * SUNCexp calls the appropriate version of cexp
 *
 * SUNRceil calls the appropriate version of ceil
 * SUNCceil calls the appropriate version of ceil on the real part of the argument
 * -----------------------------------------------------------------
 */

#ifndef SUNMIN
#define SUNMIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef SUNCMIN
#define SUNCMIN(A, B) ((SUN_REAL(A)) < (SUN_REAL(B)) ? (A) : (B))
#endif

#ifndef SUNMAX
#define SUNMAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef SUNCMAX
#define SUNCMAX(A, B) ((SUN_REAL(A)) > (SUN_REAL(B)) ? (A) : (B))
#endif

#ifndef SUNSQR
#define SUNSQR(A) ((A) * (A))
#endif

#ifndef SUNCSQR
#define SUNCSQR(A) ((A) * SUNCONJ(A))
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
 * Function : SUNCsqrt
 * -----------------------------------------------------------------
 * Usage : suncomplextype sqrt_x;
 *         sqrt_x = SUNCsqrt(x);
 * -----------------------------------------------------------------
 * SUNCsqrt(x) returns the complex square root of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNCsqrt
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCsqrt(x) (csqrt((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCsqrt(x) (csqrtf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCsqrt(x) (csqrtl((x)))
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
 * Function : SUNCabs
 * -----------------------------------------------------------------
 * Usage : suncomplextype abs_x;
 *         abs_x = SUNCabs(x);
 * -----------------------------------------------------------------
 * SUNCabs(x) returns the complex absolute value of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNCabs
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCabs(x) (cabs((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCabs(x) (cabsf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCabs(x) (cabsl((x)))
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
 * Function : SUNCexp
 * -----------------------------------------------------------------
 * Usage : suncomplextype exp_x;
 *         exp_x = SUNCexp(x);
 * -----------------------------------------------------------------
 * SUNCexp(x) returns e^x (base-e exponential function).
 * -----------------------------------------------------------------
 */

#ifndef SUNCexp
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCexp(x) (cexp((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCexp(x) (cexpf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCexp(x) (cexpl((x)))
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
 * Function : SUNCceil
 * -----------------------------------------------------------------
 * Usage : suncomplextype ceil_x;
 *         ceil_x = SUNCceil(x);
 * -----------------------------------------------------------------
 * SUNCceil(x) returns the smallest integer value not less than the
 * real part of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNCceil
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCceil(x) (ceil((SUN_REAL(x))))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCceil(x) (ceilf((SUN_REAL(x))))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCceil(x) (ceill((SUN_REAL(x))))
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
 * Function : SUNCpowerI
 * -----------------------------------------------------------------
 * Usage : int exponent;
 *         suncomplextype base, ans;
 *         ans = SUNCpowerI(base,exponent);
 * -----------------------------------------------------------------
 * SUNCpowerI returns the value of base^exponent, where base is of type
 * suncomplextype and exponent is of type int.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT suncomplextype SUNCpowerI(suncomplextype base, int exponent);

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
 * Function : SUNCpowerC
 * -----------------------------------------------------------------
 * Usage : suncomplextype base, exponent, ans;
 *         ans = SUNCpowerC(base,exponent);
 * -----------------------------------------------------------------
 * SUNCpowerR returns the value of base^exponent, where both base and
 * exponent are of type suncomplextype.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT suncomplextype SUNCpowerC(suncomplextype base, suncomplextype exponent);

/*
 * -----------------------------------------------------------------
 * Function : SUNRCompare
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         sunrealtype a, b;
 *         isNotEqual = SUNRCompare(a, b);
 * -----------------------------------------------------------------
 * SUNRCompare returns 0 if the relative difference of a and b is
 * less than or equal to 10*machine epsilon. If the relative
 * difference is greater than 10*machine epsilon, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b);

/*
 * -----------------------------------------------------------------
 * Function : SUNCCompare
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         suncomplextype a, b;
 *         isNotEqual = SUNCCompare(a, b);
 * -----------------------------------------------------------------
 * SUNCCompare returns 0 if the relative difference of a and b is
 * less than or equal to 10*machine epsilon. If the relative
 * difference is greater than 10*machine epsilon, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunbooleantype SUNCCompare(suncomplextype a, suncomplextype b);

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
 * Function : SUNCCompareTol
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         suncomplextype a, b, tol;
 *         isNotEqual = SUNCCompareTol(a, b, tol);
 * -----------------------------------------------------------------
 * SUNCCompareTol returns 0 if the relative difference of a and b is
 * less than or equal to the provided tolerance. If the relative
 * difference is greater than the tolerance, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunbooleantype SUNCCompareTol(suncomplextype a, suncomplextype b,
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
