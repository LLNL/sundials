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
 * Function : SUNCONJ
 * -----------------------------------------------------------------
 * Usage : sunscalartype sqrt_x;
 *         sqrt_x = SUNCONJ(x);
 * -----------------------------------------------------------------
 * SUNCONJ(x) returns the complex conjugate of x if x is complex;
 * if x is real then it does nothing.
 * -----------------------------------------------------------------
 */

#ifndef SUNCONJ
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNCONJ(x) (x)
#else
#define SUNCONJ(x) SUNCCONJ(x)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Macros
 * -----------------------------------------------------------------
 * SUNRMIN(A,B) returns the minimum of the real numbers A and B
 * SUNCMIN(A,B) returns the whichever of A and B has minimum real part
 * SUNMIN(A,B) uses whichever of the above is mapped to sunscalartype
 *
 * SUNRMAX(A,B) returns the maximum of the real numbers A and B
 * SUNCMAX(A,B) returns whichever of A and B has maximum real part
 * SUNMAX(A,B) uses whichever of the above is mapped to sunscalartype
 *
 * SUNRSQR(A) returns A^2
 * SUNCSQR(A) returns A*conj(A)
 * SUNSQR(A) uses whichever of the above is mapped to sunscalartype
 *
 * SUNRsqrt calls the appropriate version of sqrt (real numbers)
 * SUNCsqrt calls the appropriate version of csqrt
 * SUNsqrt uses whichever of the above is mapped to sunrealtype
 *
 * SUNRabs calls the appropriate version of abs (real numbers)
 * SUNCabs calls the appropriate version of cabs
 * SUNabs uses whichever of the above is mapped to sunrealtype
 *
 * SUNRexp calls the appropriate version of exp (real numbers)
 * SUNCexp calls the appropriate version of cexp
 * SUNexp uses whichever of the above is mapped to sunrealtype
 *
 * SUNRceil calls the appropriate version of ceil (real numbers)
 * SUNCceil calls the appropriate version of ceil on the real part of the argument
 * SUNceil uses whichever of the above is mapped to sunrealtype
 * -----------------------------------------------------------------
 */

#ifndef SUNRMIN
#define SUNRMIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef SUNCMIN
#define SUNCMIN(A, B) ((SUN_REAL(A)) < (SUN_REAL(B)) ? (A) : (B))
#endif

#ifndef SUNMIN
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNMIN(A, B) SUNRMIN(A, B)
#else
#define SUNMIN(A, B) SUNCMIN(A, B)
#endif
#endif

#ifndef SUNRMAX
#define SUNRMAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef SUNCMAX
#define SUNCMAX(A, B) ((SUN_REAL(A)) > (SUN_REAL(B)) ? (A) : (B))
#endif

#ifndef SUNMAX
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNMAX(A, B) SUNRMAX(A, B)
#else
#define SUNMAX(A, B) SUNCMAX(A, B)
#endif
#endif

#ifndef SUNRSQR
#define SUNRSQR(A) ((A) * (A))
#endif

#ifndef SUNCSQR
#define SUNCSQR(A) ((A) * SUNCCONJ(A))
#endif

#ifndef SUNSQR
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNSQR(A) SUNRSQR(A)
#else
#define SUNSQR(A) SUNCSQR(A)
#endif
#endif


/*
 * -----------------------------------------------------------------
 * Function : SUNRsqrt, SUNCsqrt, and SUNsqrt
 * -----------------------------------------------------------------
 * Usage : sunrealtype sqrt_x;
 *         sqrt_x = SUNRsqrt(x);
 *         suncomplextype sqrt_y;
 *         sqrt_y = SUNCsqrt(y);
 *         sunscalartype sqrt_z;
 *         sqrt_z = SUNsqrt(z);
 * -----------------------------------------------------------------
 * SUNRsqrt(x) returns the square root of x. If x < ZERO, then
 * SUNRsqrt returns ZERO.
 * SUNCsqrt(y) returns the complex square root of x.
 * SUNsqrt(z) uses whichever of the above is mapped to sunscalartype
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

#ifndef SUNsqrt
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNsqrt(x) SUNRsqrt(x)
#else
#define SUNsqrt(x) SUNCsqrt(x)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRabs, SUNCabs, SUNabs
 * -----------------------------------------------------------------
 * Usage : sunrealtype abs_x;
 *         abs_x = SUNRabs(x);
 *         suncomplextype abs_y;
 *         abs_y = SUNCabs(y);
 *         sunscalartype abs_z;
 *         abs_z = SUNabs(z);
 * -----------------------------------------------------------------
 * SUNRabs(x) returns the absolute value of x.
 * SUNCabs(x) returns the complex absolute value of x.
 * SUNabs(x) uses whichever of the above is mapped to sunscalartype.
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

#ifndef SUNabs
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNabs(x) SUNRabs(x)
#else
#define SUNabs(x) SUNCabs(x)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRexp, SUNCexp, SUNexp
 * -----------------------------------------------------------------
 * Usage : sunrealtype exp_x;
 *         exp_x = SUNRexp(x);
 *         suncomplextype exp_y;
 *         exp_y = SUNCexp(y);
 *         sunscalartype exp_z;
 *         exp_z = SUNexp(z);
 * -----------------------------------------------------------------
 * SUNRexp(x) returns e^x (base-e exponential function) for
 * real-valued x.
 * SUNCexp(y) returns e^y (base-e exponential function) for
 * complex-valued y.
 * SUNexp(z) uses whichever of the above is mapped to sunscalartype.
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

#ifndef SUNexp
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNexp(x) SUNRexp(x)
#else
#define SUNexp(x) SUNCexp(x)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRceil, SUNCceil, SUNceil
 * -----------------------------------------------------------------
 * Usage : sunrealtype ceil_x;
 *         ceil_x = SUNRceil(x);
 *         suncomplextype ceil_y;
 *         ceil_y = SUNCceil(y);
 *         sunscalartype ceil_z;
 *         ceil_z = SUNceil(z);
 * -----------------------------------------------------------------
 * SUNRceil(x) returns the smallest integer value not less than x.
 * SUNCceil(x) returns the smallest integer value not less than the
 * real part of x.
 * SUNceil(z) uses whichever of the above is mapped to sunscalartype.
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

#ifndef SUNceil
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNceil(x) SUNRceil(x)
#else
#define SUNceil(x) SUNCceil(x)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRpowerI, SUNCpowerI, SUNpowerI
 * -----------------------------------------------------------------
 * Usage : int exponent;
 *         sunrealtype rbase, rans;
 *         rans = SUNRpowerI(rbase,exponent);
 *         suncomplextype cbase, cans;
 *         cans = SUNCpowerI(cbase,exponent);
 *         sunscalartype base, ans;
 *         ans = SUNpowerI(base,exponent);
 * -----------------------------------------------------------------
 * SUNRpowerI returns the value of base^exponent, where base is of
 * type sunrealtype and exponent is of type int.
 * SUNCpowerI returns the value of base^exponent, where base is of
 * type suncomplextype and exponent is of type int.
 * SUNpowerI uses whichever of the above is mapped to sunscalartype.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunrealtype SUNRpowerI(sunrealtype base, int exponent);
SUNDIALS_EXPORT suncomplextype SUNCpowerI(suncomplextype base, int exponent);

#ifndef SUNpowerI
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNpowerI(base,exponent) SUNRpowerI(base,exponent)
#else
#define SUNpowerI(base,exponent) SUNCpowerI(base,exponent)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRpowerR, SUNCpowerC, SUNpower
 * -----------------------------------------------------------------
 * Usage : sunrealtype rbase, rexponent, rans;
 *         rans = SUNRpowerR(rbase,rexponent);
 *         suncomplextype cbase, cexponent, cans;
 *         cans = SUNRpowerR(cbase,cexponent);
 *         sunscalartype base, exponent, ans;
 *         ans = SUNpower(base,exponent);
 * -----------------------------------------------------------------
 * SUNRpowerR returns the value of base^exponent, where both base and
 * exponent are of type sunrealtype. If base < ZERO, then SUNRpowerR
 * returns ZERO.
 * SUNCpowerC returns the value of base^exponent, where both base and
 * exponent are of type suncomplextype.
 * SUNpower uses whichever of the above is mapped to sunscalartype.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunrealtype SUNRpowerR(sunrealtype base, sunrealtype exponent);

#ifndef SUNCpowerC
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCpowerC(base,exponent) (cpow(base, exponent))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCpowerC(base,exponent) (cpowf(base, exponent))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCpowerC(base,exponent) (cpowl(base, exponent))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

#ifndef SUNpower
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNpower(base,exponent) SUNRpowerR(base,exponent)
#else
#define SUNpower(base,exponent) SUNCpowerC(base,exponent)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRCompare, SUNCCompare, SUNCompare
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         sunrealtype ra, rb;
 *         isNotEqual = SUNRCompare(ra, rb);
 *         suncomplextype ca, cb;
 *         isNotEqual = SUNCCompare(ca, cb);
 *         sunscalartype a, b;
 *         isNotEqual = SUNCompare(a, b);
 * -----------------------------------------------------------------
 * These functions return 0 if the relative difference of a and b is
 * less than or equal to 10*machine epsilon. If the relative
 * difference is greater than 10*machine epsilon, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * SUNRCompare is designed for real-valued inputs.
 * SUNCCompare is designed for complex-valued inputs.
 * SUNCompare uses whichever of the above is mapped to sunscalartype.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b);
SUNDIALS_EXPORT sunbooleantype SUNCCompare(suncomplextype a, suncomplextype b);

#ifndef SUNCompare
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNCompare(a,b) SUNRCompare(a,b)
#else
#define SUNCompare(a,b) SUNCCompare(a,b)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRCompareTol, SUNCompareTol, SUNCompareTol
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         sunrealtype tol;
 *         sunrealtype ra, rb;
 *         isNotEqual = SUNRCompareTol(ra, rb, tol);
 *         suncomplextype ca, cb;
 *         isNotEqual = SUNRCompareTol(ca, cb, tol);
 *         sunscalartype a, b;
 *         isNotEqual = SUNCompareTol(a, b, tol);
 * -----------------------------------------------------------------
 * These return 0 if the relative difference of a and b is
 * less than or equal to the provided tolerance. If the relative
 * difference is greater than the tolerance, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * SUNRCompareTol is designed for real-valued inputs.
 * SUNCCompareTol is designed for complex-valued inputs.
 * SUNCompareTol uses whichever of the above is mapped to
 * sunscalartype.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunbooleantype SUNRCompareTol(sunrealtype a, sunrealtype b,
                                              sunrealtype tol);
SUNDIALS_EXPORT sunbooleantype SUNCCompareTol(suncomplextype a, suncomplextype b,
                                              sunrealtype tol);

#ifndef SUNCompareTol
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNCompareTol(a,b) SUNRCompareTol(a,b,tol)
#else
#define SUNCompareTol(a,b) SUNCCompareTol(a,b,tol)
#endif
#endif

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
