/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2005-01-24 22:29:17 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file for a simple C-language math library. The
 * routines listed here work with the type realtype as defined in
 * the header file shared/include/sundialstypes.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALSMATH_H
#define _SUNDIALSMATH_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Macros : MIN and MAX
 * -----------------------------------------------------------------
 * MIN(A,B) returns the minimum of A and B
 *
 * MAX(A,B) returns the maximum of A and B
 * -----------------------------------------------------------------
 */

#ifndef MIN
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef ABS
#define ABS RAbs
#endif

#ifndef SQR
#define SQR RPower2
#endif

/*
 * -----------------------------------------------------------------
 * Function : RPowerI
 * -----------------------------------------------------------------
 * Usage : int exponent;
 *         realtype base, ans;
 *         ans = RPowerI(base,exponent);
 * -----------------------------------------------------------------
 * RPowerI returns the value of base^exponent, where base is of type
 * realtype and exponent is of type int.
 * -----------------------------------------------------------------
 */

realtype RPowerI(realtype base, int exponent);

/*
 * -----------------------------------------------------------------
 * Function : RPowerR
 * -----------------------------------------------------------------
 * Usage : realtype base, exponent, ans;
 *         ans = RPowerR(base,exponent);
 * -----------------------------------------------------------------
 * RPowerR returns the value of base^exponent, where both base and
 * exponent are of type realtype. If base < ZERO, then RPowerR
 * returns ZERO.
 * -----------------------------------------------------------------
 */

realtype RPowerR(realtype base, realtype exponent);

/*
 * -----------------------------------------------------------------
 * Function : RSqrt
 * -----------------------------------------------------------------
 * Usage : realtype sqrt_x;
 *         sqrt_x = RSqrt(x);
 * -----------------------------------------------------------------
 * RSqrt(x) returns the square root of x. If x < ZERO, then RSqrt
 * returns ZERO.
 * -----------------------------------------------------------------
 */

realtype RSqrt(realtype x);

/*
 * -----------------------------------------------------------------
 * Function : RAbs (a.k.a. ABS)
 * -----------------------------------------------------------------
 * Usage : realtype abs_x;
 *         abs_x = RAbs(x);
 * -----------------------------------------------------------------
 * RAbs(x) returns the absolute value of x.
 * -----------------------------------------------------------------
 */

realtype RAbs(realtype x);

/*
 * -----------------------------------------------------------------
 * Function : RPower2 (a.k.a. SQR)
 * -----------------------------------------------------------------
 * Usage : realtype sqr_x;
 *         sqr_x = RPower2(x);
 * -----------------------------------------------------------------
 * RPower2(x) returns x^2.
 * -----------------------------------------------------------------
 */

realtype RPower2(realtype x);

#ifdef __cplusplus
}
#endif

#endif
