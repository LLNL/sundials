/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for a simple C-language math library. The
 * routines listed here work with the type realtype as defined in
 * the header file sundials_types.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALSMATH_H
#define _SUNDIALSMATH_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_types.h>

/*
 * -----------------------------------------------------------------
 * Macros : MIN and MAX
 * -----------------------------------------------------------------
 * MIN(A,B) returns the minimum of A and B
 *
 * MAX(A,B) returns the maximum of A and B
 *
 * SQR(A) returns A^2
 * -----------------------------------------------------------------
 */

#ifndef SUN_MIN
#define SUN_MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef SUN_MAX
#define SUN_MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef SUN_SQR
#define SUN_SQR(A) ((A)*(A))
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

SUNDIALS_EXPORT realtype RPowerI(realtype base, int exponent);

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

SUNDIALS_EXPORT realtype RPowerR(realtype base, realtype exponent);

/*
 * -----------------------------------------------------------------
 * Function : SUN_SQRT
 * -----------------------------------------------------------------
 * Usage : realtype sqrt_x;
 *         sqrt_x = SUN_SQRT(x);
 * -----------------------------------------------------------------
 * SUN_SQRT(x) returns the square root of x. If x < ZERO, then
 * SUN_SQRT returns ZERO.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT realtype SUN_SQRT(realtype x);

/*
 * -----------------------------------------------------------------
 * Function : SUN_ABS
 * -----------------------------------------------------------------
 * Usage : realtype abs_x;
 *         abs_x = SUN_ABS(x);
 * -----------------------------------------------------------------
 * SUN_ABS(x) returns the absolute value of x.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT realtype SUN_ABS(realtype x);

/*
 * -----------------------------------------------------------------
 * Function : SUN_EXP
 * -----------------------------------------------------------------
 * Usage : realtype exp_x;
 *         exp_x = SUN_EXP(x);
 * -----------------------------------------------------------------
 * SUN_EXP(x) returns e^x (base-e exponential function).
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT realtype SUN_EXP(realtype x);

#ifdef __cplusplus
}
#endif

#endif
