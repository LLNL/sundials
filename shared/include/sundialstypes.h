/*******************************************************************
 *                                                                 *
 * File          : sundialstypes.h                                 *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, and          *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This header file exports two types: realtype and booleantype,   *
 *as well as the constants TRUE and FALSE.                         *
 *                                                                 *
 * Users should #include "sundialstypes.h" in any file that        *
 * shhould be easily modifiable to work with a different real      *
 * type and use the exported name realtype within such a file.     *
 *                                                                 *
 * The constant SUNDIALS_SINGLE_PRECISION indicates the underlying *
 * type for realtype. It is set at the configuration stage.        *
 *                                                                 *
 * The legal types for realtype are float and double.              *
 *                                                                 *
 * The macro RCONST gives a user a convenient way to define real   *
 * constants. To use the real constant 1.0, for example, the       *
 * user should write                                               *
 *                                                                 *
 * #define ONE RCONST(1.0)                                         *
 *                                                                 *
 * If realtype is double, then RCONST(1.0) expands to 1.0.         *
 * If realtype is float, then RCONST(1.0) expands to 1.0F.         *
 * There is never a need to explicitly cast 1.0 to (realtype).     *
 *                                                                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
  
#ifndef _sundialstypes_h
#define _sundialstypes_h
  
#include <float.h>
  
#if defined(SUNDIALS_SINGLE_PRECISION)

typedef float realtype;
#define RCONST(x) x##F
#define BIG_REAL FLT_MAX
#define UNIT_ROUNDOFF FLT_EPSILON

#else

typedef double realtype;
#define RCONST(x) x
#define BIG_REAL DBL_MAX
#define UNIT_ROUNDOFF DBL_EPSILON

#endif

/******************************************************************
 *                                                                *
 * Type : booleantype                                             *
 * Constants : FALSE, TRUE                                        *
 *----------------------------------------------------------------*
 * ANSI C does not have a built-in boolean type. Below is the     *
 * definition for a new type booleantype. The advantage of using  *
 * the name booleantype (instead of int) is an increase in code   *
 * readability.                                                   *
 * It allows the programmer to make a distinction between int and *
 * boolean data. Variables of type booleantype are intended to    *
 * have only the two values FALSE and TRUE which are defined      *
 * below to be equal to 0 and 1, respectively.                    *
 *                                                                *
 ******************************************************************/

#ifndef booleantype
#define booleantype int
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#endif

#ifdef __cplusplus
}
#endif
