/*******************************************************************
 *                                                                 *
 * File          : sundialstypes.h                                 *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This header file exports three types: realtype, integertype,    *
 * and booleantype, as well as the constants TRUE and FALSE.       *
 *                                                                 *
 * Users should #include "sundialstypes.h" in any file that        *
 * shhould be easily modifiable to work with different real or     *
 * integer types and use the exported names realtype and           *
 * integertype within such a file.                                 *
 * The types for realtype and integertype below have been set to   *
 * double and long int, respectively. A user should modify these   *
 * type declarations as he/she sees fit. For example, if a user    *
 * wants the work with type float because double precision         *
 * floating point arithmetic is too expensive on the user's        *
 * machine, then the definition below should be changed to:        *
 *                                                                 *
 * typedef float realtype;                                         *
 *                                                                 *
 * Similarly, if a user does not need to work with extremely large *
 * integers (see the system header file <limits.h> for the limits  *
 * on type int and long int on your machine), then the user        *
 * should change the definition below to:                          *
 *                                                                 *
 * typedef int integertype;                                        *
 *                                                                 *
 * The constants SUNDIALS_FLOAT, SUNDIALS_DOUBLE, SUNDIALS_INT,    *
 * SUNDIALS_LONG indicate the underlying types for realtype and    *
 * integertype. They should be set as follows:                     *
 *                                                                 *
 * (1) #define SUNDIALS_FLOAT 1                                    *
 *     #define SUNDIALS_DOUBLE 0     (real is float)               *
 *                                                                 *
 * (2) #define SUNDIALS_FLOAT 0                                    *
 *     #define SUNDIALS_DOUBLE 1     (real is double)              *
 *                                                                 *
 * (3) #define SUNDIALS_INT 1                                      *
 *     #define SUNDIALS_LONG 0   (integer is int)                  *
 *                                                                 *
 * (4) #define SUNDIALS_INT 0                                      *
 *     #define SUNDIALS_LONG 1   (integer is long int)             *
 *                                                                 *
 * Thus the legal types for realtype are float and double, while   *
 * the legal types for integertype are int and long int.           *
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


/******************************************************************
 *                                                                *
 * Types : realtype, integertype                                  *
 *----------------------------------------------------------------*
 * The types realtype and integertype are currently set to double *
 * and int, respectively. See the documentation at the top for    *
 * usage details and a description of associated constants and    *
 * macros.                                                        *
 *                                                                *
 ******************************************************************/

typedef double realtype;
typedef long int integertype;

#define SUNDIALS_FLOAT  0
#define SUNDIALS_DOUBLE 1

#define SUNDIALS_LONG 1
#define SUNDIALS_INT  0

#if SUNDIALS_FLOAT

#define RCONST(x) x##F

#elif SUNDIALS_DOUBLE

#define RCONST(x) x

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
