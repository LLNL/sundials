/* -----------------------------------------------------------------
 * Programmer(s): Scott Cohen, Alan Hindmarsh, Radu Serban,
 *                Aaron Collier, and Slaven Peles @ LLNL
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
 * This header file exports three types: sunrealtype, sunindextype,
 * and sunbooleantype, as well as the constants SUNTRUE and SUNFALSE.
 *
 * Users should include the header file sundials_types.h in every
 * program file and use the exported name sunrealtype instead of
 * float, double or long double.
 *
 * The constants SUNDIALS_SINGLE_PRECISION, SUNDIALS_DOUBLE_PRECISION
 * and SUNDIALS_LONG_DOUBLE_PRECISION indicate the underlying data
 * type of sunrealtype.
 *
 * The legal types for sunrealtype are float, double and long double.
 *
 * The constants SUNDIALS_INT64_T and SUNDIALS_INT32_T indicate
 * the underlying data type of sunindextype -- the integer data type
 * used for vector and matrix indices.
 *
 * Data types are set at the configuration stage.
 *
 * The macro SUN_RCONST gives the user a convenient way to define
 * real-valued literal constants. To use the constant 1.0, for example,
 * the user should write the following:
 *
 *   #define ONE SUN_RCONST(1.0)
 *
 * If sunrealtype is defined as a double, then SUN_RCONST(1.0) expands
 * to 1.0. If sunrealtype is defined as a float, then SUN_RCONST(1.0)
 * expands to 1.0F. If sunrealtype is defined as a long double,
 * then SUN_RCONST(1.0) expands to 1.0L. There is never a need to
 * explicitly cast 1.0 to (sunrealtype). The macro can be used for
 * literal constants only. It cannot be used for expressions.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_TYPES_H
#define _SUNDIALS_TYPES_H

#include <float.h>
#include <stddef.h>
#include <stdint.h>
#include <sundials/sundials_config.h>

#if SUNDIALS_MPI_ENABLED
#include <mpi.h>
#endif

/*
 *------------------------------------------------------------------
 * Macro _SUNDIALS_STRUCT_
 * The _SUNDIALS_STRUCT_ macro is defined as a `struct` unless
 * generating the SWIG interfaces - in that case it is defined as
 * nothing. This is needed to work around a bug in SWIG which prevents
 * it from properly parsing our generic module structures.
 *------------------------------------------------------------------
 */
#ifdef SWIG
#define _SUNDIALS_STRUCT_
#else
#define _SUNDIALS_STRUCT_ struct
#endif

/*
 *------------------------------------------------------------------
 * Type sunrealtype
 * Macro SUN_RCONST
 * Constants SUN_SMALL_REAL, SUN_BIG_REAL, and SUN_UNIT_ROUNDOFF
 *------------------------------------------------------------------
 */

#if defined(SUNDIALS_SINGLE_PRECISION)

typedef float sunrealtype;
#define SUN_RCONST(x)     x##F
#define SUN_BIG_REAL      FLT_MAX
#define SUN_SMALL_REAL    FLT_MIN
#define SUN_UNIT_ROUNDOFF FLT_EPSILON

#elif defined(SUNDIALS_DOUBLE_PRECISION)

typedef double sunrealtype;
#define SUN_RCONST(x)     x
#define SUN_BIG_REAL      DBL_MAX
#define SUN_SMALL_REAL    DBL_MIN
#define SUN_UNIT_ROUNDOFF DBL_EPSILON

#elif defined(SUNDIALS_EXTENDED_PRECISION)

typedef long double sunrealtype;
#define SUN_RCONST(x)     x##L
#define SUN_BIG_REAL      LDBL_MAX
#define SUN_SMALL_REAL    LDBL_MIN
#define SUN_UNIT_ROUNDOFF LDBL_EPSILON

#endif

/*
 *------------------------------------------------------------------
 * Type suncomplextype
 * Macros SUN_CCONST, SUN_I
 *------------------------------------------------------------------
 */

/* safeguard the #include to avoid namespace conflicts with Kokkos, CUDA and HIP,
   and define "suncomplextype" based on the requested floating-point precision */
#if defined(__cplusplus)  /* C++ complex support */

  #if defined(SUNDIALS_KOKKOS_ENABLED)
    #define suncomplexlib Kokkos
    #include <Kokkos_Complex.hpp>
  #elif (defined(SUNDIALS_CUDA_ENABLED) || defined(SUNDIALS_HIP_ENABLED))
    #define suncomplexlib thrust
    #include <thrust/complex.h>
  #elif defined(SUNDIALS_EXTENDED_PRECISION)
    #include <complex.h>
  #else
    #define suncomplexlib std
    #include <complex>
  #endif

  #if defined(SUNDIALS_SINGLE_PRECISION)
typedef suncomplexlib::complex<float>  suncomplextype;
  #elif defined(SUNDIALS_DOUBLE_PRECISION)
typedef suncomplexlib::complex<double> suncomplextype;
  #elif defined(SUNDIALS_EXTENDED_PRECISION)
typedef __complex128                   suncomplextype;
  #endif

#else       /* C99 complex support */
  #include <complex.h>

  #if defined(SUNDIALS_SINGLE_PRECISION)
    #if defined(WIN32) || defined(_WIN32)
typedef _Fcomplex             suncomplextype;
    #else
typedef float _Complex        suncomplextype;
    #endif
  #elif defined(SUNDIALS_DOUBLE_PRECISION)
    #if defined(WIN32) || defined(_WIN32)
typedef _Dcomplex             suncomplextype;
    #else
typedef double _Complex       suncomplextype;
    #endif
  #elif defined(SUNDIALS_EXTENDED_PRECISION)
    #if defined(WIN32) || defined(_WIN32)
typedef _Lcomplex             suncomplextype;
    #else
typedef long double _Complex  suncomplextype;
    #endif
  #endif
#endif /* !__cplusplus */

#if defined(__cplusplus) && !defined(SUNDIALS_EXTENDED_PRECISION)
  #define SUN_I              (suncomplextype(SUN_RCONST(0.0), SUN_RCONST(1.0)))
  #define SUN_CCONST(x, y)   (suncomplextype((x), (y)))
#else
  #define SUN_I              (_Complex_I)
  #define SUN_CCONST(x, y)   (SUN_RCONST(x) + SUN_RCONST(y) * SUN_I)
#endif


/*
 *------------------------------------------------------------------
 * Type sunscalartype
 *------------------------------------------------------------------
 */

#if defined(SUNDIALS_SCALAR_TYPE_REAL)
typedef sunrealtype sunscalartype;
#elif defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
typedef suncomplextype sunscalartype;
#else
#error                                                                  \
  "SUNDIALS scalar type not defined, report to github.com/LLNL/sundials/issues"
#endif

/*
 *------------------------------------------------------------------
 * Type : sunindextype
 *------------------------------------------------------------------
 * Defines integer type to be used for vector and matrix indices.
 * User can build sundials to use 32- or 64-bit signed integers.
 * If compiler does not support portable data types, the SUNDIALS
 * CMake build system tries to find a type of the desired size.
 *------------------------------------------------------------------
 */

typedef SUNDIALS_INDEX_TYPE sunindextype;

/*
 *------------------------------------------------------------------
 * Type : sunbooleantype
 *------------------------------------------------------------------
 * Constants : SUNFALSE and SUNTRUE
 *------------------------------------------------------------------
 * ANSI C does not have a built-in boolean data type. Below is the
 * definition for a new type called sunbooleantype. The advantage of
 * using the name sunbooleantype (instead of int) is an increase in
 * code readability. It also allows the programmer to make a
 * distinction between int and boolean data. Variables of type
 * sunbooleantype are intended to have only the two values SUNFALSE and
 * SUNTRUE which are defined below to be equal to 0 and 1,
 * respectively.
 *------------------------------------------------------------------
 */

#ifndef sunbooleantype
#define sunbooleantype int
#endif

#ifndef SUNFALSE
#define SUNFALSE 0
#endif

#ifndef SUNTRUE
#define SUNTRUE 1
#endif

/*
 *------------------------------------------------------------------
 * Type : SUNOutputFormat
 *------------------------------------------------------------------
 * Constants for different output formats
 *------------------------------------------------------------------
 */

typedef enum
{
  SUN_OUTPUTFORMAT_TABLE,
  SUN_OUTPUTFORMAT_CSV
} SUNOutputFormat;

/*
 *------------------------------------------------------------------
 * Type : SUNErrCode
 *------------------------------------------------------------------
 * Error code type
 *------------------------------------------------------------------
 */

typedef int SUNErrCode;

/* -----------------------------------------------------------------------------
 * Forward declarations of SUNDIALS objects
 * ---------------------------------------------------------------------------*/

/* SUNDIALS context -- see sundials_context_impl.h */
typedef struct SUNContext_* SUNContext;

/* SUNDIALS error handler -- see sundials_errors.h */
typedef struct SUNErrHandler_* SUNErrHandler;

/* SUNDIALS profiler */
typedef struct SUNProfiler_* SUNProfiler;

/* SUNDIALS logger */
typedef struct SUNLogger_* SUNLogger;

/* -----------------------------------------------------------------------------
 * SUNDIALS function types
 * ---------------------------------------------------------------------------*/

#ifdef __cplusplus /* enable C++ linkage for SUNDIALS function pointer types */
extern "C" {
#endif

/* Error handler function */
typedef void (*SUNErrHandlerFn)(int line, const char* func, const char* file,
                                const char* msg, SUNErrCode err_code,
                                void* err_user_data, SUNContext sunctx);

#ifdef __cplusplus
}
#endif

/*
 *------------------------------------------------------------------
 * Type : SUNComm
 *------------------------------------------------------------------
 * SUNComm replaces MPI_Comm use in SUNDIALS code. It maps to
 * MPI_Comm when MPI is enabled.
 *------------------------------------------------------------------
 */

/* We don't define SUN_COMM_NULL when SWIG is processing the header
    because we manually insert the wrapper code for SUN_COMM_NULL
    (and %ignoring it in the SWIG code doesn't seem to work). */

#if SUNDIALS_MPI_ENABLED
#ifndef SWIG
#define SUN_COMM_NULL MPI_COMM_NULL
#endif
typedef MPI_Comm SUNComm;
#else
#ifndef SWIG
#define SUN_COMM_NULL 0
#endif
typedef int SUNComm;
#endif

#endif /* _SUNDIALS_TYPES_H */
