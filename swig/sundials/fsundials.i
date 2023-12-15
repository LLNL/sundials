// ---------------------------------------------------------------
// Programmer: Seth R. Johnson @ ORNL
//             Cody J. Balos @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2002-2023, Lawrence Livermore National Security
// and Southern Methodist University.
// All rights reserved.
//
// See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-3-Clause
// SUNDIALS Copyright End
// ---------------------------------------------------------------
// Swig interface file
// ---------------------------------------------------------------

// By default, wrap all constants as native fortran PARAMETERs
%fortranconst;

// Inform SWIG of the SUNDIALS_EXPORT macro
#define SUNDIALS_EXPORT
#define SUNDIALS_DEPRECATED_EXPORT
#define SUNDIALS_DEPRECATED_EXPORT_MSG(msg)
#define SUNDIALS_STATIC_INLINE

// All modules need sundials_types
%import "../sundials/fsundials_types_mod.i"

// Prefix all functions with F
// E.g. CVodeCreate -> FCVodeCreate
%rename("F%s", %$isfunction) "";

// Macro for creating the Fortran derived types for the generic SUNDIALS structures
%define %sundials_generic(TYPE)
  %fortran_struct(_generic_ ## TYPE ## _Ops);
  %typemap(ctype) _generic_ ## TYPE * "TYPE ## _Ops";
  %rename(TYPE ## _Ops) _generic_ ## TYPE ## _Ops;
  %fortran_struct(_generic_ ## TYPE);
  %typemap(ctype) _generic_ ## TYPE * "TYPE";
  %rename(TYPE) _generic_ ## TYPE;
%enddef

// Treat sundials generics as void pointers
%sundials_generic(N_Vector)
%sundials_generic(SUNLinearSolver)
%sundials_generic(SUNNonlinearSolver)
%sundials_generic(SUNMatrix)
%sundials_generic(SUNAdaptController)

// Treat FILE* as an opaque pointer
%apply void* { FILE* };

// Treat array of N_Vectors as an opaque pointer
%apply void* { N_Vector* };

// Assume sunrealtype* is an array of doubles
%apply double[] { sunrealtype* };

// Assume sunindextype* is an array of long int
%apply long int[] { sunindextype* };

// Assume int* is an array of integers
%apply int[] { int* };

// Assume long int * is an array of long int
%apply long int[] { long int* };

// Treat all ** as an opaque pointer
%apply void** { SWIGTYPE ** };

%include "fcopyright.i"
