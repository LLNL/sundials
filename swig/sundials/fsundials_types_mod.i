// ---------------------------------------------------------------
// Programmer: Seth R. Johnson @ ORNL
//             Cody J. Balos @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2002-2019, Lawrence Livermore National Security
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

%module fsundials_types_mod

%include "../sundials/fsundials.i"
%include <stdint.i>

// Inform SWIG of the configure-provided types
#define SUNDIALS_INT64_T
#define SUNDIALS_INDEX_TYPE int64_t
#define SUNDIALS_DOUBLE_PRECISION
#define booleantype bool

// Insert code into the C wrapper to check that the sizes match
%{
#include "sundials/sundials.h"

#ifndef SUNDIALS_DOUBLE_PRECISION
#error "The Fortran bindings are only targeted at double-precision"
#endif
const char _your_int_size_is_wrong = ""[sizeof(SUNDIALS_INDEX_TYPE) != 8 ? -1 : 0];
%}

// Process and wrap functions in the following files
%include "sundials/sundials_types.h"

