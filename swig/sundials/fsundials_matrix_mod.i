// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
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

%module fsundials_matrix_mod

// Load the typedefs and generate a "use fsundials_types_mod" statement in the module  
%import "../sundials/fsundials_context_mod.i"  
%import "../sundials/fsundials_types_mod.i"  
%import "../sundials/fsundials_nvector_mod.i"  

%include "../sundials/fcopyright.i"  

%{
#include "sundials/sundials_matrix.h"
%}

// Process and wrap functions in the following files
%include "sundials/sundials_matrix.h"

