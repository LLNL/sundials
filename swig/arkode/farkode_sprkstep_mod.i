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

%module farkode_sprkstep_mod

%include "../sundials/fsundials.i"

// include the header file(s) in the c wrapper that is generated
%{
#include "arkode/arkode_sprkstep.h"
%}

// Load the typedefs and generate a "use" statements in the module
%import "farkode_mod.i"

// Process definitions from these files
%include "arkode/arkode_sprkstep.h"
