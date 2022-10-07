// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2002-2022, Lawrence Livermore National Security
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

// Include shared configuration
%include "../sundials/fsundials.i"

%{
#include "sundials/sundials_linearsolver.h"
%}

// Load the typedefs and generate "use" statements
%import "../sundials/fsundials_linearsolver_mod.i"

// Macro for creating an interface to an N_Vector
%define %sunlinsol_impl(TYPE)
  %ignore _SUNLinearSolverContent_## TYPE ##;
%enddef

