// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
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

%module sundials_linearsolver

%include "../sundials/sundials.i"

%import "../sundials/sundials_nvector.i"
%import "../sundials/sundials_matrix.i"

// insert the include into the swig wrapper
%{
#include "sundials/sundials_iterative.h"
#include "sundials/sundials_linearsolver.h"
%}

// Process and wrap functions in the following files
%include "sundials/sundials_iterative.h"
%include "sundials/sundials_linearsolver.h"
