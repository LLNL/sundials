// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2002-2025, Lawrence Livermore National Security
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

// Ignore command-line processing functions since they are not supported in Fortran
%ignore SUNNonlinSolSetOptions;

// insert the include into the swig wrapper
%{
#include "sundials/sundials_nonlinearsolver.h"
%}

%ignore SUN_NLS_MSG_RESIDUAL;

// Process and wrap functions in the following files
%include "sundials/sundials_nonlinearsolver.h"

