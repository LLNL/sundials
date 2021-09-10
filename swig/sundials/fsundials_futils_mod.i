// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2002-2021, Lawrence Livermore National Security
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

%module fsundials_futils_mod

%include "../sundials/fsundials.i"

// Insert code into the C wrapper to check that the sizes match
%{
#include "sundials/sundials_futils.h"
%}

// Process and wrap functions in the following files
%include "sundials/sundials_futils.h"

