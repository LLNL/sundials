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

%module fsundials_context_mod

%include "../sundials/fsundials.i"

// insert the include into the swig wrapper
%{
#include "sundials/sundials_context.h"
%}

%apply void* { SUNContext };
%apply void* { SUNProfiler };
%apply void** { SUNProfiler* };
%apply void* { SUNLogger };
%apply void** { SUNLogger* };

// Process and wrap functions in the following files
%include "sundials/sundials_context.h"
