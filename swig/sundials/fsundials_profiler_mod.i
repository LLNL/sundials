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

%module fsundials_profiler_mod

%include "../sundials/fsundials.i"

// insert the include into the swig wrapper
%{
#include "sundials/sundials_profiler.h"
#if SUNDIALS_MPI_ENABLED
#include <mpi.h>
#endif
%}

%apply void* { SUNProfiler };
%apply void** { SUNProfiler* };

// Utility class for C++ only.
%ignore SUNProfilerMarkScope;

// Process and wrap functions in the following files
%include "sundials/sundials_profiler.h"
