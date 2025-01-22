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


// insert the include into the swig wrapper
%{
#include "sundials/sundials_memory.h"
%}

%ignore SUNMemory_;
%apply void* { SUNMemory };
%apply void** { SUNMemory* };

%fortran_struct(SUNMemoryHelper_Ops_);
%typemap(ctype) SUNMemoryHelper_Ops_* "SUNAdjointCheckpointScheme_Ops";
%rename(SUNMemoryHelper_Ops) SUNMemoryHelper_Ops_;

%fortran_struct(SUNMemoryHelper_);
%typemap(ctype) SUNMemoryHelper_* "SUNMemoryHelper";
%rename(SUNMemoryHelper) SUNMemoryHelper;

// Process and wrap functions in the following files
%include "sundials/sundials_memory.h"
