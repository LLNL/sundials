// ---------------------------------------------------------------
// Programmer: Steven B. Roberts @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2002-2024, Lawrence Livermore National Security
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

%{
#include "sundials/sundials_stepper.h"
%}

%apply void* { SUNStepper };

// Process and wrap functions in the following files
%include "sundials/sundials_stepper.h"

