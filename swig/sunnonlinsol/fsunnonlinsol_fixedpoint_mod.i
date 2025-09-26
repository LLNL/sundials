// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2025, Lawrence Livermore National Security,
// University of Maryland Baltimore County, and the SUNDIALS contributors.
// Copyright (c) 2013, Lawrence Livermore National Security
// and Southern Methodist University.
// Copyright (c) 2002, Lawrence Livermore National Security.
// All rights reserved.
//
// See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-3-Clause
// SUNDIALS Copyright End
// ---------------------------------------------------------------
// Swig interface file
// ---------------------------------------------------------------

%module fsunnonlinsol_fixedpoint_mod

// include code common to all implementations
%include "fsunnonlinsol.i"


%{
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"
%}

%sunnonlinsol_impl(FixedPoint)

// Process and wrap functions in the following files
%include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

