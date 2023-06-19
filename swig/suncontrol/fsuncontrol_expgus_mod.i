// ---------------------------------------------------------------
// Programmer: Daniel R. Reynolds @ SMU
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

%module fsuncontrol_expgus_mod

// include code common to all implementations
%include "fsuncontrol.i"

%{
#include "suncontrol/suncontrol_expgus.h"
%}

%suncontrol_impl(ExpGus)

// Process and wrap functions in the following files
%include "suncontrol/suncontrol_expgus.h"
