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

%module fsunlinsol_spgmr_mod

// include code common to all nvector implementations
%include "fsunlinsol.i"

%{
#include "sunlinsol/sunlinsol_spgmr.h"
%}

// sunlinsol_impl macro defines some ignore and inserts with the linear solver name appended
%sunlinsol_impl(SPGMR)

// Process and wrap functions in the following files
%include "sunlinsol/sunlinsol_spgmr.h"

