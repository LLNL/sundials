// ---------------------------------------------------------------
// Programmer: Mustafa Aggul @ SMU
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

%module fsundomeigest_arnoldi_mod

// include code common to all implementations
%include "fsundomeigest.i"

%{
#include "sundomeigest/sundomeigest_arnoldi.h"
%}

// sundomeigest_impl macro defines some ignore and inserts with the estimator name appended
%sundomeigest_impl(ARNOLDI)

// Process and wrap functions in the following files
%include "sundomeigest/sundomeigest_arnoldi.h"