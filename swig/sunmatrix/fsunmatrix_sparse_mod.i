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

%module fsunmatrix_sparse_mod

// include code common to all nvector implementations
%include "fsunmatrix.i"

%{
#include "sunmatrix/sunmatrix_sparse.h"
%}

// sunmatrix_impl macro defines some ignore and inserts with the matrix name appended
%sunmatrix_impl(Sparse)

// Process and wrap functions in the following files
%include "sunmatrix/sunmatrix_sparse.h"

