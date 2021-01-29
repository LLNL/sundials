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

%module fnvector_pthreads_mod

// include code common to all nvector implementations
%include "fnvector.i"

// include the header file in the swig wraper
%{
#include "nvector/nvector_pthreads.h"
%}

// nvector_impl macro defines some ignore and inserts with the vector name appended
%nvector_impl(Pthreads)

// ignore the Pthreads_Data struct  
%ignore _Pthreads_Data;
  
// Process and wrap functions in the following files
%include "nvector/nvector_pthreads.h"

