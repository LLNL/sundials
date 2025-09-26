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

%module fkinsol_mod

%include "../sundials/fsundials.i"

// Ignore command-line processing functions since they are not supported in Fortran
%ignore KINSetOptions;

%{
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_bbdpre.h"
#include "kinsol/kinsol_ls.h"
%}

%import "../sundials/fsundials_core_mod.i"

// Process definitions from these files
%include "kinsol/kinsol.h"
%include "kinsol/kinsol_bbdpre.h"
%include "kinsol/kinsol_ls.h"
