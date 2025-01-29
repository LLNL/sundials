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

%module fcvodes_mod

%include "../sundials/fsundials.i"

%{
#include "cvodes/cvodes.h"
#include "cvodes/cvodes_bandpre.h"
#include "cvodes/cvodes_bbdpre.h"
#include "cvodes/cvodes_diag.h"
#include "cvodes/cvodes_ls.h"
%}

%import "../sundials/fsundials_core_mod.i"

// Process definitions from these files
%include "cvodes/cvodes.h"
%include "cvodes/cvodes_bandpre.h"
%include "cvodes/cvodes_bbdpre.h"
%include "cvodes/cvodes_diag.h"
%include "cvodes/cvodes_ls.h"

