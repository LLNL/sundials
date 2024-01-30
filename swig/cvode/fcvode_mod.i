// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
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

%module fcvode_mod

%include "../sundials/fsundials.i"

%{
#include "cvode/cvode.h"
#include "cvode/cvode_bandpre.h"
#include "cvode/cvode_bbdpre.h"
#include "cvode/cvode_diag.h"
#include "cvode/cvode_ls.h"
#include "cvode/cvode_proj.h"
%}

%import "../sundials/fsundials_core_mod.i"

// Process definitions from these files
%include "cvode/cvode.h"
%include "cvode/cvode_bandpre.h"
%include "cvode/cvode_bbdpre.h"
%include "cvode/cvode_diag.h"
%include "cvode/cvode_ls.h"
%include "cvode/cvode_proj.h"
