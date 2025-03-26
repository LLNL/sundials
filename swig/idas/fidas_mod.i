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

%module fidas_mod

%include "../sundials/fsundials.i"

%{
#include "idas/idas.h"
#include "idas/idas_bbdpre.h"
#include "idas/idas_ls.h"
%}

%import "../sundials/fsundials_core_mod.i"

// Process definitions from these files
%include "idas/idas.h"
%include "idas/idas_bbdpre.h"
%include "idas/idas_ls.h"

