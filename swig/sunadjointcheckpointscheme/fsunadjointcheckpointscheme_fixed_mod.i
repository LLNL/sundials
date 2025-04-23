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

%module fsunadjointcheckpointscheme_fixed_mod

// Include shared configuration
%include "../sundials/fsundials.i"

%include <stdint.i>

%{
#include "sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h"
%}

%import "../sundials/fsundials_core_mod.i"

// Process and wrap functions in the following files
%include  "sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h"
