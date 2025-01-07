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

%include <stdint.i>

%{
#include "sunadjoint/sunadjoint_stepper.h"
#include "sunadjoint/sunadjoint_checkpointscheme.h"
%}


%ignore SUNAdjointStepper_;
%apply void* { SUNAdjointStepper };
%apply void** { SUNAdjointStepper* };

%fortran_struct(SUNAdjointCheckpointScheme_Ops_);
%typemap(ctype) SUNAdjointCheckpointScheme_Ops_* "SUNAdjointCheckpointScheme_Ops";
%rename(SUNAdjointCheckpointScheme_Ops) SUNAdjointCheckpointScheme_Ops_;

%fortran_struct(SUNAdjointCheckpointScheme_);
%typemap(ctype) SUNAdjointCheckpointScheme_* "SUNAdjointCheckpointScheme";
%rename(SUNAdjointCheckpointScheme) SUNAdjointCheckpointScheme_;


// Process and wrap functions in the following files
%include "sunadjoint/sunadjoint_checkpointscheme.h"
%include "sunadjoint/sunadjoint_stepper.h"

