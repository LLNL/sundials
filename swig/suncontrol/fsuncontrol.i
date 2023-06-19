// ---------------------------------------------------------------
// Programmer: Daniel R. Reynolds @ SMU
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2002-2023, Lawrence Livermore National Security
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

// Include shared configuration
%include "../sundials/fsundials.i"

%{
#include "sundials/sundials_control.h"
%}

// Load the typedefs and generate "use" statements
%import "../sundials/fsundials_control_mod.i"

// Macro for creating an interface to a SUNControl
%define %suncontrol_impl(TYPE)
  %ignore _SUNControlContent_## TYPE ##;
%enddef
