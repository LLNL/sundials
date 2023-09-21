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
#include "sundials/sundials_adaptcontroller.h"
%}

// Load the typedefs and generate "use" statements
%import "../sundials/fsundials_adaptcontroller_mod.i"

// Macro for creating an interface to a SUNAdaptController
%define %sunadaptcontroller_impl(TYPE)
  %ignore _SUNAdaptControllerContent_## TYPE ##;
%enddef
