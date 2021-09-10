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

// Include shared configuration
%include "../sundials/fsundials.i"

%{
#include "sundials/sundials_nvector.h"
%}

// Load the typedefs and generate "use" statements
%import "../sundials/fsundials_nvector_mod.i"

// Macro for creating an interface to an N_Vector
%define %nvector_impl(TYPE)
  %ignore _N_VectorContent_## TYPE ##;
  // Ignore functions with arrays of vector arrays since they are not supported
  %ignore N_VScaleAddMultiVectorArray_## TYPE ##;
  %ignore N_VLinearCombinationVectorArray_## TYPE ##;
  %ignore N_VEnableScaleAddMultiVectorArray_## TYPE ##;
  %ignore N_VEnableLinearCombinationVectorArray_## TYPE ##;
%enddef

