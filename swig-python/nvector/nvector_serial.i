// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2002-2019, Lawrence Livermore National Security
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

//%module nvector_serial
//%import "../sundials/sundials.i"

//%include "../nvector/nvector.i"
%nvector_impl(Serial)

// include the header file in the swig wrapper
%{
#include "nvector/nvector_serial.h"
%}

// Process and wrap functions in the following files
%include "nvector/nvector_serial.h"
