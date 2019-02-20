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

%module fnvector_serial_mod

// include code common to all nvector implementations
%include "fnvector.i"

// nvector_impl macro defines some ignore and inserts with the vector name appended
%nvector_impl(Serial)

// Process and wrap functions in the following files
%include "nvector/nvector_serial.h"

