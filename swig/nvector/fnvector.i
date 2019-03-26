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

// Include shared configuration
%include "../sundials/fsundials.i"

// Load the typedefs and generate "use" statements
%import "../sundials/fnvector_mod.i"

// Assume double* is an array of doubles
// NOTE: F2003 does not allow assumed-shape, but does allow assumed-size
%typemap(bindc, in="real(C_DOUBLE), dimension(*)") realtype* v_data "type(C_PTR)"

// Macro for creating an interface to an N_Vector
%define %nvector_impl(TYPE)
  %ignore _N_VectorContent_## TYPE ##;

  %insert(fdecl) %{
  public :: FN_VGetData_## TYPE ##
  %}

  %insert(fsubprograms) %{
  subroutine FN_VGetData_## TYPE ##(vec, vdata)

      use, intrinsic :: iso_c_binding
      implicit none

      type(C_PTR)        :: vec
      integer(C_INT64_T) :: len
      type(C_PTR)        :: cptr
      real(C_DOUBLE), dimension(:), pointer :: vdata

      len = FN_VGetLength_## TYPE ##(vec)
      cptr = FN_VGetArrayPointer_## TYPE ##(vec)

      call c_f_pointer(cptr, vdata, (/len/))

  end subroutine FN_VGetData_## TYPE ##
  %}
%enddef

