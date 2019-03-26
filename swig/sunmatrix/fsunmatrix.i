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
%import "../sundials/fsunmatrix_mod.i"

// Macro for creating an interface to an N_Vector
%define %sunmatrix_impl(TYPE)
  %ignore _SUNMatrixContent_## TYPE ##;

  %insert(fdecl) %{
  public :: FSUNMatGetData_## TYPE ##
  %}

  %insert(fsubprograms) %{
  subroutine FSUNMatGetData_## TYPE ##(mat, mdata)

      use, intrinsic :: iso_c_binding
      implicit none

      type(C_PTR)        :: mat
      integer(C_INT64_T) :: M, N
      type(C_PTR)        :: cptr
      real(C_DOUBLE), dimension(:,:), pointer :: mdata

      M = FSUN## TYPE ##Matrix_Rows(mat)
      N = FSUN## TYPE ##Matrix_Columns(mat)
      cptr = FSUN## TYPE ##Matrix_Data(mat)

      call c_f_pointer(cptr, mdata, (/M,N/))

  end subroutine FSUNMatGetData_## TYPE ##
  %}
%enddef

