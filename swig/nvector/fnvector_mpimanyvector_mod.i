// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
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

%module fnvector_mpimanyvector_mod

// include code common to all nvector implementations
%include "fnvector.i"

// include the header file in the swig wrapper
%{
#include "nvector/nvector_mpimanyvector.h"
%}

// nvector_impl macro defines some ignore and inserts with the vector name appended
%nvector_impl(MPIManyVector)

// Process and wrap functions in the following files
%include "nvector/nvector_mpimanyvector.h"

%insert("wrapper") %{
SWIGEXPORT double * _wrap_FN_VGetSubvectorArrayPointer_MPIManyVector(N_Vector farg1, int64_t const *farg2) {
  double * fresult ;
  N_Vector arg1 = (N_Vector) 0 ;
  sunindextype arg2 ;
  sunrealtype *result = 0 ;
  
  arg1 = (N_Vector)(farg1);
  arg2 = (sunindextype)(*farg2);
  result = (sunrealtype *)N_VGetSubvectorArrayPointer_MPIManyVector(arg1,arg2);
  fresult = result;
  return fresult;
}
%}

%insert("fdecl") %{
 public :: FN_VGetSubvectorArrayPointer_MPIManyVector
%}

%insert("finterfaces") %{
function swigc_FN_VGetSubvectorArrayPointer_MPIManyVector(farg1, farg2) &
bind(C, name="_wrap_FN_VGetSubvectorArrayPointer_MPIManyVector") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT64_T), intent(in) :: farg2
type(C_PTR) :: fresult
end function
%}

%insert("fsubprograms") %{
function FN_VGetSubvectorArrayPointer_MPIManyVector(v, vec_num) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(N_Vector), target, intent(inout) :: v
integer(C_INT64_T), intent(in) :: vec_num
type(C_PTR) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT64_T) :: farg2 

farg1 = c_loc(v)
farg2 = vec_num
fresult = swigc_FN_VGetSubvectorArrayPointer_MPIManyVector(farg1, farg2)
call c_f_pointer(fresult, swig_result, [FN_VGetSubvectorLocalLength_MPIManyVector(v, vec_num)])
end function
%}
