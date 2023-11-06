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

%module fsunmatrix_dense_mod

// include code common to all nvector implementations
%include "fsunmatrix.i"

%{
#include "sunmatrix/sunmatrix_dense.h"
%}

// sunmatrix_impl macro defines some ignore and inserts with the matrix name appended
%sunmatrix_impl(Dense)

// we manually insert so that the correct shape array is returned
%ignore SUNDenseMatrix_Data;
%ignore SUNDenseMatrix_Column;

// Process and wrap functions in the following files
%include "sunmatrix/sunmatrix_dense.h"

%insert("wrapper") %{
SWIGEXPORT double * _wrap_FSUNDenseMatrix_Data(SUNMatrix farg1) {
  double * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunrealtype *result = 0 ;
  
  arg1 = (SUNMatrix)(farg1);
  result = (sunrealtype *)SUNDenseMatrix_Data(arg1);
  fresult = result;
  return fresult;
}

SWIGEXPORT double * _wrap_FSUNDenseMatrix_Column(SUNMatrix farg1, int64_t const *farg2) {
  double * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunindextype arg2 ;
  sunrealtype *result = 0 ;
  
  arg1 = (SUNMatrix)(farg1);
  arg2 = (sunindextype)(*farg2);
  result = (sunrealtype *)SUNDenseMatrix_Column(arg1,arg2);
  fresult = result;
  return fresult;
}
%}

%insert("fdecl") %{
 public :: FSUNDenseMatrix_Data
 public :: FSUNDenseMatrix_Column
%}

%insert("finterfaces") %{
function swigc_FSUNDenseMatrix_Data(farg1) &
bind(C, name="_wrap_FSUNDenseMatrix_Data") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function

function swigc_FSUNDenseMatrix_Column(farg1, farg2) &
bind(C, name="_wrap_FSUNDenseMatrix_Column") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT64_T), intent(in) :: farg2
type(C_PTR) :: fresult
end function
%}

%insert("fsubprograms") %{
function FSUNDenseMatrix_Data(a) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(SUNMatrix), target, intent(inout) :: a
type(C_PTR) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(a)
fresult = swigc_FSUNDenseMatrix_Data(farg1)
call c_f_pointer(fresult, swig_result, [FSUNDenseMatrix_LData(a)])
end function

function FSUNDenseMatrix_Column(a, j) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(SUNMatrix), target, intent(inout) :: a
integer(C_INT64_T), intent(in) :: j
type(C_PTR) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT64_T) :: farg2 

farg1 = c_loc(a)
farg2 = j
fresult = swigc_FSUNDenseMatrix_Column(farg1, farg2)
call c_f_pointer(fresult, swig_result, [FSUNDenseMatrix_Rows(a)])
end function
%}
