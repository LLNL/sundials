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

// handle MPI comm
%include <typemaps.i>

%apply int { MPI_Comm };
%typemap(ftype) MPI_Comm
   "integer"
%typemap(fin, noblock=1) MPI_Comm {
    $1 = int($input, C_INT)
}
%typemap(fout, noblock=1) MPI_Comm {
    $result = int($1)
}

%typemap(in, noblock=1) MPI_Comm {
%#ifdef SUNDIALS_MPI_ENABLED
    $1 = MPI_Comm_f2c(%static_cast(*$input, MPI_Fint));
%#else
    $1 = *$input;
%#endif
}
%typemap(out, noblock=1) MPI_Comm {
%#ifdef SUNDIALS_MPI_ENABLED
    $result = %static_cast(MPI_Comm_c2f($1), int);
%#else
    $result = $1;
%#endif
}

// Process and wrap functions in the following files
%include "nvector/nvector_mpimanyvector.h"

%insert("wrapper") %{
SWIGEXPORT double * _wrap_FN_VGetSubvectorArrayPointer_MPIManyVector(N_Vector farg1, int64_t const *farg2) {
  double * fresult ;
  N_Vector arg1 = (N_Vector) 0 ;
  sunindextype arg2 ;
  realtype *result = 0 ;
  
  arg1 = (N_Vector)(farg1);
  arg2 = (sunindextype)(*farg2);
  result = (realtype *)N_VGetSubvectorArrayPointer_MPIManyVector(arg1,arg2);
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
