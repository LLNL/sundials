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

%module fnvector_parallel_mod

// include code common to all nvector implementations
%include "fnvector.i"

// include the header file in the swig wrapper
%{
#include "nvector/nvector_parallel.h"
%}

// nvector_impl macro defines some ignore and inserts with the vector name appended
%nvector_impl(Parallel)

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
%#if SUNDIALS_MPI_ENABLED
    $1 = MPI_Comm_f2c(%static_cast(*$input, MPI_Fint));
%#else
    $1 = *$input;
%#endif
}
%typemap(out, noblock=1) MPI_Comm {
%#if SUNDIALS_MPI_ENABLED
    $result = %static_cast(MPI_Comm_c2f($1), int);
%#else
    $result = $1;
%#endif
}

// Process and wrap functions in the following files
%include "nvector/nvector_parallel.h"

%insert("wrapper") %{
SWIGEXPORT double * _wrap_FN_VGetArrayPointer_Parallel(N_Vector farg1) {
  double * fresult ;
  N_Vector arg1 = (N_Vector) 0 ;
  sunrealtype *result = 0 ;
  
  arg1 = (N_Vector)(farg1);
  result = (sunrealtype *)N_VGetArrayPointer_Parallel(arg1);
  fresult = result;
  return fresult;
}
%}

%insert("fdecl") %{
 public :: FN_VGetArrayPointer_Parallel
%}

%insert("finterfaces") %{
function swigc_FN_VGetArrayPointer_Parallel(farg1) &
bind(C, name="_wrap_FN_VGetArrayPointer_Parallel") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function
%}

%insert("fsubprograms") %{
function FN_VGetArrayPointer_Parallel(v) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(N_Vector), target, intent(inout) :: v
type(C_PTR) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(v)
fresult = swigc_FN_VGetArrayPointer_Parallel(farg1)
call c_f_pointer(fresult, swig_result, [FN_VGetLocalLength_Parallel(v)])
end function
%}
