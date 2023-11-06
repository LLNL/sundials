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

%module fsundials_logger_mod

%include "../sundials/fsundials.i"

// insert the include into the swig wrapper
%{
#include "sundials/sundials_logger.h"
#if SUNDIALS_MPI_ENABLED
#include <mpi.h>
#endif
%}

%apply void* { SUNLogger };
%apply void** { SUNLogger* };


// We have to manually insert the wrapper code for SUNLogger_Create
// to handle the Fortran to MPI MPI_Comm translation.
//%ignore SUNLogger_Create;
// %ignore SUNLogger_Destroy;

// Process and wrap functions in the following files
%include "sundials/sundials_logger.h"

// %insert("wrapper") %{
// SWIGEXPORT int _wrap_FSUNLogger_Create(void *farg1, int const *farg2, void *farg3) {
//   int fresult ;
//   void *arg1 = (void *) 0 ;
//   int arg2 ;
//   SUNLogger *arg3 = (SUNLogger *) 0 ;
//   int result;
// #if SUNDIALS_MPI_ENABLED
//   MPI_Comm comm;
// #endif

//   arg1 = (void *)(farg1);
//   arg2 = (int)(*farg2);
//   arg3 = (SUNLogger *)(farg3);
// #if SUNDIALS_MPI_ENABLED
//   if (arg1 != NULL) {
//     comm = MPI_Comm_f2c(*((MPI_Fint *) arg1));
//     result = (int)SUNLogger_Create((void*)&comm,arg2,arg3);
//   }
//   else {
//     result = (int)SUNLogger_Create(arg1,arg2,arg3);
//   }
// #else
//   result = (int)SUNLogger_Create(arg1,arg2,arg3);
// #endif
//   fresult = (int)(result);
//   return fresult;
// }
// %}

// %insert("fdecl") %{
//   public :: FSUNLogger_Create
// %}

// %insert("finterfaces") %{
// function swigc_FSUNLogger_Create(farg1, farg2, farg3) &
// bind(C, name="_wrap_FSUNLogger_Create") &
// result(fresult)
// use, intrinsic :: ISO_C_BINDING
// type(C_PTR), value :: farg1
// integer(C_INT), intent(in) :: farg2
// type(C_PTR), value :: farg3
// integer(C_INT) :: fresult
// end function
// %}

// %insert("fsubprograms") %{
// function FSUNLogger_Create(comm, output_rank, logger) &
// result(swig_result)
// use, intrinsic :: ISO_C_BINDING
// integer(C_INT) :: swig_result
// type(C_PTR) :: comm
// integer(C_INT), intent(in) :: output_rank
// type(C_PTR), target, intent(inout) :: logger
// integer(C_INT) :: fresult
// type(C_PTR) :: farg1
// integer(C_INT) :: farg2
// type(C_PTR) :: farg3

// farg1 = comm
// farg2 = output_rank
// farg3 = c_loc(logger)
// fresult = swigc_FSUNLogger_Create(farg1, farg2, farg3)
// swig_result = fresult
// end function
// %}
