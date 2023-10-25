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

%module fsundials_profiler_mod

// Load the typedefs and generate a "use fsundials_types_mod" statement in the module
// %import "../sundials/fsundials_context_mod.i"
%import "../sundials/fsundials_types_mod.i"

%include "../sundials/fcopyright.i"

// insert the include into the swig wrapper
%{
#include "sundials/sundials_profiler.h"
#if SUNDIALS_MPI_ENABLED
#include <mpi.h>
#endif
%}

%apply void* { SUNProfiler };
%apply void** { SUNProfiler* };

// Utility class for C++ only.
%ignore SUNProfilerMarkScope;

// We have to manually insert the wrapper code for SUNProfiler_Create
// to handle the Fortran to MPI MPI_Comm translation.
%ignore SUNProfiler_Create;

// Process and wrap functions in the following files
%include "sundials/sundials_profiler.h"

%insert("wrapper") %{
SWIGEXPORT int _wrap_FSUNProfiler_Create(void *farg1, SwigArrayWrapper *farg2, void *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  char *arg2 = (char *) 0 ;
  SUNProfiler *arg3 = (SUNProfiler *) 0 ;
  int result;
#if SUNDIALS_MPI_ENABLED
  MPI_Comm comm;
#endif

  arg1 = (void *)(farg1);
  arg2 = (char *)(farg2->data);
  arg3 = (SUNProfiler *)(farg3);
#if SUNDIALS_MPI_ENABLED
  if (arg1 != NULL) {
    comm = MPI_Comm_f2c(*((MPI_Fint *) arg1));
    result = (int)SUNProfiler_Create((void*)&comm,(char const *)arg2,arg3);
  }
  else {
    result = (int)SUNProfiler_Create(arg1,(char const *)arg2,arg3);
  }
#else
  result = (int)SUNProfiler_Create(arg1,(char const *)arg2,arg3);
#endif
  fresult = (int)(result);
  return fresult;
}
%}

%insert("fdecl") %{
  public :: FSUNProfiler_Create
%}

%insert("finterfaces") %{
function swigc_FSUNProfiler_Create(farg1, farg2, farg3) &
bind(C, name="_wrap_FSUNProfiler_Create") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
import :: swigarraywrapper
type(C_PTR), value :: farg1
type(SwigArrayWrapper) :: farg2
type(C_PTR), value :: farg3
integer(C_INT) :: fresult
end function
%}

%insert("fsubprograms") %{
function FSUNProfiler_Create(comm, title, p) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: comm
character(kind=C_CHAR, len=*), target :: title
character(kind=C_CHAR), dimension(:), allocatable, target :: farg2_chars
type(C_PTR), target, intent(inout) :: p
integer(C_INT) :: fresult
type(C_PTR) :: farg1
type(SwigArrayWrapper) :: farg2
type(C_PTR) :: farg3

farg1 = comm
call SWIG_string_to_chararray(title, farg2_chars, farg2)
farg3 = c_loc(p)
fresult = swigc_FSUNProfiler_Create(farg1, farg2, farg3)
swig_result = fresult
end function
%}
