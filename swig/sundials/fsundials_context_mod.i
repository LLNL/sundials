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

%module fsundials_context_mod

// Load the typedefs and generate a "use fsundials_types_mod" statement in the module
%import "../sundials/fsundials_types_mod.i"

%include "../sundials/fcopyright.i"

// insert the include into the swig wrapper
%{
#include "sundials/sundials_context.h"
#include "sundials/sundials_profiler.h"
%}

%apply void* { SUNContext };
%apply void* { SUNProfiler };
%apply void** { SUNProfiler* };
%apply void* { SUNLogger };
%apply void** { SUNLogger* };

// We have to manually insert the wrapper code for SUNContext_Create,
// and SUNContext_Free to handle the Fortran to MPI MPI_Comm translation.
%ignore SUNContext_Create;
%ignore SUNContext_Free;

// Process and wrap functions in the following files
%include "sundials/sundials_context.h"

%insert("wrapper") %{
SWIGEXPORT int _wrap_FSUNContext_Free(void *farg1) {
  int fresult ;
  SUNContext *arg1 = (SUNContext *) 0 ;
  int result;
#ifdef SUNDIALS_BUILD_WITH_PROFILING
  SUNProfiler profiler;
#endif

  arg1 = (SUNContext *)(farg1);
#ifdef SUNDIALS_BUILD_WITH_PROFILING
  result = (int)SUNContext_GetProfiler(*arg1,&profiler);
  result = (int)SUNContext_Free(arg1);
  result = (int)SUNProfiler_Free(&profiler);
#else
  result = (int)SUNContext_Free(arg1);
#endif
  fresult = (int)(result);
  return fresult;
}

SWIGEXPORT int _wrap_FSUNContext_Create(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  SUNContext *arg2 = (SUNContext *) 0 ;
  int result;

  arg1 = (void *)(farg1);
  arg2 = (SUNContext *)(farg2);
  result = (int)SUNContext_Create(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}
%}

%insert("fbegin") %{
#include "sundials/sundials_config.h"
%}

%insert("fdecl") %{
public :: FSUNContext_Free
public :: FSUNContext_Create
%}

%insert("finterfaces") %{
function swigc_FSUNContext_Free(farg1) &
bind(C, name="_wrap_FSUNContext_Free") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

function swigc_FSUNContext_Create(farg1, farg2) &
bind(C, name="_wrap_FSUNContext_Create") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
integer(C_INT) :: fresult
end function

%}

%insert("fsubprograms") %{
function FSUNContext_Free(ctx) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR), target, intent(inout) :: ctx
integer(C_INT) :: fresult
type(C_PTR) :: farg1

farg1 = c_loc(ctx)
fresult = swigc_FSUNContext_Free(farg1)
swig_result = fresult
end function

function FSUNContext_Create(comm, ctx) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
#ifdef SUNDIALS_BUILD_WITH_PROFILING
use fsundials_profiler_mod, only : FSUNProfiler_Create
#endif
integer(C_INT) :: swig_result
type(C_PTR) :: comm
type(C_PTR), target, intent(inout) :: ctx
integer(C_INT) :: fresult
type(C_PTR) :: farg1
type(C_PTR) :: farg2
#ifdef SUNDIALS_BUILD_WITH_PROFILING
type(C_PTR) :: profiler
#endif

farg1 = comm
farg2 = c_loc(ctx)
fresult = swigc_FSUNContext_Create(c_null_ptr, farg2)
#ifdef SUNDIALS_BUILD_WITH_PROFILING
fresult = FSUNProfiler_Create(farg1, "FSUNContext Default", profiler)
fresult = swigc_FSUNContext_SetProfiler(ctx, profiler)
#endif
swig_result = fresult
end function
%}
