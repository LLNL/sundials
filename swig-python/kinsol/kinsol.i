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

%module(directors="1") kinsol

%include "stdint.i"

// ----------------------------------
// numpy swig setup stuff.
// This must only happen in one file.
// ----------------------------------

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}

// -------------------------------
// Bring in shared sundials stuff.
// -------------------------------

%include "../sundials/sundials.i"

// ---------------------
// KINSOL specific stuff
// --------------------- 

// KINInit cannot be called from Python.
// Instead, users should call KINInitPy.
%ignore KINInit;

// We hijack KINSetUserData to pass out director class
// objects. So, hide the function from users.
%ignore KINSetUserData;

%{
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_bbdpre.h"
#include "kinsol/kinsol_ls.h"
#include "callbacks.h"
%}

// KINSysPyFn is a 'director' class
%feature("director") KINSysPyFn;

// Apply typemap for Get functions that use an argout variable
%typemap(in, numinputs=0) realtype* (realtype temp) {
  $1 = &temp;
}
%typemap(argout) realtype* {
  $result = SWIG_Python_AppendOutput($result, PyFloat_FromDouble(*$1));
}
%typemap(in, numinputs=0) long* (long temp) {
  $1 = &temp;
}
%typemap(argout) long* {
  $result = SWIG_Python_AppendOutput($result, PyLong_FromLong(*$1));
}

// Process definitions from these files
%include "kinsol/kinsol.h"
%include "kinsol/kinsol_bbdpre.h"
%include "kinsol/kinsol_ls.h"
%include "callbacks.h"

// Insert helper code for setting sysfn
%pythoncode
%{
def WrapPythonSysFn(user_sysfun):
  caller = KINSysFnCaller()
  caller.setFn(KINSysPyFnPyChild(user_sysfun).__disown__())
  return caller

# inherits from the C++ KinSysPyFn class
class KINSysPyFnPyChild(KINSysPyFn):
  def __init__(self, user_sysfun):
    KINSysPyFn.__init__(self)
    self.user_sysfun = user_sysfun

  def actual_sysfun(self, y, g, user_data):
    self.user_sysfun(y, g, user_data)
    return 0
%}
