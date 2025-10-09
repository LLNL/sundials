/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2025, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file is the entrypoint for the Python binding code for the
 * SUNDIALS SUNNonlinearSolver class. It contains hand-written code 
 * for functions that require special treatment, and includes the
 * generated code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>

#include <sundials/sundials_nonlinearsolver.hpp>

#include "sundials_nonlinearsolver_usersupplied.hpp"

namespace nb = nanobind;

using namespace sundials::experimental;

void bind_sunnonlinearsolver(nb::module_& m)
{
#include "sundials_nonlinearsolver_generated.hpp"

  nb::class_<SUNNonlinearSolverView>(m, "SUNNonlinearSolverView")
    .def_static("Create", &SUNNonlinearSolverView::Create<SUNNonlinearSolver>)
    .def("get", nb::overload_cast<>(&SUNNonlinearSolverView::get, nb::const_),
         nb::rv_policy::reference);

  m.def(
    "SUNNonlinSolSetup",
    [](SUNNonlinearSolver NLS, N_Vector y)
    {
      if (!NLS->python)
      {
        NLS->python = SUNNonlinearSolverFunctionTable_Alloc();
      }
      return SUNNonlinSolSetup(NLS, y, NLS->python);
    },
    nb::arg("NLS"), nb::arg("y"));

  m.def(
    "SUNNonlinSolSolve",
    [](SUNNonlinearSolver NLS, N_Vector y0, N_Vector y, N_Vector w,
       sunrealtype tol, sunbooleantype callLSetup)
    {
      if (!NLS->python)
      {
        NLS->python = SUNNonlinearSolverFunctionTable_Alloc();
      }
      return SUNNonlinSolSolve(NLS, y0, y, w, tol, callLSetup, NLS->python);
    },
    nb::arg("NLS"), nb::arg("y0"), nb::arg("y"), nb::arg("w"), nb::arg("tol"),
    nb::arg("callLSetup"));

  m.def(
    "SUNNonlinSolSetSysFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolSysFn>> SysFn) -> SUNErrCode
    {
      if (!NLS->python)
      {
        NLS->python = SUNNonlinearSolverFunctionTable_Alloc();
      }
      auto fntable = static_cast<SUNNonlinearSolverFunctionTable*>(NLS->python);
      fntable->sysfn = nb::cast(SysFn);
      if (SysFn)
      {
        return SUNNonlinSolSetSysFn(NLS, sunnonlinearsolver_sysfn_wrapper);
      }
      else { return SUNNonlinSolSetSysFn(NLS, nullptr); }
    },
    nb::arg("NLS"), nb::arg("SysFn").none());

  m.def(
    "SUNNonlinSolSetLSetupFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolLSetupFn>> SetupFn) -> SUNErrCode
    {
      if (!NLS->python)
      {
        NLS->python = SUNNonlinearSolverFunctionTable_Alloc();
      }
      auto fntable = static_cast<SUNNonlinearSolverFunctionTable*>(NLS->python);
      fntable->lsetupfn = nb::cast(SetupFn);
      if (SetupFn)
      {
        return SUNNonlinSolSetLSetupFn(NLS, sunnonlinearsolver_lsetupfn_wrapper);
      }
      else { return SUNNonlinSolSetLSetupFn(NLS, nullptr); }
    },
    nb::arg("NLS"), nb::arg("SetupFn").none());

  m.def(
    "SUNNonlinSolSetLSolveFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolLSolveFn>> SolveFn) -> SUNErrCode
    {
      if (!NLS->python)
      {
        NLS->python = SUNNonlinearSolverFunctionTable_Alloc();
      }
      auto fntable = static_cast<SUNNonlinearSolverFunctionTable*>(NLS->python);
      fntable->lsolvefn = nb::cast(SolveFn);
      if (SolveFn)
      {
        return SUNNonlinSolSetLSolveFn(NLS, sunnonlinearsolver_lsolvefn_wrapper);
      }
      else { return SUNNonlinSolSetLSolveFn(NLS, nullptr); }
    },
    nb::arg("NLS"), nb::arg("SolveFn").none());

  m.def(
    "SUNNonlinSolSetConvTestFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolConvTestFn>> CTestFn) -> SUNErrCode
    {
      if (!NLS->python)
      {
        NLS->python = SUNNonlinearSolverFunctionTable_Alloc();
      }
      auto fntable = static_cast<SUNNonlinearSolverFunctionTable*>(NLS->python);
      fntable->convtestfn = nb::cast(CTestFn);
      if (CTestFn)
      {
        return SUNNonlinSolSetConvTestFn(NLS,
                                         sunnonlinearsolver_convtestfn_wrapper,
                                         NLS->python);
      }
      else { return SUNNonlinSolSetConvTestFn(NLS, nullptr, nullptr); }
    },
    nb::arg("NLS"), nb::arg("CTestFn").none());
}
