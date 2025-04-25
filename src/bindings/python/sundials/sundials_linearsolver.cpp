/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file is the entrypoint for the Python binding code for the
 * SUNDIALS N_Vector class. It contains hand-written code for functions
 * that require special treatment, and includes the generated code
 * produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_linearsolver.hpp>

namespace nb = nanobind;

void bind_linearsolver(nb::module_& m)
{
  nb::class_<sundials::experimental::SUNLinearSolverView>(m,
                                                          "SUNLinearSolverView")
    .def(nb::init<>())
    .def(nb::init<_generic_SUNLinearSolver*>())
    .def("get",
         nb::overload_cast<>(&sundials::experimental::SUNLinearSolverView::get,
                             nb::const_),
         nb::rv_policy::reference);

  nb::class_<_generic_SUNLinearSolver_Ops>(m,
                                            "__generic_SUNLinearSolver_Ops", "Structure containing function pointers to linear solver operations")
    .def(nb::init<>());

  nb::class_<_generic_SUNLinearSolver>(m, "_generic_SUNLinearSolver")
    .def(
      "__init__",
      [](_generic_SUNLinearSolver* self, SUNLinearSolver_Ops ops = SUNLinearSolver_Ops(),
         SUNContext sunctx = SUNContext())
      {
        new (self) _generic_SUNLinearSolver(); // placement new
        auto r    = self;
        r->ops    = ops;
        r->sunctx = sunctx;
      },
      nb::arg("ops") = N_Vector_Ops(), nb::arg("sunctx") = SUNContext())
    .def_rw("content", &_generic_SUNLinearSolver::content, "")
    .def_rw("ops", &_generic_SUNLinearSolver::ops, "")
    .def_rw("sunctx", &_generic_SUNLinearSolver::sunctx, "");

  m.def("SUNLinSolNewEmpty", &SUNLinSolNewEmpty,
        nb::rv_policy::reference);
  m.def("SUNLinSolFreeEmpty", &SUNLinSolFreeEmpty);
  m.def("SUNLinSolGetType", &SUNLinSolGetType);
  m.def("SUNLinSolGetID", &SUNLinSolGetID);
  m.def("SUNLinSolSetATimes", &SUNLinSolSetATimes);
  m.def("SUNLinSolSetPreconditioner", &SUNLinSolSetPreconditioner);
  m.def("SUNLinSolSetScalingVectors", &SUNLinSolSetScalingVectors);
  m.def("SUNLinSolSetZeroGuess", &SUNLinSolSetZeroGuess);
  m.def("SUNLinSolInitialize", &SUNLinSolInitialize);
  m.def("SUNLinSolSetup", &SUNLinSolSetup);
  m.def("SUNLinSolSolve", &SUNLinSolSolve);
  m.def("SUNLinSolNumIters", &SUNLinSolNumIters);
  m.def("SUNLinSolResNorm", &SUNLinSolResNorm);
  m.def("SUNLinSolResid", &SUNLinSolResid);
  m.def("SUNLinSolLastFlag", &SUNLinSolLastFlag);
}
