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
 * SUNDIALS SUNNonlinearSolver class. It contains hand-written code 
 * for functions that require special treatment, and includes the
 * generated code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_nonlinearsolver.hpp>

namespace nb = nanobind;

void bind_sunnonlinearsolver(nb::module_& m)
{
#include "sundials_nonlinearsolver_generated.hpp"

  nb::class_<sundials::experimental::SUNNonlinearSolverView>(m,
                                                             "SUNNonlinearSolverView")
    .def(nb::init<>())
    .def(nb::init<_generic_SUNNonlinearSolver*>())
    .def("get",
         nb::overload_cast<>(&sundials::experimental::SUNNonlinearSolverView::get,
                             nb::const_),
         nb::rv_policy::reference);
}
