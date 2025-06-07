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
 * SUNDIALS SUNStepper class. It contains hand-written code 
 * for functions that require special treatment, and includes the
 * generated code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_stepper.hpp>

#include "sundials_stepper_impl.h"

namespace nb = nanobind;

using SUNStepperView = sundials::experimental::SUNStepperView;

void bind_sunstepper(nb::module_& m)
{
#include "sundials_stepper_generated.hpp"

  nb::class_<SUNStepper_>(m, "SUNStepper_");

  nb::class_<SUNStepperView>(m, "SUNStepperView")
    .def(nb::init<SUNStepper>())
    .def("get",
         nb::overload_cast<>(&SUNStepperView::get,
                             nb::const_),
         nb::rv_policy::reference)
    .def_static("Create", SUNStepperView::Create<SUNContext>);
}
