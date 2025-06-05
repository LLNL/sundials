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
 * SUNDIALS SUNAdjointStepper class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>

#include <sundials/sundials_adjointstepper.hpp>

#include "sundials_adjointstepper_impl.h"
#include "sundials_adjointcheckpointscheme_impl.h"
#include "sundials_stepper_impl.h"

namespace nb = nanobind;

using SUNAdjointStepperView = sundials::experimental::SUNAdjointStepperView;

void bind_sunadjointstepper(nb::module_& m)
{
#include "sundials_adjointstepper_generated.hpp"

  nb::class_<SUNAdjointStepperView>(m, "SUNAdjointStepperView")
    .def(nb::init<>())
    .def(nb::init<SUNAdjointStepper_*>())
    .def("get",
         nb::overload_cast<>(&SUNAdjointStepperView::get, nb::const_),
         nb::rv_policy::reference)
    .def_static("make_view", &SUNAdjointStepperView::make_view<SUNStepper, sunbooleantype, SUNStepper, sunbooleantype, suncountertype, sunrealtype, N_Vector, SUNAdjointCheckpointScheme, SUNContext>);
}
