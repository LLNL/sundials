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
 * SUNDIALS SUNAdjointStepper class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_adjointstepper.hpp>

#include "sundials_adjointcheckpointscheme_impl.h"
#include "sundials_adjointstepper_impl.h"
#include "sundials_stepper_impl.h"

namespace nb = nanobind;

using SUNAdjointStepperView = sundials::experimental::SUNAdjointStepperView;

void bind_sunadjointstepper(nb::module_& m)
{
#include "pysundials_adjointstepper_generated.hpp"

  nb::class_<SUNAdjointStepperView>(m, "SUNAdjointStepperView")
    .def("get", nb::overload_cast<>(&SUNAdjointStepperView::get, nb::const_),
         nb::rv_policy::reference)
    .def_static("Create",
                &SUNAdjointStepperView::Create<
                  SUNStepper, sunbooleantype, SUNStepper, sunbooleantype, suncountertype,
                  sunrealtype, N_Vector, SUNAdjointCheckpointScheme, SUNContext>);
}
