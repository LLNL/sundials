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
 * SUNDIALS SUNAdaptController class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>

#include <sundials/sundials_adaptcontroller.hpp>

namespace nb = nanobind;

using namespace sundials::experimental;

void bind_sunadaptcontroller(nb::module_& m)
{
#include "pysundials_adaptcontroller_generated.hpp"

  nb::class_<SUNAdaptControllerView>(m, "SUNAdaptControllerView")
    .def_static("Create", &SUNAdaptControllerView::Create<SUNAdaptController>)
    .def("get", nb::overload_cast<>(&SUNAdaptControllerView::get, nb::const_),
         nb::rv_policy::reference);
}
