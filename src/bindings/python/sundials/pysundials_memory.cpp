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
 * SUNDIALS SUNMemoryHelper class. It contains hand-written code for 
 * functions that require special treatment, and includes the 
 * generated code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_memory.hpp>

namespace nb = nanobind;

using namespace sundials::experimental;

void bind_sunmemory(nb::module_& m)
{
#include "pysundials_memory_generated.hpp"

  nb::class_<SUNMemoryHelperView>(m, "SUNMemoryHelperView")
    .def_static("Create", &SUNMemoryHelperView::Create<SUNMemoryHelper>)
    .def("get", nb::overload_cast<>(&SUNMemoryHelperView::get, nb::const_),
         nb::rv_policy::reference);
}
