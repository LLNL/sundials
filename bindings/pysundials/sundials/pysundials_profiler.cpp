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
 * SUNDIALS SUNProfiler class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_profiler.hpp>

#include "sundials/sundials_types.h"
#include "sundials_profiler_impl.h"

namespace nb = nanobind;

using SUNProfilerView = sundials::experimental::SUNProfilerView;

void bind_sunprofiler(nb::module_& m)
{
#include "pysundials_profiler_generated.hpp"

  nb::class_<SUNProfiler_>(m, "SUNProfiler_");

  nb::class_<SUNProfilerView>(m, "SUNProfilerView")
    .def("get", nb::overload_cast<>(&SUNProfilerView::get, nb::const_),
         nb::rv_policy::reference)
    .def_static("Create", &SUNProfilerView::Create<SUNComm, const char*>);
}
