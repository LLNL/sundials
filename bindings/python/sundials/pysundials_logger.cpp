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
 * SUNDIALS SUNLogger class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>

#include <sundials/sundials_logger.hpp>

#include "sundials/sundials_types.h"
#include "sundials_logger_impl.h"

namespace nb = nanobind;

using SUNLoggerView = sundials::experimental::SUNLoggerView;

void bind_sunlogger(nb::module_& m)
{
#include "pysundials_logger_generated.hpp"
  nb::class_<SUNLogger_>(m, "SUNLogger_");

  nb::class_<SUNLoggerView>(m, "SUNLoggerView")
    .def("get", nb::overload_cast<>(&SUNLoggerView::get, nb::const_),
         nb::rv_policy::reference)
    .def_static("Create", &SUNLoggerView::Create<SUNComm>)
    .def_static("Create", &SUNLoggerView::Create<SUNComm, int>);
}
