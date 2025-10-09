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
 * SUNDIALS N_Vector class. It contains hand-written code for functions
 * that require special treatment, and includes the generated code
 * produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/tuple.h>

#include <sundials/sundials_domeigestimator.hpp>

#include "sundials_domeigestimator_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace pysundials {

void bind_sundomeigestimator(nb::module_& m)
{
#include "sundials_domeigestimator_generated.hpp"

  nb::class_<SUNDomEigEstimatorView>(m, "SUNDomEigEstimatorView")
    .def_static("Create", &SUNDomEigEstimatorView::Create<SUNDomEigEstimator>)
    .def("get", nb::overload_cast<>(&SUNDomEigEstimatorView::get, nb::const_),
         nb::rv_policy::reference);

  m.def(
    "SUNDomEigEstimator_SetATimes",
    [](SUNDomEigEstimator dee,
       std::function<std::remove_pointer_t<SUNATimesFn>> ATimes) -> SUNErrCode
    {
      if (!dee->python)
      {
        dee->python = SUNDomEigEstimatorFunctionTable_Alloc();
      }

      auto fntable = static_cast<SUNDomEigEstimatorFunctionTable*>(dee->python);

      fntable->atimes = nb::cast(ATimes);

      if (ATimes)
      {
        return SUNDomEigEstimator_SetATimes(dee, fntable,
                                            sundomeigestimator_atimes_wrapper);
      }
      else { return SUNDomEigEstimator_SetATimes(dee, fntable, nullptr); }
    },
    nb::arg("DEE"), nb::arg("ATimes").none());
}

} // namespace pysundials
