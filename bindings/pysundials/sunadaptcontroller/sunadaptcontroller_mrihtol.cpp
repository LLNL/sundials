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
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <sunadaptcontroller/sunadaptcontroller_mrihtol.h>

namespace nb = nanobind;

namespace pysundials {

void bind_sunadaptcontroller_mrihtol(nb::module_& m)
{
  m.def("SUNAdaptController_MRIHTol", &SUNAdaptController_MRIHTol,
        nb::rv_policy::reference);
}

} // namespace pysundials
