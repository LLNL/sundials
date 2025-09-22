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
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

namespace nb = nanobind;

void bind_sunnonlinsol_fixedpoint(nb::module_& m)
{
  m.def("SUNNonlinSol_FixedPoint", &SUNNonlinSol_FixedPoint,
        nb::rv_policy::reference);
}
