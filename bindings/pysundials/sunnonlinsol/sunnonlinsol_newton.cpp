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
#include <sunnonlinsol/sunnonlinsol_newton.h>

namespace nb = nanobind;

void bind_sunnonlinsol_newton(nb::module_& m)
{
  m.def("SUNNonlinSol_Newton", &SUNNonlinSol_Newton, nb::rv_policy::reference);
}
