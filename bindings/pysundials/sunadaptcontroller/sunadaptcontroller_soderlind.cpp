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
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

namespace nb = nanobind;

void bind_sunadaptcontroller_soderlind(nb::module_& m)
{
  m.def("SUNAdaptController_Soderlind", &SUNAdaptController_Soderlind,
        nb::rv_policy::reference);
}
