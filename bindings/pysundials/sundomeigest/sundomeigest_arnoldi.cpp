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
#include <sundomeigest/sundomeigest_arnoldi.h>

namespace nb = nanobind;

void bind_sundomeigest_arnoldi(nb::module_& m)
{
  m.def("SUNDomEigEstimator_Arnoldi", &SUNDomEigEstimator_Arnoldi,
        nb::rv_policy::reference);
}
