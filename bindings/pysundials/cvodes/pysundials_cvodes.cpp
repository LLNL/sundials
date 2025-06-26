/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

#include <sundials/sundials_core.hpp>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes.hpp>
#include <cvodes/cvodes_ls.h>

#include "sundials_adjointcheckpointscheme_impl.h"

// #include "pysundials_cvodes_usersupplied.hpp"

namespace nb = nanobind;

using namespace sundials::experimental;

// Forward declarations of functions defined in other translation units


void bind_cvodes(nb::module_& m)
{
#include "pysundials_cvodes_generated.hpp"

  nb::class_<CVodeView>(m, "CVodeView")
    .def_static("Create", &CVodeView::Create<void*>)
    .def("get", nb::overload_cast<>(&CVodeView::get, nb::const_),
         nb::rv_policy::reference);

}
