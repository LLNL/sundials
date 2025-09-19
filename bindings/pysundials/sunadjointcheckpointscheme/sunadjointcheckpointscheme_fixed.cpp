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
#include <nanobind/stl/tuple.h>

#include <sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h>
#include <sundials/sundials_core.hpp>

#include "sundials_adjointcheckpointscheme_impl.h"

namespace nb = nanobind;

void bind_sunadjointcheckpointscheme_fixed(nb::module_& m)
{
  m.def(
    "SUNAdjointCheckpointScheme_Create_Fixed",
    [](SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, suncountertype interval,
       suncountertype estimate, sunbooleantype keep, SUNContext sunctx)
    {
      SUNAdjointCheckpointScheme check_scheme = nullptr;

      int status = SUNAdjointCheckpointScheme_Create_Fixed(io_mode, mem_helper,
                                                           interval, estimate,
                                                           keep, sunctx,
                                                           &check_scheme);

      return std::make_tuple(status, check_scheme);
    },
    nb::rv_policy::reference);
}
