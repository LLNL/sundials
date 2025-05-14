/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * C++ view of SUNDIALS SUNAdaptController
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ADAPTCONTROLLER_HPP
#define _SUNDIALS_ADAPTCONTROLLER_HPP

#include <sundials/sundials_base.hpp>
#include <sundials/sundials_adaptcontroller.h>
#include "sundials_nvector.h"

namespace sundials {
namespace experimental {

struct SUNAdaptControllerDeleter
{
  void operator()(SUNAdaptController C)
  {
    if (C) { SUNAdaptController_Destroy(C); }
  }
};

using SUNAdaptControllerView = ClassView<SUNAdaptController, SUNAdaptControllerDeleter>;

} // namespace experimental
} // namespace sundials

#endif
