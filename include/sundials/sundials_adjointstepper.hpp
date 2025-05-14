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
 * C++ view of SUNDIALS SUNAdjointStepper
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ADJOINTSTEPPER_HPP
#define _SUNDIALS_ADJOINTSTEPPER_HPP

#include <sundials/sundials_base.hpp>
#include <sundials/sundials_adjointstepper.h>

namespace sundials {
namespace experimental {

struct SUNAdjointStepperDeleter
{
  void operator()(SUNAdjointStepper self)
  {
    if (self) { SUNAdjointStepper_Destroy(&self); }
  }
};

using SUNAdjointStepperView = ClassView<SUNAdjointStepper, SUNAdjointStepperDeleter>;

} // namespace experimental
} // namespace sundials

#endif
