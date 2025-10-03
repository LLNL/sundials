/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * -----------------------------------------------------------------------------
 * C++ view of SUNDIALS SUNAdjointStepper
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ADJOINTSTEPPER_HPP
#define _SUNDIALS_ADJOINTSTEPPER_HPP

#include <utility>

#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_classview.hpp>

namespace sundials {
namespace experimental {

struct SUNAdjointStepperDeleter
{
  void operator()(SUNAdjointStepper self) { SUNAdjointStepper_Destroy(&self); }
};

class SUNAdjointStepperView
  : public ClassView<SUNAdjointStepper, SUNAdjointStepperDeleter>
{
public:
  using ClassView<SUNAdjointStepper, SUNAdjointStepperDeleter>::ClassView;

  template<typename... Args>
  static SUNAdjointStepperView Create(Args&&... args)
  {
    SUNAdjointStepper stepper;
    SUNAdjointStepper_Create(std::forward<Args>(args)..., &stepper);
    return SUNAdjointStepperView(stepper);
  }
};

template<>
SUNAdjointStepperView SUNAdjointStepperView::Create(SUNAdjointStepper& stepper)
{
  return SUNAdjointStepperView(std::forward<SUNAdjointStepper>(stepper));
}

template<>
SUNAdjointStepperView SUNAdjointStepperView::Create(SUNAdjointStepper&& stepper)
{
  return SUNAdjointStepperView(std::forward<SUNAdjointStepper>(stepper));
}

} // namespace experimental
} // namespace sundials

#endif
