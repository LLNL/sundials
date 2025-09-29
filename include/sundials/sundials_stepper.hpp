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
 * C++ view of SUNDIALS SUNStepper
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_STEPPER_HPP
#define _SUNDIALS_STEPPER_HPP

#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_stepper.h>

namespace sundials {

namespace experimental {

struct SUNStepperDeleter
{
  void operator()(SUNStepper self)
  {
    if (self) { SUNStepper_Destroy(&self); }
  }
};

class SUNStepperView : public ClassView<SUNStepper, SUNStepperDeleter>
{
public:
  using ClassView<SUNStepper, SUNStepperDeleter>::ClassView;

  template<typename... Args>
  static SUNStepperView Create(Args&&... args);
};

template<>
SUNStepperView SUNStepperView::Create(SUNStepper&& stepper)
{
  return SUNStepperView(std::forward<SUNStepper>(stepper));
}

template<typename... Args>
SUNStepperView SUNStepperView::Create(Args&&... args)
{
  SUNStepper stepper;
  SUNStepper_Create(std::forward<Args>(args)..., &stepper);
  return SUNStepperView(std::move(stepper));
}

} // namespace experimental
} // namespace sundials

#endif
