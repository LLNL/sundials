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
 * C++ specific ARKODE definitions.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ARKODE_MRISTEP_HPP
#define _SUNDIALS_ARKODE_MRISTEP_HPP

#include <arkode/arkode_mristep.h>
#include <sundials/sundials_classview.hpp>

namespace sundials {
namespace experimental {

struct MRIStepCouplingDeleter
{
  void operator()(MRIStepCoupling t)
  {
    if (t) { MRIStepCoupling_Free(t); }
  }
};

class MRIStepCouplingView
  : public ClassView<MRIStepCoupling, MRIStepCouplingDeleter>
{
public:
  using ClassView<MRIStepCoupling, MRIStepCouplingDeleter>::ClassView;
  template<typename... Args>
  static MRIStepCouplingView Create(Args&&... args);
};

template<typename... Args>
MRIStepCouplingView MRIStepCouplingView::Create(Args&&... args)
{
  MRIStepCoupling table = MRIStepCoupling_Create(std::forward<Args>(args)...);
  return MRIStepCouplingView(table);
}

struct MRIStepInnerStepperDeleter
{
  void operator()(MRIStepInnerStepper s)
  {
    if (s) { MRIStepInnerStepper_Free(&s); }
  }
};

class MRIStepInnerStepperView
  : public ClassView<MRIStepInnerStepper, MRIStepInnerStepperDeleter>
{
public:
  using ClassView<MRIStepInnerStepper, MRIStepInnerStepperDeleter>::ClassView;

  template<typename... Args>
  static MRIStepInnerStepperView Create(Args&&... args);
};

template<typename... Args>
MRIStepInnerStepperView MRIStepInnerStepperView::Create(Args&&... args)
{
  return MRIStepInnerStepperView(std::forward<Args>(args)...);
}

} // namespace experimental
} // namespace sundials

#endif
