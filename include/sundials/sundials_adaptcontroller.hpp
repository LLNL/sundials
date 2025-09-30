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
 * C++ view of SUNDIALS SUNAdaptController
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ADAPTCONTROLLER_HPP
#define _SUNDIALS_ADAPTCONTROLLER_HPP

#include <sundials/sundials_adaptcontroller.h>
#include <sundials/sundials_classview.hpp>
#include "sundials_nvector.h"

namespace sundials {
namespace experimental {

struct SUNAdaptControllerDeleter
{
  void operator()(SUNAdaptController C) { SUNAdaptController_Destroy(C); }
};

class SUNAdaptControllerView
  : public ClassView<SUNAdaptController, SUNAdaptControllerDeleter>
{
public:
  using ClassView<SUNAdaptController, SUNAdaptControllerDeleter>::ClassView;
  template<typename... Args>
  static SUNAdaptControllerView Create(Args&&... args);
};

template<typename... Args>
SUNAdaptControllerView SUNAdaptControllerView::Create(Args&&... args)
{
  return SUNAdaptControllerView(std::forward<Args>(args)...);
}

} // namespace experimental
} // namespace sundials

#endif
