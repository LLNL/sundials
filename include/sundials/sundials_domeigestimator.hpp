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
 * C++ view of SUNDIALS SUNDomEigEstimator
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_DOMEIGESTIMATOR_HPP
#define _SUNDIALS_DOMEIGESTIMATOR_HPP

#include <sundials/sundials_base.hpp>
#include <sundials/sundials_domeigestimator.h>

namespace sundials {
namespace impl {
using BaseDomEigEstimator =
  BaseObject<SUNDomEigEstimator_, SUNDomEigEstimator_Ops_>;
} // namespace impl

namespace experimental {
struct SUNDomEigEstimatorDeleter
{
  void operator()(SUNDomEigEstimator DEE)
  {
    if (DEE) { SUNDomEigEstimator_Destroy(&DEE); }
  }
};

class SUNDomEigEstimatorView
  : public ClassView<SUNDomEigEstimator, SUNDomEigEstimatorDeleter>
{
public:
  using ClassView<SUNDomEigEstimator, SUNDomEigEstimatorDeleter>::ClassView;
  template<typename... Args>
  static SUNDomEigEstimatorView Create(Args&&... args);
};

template<typename... Args>
SUNDomEigEstimatorView SUNDomEigEstimatorView::Create(Args&&... args)
{
  return SUNDomEigEstimatorView(std::forward<Args>(args)...);
}

} // namespace experimental
} // namespace sundials

#endif
