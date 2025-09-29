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
 * C++ specific CVODES definitions.
 * ---------------------------------------------------------------------------*/

#ifndef _CVODES_HPP
#define _CVODES_HPP

extern "C" {
#include <cvodes/cvodes.h>
}

#include <sundials/sundials_base.hpp>

namespace sundials {
namespace experimental {

struct CVodeDeleter
{
  void operator()(void* v)
  {
    if (v) { CVodeFree(&v); }
  }
};

class CVodeView : public ClassView<void*, CVodeDeleter>
{
public:
  using ClassView<void*, CVodeDeleter>::ClassView;
  template<typename... Args>
  static CVodeView Create(Args&&... args);
};

template<typename... Args>
CVodeView CVodeView::Create(Args&&... args)
{
  return CVodeView(std::forward<Args>(args)...);
}

} // namespace experimental
} // namespace sundials

#endif
