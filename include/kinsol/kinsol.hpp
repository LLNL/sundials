/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * C++ specific KINSOL definitions.
 * ---------------------------------------------------------------------------*/

#ifndef _KINSOL_HPP
#define _KINSOL_HPP

extern "C" {
#include <kinsol/kinsol.h>
}

#include <sundials/sundials_base.hpp>

namespace sundials {
namespace experimental {

struct KINSOLDeleter {
  void operator()(void* v) {
    if (v) { KINFree(&v); }
  }
};

class KINSOLView : public ClassView<void*, KINSOLDeleter> {
public:
  using ClassView<void*, KINSOLDeleter>::ClassView;
  template<typename... Args>
  static KINSOLView Create(Args&&... args);
};

template<typename... Args>
KINSOLView KINSOLView::Create(Args&&... args) {
  return KINSOLView(std::forward<Args>(args)...);
}

} // namespace experimental
} // namespace sundials

#endif
