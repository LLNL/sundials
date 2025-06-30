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
 * C++ specific IDAS definitions.
 * ---------------------------------------------------------------------------*/

#ifndef _IDAS_HPP
#define _IDAS_HPP

extern "C" {
#include <idas/idas.h>
}

#include <sundials/sundials_base.hpp>

namespace sundials {
namespace experimental {

struct IDADeleter
{
  void operator()(void* v)
  {
    if (v) { IDAFree(&v); }
  }
};

class IDAView : public ClassView<void*, IDADeleter>
{
public:
  using ClassView<void*, IDADeleter>::ClassView;
  template<typename... Args>
  static IDAView Create(Args&&... args);
};

template<typename... Args>
IDAView IDAView::Create(Args&&... args)
{
  return IDAView(std::forward<Args>(args)...);
}

} // namespace experimental
} // namespace sundials

#endif
