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
 * C++ view of ARKODE memory
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ARKODE_HPP
#define _SUNDIALS_ARKODE_HPP

#include <sundials/sundials_base.hpp>
#include <arkode/arkode.h>
#include <arkode/arkode_erkstep.h>

namespace sundials {
namespace experimental {

struct ARKodeDeleter
{
  void operator()(void* v)
  {
    if (v) { ARKodeFree(&v); }
  }
};

class ARKodeView : public ClassView<void*, ARKodeDeleter>
{ 
public:
  using ClassView<void*, ARKodeDeleter>::ClassView;
  template<typename... Args>
  static ARKodeView Create(Args&&... args);
};

template<typename... Args>
ARKodeView ARKodeView::Create(Args&&... args)
{
  return ARKodeView(std::forward<Args>(args)...);
}

} // namespace experimental
} // namespace sundials

#endif
