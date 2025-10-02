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
 * C++ specific IDA definitions.
 * ---------------------------------------------------------------------------*/

#ifndef _IDA_HPP
#define _IDA_HPP

#include <sundials/sundials_classview.hpp>

#include <ida/ida.h>

namespace sundials {
namespace experimental {

struct IDADeleter
{
  void operator()(void* v) { IDAFree(&v); }
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
