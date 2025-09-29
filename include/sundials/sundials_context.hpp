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
 * C++ interface to the SUNDIALS context object
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_CONTEXT_HPP
#define _SUNDIALS_CONTEXT_HPP

#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_context.h>
#include <sundials/sundials_convertibleto.hpp>
#include <sundials/sundials_types.h>

namespace sundials {

struct SUNContextDeleter
{
  void operator()(SUNContext sunctx) { SUNContext_Free(&sunctx); }
};

class SUNContextView
  : public sundials::experimental::ClassView<SUNContext, SUNContextDeleter>
{
public:
  using ClassView<SUNContext, SUNContextDeleter>::ClassView;

  SUNContextView(SUNComm comm = SUN_COMM_NULL)
  {
    SUNContext_Create(comm, &object_);
  }

  template<typename... Args>
  static SUNContextView Create(Args&&... args);
};

template<typename... Args>
SUNContextView SUNContextView::Create(Args&&... args)
{
  return SUNContextView(std::forward<Args>(args)...);
}

// We provide this for backwards compatibility
using Context = SUNContextView;

} // namespace sundials

#endif // _SUNDIALS_CONTEXT_HPP
