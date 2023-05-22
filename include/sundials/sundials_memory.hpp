/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * C++ view of SUNDIALS SUNMemoryHelper
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_MEMORY_HPP
#define _SUNDIALS_MEMORY_HPP

#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_memory.h>

namespace sundials {
namespace impl {
using BaseMemoryHelper = BaseObject<_SUNMemoryHelper, _SUNMemoryHelper_Ops>;
} // namespace impl

namespace experimental {
struct SUNMemoryHelperDeleter
{
  void operator()(SUNMemoryHelper helper)
  {
    if (helper) SUNMemoryHelper_Destroy(helper);
  }
};
using SUNMemoryHelperView = ClassView<SUNMemoryHelper, SUNMemoryHelperDeleter>;
} // namespace experimental
} // namespace sundials

#endif
