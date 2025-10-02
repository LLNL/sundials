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
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ARKODE_SPRKSTEP_HPP
#define _SUNDIALS_ARKODE_SPRKSTEP_HPP

#include <utility>

#include <arkode/arkode_sprk.h>
#include <sundials/sundials_classview.hpp>

namespace sundials {
namespace experimental {

struct ARKodeSPRKTableDeleter
{
  void operator()(ARKodeSPRKTable t) { ARKodeSPRKTable_Free(t); }
};

class ARKodeSPRKTableView
  : public ClassView<ARKodeSPRKTable, ARKodeSPRKTableDeleter>
{
public:
  using ClassView<ARKodeSPRKTable, ARKodeSPRKTableDeleter>::ClassView;
  
  template<typename... Args>
  static ARKodeSPRKTableView Create(Args&&... args)
  {
    return ARKodeSPRKTableView(std::forward<Args>(args)...);
  }
};

} // namespace experimental
} // namespace sundials

#endif
