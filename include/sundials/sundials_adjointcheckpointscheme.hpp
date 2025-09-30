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
 * C++ view of SUNDIALS SUNAdjointCheckpointScheme
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ADJOINTCHECKPOINTSCHEME_HPP
#define _SUNDIALS_ADJOINTCHECKPOINTSCHEME_HPP

#include <utility>

#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_classview.hpp>

namespace sundials {
namespace experimental {

struct SUNAdjointCheckpointSchemeDeleter
{
  void operator()(SUNAdjointCheckpointScheme self)
  {
    SUNAdjointCheckpointScheme_Destroy(&self);
  }
};

class SUNAdjointCheckpointSchemeView
  : public ClassView<SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeDeleter>
{
public:
  using ClassView<SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeDeleter>::ClassView;

  template<typename... Args>
  static SUNAdjointCheckpointSchemeView Create(Args&&... args)
  {
    return SUNAdjointCheckpointSchemeView(std::forward<Args>(args)...);
  }
};

} // namespace experimental
} // namespace sundials

#endif
