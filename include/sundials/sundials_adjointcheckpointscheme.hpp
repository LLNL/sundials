/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
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

#include <sundials/sundials_base.hpp>
#include <sundials/sundials_adjointcheckpointscheme.h>

namespace sundials {
namespace experimental {

struct SUNAdjointCheckpointSchemeDeleter
{
  void operator()(SUNAdjointCheckpointScheme self)
  {
    if (self) { SUNAdjointCheckpointScheme_Destroy(&self); }
  }
};

using SUNAdjointCheckpointSchemeView = ClassView<SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeDeleter>;

} // namespace experimental
} // namespace sundials

#endif
