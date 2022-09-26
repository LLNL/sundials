/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * C++ view of SUNDIALS SUNNonlinearSolver
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_NONLINEARSOLVER_HPP
#define _SUNDIALS_NONLINEARSOLVER_HPP

#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_nonlinearsolver.h>

namespace sundials {
namespace impl {

struct SUNNonlinearSolverDeleter
{
  void operator()(SUNNonlinearSolver NLS)
  {
    if (NLS) SUNNonlinSolFree(NLS);
  }
};

using BaseNonlinearSolver    = BaseObject<_generic_SUNNonlinearSolver, _generic_SUNNonlinearSolver_Ops>;
using SUNNonlinearSolverView = ClassView<SUNNonlinearSolver, SUNNonlinearSolverDeleter>;

} // namespace impl
} // namespace sundials

#endif
