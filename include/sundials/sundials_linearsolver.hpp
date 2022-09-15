/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_LINEARSOLVER_HPP
#define _SUNDIALS_LINEARSOLVER_HPP

#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_linearsolver.h>

namespace sundials {
namespace impl {

struct SUNLinearSolverDeleter
{
  constexpr void operator()(SUNLinearSolver LS)
  {
    if (LS) SUNLinSolFree(LS);
  }
};

using BaseLinearSolver    = BaseObject<_generic_SUNLinearSolver, _generic_SUNLinearSolver_Ops>;
using SUNLinearSolverView = ClassView<SUNLinearSolver, SUNLinearSolverDeleter>;

} // namespace impl
} // namespace sundials

#endif
