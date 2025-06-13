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
 * C++ view of SUNDIALS SUNLinaerSolver
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_LINEARSOLVER_HPP
#define _SUNDIALS_LINEARSOLVER_HPP

#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_linearsolver.h>

namespace sundials {
namespace impl {
using BaseLinearSolver =
  BaseObject<_generic_SUNLinearSolver, _generic_SUNLinearSolver_Ops>;
} // namespace impl

namespace experimental {
struct SUNLinearSolverDeleter
{
  void operator()(SUNLinearSolver LS)
  {
    if (LS) { SUNLinSolFree(LS); }
  }
};

class SUNLinearSolverView : public ClassView<SUNLinearSolver, SUNLinearSolverDeleter>
{ 
public:
  using ClassView<SUNLinearSolver, SUNLinearSolverDeleter>::ClassView;
  template<typename... Args>
  static SUNLinearSolverView Create(Args&&... args);
};

template<typename... Args>
SUNLinearSolverView SUNLinearSolverView::Create(Args&&... args)
{
  return SUNLinearSolverView(std::forward<Args>(args)...);
}

} // namespace experimental
} // namespace sundials

#endif
