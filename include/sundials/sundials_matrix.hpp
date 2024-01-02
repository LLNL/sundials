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
 * C++ view of SUNDIALS SUNMatrix
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_MATRIX_HPP
#define _SUNDIALS_MATRIX_HPP

#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_matrix.h>

namespace sundials {
namespace impl {
using BaseMatrix = BaseObject<_generic_SUNMatrix, _generic_SUNMatrix_Ops>;
} // namespace impl

namespace experimental {
struct SUNMatrixDeleter
{
  void operator()(SUNMatrix A)
  {
    if (A) { SUNMatDestroy(A); }
  }
};

using SUNMatrixView = ClassView<SUNMatrix, SUNMatrixDeleter>;
} // namespace experimental
} // namespace sundials

#endif
