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
 * C++ view of SUNDIALS NVector
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_NVECTOR_HPP
#define _SUNDIALS_NVECTOR_HPP

#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_nvector.h>

namespace sundials {
namespace impl {
using BaseNVector = BaseObject<_generic_N_Vector, _generic_N_Vector_Ops>;
} // namespace impl

namespace experimental {
struct NVectorDeleter
{
  void operator()(N_Vector v)
  {
    if (v) { N_VDestroy(v); }
  }
};

using NVectorView = ClassView<N_Vector, NVectorDeleter>;
} // namespace experimental
} // namespace sundials

#endif
