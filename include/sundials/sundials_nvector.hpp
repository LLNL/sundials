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
 * C++ view of SUNDIALS NVector
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_NVECTOR_HPP
#define _SUNDIALS_NVECTOR_HPP

#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_nvector.h>

namespace sundials {
namespace impl {

struct NVectorDeleter
{
  void operator()(N_Vector v)
  {
    if (v) N_VDestroy(v);
  }
};

using BaseNvector = BaseObject<_generic_N_Vector, _generic_N_Vector_Ops>;
using NVectorView = ClassView<N_Vector, NVectorDeleter>;

} // namespace impl
} // namespace sundials

#endif
