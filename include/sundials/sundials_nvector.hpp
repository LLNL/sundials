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
 * C++ view of SUNDIALS NVector
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_NVECTOR_HPP
#define _SUNDIALS_NVECTOR_HPP

#include <utility>

#include <sundials/sundials_base.hpp>
#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_nvector.h>

namespace sundials {
namespace impl {
using BaseNVector = BaseObject<_generic_N_Vector, _generic_N_Vector_Ops>;
} // namespace impl

namespace experimental {

struct NVectorDeleter
{
  void operator()(N_Vector v) { N_VDestroy(v); }
};

class NVectorView : public ClassView<N_Vector, NVectorDeleter>
{
public:
  using ClassView<N_Vector, NVectorDeleter>::ClassView;

  template<typename... Args>
  static NVectorView Create(Args&&... args)
  {
    return NVectorView(std::forward<Args>(args)...);
  }
};

} // namespace experimental
} // namespace sundials

#endif
