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
 * ----------------------------------------------------------------*/

#ifndef _SUNMATRIX_HPP
#define _SUNMATRIX_HPP

#include <sundials/core/sundials_base.hpp>
#include <sundials/sundials_matrix.h>

namespace sundials {
namespace impl {

using BaseMatrix = BaseObject<_generic_SUNMatrix, _generic_SUNMatrix_Ops>;

} // namespace impl
} // namespace sundials

#endif // _SUNMATRIX_HPP