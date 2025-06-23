/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#ifndef _PYSUNDIALS_HELPERS_HPP
#define _PYSUNDIALS_HELPERS_HPP

#include <nanobind/nanobind.h>

namespace nb = nanobind;

namespace pysundials {

template<typename FnType, typename CallbackType>
int user_supplied_fn_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                             void* user_data,
                             nb::object CallbackType::* fn_member)
{
  auto fn_table = static_cast<CallbackType*>(user_data);
  auto fn       = nb::cast<std::function<FnType>>(fn_table->*fn_member);
  return fn(t, y, ydot, user_data);
}

} // namespace pysundials

#endif
