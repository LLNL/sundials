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

template<typename FnType, typename FnTableType, std::size_t UserDataArg, typename... Args>
int user_supplied_fn_caller(nb::object FnTableType::*fn_member, Args... args)
{
  constexpr size_t N            = sizeof...(Args);
  constexpr int user_data_index = N - UserDataArg;
  auto args_tuple               = std::tuple<Args...>(args...);

  // Extract user_data from the specified position
  void* user_data = std::get<user_data_index>(args_tuple);

  // Cast user_data to FnTableType*
  auto fn_table = static_cast<FnTableType*>(user_data);
  auto fn       = nb::cast<std::function<FnType>>(fn_table->*fn_member);

  // Pass nullptr as user_data since we do not want the user to mess with user_data (which holds our function table)
  std::get<user_data_index>(args_tuple) = nullptr;
  return std::apply([&](auto&&... call_args)
                    { return fn(call_args...); }, args_tuple);
}

} // namespace pysundials

#endif
