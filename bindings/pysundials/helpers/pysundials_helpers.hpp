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

// Helper for tuple slicing
// The template lambda used requires C++20
template<std::size_t Start, std::size_t End, typename Tuple>
constexpr auto tuple_slice(const Tuple& t)
{
  return [&]<std::size_t... I>(std::index_sequence<I...>)
  {
    return std::make_tuple(std::get<Start + I>(t)...);
  }(std::make_index_sequence<End - Start>{});
}

template<typename FnType, typename FnTableType, typename... Args>
int user_supplied_fn_caller(nb::object FnTableType::* fn_member, Args... args)
{
  constexpr size_t N = sizeof...(Args);
  auto args_tuple    = std::tuple<Args...>(args...);

  // Extract user_data (assumed to be the last argument)
  void* user_data = std::get<N - 1>(args_tuple);

  // Cast user_data to FnTableType*
  auto fn_table = static_cast<FnTableType*>(user_data);
  auto fn       = nb::cast<std::function<FnType>>(fn_table->*fn_member);

  // Call the function with all arguments except the last (user_data), passing nullptr as user_data
  // to the user-supplied Python function. This ensures that the user cannot mess with user_data,
  // which we are using to store the function table.
  auto sliced_args = tuple_slice<0, N - 1>(args_tuple);
  return std::apply([&](auto&&... call_args)
                    { return fn(call_args..., nullptr); }, sliced_args);
}

} // namespace pysundials

#endif
