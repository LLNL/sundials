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

/// This function will call a user-supplied Python function through C++ side wrappers
/// \tparam FnType is the function signature, e.g., std::remove_pointer_t<CVRhsFn>
/// \tparam FnTableType is the struct function table that holds the user-supplied Python functions as std::function
/// \tparam UserDataArg is the index of the void* user_data argument of the C function. We are counting from the last arg to the first arg, so if user_data is the last arg then this should be 1.
/// \tparam Args is the template parameter pack that contains all of the types of the function arguments to the C function
///
/// \param fn_member is the name of the function in the FnTableType to call
/// \param args is the arguments to the C function, which will be forwarded to the user-supplied Python function, except user_data, which is intercepted and passed as a nullptr.
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
  return std::apply([&](auto&&... call_args) { return fn(call_args...); },
                    args_tuple);
}

/// This function will call a user-supplied Python function through C++ side wrappers
/// \tparam FnType is the function signature, e.g., std::remove_pointer_t<CVRhsFn>
/// \tparam FnTableType is the struct function table that holds the user-supplied Python functions as std::function
/// \tparam Args is the template parameter pack that contains all of the types of the function arguments to the C function
///
/// \param fn_member is the name of the function in the FnTableType to call
/// \param args is the arguments to the C function, which will be forwarded to the user-supplied Python function.
template<typename FnType, typename FnTableType, typename... Args>
int user_supplied_fn_caller(nb::object FnTableType::*fn_member, void* user_data,
                            Args... args)
{
  auto args_tuple = std::tuple<Args...>(args...);

  // Cast user_data to FnTableType*
  auto fn_table = static_cast<FnTableType*>(user_data);
  auto fn       = nb::cast<std::function<FnType>>(fn_table->*fn_member);

  return std::apply([&](auto&&... call_args) { return fn(call_args...); },
                    args_tuple);
}

/// This function will call a user-supplied Python function through C++ side wrappers
/// \tparam FnType is the function signature, e.g., std::remove_pointer_t<CVRhsFn>
/// \tparam FnTableType is the struct function table that holds the user-supplied Python functions as std::function
/// \tparam Args is the template parameter pack that contains all of the types of the function arguments to the C function
///
/// \param fn_member is the name of the function in the FnTableType to call
/// \param args is the arguments to the C function, which will be forwarded to the user-supplied Python function.
template<typename FnType, typename FnTableType, typename T, typename... Args>
int user_supplied_fn_caller(nb::object FnTableType::*fn_member, Args... args)
{
  auto args_tuple = std::tuple<Args...>(args...);

  // Cast object->python to FnTableType*
  auto object   = static_cast<T>(std::get<0>(args_tuple));
  auto fn_table = static_cast<FnTableType*>(object->python);
  auto fn       = nb::cast<std::function<FnType>>(fn_table->*fn_member);

  return std::apply([&](auto&&... call_args) { return fn(call_args...); },
                    args_tuple);
}

// TODO(CJB): implement custom exception type so all messages follow a specific format

} // namespace pysundials

#endif
