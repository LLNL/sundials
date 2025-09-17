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

#ifndef _PYSUNDIALS_KINSOL_USERSUPPLIED_HPP
#define _PYSUNDIALS_KINSOL_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>

#include <nanobind/nanobind.h>

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_ls.h>

#include <sundials/sundials_core.hpp>

#include "pysundials_helpers.hpp"

namespace nb = nanobind;

///////////////////////////////////////////////////////////////////////////////
// KINSOL user-supplied function table
// Every package-level user-supplied function must be in this table.
// The user-supplied function table is passed to KINSOL as user_data.
///////////////////////////////////////////////////////////////////////////////

struct kinsol_user_supplied_fn_table
{
  // KINSOL user-supplied function pointers
  nb::object sysfn, dampingfn, depthfn;

  // KINSOL LS user-supplied function pointers
  nb::object lsjacfn, lsjactimesvecfn, lsjtvsysfn, lsprecsetupfn, lsprecsolvefn;
};

///////////////////////////////////////////////////////////////////////////////
// KINSOL user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline kinsol_user_supplied_fn_table* kinsol_user_supplied_fn_table_alloc()
{
  // We must use malloc since KINFree calls free
  auto fn_table = static_cast<kinsol_user_supplied_fn_table*>(
    std::malloc(sizeof(kinsol_user_supplied_fn_table)));

  // Zero out the memory
  std::memset(fn_table, 0, sizeof(kinsol_user_supplied_fn_table));

  return fn_table;
}

template<typename... Args>
inline int kinsol_sysfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<KINSysFn>, kinsol_user_supplied_fn_table,
    1>(&kinsol_user_supplied_fn_table::sysfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int kinsol_dampingfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<KINDampingFn>, kinsol_user_supplied_fn_table,
    2>(&kinsol_user_supplied_fn_table::dampingfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int kinsol_depthfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<KINDepthFn>, kinsol_user_supplied_fn_table,
    3>(&kinsol_user_supplied_fn_table::depthfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int kinsol_lsjacfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<KINLsJacFn>, kinsol_user_supplied_fn_table,
    3>(&kinsol_user_supplied_fn_table::lsjacfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int kinsol_lsprecsetupfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<KINLsPrecSetupFn>, kinsol_user_supplied_fn_table,
    1>(&kinsol_user_supplied_fn_table::lsprecsetupfn,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int kinsol_lsprecsolvefn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<KINLsPrecSolveFn>, kinsol_user_supplied_fn_table,
    1>(&kinsol_user_supplied_fn_table::lsprecsolvefn,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int kinsol_lsjactimesvecfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<KINLsJacTimesVecFn>, kinsol_user_supplied_fn_table,
    1>(&kinsol_user_supplied_fn_table::lsjactimesvecfn,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int kinsol_lsjtvsysfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<KINSysFn>, kinsol_user_supplied_fn_table,
    1>(&kinsol_user_supplied_fn_table::lsjtvsysfn, std::forward<Args>(args)...);
}

#endif // _PYSUNDIALS_KINSOL_USERSUPPLIED_HPP
