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

#ifndef _PYSUNDIALS_IDAS_USERSUPPLIED_HPP
#define _PYSUNDIALS_IDAS_USERSUPPLIED_HPP

#include <idas/idas.h>
#include <idas/idas_ls.h>

#include <sundials/sundials_core.hpp>

#include "pysundials_helpers.hpp"

///////////////////////////////////////////////////////////////////////////////
// IDAS user-supplied function table
// Every integrator-level user-supplied function must be in this table.
// The user-supplied function table is passed to IDAS as user_data.
///////////////////////////////////////////////////////////////////////////////

struct idas_user_supplied_fn_table
{
  // user-supplied function pointers
  nb::object res, rootfn, ewtn, rwtn, resNLS;

  // idas_ls user-supplied function pointers
  nb::object lsjacfn, lsprecsetupfn, lsprecsolvefn, lsjactimessetupfn,
    lsjactimesvecfn, lsjacresfn;

  // idas quadrature user-supplied function pointers
  nanobind::object resQ, resQS;

  // idas FSA user-supplied function pointers
  nanobind::object resS;
};

///////////////////////////////////////////////////////////////////////////////
// IDAS user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline idas_user_supplied_fn_table* idas_user_supplied_fn_table_alloc()
{
  // We must use malloc since IDASFree calls free
  auto fn_table = static_cast<idas_user_supplied_fn_table*>(
    std::malloc(sizeof(idas_user_supplied_fn_table)));

  // Zero out the memory
  std::memset(fn_table, 0, sizeof(idas_user_supplied_fn_table));

  return fn_table;
}

template<typename... Args>
inline int idas_res_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<std::remove_pointer_t<IDAResFn>,
                                             idas_user_supplied_fn_table,
                                             1>(&idas_user_supplied_fn_table::res,
                                                std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_rootfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDARootFn>, idas_user_supplied_fn_table,
    1>(&idas_user_supplied_fn_table::rootfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_ewtfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<std::remove_pointer_t<IDAEwtFn>,
                                             idas_user_supplied_fn_table,
                                             1>(&idas_user_supplied_fn_table::ewtn,
                                                std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_nlsresfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDAResFn>, idas_user_supplied_fn_table,
    1>(&idas_user_supplied_fn_table::resNLS, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjacfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacFn>, idas_user_supplied_fn_table,
    4>(&idas_user_supplied_fn_table::lsjacfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsprecsetupfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSetupFn>, idas_user_supplied_fn_table,
    1>(&idas_user_supplied_fn_table::lsprecsetupfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsprecsolvefn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSolveFn>, idas_user_supplied_fn_table,
    1>(&idas_user_supplied_fn_table::lsprecsolvefn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjactimessetupfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesSetupFn>, idas_user_supplied_fn_table,
    1>(&idas_user_supplied_fn_table::lsjactimessetupfn,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjactimesvecfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesVecFn>, idas_user_supplied_fn_table,
    3>(&idas_user_supplied_fn_table::lsjactimesvecfn,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjacresfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDAResFn>, idas_user_supplied_fn_table,
    1>(&idas_user_supplied_fn_table::lsjacresfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_resQ_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDAQuadRhsFn>, idas_user_supplied_fn_table,
    1>(&idas_user_supplied_fn_table::resQ, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_resQS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDAQuadSensRhsFn>, idas_user_supplied_fn_table,
    4>(&idas_user_supplied_fn_table::resQS, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_resS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDASensResFn>, idas_user_supplied_fn_table,
    4>(&idas_user_supplied_fn_table::resS, std::forward<Args>(args)...);
}

///////////////////////////////////////////////////////////////////////////////
// IDAS Adjoint user-supplied functions
///////////////////////////////////////////////////////////////////////////////

struct idasa_user_supplied_fn_table
{
  // idas adjoint user-supplied function pointers
  nb::object resB, resQB, resQBs;

  // idas_ls adjoint user-supplied function pointers
  nb::object lsjacfnB, lsjacfnBS, lsprecsetupfnB, lsprecsetupfnBS,
    lsprecsolvefnB, lsprecsolvefnBS, lsjactimessetupfnB, lsjactimessetupfnBS,
    lsjactimesvecfnB, lsjactimesvecfnBS;
};

inline idasa_user_supplied_fn_table* idasa_user_supplied_fn_table_alloc()
{
  // We must use malloc since IDASFree calls free
  auto fn_table = static_cast<idasa_user_supplied_fn_table*>(
    std::malloc(sizeof(idasa_user_supplied_fn_table)));

  // Zero out the memory
  std::memset(fn_table, 0, sizeof(idasa_user_supplied_fn_table));

  return fn_table;
}

template<typename... Args>
inline int idas_resB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<std::remove_pointer_t<IDAResFnB>,
                                             idasa_user_supplied_fn_table,
                                             1>(&idasa_user_supplied_fn_table::resB,
                                                std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_resQB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDAQuadRhsFnB>, idasa_user_supplied_fn_table,
    1>(&idasa_user_supplied_fn_table::resQB, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_resQBs_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDAQuadRhsFnBS>, idasa_user_supplied_fn_table,
    1>(&idasa_user_supplied_fn_table::resQBs, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjacfnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacFnB>, idasa_user_supplied_fn_table,
    4>(&idasa_user_supplied_fn_table::lsjacfnB, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjacfnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacFnBS>, idasa_user_supplied_fn_table,
    4>(&idasa_user_supplied_fn_table::lsjacfnBS, std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsprecsetupfnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSetupFnB>, idasa_user_supplied_fn_table,
    1>(&idasa_user_supplied_fn_table::lsprecsetupfnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsprecsetupfnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSetupFnBS>, idasa_user_supplied_fn_table,
    1>(&idasa_user_supplied_fn_table::lsprecsetupfnBS,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsprecsolvefnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSolveFnB>, idasa_user_supplied_fn_table,
    1>(&idasa_user_supplied_fn_table::lsprecsolvefnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsprecsolvefnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSolveFnBS>, idasa_user_supplied_fn_table,
    1>(&idasa_user_supplied_fn_table::lsprecsolvefnBS,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjactimessetupfnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesSetupFnB>, idasa_user_supplied_fn_table,
    1>(&idasa_user_supplied_fn_table::lsjactimessetupfnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjactimessetupfnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesSetupFnBS>, idasa_user_supplied_fn_table,
    1>(&idasa_user_supplied_fn_table::lsjactimessetupfnBS,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjactimesvecfnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesVecFnB>, idasa_user_supplied_fn_table,
    3>(&idasa_user_supplied_fn_table::lsjactimesvecfnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int idas_lsjactimesvecfnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesVecFnBS>, idasa_user_supplied_fn_table,
    3>(&idasa_user_supplied_fn_table::lsjactimesvecfnBS,
       std::forward<Args>(args)...);
}

#endif
