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

#ifndef _PYSUNDIALS_CVODE_USERSUPPLIED_HPP
#define _PYSUNDIALS_CVODE_USERSUPPLIED_HPP

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_ls.h>

#include <sundials/sundials_core.hpp>

#include "pysundials_helpers.hpp"

///////////////////////////////////////////////////////////////////////////////
// CVODE user-supplied function table
// Every integrator-level user-supplied function must be in this table.
// The user-supplied function table is passed to CVODE as user_data.
///////////////////////////////////////////////////////////////////////////////

struct cvode_user_supplied_fn_table
{
  // user-supplied function pointers
  nb::object f, rootfn, ewtn, rwtn, nlsfi, projfn;

  // cvode_ls user-supplied function pointers
  nb::object lsjacfn, lsprecsetupfn, lsprecsolvefn, lsjactimessetupfn,
    lsjactimesvecfn, lslinsysfn, lsjacrhsfn;

  // cvode quadrature user-supplied function pointers
  nanobind::object fQS;

  // cvode FSA user-supplied function pointers
  nanobind::object fS, fS1;
};

///////////////////////////////////////////////////////////////////////////////
// CVODE user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline cvode_user_supplied_fn_table* cvode_user_supplied_fn_table_alloc()
{
  // We must use malloc since CVODEFree calls free
  auto fn_table = static_cast<cvode_user_supplied_fn_table*>(
    std::malloc(sizeof(cvode_user_supplied_fn_table)));

  // Zero out the memory
  std::memset(fn_table, 0, sizeof(cvode_user_supplied_fn_table));

  return fn_table;
}

template<typename... Args>
inline int cvode_f_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<std::remove_pointer_t<CVRhsFn>,
                                             cvode_user_supplied_fn_table,
                                             1>(&cvode_user_supplied_fn_table::f,
                                                std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_rootfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVRootFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::rootfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_ewtfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVEwtFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::ewtn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_nlsrhsfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::nlsfi, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjacfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacFn>, cvode_user_supplied_fn_table,
    4>(&cvode_user_supplied_fn_table::lsjacfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsprecsetupfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSetupFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsprecsetupfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsprecsolvefn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSolveFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsprecsolvefn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimessetupfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesSetupFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsjactimessetupfn,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimesvecfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesVecFn>, cvode_user_supplied_fn_table,
    2>(&cvode_user_supplied_fn_table::lsjactimesvecfn,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lslinsysfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsLinSysFn>, cvode_user_supplied_fn_table,
    4>(&cvode_user_supplied_fn_table::lslinsysfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjacrhsfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsjacrhsfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_fqs_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVQuadSensRhsFn>, cvode_user_supplied_fn_table,
    3>(&cvode_user_supplied_fn_table::fQS, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_fs_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<std::remove_pointer_t<CVSensRhsFn>,
                                             cvode_user_supplied_fn_table,
                                             3>(&cvode_user_supplied_fn_table::fS,
                                                std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_fs1_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVSensRhs1Fn>, cvode_user_supplied_fn_table,
    3>(&cvode_user_supplied_fn_table::fS1, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_projfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVProjFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::projfn, std::forward<Args>(args)...);
}

///////////////////////////////////////////////////////////////////////////////
// CVODE Adjoint user-supplied functions
///////////////////////////////////////////////////////////////////////////////

struct cvodea_user_supplied_fn_table
{
  // cvode adjoint user-supplied function pointers
  nb::object fB, fQB, fQBs;

  // cvode_ls adjoint user-supplied function pointers
  nb::object lsjacfnB, lsjacfnBS, lsprecsetupfnB, lsprecsetupfnBS,
    lsprecsolvefnB, lsprecsolvefnBS, lsjactimessetupfnB, lsjactimessetupfnBS,
    lsjactimesvecfnB, lsjactimesvecfnBS, lslinsysfnB, lslinsysfnBS;
};

inline cvodea_user_supplied_fn_table* cvodea_user_supplied_fn_table_alloc()
{
  // We must use malloc since CVODEFree calls free
  auto fn_table = static_cast<cvodea_user_supplied_fn_table*>(
    std::malloc(sizeof(cvodea_user_supplied_fn_table)));

  // Zero out the memory
  std::memset(fn_table, 0, sizeof(cvodea_user_supplied_fn_table));

  return fn_table;
}

template<typename... Args>
inline int cvode_fB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<std::remove_pointer_t<CVRhsFnB>,
                                             cvodea_user_supplied_fn_table,
                                             1>(&cvodea_user_supplied_fn_table::fB,
                                                std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_fQB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVQuadRhsFnB>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::fQB, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_fQBs_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVQuadRhsFnBS>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::fQBs, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjacfnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacFnB>, cvodea_user_supplied_fn_table,
    4>(&cvodea_user_supplied_fn_table::lsjacfnB, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjacfnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacFnBS>, cvodea_user_supplied_fn_table,
    4>(&cvodea_user_supplied_fn_table::lsjacfnBS, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsprecsetupfnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSetupFnB>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::lsprecsetupfnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsprecsetupfnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSetupFnBS>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::lsprecsetupfnBS,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsprecsolvefnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSolveFnB>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::lsprecsolvefnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsprecsolvefnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSolveFnBS>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::lsprecsolvefnBS,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimessetupfnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesSetupFnB>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::lsjactimessetupfnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimessetupfnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesSetupFnBS>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::lsjactimessetupfnBS,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimesvecfnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesVecFnB>, cvodea_user_supplied_fn_table,
    2>(&cvodea_user_supplied_fn_table::lsjactimesvecfnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimesvecfnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesVecFnBS>, cvodea_user_supplied_fn_table,
    2>(&cvodea_user_supplied_fn_table::lsjactimesvecfnBS,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lslinsysfnB_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsLinSysFnB>, cvodea_user_supplied_fn_table,
    4>(&cvodea_user_supplied_fn_table::lslinsysfnB, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lslinsysfnBS_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsLinSysFnBS>, cvodea_user_supplied_fn_table,
    4>(&cvodea_user_supplied_fn_table::lslinsysfnBS, std::forward<Args>(args)...);
}

#endif
