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
  nb::object rootfn;
  nb::object ewtn;
  nb::object rwtn;
  nb::object nlsfi;
  nb::object f;

  // cvode_ls user-supplied function pointers
  nb::object lsjacfn;
  nb::object lsprecsetupfn;
  nb::object lsprecsolvefn;
  nb::object lsjactimessetupfn;
  nb::object lsjactimesvecfn;
  nb::object lslinsysfn;
  nb::object lsjacrhsfn;

  // cvode quadrature user-supplied function pointers
  nanobind::object fQS;

  // cvode FSA user-supplied function pointers
  nanobind::object fS;
  nanobind::object fS1;
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

inline int cvode_rootfn_wrapper(sunrealtype t, N_Vector y, sunrealtype* gout,
                                void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVRootFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::rootfn, t, y, gout, user_data);
}

inline int cvode_ewtfn_wrapper(N_Vector y, N_Vector ewt, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVEwtFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::ewtn, y, ewt, user_data);
}

inline int cvode_nlsrhsfn_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                                  void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::nlsfi, t, y, ydot, user_data);
}

inline int cvode_lsjacfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                 SUNMatrix Jac, void* user_data, N_Vector tmp1,
                                 N_Vector tmp2, N_Vector tmp3)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacFn>, cvode_user_supplied_fn_table,
    4>(&cvode_user_supplied_fn_table::lsjacfn, t, y, fy, Jac, user_data, tmp1,
       tmp2, tmp3);
}

inline int cvode_lsprecsetupfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                       sunbooleantype jok,
                                       sunbooleantype* jcurPtr,
                                       sunrealtype gamma, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSetupFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsprecsetupfn, t, y, fy, jok, jcurPtr,
       gamma, user_data);
}

inline int cvode_lsprecsolvefn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                       N_Vector r, N_Vector z, sunrealtype gamma,
                                       sunrealtype delta, int lr, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSolveFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsprecsolvefn, t, y, fy, r, z, gamma,
       delta, lr, user_data);
}

inline int cvode_lsjactimessetupfn_wrapper(sunrealtype t, N_Vector y,
                                           N_Vector fy, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesSetupFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsjactimessetupfn, t, y, fy, user_data);
}

inline int cvode_lsjactimesvecfn_wrapper(N_Vector v, N_Vector Jv, sunrealtype t,
                                         N_Vector y, N_Vector fy,
                                         void* user_data, N_Vector tmp)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesVecFn>, cvode_user_supplied_fn_table,
    2>(&cvode_user_supplied_fn_table::lsjactimesvecfn, v, Jv, t, y, fy,
       user_data, tmp);
}

inline int cvode_lslinsysfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                    SUNMatrix A, sunbooleantype jok,
                                    sunbooleantype* jcur, sunrealtype gamma,
                                    void* user_data, N_Vector tmp1,
                                    N_Vector tmp2, N_Vector tmp3)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsLinSysFn>, cvode_user_supplied_fn_table,
    4>(&cvode_user_supplied_fn_table::lslinsysfn, t, y, fy, A, jok, jcur, gamma,
       user_data, tmp1, tmp2, tmp3);
}

inline int cvode_lsjacrhsfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                    void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsjacrhsfn, t, y, fy, user_data);
}

inline int cvode_f_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                           void* user_data)
{
  return pysundials::user_supplied_fn_caller<std::remove_pointer_t<CVRhsFn>,
                                             cvode_user_supplied_fn_table,
                                             1>(&cvode_user_supplied_fn_table::f,
                                                t, y, ydot, user_data);
}

inline int cvode_fqs_wrapper(int Ns, sunrealtype t, N_Vector y, N_Vector* yS,
                             N_Vector yQdot, N_Vector* yQSdot, void* user_data,
                             N_Vector tmp, N_Vector tmpQ)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<CVQuadSensRhsFn>, cvode_user_supplied_fn_table,
    3>(&cvode_user_supplied_fn_table::fQS, Ns, t, y, yS, yQdot, yQSdot,
       user_data, tmp, tmpQ);
}

inline int cvode_fs_wrapper(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                            N_Vector* yS, N_Vector* ySdot, void* user_data,
                            N_Vector tmp1, N_Vector tmp2)
{
  return pysundials::user_supplied_fn_caller<std::remove_pointer_t<CVSensRhsFn>,
                                             cvode_user_supplied_fn_table,
                                             3>(&cvode_user_supplied_fn_table::fS,
                                                Ns, t, y, ydot, yS, ySdot,
                                                user_data, tmp1, tmp2);
}

inline int cvode_fs1_wrapper(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                             int iS, N_Vector yS, N_Vector ySdot,
                             void* user_data, N_Vector tmp1, N_Vector tmp2)
{
  return pysundials::user_supplied_fn_caller<std::remove_pointer_t<CVSensRhs1Fn>,
                                             cvode_user_supplied_fn_table,
                                             3>(&cvode_user_supplied_fn_table::fS1,
                                                Ns, t, y, ydot, iS, yS, ySdot,
                                                user_data, tmp1, tmp2);
}

#endif
