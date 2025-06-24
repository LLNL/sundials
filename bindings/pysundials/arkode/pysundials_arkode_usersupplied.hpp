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

#ifndef _PYSUNDIALS_ARKODE_USERSUPPLIED_HPP
#define _PYSUNDIALS_ARKODE_USERSUPPLIED_HPP

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_lsrkstep.h>
#include <arkode/arkode_mristep.h>

#include <sundials/sundials_core.hpp>

#include "pysundials_helpers.hpp"

///////////////////////////////////////////////////////////////////////////////
// ARKODE user-supplied function table
// Every integrator-level user-supplied function must be in this table.
// The user-supplied function table is passed to ARKODE as user_data.
///////////////////////////////////////////////////////////////////////////////

struct arkode_user_supplied_fn_table
{
  // common user-supplied function pointers
  nb::object rootfn;
  nb::object ewtn;
  nb::object rwtn;
  nb::object adaptfn;
  nb::object expstabfn;
  nb::object vecresizefn;
  nb::object postprocessfn;
  nb::object stagepredictfn;
  nb::object relaxfn;
  nb::object relaxjacfn;
  nb::object nlsfi;

  // akrode_ls user-supplied function pointers
  nb::object lsjacfn;
  nb::object lsmassfn;
  nb::object lsprecsetupfn;
  nb::object lsprecsolvefn;
  nb::object lsjactimessetupfn;
  nb::object lsjactimesvecfn;
  nb::object lslinsysfn;
  nb::object lsmass_timessetupfn;
  nb::object lsmass_timesvecfn;
  nb::object lsmassprecsetupfn;
  nb::object lsmassprecsolvefn;
  nb::object lsjacrhsfn;

  // erkstep-specific user-supplied function pointers
  nb::object erkstep_f;

  // arkstep-specific user-supplied function pointers
  nb::object arkstep_fe;
  nb::object arkstep_fi;
};

///////////////////////////////////////////////////////////////////////////////
// ARKODE user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline arkode_user_supplied_fn_table* arkode_user_supplied_fn_table_alloc()
{
  // We must use malloc since ARKodeFree calls free
  return static_cast<arkode_user_supplied_fn_table*>(
    std::malloc(sizeof(arkode_user_supplied_fn_table)));
}

inline int arkode_rootfn_wrapper(sunrealtype t, N_Vector y, sunrealtype* gout,
                                 void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRootFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::rootfn, t, y,
                                   gout, user_data);
}

inline int arkode_ewtfn_wrapper(N_Vector y, N_Vector ewt, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKEwtFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::ewtn, y, ewt,
                                   user_data);
}

inline int arkode_rwtn_wrapper(N_Vector y, N_Vector rwt, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRwtFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::rwtn, y, rwt,
                                   user_data);
}

inline int arkode_adaptfn_wrapper(N_Vector y, sunrealtype t, sunrealtype h1,
                                  sunrealtype h2, sunrealtype h3, sunrealtype e1,
                                  sunrealtype e2, sunrealtype e3, int q, int p,
                                  sunrealtype* hnew, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKAdaptFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::adaptfn, y, t,
                                   h1, h2, h3, e1, e2, e3, q, p, hnew, user_data);
}

inline int arkode_expstabfn_wrapper(N_Vector y, sunrealtype t,
                                    sunrealtype* hstab, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKExpStabFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::expstabfn, y,
                                   t, hstab, user_data);
}

inline int arkode_vecresizefn_wrapper(N_Vector y, N_Vector ytemplate,
                                      void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKVecResizeFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::vecresizefn,
                                   y, ytemplate, user_data);
}

inline int arkode_postprocessfn_wrapper(sunrealtype t, N_Vector y, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKPostProcessFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::postprocessfn,
                                   t, y, user_data);
}

inline int arkode_stagepredictfn_wrapper(sunrealtype t, N_Vector zpred,
                                         void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKStagePredictFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::stagepredictfn,
                                   t, zpred, user_data);
}

inline int arkode_nlsrhsfn_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                                   void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::nlsfi, t, y,
                                   ydot, user_data);
}

inline int arkode_relaxfn_wrapper(N_Vector y, sunrealtype* r, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRelaxFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::relaxfn, y,
                                   r, user_data);
}

inline int arkode_relaxjacfn_wrapper(N_Vector y, N_Vector J, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRelaxJacFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::relaxjacfn,
                                   y, J, user_data);
}

inline int arkode_lsjacfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                  SUNMatrix Jac, void* user_data, N_Vector tmp1,
                                  N_Vector tmp2, N_Vector tmp3)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsJacFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsjacfn, t,
                                   y, fy, Jac, user_data, tmp1, tmp2, tmp3);
}

inline int arkode_lsmassfn_wrapper(sunrealtype t, SUNMatrix M, void* user_data,
                                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsmassfn, t,
                                   M, user_data, tmp1, tmp2, tmp3);
}

inline int arkode_lsprecsetupfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                        sunbooleantype jok,
                                        sunbooleantype* jcurPtr,
                                        sunrealtype gamma, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsPrecSetupFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsprecsetupfn,
                                   t, y, fy, jok, jcurPtr, gamma, user_data);
}

inline int arkode_lsprecsolvefn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                        N_Vector r, N_Vector z,
                                        sunrealtype gamma, sunrealtype delta,
                                        int lr, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsPrecSolveFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsprecsolvefn,
                                   t, y, fy, r, z, gamma, delta, lr, user_data);
}

inline int arkode_lsjactimessetupfn_wrapper(sunrealtype t, N_Vector y,
                                            N_Vector fy, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsJacTimesSetupFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsjactimessetupfn,
                                   t, y, fy, user_data);
}

inline int arkode_lsjactimesvecfn_wrapper(N_Vector v, N_Vector Jv,
                                          sunrealtype t, N_Vector y, N_Vector fy,
                                          void* user_data, N_Vector tmp)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsJacTimesVecFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsjactimesvecfn,
                                   v, Jv, t, y, fy, user_data, tmp);
}

inline int arkode_lslinsysfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                     SUNMatrix A, SUNMatrix M,
                                     sunbooleantype jok, sunbooleantype* jcur,
                                     sunrealtype gamma, void* user_data,
                                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsLinSysFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lslinsysfn,
                                   t, y, fy, A, M, jok, jcur, gamma, user_data,
                                   tmp1, tmp2, tmp3);
}

inline int arkode_lsmass_timessetupfn_wrapper(sunrealtype t, void* mtimes_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassTimesSetupFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsmass_timessetupfn,
                                   t, mtimes_data);
}

inline int arkode_lsmass_timesvecfn_wrapper(N_Vector v, N_Vector Mv,
                                            sunrealtype t, void* mtimes_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassTimesVecFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsmass_timesvecfn,
                                   v, Mv, t, mtimes_data);
}

inline int arkode_lsmassprecsetupfn_wrapper(sunrealtype t, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassPrecSetupFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsmassprecsetupfn,
                                   t, user_data);
}

inline int arkode_lsmassprecsolvefn_wrapper(sunrealtype t, N_Vector r,
                                            N_Vector z, sunrealtype delta,
                                            int lr, void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassPrecSolveFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsmassprecsolvefn,
                                   t, r, z, delta, lr, user_data);
}

inline int arkode_lsjacrhsfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                     void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::lsjacrhsfn,
                                   t, y, fy, user_data);
}

///////////////////////////////////////////////////////////////////////////////
// ERKStep user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline int erkstep_f_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                             void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::erkstep_f, t,
                                   y, ydot, user_data);
}

///////////////////////////////////////////////////////////////////////////////
// ARKStep user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline int arkstep_fe_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                              void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::arkstep_fe,
                                   t, y, ydot, user_data);
}

inline int arkstep_fi_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                              void* user_data)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>,
    arkode_user_supplied_fn_table>(&arkode_user_supplied_fn_table::arkstep_fi,
                                   t, y, ydot, user_data);
}

#endif
