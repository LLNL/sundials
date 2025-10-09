/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2025, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#include <optional>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <sundials/sundials_core.hpp>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes.hpp>
#include <cvodes/cvodes_ls.h>
#include <cvodes/cvodes_proj.h>

#include "cvodes/cvodes_impl.h"
#include "sundials_adjointcheckpointscheme_impl.h"

#include "cvodes_usersupplied.hpp"

namespace nb = nanobind;

using namespace sundials::experimental;

#define BIND_CVODE_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)               \
  m.def(                                                                       \
    #NAME,                                                                     \
    [](void* cv_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn)         \
    {                                                                          \
      void* user_data = nullptr;                                               \
      CVodeGetUserData(cv_mem, &user_data);                                    \
      if (!user_data)                                                          \
        throw std::runtime_error(                                              \
          "Failed to get Python function table from CVODE memory");            \
      auto fntable    = static_cast<cvode_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER = nb::cast(fn);                                          \
      if (fn) { return NAME(cv_mem, &WRAPPER); }                               \
      else { return NAME(cv_mem, nullptr); }                                   \
    },                                                                         \
    __VA_ARGS__)

#define BIND_CVODE_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,       \
                             MEMBER2, WRAPPER2, ...)                            \
  m.def(                                                                        \
    #NAME,                                                                      \
    [](void* cv_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1,        \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                      \
    {                                                                           \
      void* user_data = nullptr;                                                \
      CVodeGetUserData(cv_mem, &user_data);                                     \
      if (!user_data)                                                           \
        throw std::runtime_error(                                               \
          "Failed to get Python function table from CVODE memory");             \
      auto fntable     = static_cast<cvode_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER1 = nb::cast(fn1);                                         \
      fntable->MEMBER2 = nb::cast(fn2);                                         \
      if (fn1) { return NAME(cv_mem, WRAPPER1, WRAPPER2); }                     \
      else { return NAME(cv_mem, nullptr, WRAPPER2); }                          \
    },                                                                          \
    __VA_ARGS__)

#define BIND_CVODEB_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)                 \
  m.def(                                                                          \
    #NAME,                                                                        \
    [](void* cv_mem, int which, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                             \
      void* user_data = nullptr;                                                  \
      CVodeGetUserDataB(cv_mem, which, &user_data);                               \
      if (!user_data)                                                             \
        throw std::runtime_error(                                                 \
          "Failed to get Python function table from CVODE memory");               \
      auto fntable    = static_cast<cvodea_user_supplied_fn_table*>(user_data);   \
      fntable->MEMBER = nb::cast(fn);                                             \
      if (fn) { return NAME(cv_mem, which, &WRAPPER); }                           \
      else { return NAME(cv_mem, which, nullptr); }                               \
    },                                                                            \
    __VA_ARGS__)

#define BIND_CVODEB_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,       \
                              MEMBER2, WRAPPER2, ...)                            \
  m.def(                                                                         \
    #NAME,                                                                       \
    [](void* cv_mem, int which,                                                  \
       std::function<std::remove_pointer_t<FN_TYPE1>> fn1,                       \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                       \
    {                                                                            \
      void* user_data = nullptr;                                                 \
      CVodeGetUserDataB(cv_mem, which, &user_data);                              \
      if (!user_data)                                                            \
        throw std::runtime_error(                                                \
          "Failed to get Python function table from CVODE memory");              \
      auto fntable     = static_cast<cvodea_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER1 = nb::cast(fn1);                                          \
      fntable->MEMBER2 = nb::cast(fn2);                                          \
      if (fn1) { return NAME(cv_mem, which, WRAPPER1, WRAPPER2); }               \
      else { return NAME(cv_mem, which, nullptr, WRAPPER2); }                    \
    },                                                                           \
    __VA_ARGS__)

void bind_cvodes(nb::module_& m)
{
#include "cvodes_generated.hpp"

  nb::class_<CVodeView>(m, "CVodeView")
    .def_static("Create", &CVodeView::Create<void*>)
    .def("get", nb::overload_cast<>(&CVodeView::get, nb::const_),
         nb::rv_policy::reference);

  m.def("CVodeInit",
        [](void* cv_mem, std::function<std::remove_pointer_t<CVRhsFn>> rhs,
           sunrealtype t0, N_Vector y0)
        {
          int cv_status = CVodeInit(cv_mem, cvode_f_wrapper, t0, y0);

          // Create the user-supplied function table to store the Python user functions
          auto cb_fns = cvode_user_supplied_fn_table_alloc();

          // Smuggle the user-supplied function table into callback wrappers through the user_data pointer
          cv_status = CVodeSetUserData(cv_mem, static_cast<void*>(cb_fns));
          if (cv_status != CV_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error("Failed to set user data in CVODE memory");
          }

          // Ensure CVodeFree will free the user-supplied function table
          cv_status = cvSetOwnUserData(cv_mem, SUNTRUE);
          if (cv_status != CV_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error(
              "Failed to set user data ownership in CVODE memory");
          }

          // Finally, set the RHS function
          cb_fns->f = nb::cast(rhs);

          return cv_status;
        });

  m.def("CVodeRootInit",
        [](void* cv_mem, int nrtfn,
           std::function<std::remove_pointer_t<CVRootFn>> fn)
        {
          void* user_data = nullptr;
          CVodeGetUserData(cv_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvode_user_supplied_fn_table*>(user_data);
          fntable->rootfn = nb::cast(fn);
          return CVodeRootInit(cv_mem, nrtfn, &cvode_rootfn_wrapper);
        });

  m.def("CVodeQuadInit",
        [](void* cv_mem, std::function<std::remove_pointer_t<CVQuadRhsFn>> fQ,
           N_Vector yQ0)
        {
          void* user_data = nullptr;
          CVodeGetUserData(cv_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvode_user_supplied_fn_table*>(user_data);
          fntable->fQ  = nb::cast(fQ);
          return CVodeQuadInit(cv_mem, &cvode_fQ_wrapper, yQ0);
        });

  BIND_CVODE_CALLBACK(CVodeWFtolerances, CVEwtFn, ewtn, cvode_ewtfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("efun").none());

  BIND_CVODE_CALLBACK(CVodeSetNlsRhsFn, CVRhsFn, fNLS, cvode_nlsrhsfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("f").none());

  BIND_CVODE_CALLBACK(CVodeSetJacFn, CVLsJacFn, lsjacfn, cvode_lsjacfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("jac").none());

  BIND_CVODE_CALLBACK2(CVodeSetPreconditioner, CVLsPrecSetupFn, lsprecsetupfn,
                       cvode_lsprecsetupfn_wrapper, CVLsPrecSolveFn,
                       lsprecsolvefn, cvode_lsprecsolvefn_wrapper,
                       nb::arg("cvode_mem"), nb::arg("pset").none(),
                       nb::arg("psolve").none());

  BIND_CVODE_CALLBACK2(CVodeSetJacTimes, CVLsJacTimesSetupFn, lsjactimessetupfn,
                       cvode_lsjactimessetupfn_wrapper, CVLsJacTimesVecFn,
                       lsjactimesvecfn, cvode_lsjactimesvecfn_wrapper,
                       nb::arg("cvode_mem"), nb::arg("jtsetup").none(),
                       nb::arg("jtimes").none());

  BIND_CVODE_CALLBACK(CVodeSetLinSysFn, CVLsLinSysFn, lslinsysfn,
                      cvode_lslinsysfn_wrapper, nb::arg("cvode_mem"),
                      nb::arg("linsys").none());

  BIND_CVODE_CALLBACK(CVodeSetJacTimesRhsFn, CVLsJacTimesVecFn, lsjacrhsfn,
                      cvode_lsjacrhsfn_wrapper, nb::arg("cvode_mem"),
                      nb::arg("jtimesRhsFn").none());

  BIND_CVODE_CALLBACK(CVodeSetProjFn, CVProjFn, projfn, cvode_projfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("pfun").none());

  // We cant use std::function<std::remove_pointer_t<CVQuadSensRhsFn>>
  // because we need to replace the N_Vector* arguments with std::vector<N_Vector>.
  using CVQuadSensRhsStdFn = int(sunrealtype t, N_Vector y,
                                 std::vector<N_Vector> yQS, N_Vector ydot,
                                 std::vector<N_Vector> yQSdot, void* user_data);
  m.def("CVodeQuadSensInit",
        [](void* cv_mem, std::function<CVQuadSensRhsStdFn> fQS,
           std::vector<N_Vector> yQS0)
        {
          void* user_data = nullptr;
          CVodeGetUserData(cv_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvode_user_supplied_fn_table*>(user_data);
          fntable->fQS = nb::cast(fQS);
          return CVodeQuadSensInit(cv_mem, cvode_fQS_wrapper, yQS0.data());
        });

  using CVSensRhsStdFn = int(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                             int iS, std::vector<N_Vector> yS,
                             std::vector<N_Vector> ySdot, void* user_data);
  m.def("CVodeSensInit",
        [](void* cv_mem, int Ns, int ism, std::function<CVSensRhsStdFn> fS,
           std::vector<N_Vector> yS0)
        {
          void* user_data = nullptr;
          CVodeGetUserData(cv_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvode_user_supplied_fn_table*>(user_data);
          fntable->fS  = nb::cast(fS);
          return CVodeSensInit(cv_mem, Ns, ism, cvode_fS_wrapper, yS0.data());
        });

  using CVSensRhs1StdFn = int(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                              int iS, N_Vector yS, N_Vector ySdot,
                              void* user_data, N_Vector tmp1, N_Vector tmp2);
  m.def("CVodeSensInit1",
        [](void* cv_mem, int Ns, std::function<CVSensRhs1StdFn> fS1, int ism,
           std::vector<N_Vector> yS0)
        {
          void* user_data = nullptr;
          CVodeGetUserData(cv_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvode_user_supplied_fn_table*>(user_data);
          fntable->fS1 = nb::cast(fS1);
          return CVodeSensInit1(cv_mem, Ns, ism, cvode_fS1_wrapper, yS0.data());
        });

  ///
  // CVODES Adjoint Bindings
  ///

  m.def("CVodeInitB",
        [](void* cv_mem, int which,
           std::function<std::remove_pointer_t<CVRhsFnB>> fB, sunrealtype tB0,
           N_Vector yB0)
        {
          int cv_status = CVodeInitB(cv_mem, which, cvode_fB_wrapper, tB0, yB0);

          // Create the user-supplied function table to store the Python user functions
          auto cb_fns = cvodea_user_supplied_fn_table_alloc();

          // Smuggle the user-supplied function table into callback wrappers through the user_data pointer
          cv_status = CVodeSetUserDataB(cv_mem, which,
                                        static_cast<void*>(cb_fns));
          if (cv_status != CV_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error("Failed to set user data in CVODE memory");
          }

          // Ensure CVodeFree will free the user-supplied function table
          cv_status = cvSetOwnUserDataB(cv_mem, which, SUNTRUE);
          if (cv_status != CV_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error(
              "Failed to set user data ownership in CVODE memory");
          }

          // Finally, set the RHS function
          cb_fns->fB = nb::cast(fB);

          return cv_status;
        });

  m.def("CVodeQuadInitB",
        [](void* cv_mem, int which,
           std::function<std::remove_pointer_t<CVQuadRhsFnB>> fQB, N_Vector yQBO)
        {
          void* user_data = nullptr;
          CVodeGetUserDataB(cv_mem, which, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvodea_user_supplied_fn_table*>(user_data);
          fntable->fQB = nb::cast(fQB);
          return CVodeQuadInitB(cv_mem, which, cvode_fQB_wrapper, yQBO);
        });

  using CVQuadRhsStdFnBS = int(sunrealtype t, N_Vector y,
                               std::vector<N_Vector> yS, N_Vector yB,
                               N_Vector qBdot, void* user_dataB);
  m.def("CVodeQuadInitBS",
        [](void* cv_mem, int which, std::function<CVQuadRhsStdFnBS> fQBs,
           N_Vector yQBO)
        {
          void* user_data = nullptr;
          CVodeGetUserDataB(cv_mem, which, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvodea_user_supplied_fn_table*>(user_data);
          fntable->fQBs = nb::cast(fQBs);
          return CVodeQuadInitBS(cv_mem, which, cvode_fQBs_wrapper, yQBO);
        });

  BIND_CVODEB_CALLBACK(CVodeSetJacFnB, CVLsJacFnB, lsjacfnB,
                       cvode_lsjacfnB_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("jacB").none());

  BIND_CVODEB_CALLBACK2(CVodeSetPreconditionerB, CVLsPrecSetupFnB, lsprecsetupfnB,
                        cvode_lsprecsetupfnB_wrapper, CVLsPrecSolveFnB,
                        lsprecsolvefnB, cvode_lsprecsolvefnB_wrapper,
                        nb::arg("cv_mem"), nb::arg("which"),
                        nb::arg("psetB").none(), nb::arg("psolveB").none());

  BIND_CVODEB_CALLBACK2(CVodeSetJacTimesB, CVLsJacTimesSetupFnB,
                        lsjactimessetupfnB, cvode_lsjactimessetupfnB_wrapper,
                        CVLsJacTimesVecFnB, lsjactimesvecfnB,
                        cvode_lsjactimesvecfnB_wrapper, nb::arg("cv_mem"),
                        nb::arg("which"), nb::arg("jsetupB").none(),
                        nb::arg("jtimesB").none());

  BIND_CVODEB_CALLBACK(CVodeSetLinSysFnB, CVLsLinSysFnB, lslinsysfnB,
                       cvode_lslinsysfnB_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("linsysB").none());

  using CVLsJacStdFnBS = int(sunrealtype t, N_Vector y,
                             std::vector<N_Vector> yS, N_Vector yB,
                             N_Vector fyB, SUNMatrix JB, void* user_dataB,
                             N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
  BIND_CVODEB_CALLBACK(CVodeSetJacFnBS, CVLsJacStdFnBS, lsjacfnBS,
                       cvode_lsjacfnBS_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("jacBS").none());

  using CVLsPrecSetupStdFnBS = int(sunrealtype t, N_Vector y, N_Vector yB,
                                   N_Vector fyB, sunbooleantype jokB,
                                   sunbooleantype * jcurPtrB,
                                   sunrealtype gammaB, void* user_dataB);
  using CVLsPrecSolveStdFnBS =
    int(sunrealtype t, N_Vector y, std::vector<N_Vector> yS, N_Vector yB,
        N_Vector fyB, N_Vector rB, N_Vector zB, sunrealtype gammaB,
        sunrealtype deltaB, int lrB, void* user_dataB);
  BIND_CVODEB_CALLBACK2(CVodeSetPreconditionerBS, CVLsPrecSetupStdFnBS,
                        lsprecsetupfnBS, cvode_lsprecsetupfnBS_wrapper,
                        CVLsPrecSolveStdFnBS, lsprecsolvefnBS,
                        cvode_lsprecsolvefnBS_wrapper, nb::arg("cv_mem"),
                        nb::arg("which"), nb::arg("psetBS").none(),
                        nb::arg("psolveBS").none());

  using CVLsJacTimesSetupStdFnBS = int(sunrealtype t, N_Vector y,
                                       std::vector<N_Vector> yS, N_Vector yB,
                                       N_Vector fyB, void* jac_dataB);
  using CVLsJacTimesVecStdFnBS   = int(N_Vector vB, N_Vector JvB, sunrealtype t,
                                     N_Vector y, std::vector<N_Vector> yS,
                                     N_Vector yB, N_Vector fyB, void* jac_dataB,
                                     N_Vector tmpB);
  BIND_CVODEB_CALLBACK2(CVodeSetJacTimesBS, CVLsJacTimesSetupStdFnBS,
                        lsjactimessetupfnBS, cvode_lsjactimessetupfnBS_wrapper,
                        CVLsJacTimesVecStdFnBS, lsjactimesvecfnBS,
                        cvode_lsjactimesvecfnBS_wrapper, nb::arg("cv_mem"),
                        nb::arg("which"), nb::arg("jsetupBS").none(),
                        nb::arg("jtimesBS").none());

  using CVLsLinSysStdFnBS =
    int(sunrealtype t, N_Vector y, std::vector<N_Vector> yS, N_Vector yB,
        N_Vector fyB, SUNMatrix AB, sunbooleantype jokB, sunbooleantype * jcurB,
        sunrealtype gammaB, void* user_dataB, N_Vector tmp1B, N_Vector tmp2B,
        N_Vector tmp3B);
  BIND_CVODEB_CALLBACK(CVodeSetLinSysFnBS, CVLsLinSysStdFnBS, lslinsysfnBS,
                       cvode_lslinsysfnBS_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("linsysBS").none());
}
