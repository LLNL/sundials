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

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <sundials/sundials_core.hpp>

#include <kinsol/kinsol.h>
#include <kinsol/kinsol.hpp>
#include <kinsol/kinsol_ls.h>

#include "kinsol/kinsol_impl.h"

#include "pysundials_types.hpp"

namespace nb = nanobind;

using namespace sundials::experimental;

#include "pysundials_kinsol_usersupplied.hpp"

#define BIND_KINSOL_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)               \
  m.def(                                                                        \
    #NAME,                                                                      \
    [](void* kin_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn)         \
    {                                                                           \
      void* user_data = nullptr;                                                \
      KINGetUserData(kin_mem, &user_data);                                      \
      if (!user_data)                                                           \
        throw std::runtime_error(                                               \
          "Failed to get Python function table from KINSOL memory");            \
      auto fntable    = static_cast<kinsol_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER = nb::cast(fn);                                           \
      if (fn) { return NAME(kin_mem, WRAPPER); }                                \
      else { return NAME(kin_mem, nullptr); }                                   \
    },                                                                          \
    __VA_ARGS__)

#define BIND_KINSOL_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,       \
                              MEMBER2, WRAPPER2, ...)                            \
  m.def(                                                                         \
    #NAME,                                                                       \
    [](void* kin_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1,        \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                       \
    {                                                                            \
      void* user_data = nullptr;                                                 \
      KINGetUserData(kin_mem, &user_data);                                       \
      if (!user_data)                                                            \
        throw std::runtime_error(                                                \
          "Failed to get Python function table from KINSOL memory");             \
      auto fntable     = static_cast<kinsol_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER1 = nb::cast(fn1);                                          \
      fntable->MEMBER2 = nb::cast(fn2);                                          \
      if (fn1) { return NAME(kin_mem, WRAPPER1, WRAPPER2); }                     \
      else { return NAME(kin_mem, nullptr, WRAPPER2); }                          \
    },                                                                           \
    __VA_ARGS__)

void bind_kinsol(nb::module_& m)
{
#include "pysundials_kinsol_generated.hpp"

  nb::class_<KINView>(m, "KINView")
    .def_static("Create", &KINView::Create<void*>)
    .def("get", nb::overload_cast<>(&KINView::get, nb::const_),
         nb::rv_policy::reference);

  m.def("KINInit",
        [](void* kin_mem, std::function<std::remove_pointer_t<KINSysFn>> sysfn,
           N_Vector tmpl)
        {
          int kin_status = KINInit(kin_mem, kinsol_sysfn_wrapper, tmpl);

          auto cb_fns = kinsol_user_supplied_fn_table_alloc();
          kin_status  = KINSetUserData(kin_mem, static_cast<void*>(cb_fns));
          if (kin_status != KIN_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error(
              "Failed to set user data in KINSOL memory");
          }

          // Ensure KINFree will free the user-supplied function table
          kin_status = kinSetOwnUserData(kin_mem, SUNTRUE);
          if (kin_status != KIN_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error(
              "Failed to set user data ownership in KINSOL memory");
          }

          cb_fns->sysfn = nb::cast(sysfn);
          return kin_status;
        });

  BIND_KINSOL_CALLBACK(KINSetSysFunc, KINSysFn, sysfn, kinsol_sysfn_wrapper,
                       nb::arg("kin_mem"), nb::arg("sysfn"));

  BIND_KINSOL_CALLBACK(KINSetDampingFn, KINDampingFn, dampingfn,
                       kinsol_dampingfn_wrapper, nb::arg("kin_mem"),
                       nb::arg("damping_fn").none());

  using KINDepthStdFn =
    int(long int iter, N_Vector u_val, N_Vector g_val, N_Vector f_val,
        std::vector<N_Vector> df, pysundials::array1d R_mat, long int depth,
        void* user_data, long int* new_depth, sunbooleantype* remove_indices);
  m.def(
    "KINSetDepthFn",
    [](void* kin_mem, std::function<KINDepthStdFn> depth_fn)
    {
      void* user_data = nullptr;
      KINGetUserData(kin_mem, &user_data);
      if (!user_data)
        throw std::runtime_error(
          "Failed to get Python function table from KINSOL memory");
      auto fntable     = static_cast<kinsol_user_supplied_fn_table*>(user_data);
      fntable->depthfn = nb::cast(depth_fn);
      return KINSetDepthFn(kin_mem, kinsol_depthfn_wrapper);
    },
    nb::arg("kin_mem"), nb::arg("depth_fn").none());

  BIND_KINSOL_CALLBACK2(KINSetPreconditioner, KINLsPrecSetupFn, lsprecsetupfn,
                        kinsol_lsprecsetupfn_wrapper, KINLsPrecSolveFn,
                        lsprecsolvefn, kinsol_lsprecsolvefn_wrapper,
                        nb::arg("kin_mem"), nb::arg("psetup").none(),
                        nb::arg("psolve").none());

  BIND_KINSOL_CALLBACK(KINSetJacFn, KINSysFn, lsjacfn, kinsol_lsjacfn_wrapper,
                       nb::arg("kin_mem"), nb::arg("jac").none());

  BIND_KINSOL_CALLBACK(KINSetJacTimesVecFn, KINLsJacTimesVecFn, lsjactimesvecfn,
                       kinsol_lsjactimesvecfn_wrapper, nb::arg("kin_mem"),
                       nb::arg("jtimes").none());

  BIND_KINSOL_CALLBACK(KINSetJacTimesVecSysFn, KINSysFn, lsjtvsysfn,
                       kinsol_lsjtvsysfn_wrapper, nb::arg("kin_mem"),
                       nb::arg("jtvSysFn").none());
}
