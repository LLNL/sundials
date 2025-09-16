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

#include <optional>

#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <sundials/sundials_core.hpp>

#include <idas/idas.h>
#include <idas/idas.hpp>
#include <idas/idas_ls.h>

#include "pysundials_idas_usersupplied.hpp"

namespace nb = nanobind;

using namespace sundials::experimental;

#define BIND_IDA_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER)                      \
  m.def(#NAME,                                                                 \
        [](void* ida_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn)    \
        {                                                                      \
          void* user_data = nullptr;                                           \
          IDAGetUserData(ida_mem, &user_data);                                 \
          if (!user_data)                                                      \
            throw std::runtime_error(                                          \
              "Failed to get Python function table from IDAS memory");         \
          auto fntable = static_cast<idas_user_supplied_fn_table*>(user_data); \
          fntable->MEMBER = nb::cast(fn);                                      \
          return NAME(ida_mem, &WRAPPER);                                      \
        })

#define BIND_IDA_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,        \
                           MEMBER2, WRAPPER2, ...)                             \
  m.def(                                                                       \
    #NAME,                                                                     \
    [](void* ida_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1,      \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                     \
    {                                                                          \
      void* user_data = nullptr;                                               \
      IDAGetUserData(ida_mem, &user_data);                                     \
      if (!user_data)                                                          \
        throw std::runtime_error(                                              \
          "Failed to get Python function table from IDAS memory");             \
      auto fntable     = static_cast<idas_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER1 = nb::cast(fn1);                                        \
      fntable->MEMBER2 = nb::cast(fn2);                                        \
      if (fn1) { return NAME(ida_mem, WRAPPER1, WRAPPER2); }                   \
      else { return NAME(ida_mem, nullptr, WRAPPER2); }                        \
    },                                                                         \
    __VA_ARGS__)

#define BIND_IDAB_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER)                      \
  m.def(#NAME,                                                                  \
        [](void* ida_mem, int which,                                            \
           std::function<std::remove_pointer_t<FN_TYPE>> fn)                    \
        {                                                                       \
          void* user_data = nullptr;                                            \
          IDAGetUserDataB(ida_mem, which, &user_data);                          \
          if (!user_data)                                                       \
            throw std::runtime_error(                                           \
              "Failed to get Python function table from IDAS memory");          \
          auto fntable = static_cast<idasa_user_supplied_fn_table*>(user_data); \
          fntable->MEMBER = nb::cast(fn);                                       \
          return NAME(ida_mem, which, &WRAPPER);                                \
        })

#define BIND_IDAB_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,        \
                            MEMBER2, WRAPPER2, ...)                             \
  m.def(                                                                        \
    #NAME,                                                                      \
    [](void* ida_mem, int which,                                                \
       std::function<std::remove_pointer_t<FN_TYPE1>> fn1,                      \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                      \
    {                                                                           \
      void* user_data = nullptr;                                                \
      IDAGetUserDataB(ida_mem, which, &user_data);                              \
      if (!user_data)                                                           \
        throw std::runtime_error(                                               \
          "Failed to get Python function table from IDAS memory");              \
      auto fntable     = static_cast<idasa_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER1 = nb::cast(fn1);                                         \
      fntable->MEMBER2 = nb::cast(fn2);                                         \
      return NAME(ida_mem, which, &WRAPPER1, &WRAPPER2);                        \
    },                                                                          \
    __VA_ARGS__)

void bind_idas(nb::module_& m)
{
#include "pysundials_idas_generated.hpp"

  nb::class_<IDAView>(m, "IDAView")
    .def_static("Create", &IDAView::Create<void*>)
    .def("get", nb::overload_cast<>(&IDAView::get, nb::const_),
         nb::rv_policy::reference);

  m.def("IDAInit",
        [](void* ida_mem, std::function<std::remove_pointer_t<IDAResFn>> res,
           sunrealtype t0, N_Vector yy0, N_Vector yp0)
        {
          int ida_status = IDAInit(ida_mem, idas_res_wrapper, t0, yy0, yp0);

          auto cb_fns = idas_user_supplied_fn_table_alloc();
          ida_status  = IDASetUserData(ida_mem, static_cast<void*>(cb_fns));
          if (ida_status != IDA_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error("Failed to set user data in IDAS memory");
          }
          cb_fns->res = nb::cast(res);
          return ida_status;
        });

  m.def("IDARootInit",
        [](void* ida_mem, int nrtfn,
           std::function<std::remove_pointer_t<IDARootFn>> fn)
        {
          void* user_data = nullptr;
          IDAGetUserData(ida_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from IDAS memory");
          auto fntable = static_cast<idas_user_supplied_fn_table*>(user_data);
          fntable->rootfn = nb::cast(fn);
          return IDARootInit(ida_mem, nrtfn, &idas_rootfn_wrapper);
        });

  m.def("IDAQuadInit",
        [](void* ida_mem,
           std::function<std::remove_pointer_t<IDAQuadRhsFn>> resQ, N_Vector yQ0)
        {
          void* user_data = nullptr;
          IDAGetUserData(ida_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from IDAS memory");
          auto fntable  = static_cast<idas_user_supplied_fn_table*>(user_data);
          fntable->resQ = nb::cast(resQ);
          return IDAQuadInit(ida_mem, &idas_resQ_wrapper, yQ0);
        });

  BIND_IDA_CALLBACK(IDAWFtolerances, IDAEwtFn, ewtn, idas_ewtfn_wrapper);

  BIND_IDA_CALLBACK(IDASetNlsResFn, IDAResFn, resNLS, idas_nlsresfn_wrapper);

  BIND_IDA_CALLBACK(IDASetJacFn, IDALsJacFn, lsjacfn, idas_lsjacfn_wrapper);

  BIND_IDA_CALLBACK2(IDASetPreconditioner, IDALsPrecSetupFn, lsprecsetupfn,
                     idas_lsprecsetupfn_wrapper, IDALsPrecSolveFn, lsprecsolvefn,
                     idas_lsprecsolvefn_wrapper, nb::arg("ida_mem"),
                     nb::arg("psetup").none(), nb::arg("psolve"));

  BIND_IDA_CALLBACK2(IDASetJacTimes, IDALsJacTimesSetupFn, lsjactimessetupfn,
                     idas_lsjactimessetupfn_wrapper, IDALsJacTimesVecFn,
                     lsjactimesvecfn, idas_lsjactimesvecfn_wrapper,
                     nb::arg("ida_mem"), nb::arg("jsetup").none(),
                     nb::arg("psolve"));

  BIND_IDA_CALLBACK(IDASetJacTimesResFn, IDALsJacTimesVecFn, lsjacresfn,
                    idas_lsjacresfn_wrapper);

  // Sensitivity and quadrature sensitivity bindings (example, adapt as needed)
  using IDAQuadSensRhsStdFn = int(sunrealtype t, N_Vector yy,
                                  std::vector<N_Vector> yQS, N_Vector yp,
                                  std::vector<N_Vector> yQSdot, void* user_data);
  m.def("IDAQuadSensInit",
        [](void* ida_mem, std::function<IDAQuadSensRhsStdFn> resQS,
           std::vector<N_Vector> yQS0)
        {
          void* user_data = nullptr;
          IDAGetUserData(ida_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from IDAS memory");
          auto fntable   = static_cast<idas_user_supplied_fn_table*>(user_data);
          fntable->resQS = nb::cast(resQS);
          return IDAQuadSensInit(ida_mem, idas_resQS_wrapper, yQS0.data());
        });

  using IDASensResStdFn = int(int Ns, sunrealtype t, N_Vector yy, N_Vector yp,
                              int iS, std::vector<N_Vector> yS,
                              std::vector<N_Vector> ySdot, void* user_data);
  m.def("IDASensInit",
        [](void* ida_mem, int Ns, int ism, std::function<IDASensResStdFn> resS,
           std::vector<N_Vector> yS0, std::vector<N_Vector> ypS0)
        {
          void* user_data = nullptr;
          IDAGetUserData(ida_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from IDAS memory");
          auto fntable  = static_cast<idas_user_supplied_fn_table*>(user_data);
          fntable->resS = nb::cast(resS);
          return IDASensInit(ida_mem, Ns, ism, idas_resS_wrapper, yS0.data(),
                             ypS0.data());
        });

  ///
  // IDAS Adjoint Bindings
  ///

  m.def("IDAInitB",
        [](void* ida_mem, int which,
           std::function<std::remove_pointer_t<IDAResFnB>> resB,
           sunrealtype tB0, N_Vector yyB0, N_Vector ypB0)
        {
          int ida_status = IDAInitB(ida_mem, which, idas_resB_wrapper, tB0,
                                    yyB0, ypB0);

          auto cb_fns = idasa_user_supplied_fn_table_alloc();
          ida_status  = IDASetUserDataB(ida_mem, which,
                                        static_cast<void*>(cb_fns));
          if (ida_status != IDA_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error("Failed to set user data in IDAS memory");
          }
          cb_fns->resB = nb::cast(resB);
          return ida_status;
        });

  m.def("IDAQuadInitB",
        [](void* ida_mem, int which,
           std::function<std::remove_pointer_t<IDAQuadRhsFnB>> resQB,
           N_Vector yQBO)
        {
          void* user_data = nullptr;
          IDAGetUserDataB(ida_mem, which, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from IDAS memory");
          auto fntable = static_cast<idasa_user_supplied_fn_table*>(user_data);
          fntable->resQB = nb::cast(resQB);
          return IDAQuadInitB(ida_mem, which, idas_resQB_wrapper, yQBO);
        });

  using IDAQuadRhsStdFnBS = int(sunrealtype t, N_Vector yy,
                                std::vector<N_Vector> yS, N_Vector yyB,
                                N_Vector qBdot, void* user_dataB);
  m.def("IDAQuadInitBS",
        [](void* ida_mem, int which, std::function<IDAQuadRhsStdFnBS> resQBs,
           N_Vector yQBO)
        {
          void* user_data = nullptr;
          IDAGetUserDataB(ida_mem, which, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from IDAS memory");
          auto fntable = static_cast<idasa_user_supplied_fn_table*>(user_data);
          fntable->resQBs = nb::cast(resQBs);
          return IDAQuadInitBS(ida_mem, which, idas_resQBs_wrapper, yQBO);
        });

  BIND_IDAB_CALLBACK(IDASetJacFnB, IDALsJacFnB, lsjacfnB, idas_lsjacfnB_wrapper);

  BIND_IDAB_CALLBACK2(IDASetPreconditionerB, IDALsPrecSetupFnB, lsprecsetupfnB,
                      idas_lsprecsetupfnB_wrapper, IDALsPrecSolveFnB,
                      lsprecsolvefnB, idas_lsprecsolvefnB_wrapper,
                      nb::arg("ida_mem"), nb::arg("which"),
                      nb::arg("psetupB").none(), nb::arg("psolveB"));

  BIND_IDAB_CALLBACK2(IDASetJacTimesB, IDALsJacTimesSetupFnB, lsjactimessetupfnB,
                      idas_lsjactimessetupfnB_wrapper, IDALsJacTimesVecFnB,
                      lsjactimesvecfnB, idas_lsjactimesvecfnB_wrapper,
                      nb::arg("ida_mem"), nb::arg("which"),
                      nb::arg("jsetupB").none(), nb::arg("jtimesB"));

  using IDALsJacStdFnBS = int(sunrealtype t, N_Vector yy,
                              std::vector<N_Vector> yS, N_Vector yyB,
                              N_Vector fyB, SUNMatrix JB, void* user_dataB,
                              N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
  BIND_IDAB_CALLBACK(IDASetJacFnBS, IDALsJacStdFnBS, lsjacfnBS,
                     idas_lsjacfnBS_wrapper);

  using IDALsPrecSetupStdFnBS = int(sunrealtype t, N_Vector yy, N_Vector yyB,
                                    N_Vector fyB, sunbooleantype jokB,
                                    sunbooleantype * jcurPtrB,
                                    sunrealtype gammaB, void* user_dataB);
  using IDALsPrecSolveStdFnBS =
    int(sunrealtype t, N_Vector yy, std::vector<N_Vector> yS, N_Vector yyB,
        N_Vector fyB, N_Vector rB, N_Vector zB, sunrealtype gammaB,
        sunrealtype deltaB, int lrB, void* user_dataB);
  BIND_IDAB_CALLBACK2(IDASetPreconditionerBS, IDALsPrecSetupStdFnBS,
                      lsprecsetupfnBS, idas_lsprecsetupfnBS_wrapper,
                      IDALsPrecSolveStdFnBS, lsprecsolvefnBS,
                      idas_lsprecsolvefnBS_wrapper, nb::arg("ida_mem"),
                      nb::arg("which"), nb::arg("psetupBS").none(),
                      nb::arg("psolveBS"));

  using IDALsJacTimesSetupStdFnBS = int(sunrealtype t, N_Vector yy,
                                        std::vector<N_Vector> yS, N_Vector yyB,
                                        N_Vector fyB, void* jac_dataB);
  using IDALsJacTimesVecStdFnBS = int(N_Vector vB, N_Vector JvB, sunrealtype t,
                                      N_Vector yy, std::vector<N_Vector> yS,
                                      N_Vector yyB, N_Vector fyB,
                                      void* jac_dataB, N_Vector tmpB);
  BIND_IDAB_CALLBACK2(IDASetJacTimesBS, IDALsJacTimesSetupStdFnBS,
                      lsjactimessetupfnBS, idas_lsjactimessetupfnBS_wrapper,
                      IDALsJacTimesVecStdFnBS, lsjactimesvecfnBS,
                      idas_lsjactimesvecfnBS_wrapper, nb::arg("ida_mem"),
                      nb::arg("which"), nb::arg("jsetupBS").none(),
                      nb::arg("jtimesBS"));
}
