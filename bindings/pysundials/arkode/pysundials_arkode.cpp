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
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>

#include <sundials/sundials_core.hpp>

class SUNStepperView;

#include <arkode/arkode.h>
#include <arkode/arkode.hpp>
#include <arkode/arkode_butcher.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_sprk.hpp>

#include "arkode/arkode_impl.h"

#include "sundials_adjointcheckpointscheme_impl.h"

#include "pysundials_arkode_usersupplied.hpp"

namespace nb = nanobind;

using namespace sundials::experimental;

// Forward declarations of functions defined in other translation units
void bind_arkode_erkstep(nb::module_& m);
void bind_arkode_arkstep(nb::module_& m);
void bind_arkode_sprkstep(nb::module_& m);
void bind_arkode_lsrkstep(nb::module_& m);
void bind_arkode_mristep(nb::module_& m);
void bind_arkode_forcingstep(nb::module_& m);
void bind_arkode_splittingstep(nb::module_& m);

// ARKODE callback binding macros
#define BIND_ARKODE_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)               \
  m.def(                                                                        \
    #NAME,                                                                      \
    [](void* ark_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn)         \
    {                                                                           \
      void* user_data = nullptr;                                                \
      ARKodeGetUserData(ark_mem, &user_data);                                   \
      if (!user_data)                                                           \
        throw std::runtime_error(                                               \
          "Failed to get Python function table from ARKODE memory");            \
      auto fntable    = static_cast<arkode_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER = nb::cast(fn);                                           \
      if (fn) { return NAME(ark_mem, &WRAPPER); }                               \
      else { return NAME(ark_mem, nullptr); }                                   \
    },                                                                          \
    __VA_ARGS__)

#define BIND_ARKODE_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,       \
                              MEMBER2, WRAPPER2, ...)                            \
  m.def(                                                                         \
    #NAME,                                                                       \
    [](void* ark_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1,        \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                       \
    {                                                                            \
      void* user_data = nullptr;                                                 \
      ARKodeGetUserData(ark_mem, &user_data);                                    \
      if (!user_data)                                                            \
        throw std::runtime_error(                                                \
          "Failed to get Python function table from ARKODE memory");             \
      auto fntable     = static_cast<arkode_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER1 = nb::cast(fn1);                                          \
      fntable->MEMBER2 = nb::cast(fn2);                                          \
      if (fn1) { return NAME(ark_mem, &WRAPPER1, &WRAPPER2); }                   \
      else { return NAME(ark_mem, nullptr, &WRAPPER2); }                         \
    },                                                                           \
    __VA_ARGS__)

void bind_arkode(nb::module_& m)
{
#include "pysundials_arkode_generated.hpp"

  /////////////////////////////////////////////////////////////////////////////
  // Interface view classes for ARKODE level objects
  /////////////////////////////////////////////////////////////////////////////

  nb::class_<ARKodeView>(m, "ARKodeView")
    .def_static("Create", &ARKodeView::Create<void*>)
    .def("get", nb::overload_cast<>(&ARKodeView::get, nb::const_),
         nb::rv_policy::reference);

  nb::class_<ARKodeButcherTableView>(m, "ARKodeButcherTableView")
    .def_static("Create",
                [](int s, int q, int p, std::vector<sunrealtype> c,
                   std::vector<sunrealtype> A, std::vector<sunrealtype> b,
                   std::vector<sunrealtype> d)
                {
                  ARKodeButcherTableView::Create(s, q, p, c.data(), A.data(),
                                                 b.data(), d.data());
                })
    .def("get", nb::overload_cast<>(&ARKodeButcherTableView::get, nb::const_),
         nb::rv_policy::reference);

  nb::class_<ARKodeSPRKTableView>(m, "ARKodeSPRKTableView")
    .def_static("Create", [](int s, int q, std::vector<sunrealtype> a,
                             std::vector<sunrealtype> ahat)
                { ARKodeSPRKTableView::Create(s, q, a.data(), ahat.data()); })
    .def("get", nb::overload_cast<>(&ARKodeSPRKTableView::get, nb::const_),
         nb::rv_policy::reference);

  /////////////////////////////////////////////////////////////////////////////
  // ARKODE user-supplied function setters
  /////////////////////////////////////////////////////////////////////////////

  BIND_ARKODE_CALLBACK(ARKodeSetPostprocessStepFn, ARKPostProcessFn,
                       postprocessstepfn, arkode_postprocessstepfn_wrapper,
                       nb::arg("arkode_mem"), nb::arg("postprocessstep").none());

  BIND_ARKODE_CALLBACK(ARKodeSetPostprocessStageFn, ARKPostProcessFn,
                       postprocessstagefn, arkode_postprocessstagefn_wrapper,
                       nb::arg("arkode_mem"), nb::arg("postprocessstage").none());

  BIND_ARKODE_CALLBACK(ARKodeSetStagePredictFn, ARKStagePredictFn,
                       stagepredictfn, arkode_stagepredictfn_wrapper,
                       nb::arg("arkode_mem"), nb::arg("stagepredict").none());

  BIND_ARKODE_CALLBACK(ARKodeWFtolerances, ARKEwtFn, ewtn, arkode_ewtfn_wrapper,
                       nb::arg("arkode_mem"), nb::arg("efun").none());

  BIND_ARKODE_CALLBACK(ARKodeSetNlsRhsFn, ARKRhsFn, nlsfi,
                       arkode_nlsrhsfn_wrapper, nb::arg("arkode_mem"),
                       nb::arg("nls_fi").none());

  BIND_ARKODE_CALLBACK(ARKodeSetJacFn, ARKLsJacFn, lsjacfn,
                       arkode_lsjacfn_wrapper, nb::arg("arkode_mem"),
                       nb::arg("jac").none());

  BIND_ARKODE_CALLBACK(ARKodeSetMassFn, ARKLsMassFn, lsmassfn,
                       arkode_lsmassfn_wrapper, nb::arg("arkode_mem"),
                       nb::arg("mass").none());

  BIND_ARKODE_CALLBACK2(ARKodeSetPreconditioner, ARKLsPrecSetupFn,
                        lsprecsetupfn, arkode_lsprecsetupfn_wrapper,
                        ARKLsPrecSolveFn, lsprecsolvefn,
                        arkode_lsprecsolvefn_wrapper, nb::arg("arkode_mem"),
                        nb::arg("psetup").none(), nb::arg("psolve").none());

  BIND_ARKODE_CALLBACK2(ARKodeSetMassPreconditioner, ARKLsMassPrecSetupFn,
                        lsmassprecsetupfn, arkode_lsmassprecsetupfn_wrapper,
                        ARKLsMassPrecSolveFn, lsmassprecsolvefn,
                        arkode_lsmassprecsolvefn_wrapper, nb::arg("arkode_mem"),
                        nb::arg("psetup").none(), nb::arg("psolve").none());

  BIND_ARKODE_CALLBACK2(ARKodeSetJacTimes, ARKLsJacTimesSetupFn,
                        lsjactimessetupfn, arkode_lsjactimessetupfn_wrapper,
                        ARKLsJacTimesVecFn, lsjactimesvecfn,
                        arkode_lsjactimesvecfn_wrapper, nb::arg("arkode_mem"),
                        nb::arg("jtsetup").none(), nb::arg("jtimes").none());

  BIND_ARKODE_CALLBACK(ARKodeSetJacTimesRhsFn, ARKRhsFn, lsjacrhsfn,
                       arkode_lsjacrhsfn_wrapper, nb::arg("arkode_mem"),
                       nb::arg("jtimesRhsFn").none());

  m.def(
    "ARKodeSetMassTimes",
    [](void* ark_mem,
       std::function<std::remove_pointer_t<ARKLsMassTimesSetupFn>> msetup,
       std::function<std::remove_pointer_t<ARKLsMassTimesVecFn>> mtimes)
    {
      void* user_data = nullptr;
      ARKodeGetUserData(ark_mem, &user_data);
      if (!user_data)
        throw std::runtime_error(
          "Failed to get Python function table from ARKODE memory");
      auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
      fntable->lsmasstimessetupfn = nb::cast(msetup);
      fntable->lsmasstimesvecfn   = nb::cast(mtimes);
      return ARKodeSetMassTimes(ark_mem, &arkode_lsmasstimessetupfn_wrapper,
                                &arkode_lsmasstimesvecfn_wrapper, nullptr);
    },
    nb::arg("arkode_mem"), nb::arg("msetup").none(), nb::arg("mtimes").none());

  BIND_ARKODE_CALLBACK(ARKodeSetLinSysFn, ARKLsLinSysFn, lslinsysfn,
                       arkode_lslinsysfn_wrapper, nb::arg("arkode_mem"),
                       nb::arg("linsys").none());

  // ARKodeSetMassTimes doesnt fit the BIND_ARKODE_CALLBACK macro pattern(s)
  // due to the 4th argument for user data, so we just write it out explicitly.
  m.def(
    "ARKodeSetMassTimes",
    [](void* ark_mem,
       std::function<std::remove_pointer_t<ARKLsMassTimesSetupFn>> msetup,
       std::function<std::remove_pointer_t<ARKLsMassTimesVecFn>> mtimes)
    {
      void* user_data = nullptr;
      ARKodeGetUserData(ark_mem, &user_data);
      if (!user_data)
        throw std::runtime_error(
          "Failed to get Python function table from ARKODE memory");
      auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
      fntable->lsmasstimessetupfn = nb::cast(msetup);
      fntable->lsmasstimesvecfn   = nb::cast(mtimes);
      return ARKodeSetMassTimes(ark_mem, &arkode_lsmasstimessetupfn_wrapper,
                                &arkode_lsmasstimesvecfn_wrapper, nullptr);
    },
    nb::arg("ark_mem"), nb::arg("msetup").none(), nb::arg("mtimes").none());

  /////////////////////////////////////////////////////////////////////////////
  // Additional functions that litgen cannot generate
  /////////////////////////////////////////////////////////////////////////////

  // This function has optional arguments which litgen cannot deal with because they are followed by non-optional arguments.
  m.def("ARKodeSetMassLinearSolver", ARKodeSetMassLinearSolver,
        nb::arg("arkode_mem"), nb::arg("LS"), nb::arg("M").none(),
        nb::arg("time_dep"));

  m.def(
    "ARKodeCreateSUNStepper",
    [](void* ark_mem)
    {
      SUNStepper stepper = nullptr;

      int status = ARKodeCreateSUNStepper(ark_mem, &stepper);

      return std::make_tuple(status, stepper);
    },
    nb::rv_policy::reference);

  m.def(
    "ARKodeGetJac",
    [](void* ark_mem)
    {
      SUNMatrix J = nullptr;

      int status = ARKodeGetJac(ark_mem, &J);

      return std::make_tuple(status, J);
    },
    "WARNING: this function returns a SUNMatrix reference, DO NOT WRAP IT IN A "
    "`SUNMatrixView`. Doing so will result in a double free or worse.",
    nb::rv_policy::reference);

  m.def(
    "ARKodeGetCurrentMassMatrix",
    [](void* ark_mem)
    {
      SUNMatrix M = nullptr;

      int status = ARKodeGetJac(ark_mem, &M);

      return std::make_tuple(status, M);
    },
    "WARNING: this function returns a SUNMatrix reference, DO NOT WRAP IT IN A "
    "`SUNMatrixView`. Doing so will result in a double free or worse.",
    nb::rv_policy::reference);

  m.def(
    "ARKodeGetCurrentState",
    [](void* ark_mem)
    {
      N_Vector state = nullptr;

      int status = ARKodeGetCurrentState(ark_mem, &state);

      return std::make_tuple(status, state);
    },
    "WARNING: this function returns a N_Vector reference, DO NOT WRAP IT IN A "
    "`NVectorView`. Doing so will result in a double free or worse.",
    nb::rv_policy::reference);

  m.def(
    "ARKodeSPRKTable_ToButcher",
    [](ARKodeSPRKTable sprk_table)
    {
      ARKodeButcherTable a = nullptr;
      ARKodeButcherTable b = nullptr;

      int status = ARKodeSPRKTable_ToButcher(sprk_table, &a, &b);

      return std::make_tuple(status, a, b);
    },
    nb::rv_policy::reference);

  bind_arkode_arkstep(m);
  bind_arkode_erkstep(m);
  bind_arkode_sprkstep(m);
  bind_arkode_lsrkstep(m);
  bind_arkode_mristep(m);
  bind_arkode_forcingstep(m);
  bind_arkode_splittingstep(m);
}
