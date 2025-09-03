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
#include <stdexcept>

#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>

#include <sundials/sundials_core.hpp>

class SUNStepperView;

#include <arkode/arkode.h>
#include <arkode/arkode.hpp>
#include <arkode/arkode_ls.h>

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

// ARKODE callback binding macros (mirroring CVODES)
#define BIND_ARKODE_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER)                     \
  m.def(#NAME,                                                                   \
        [](void* ark_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn)      \
        {                                                                        \
          void* user_data = nullptr;                                             \
          ARKodeGetUserData(ark_mem, &user_data);                                \
          if (!user_data)                                                        \
            throw std::runtime_error(                                            \
              "Failed to get Python function table from ARKODE memory");         \
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data); \
          fntable->MEMBER = nb::cast(fn);                                        \
          return NAME(ark_mem, &WRAPPER);                                        \
        })

#define BIND_ARKODE_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,       \
                              MEMBER2, WRAPPER2)                                 \
  m.def(#NAME,                                                                   \
        [](void* ark_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1,    \
           std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                   \
        {                                                                        \
          void* user_data = nullptr;                                             \
          ARKodeGetUserData(ark_mem, &user_data);                                \
          if (!user_data)                                                        \
            throw std::runtime_error(                                            \
              "Failed to get Python function table from ARKODE memory");         \
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data); \
          fntable->MEMBER1 = nb::cast(fn1);                                      \
          fntable->MEMBER2 = nb::cast(fn2);                                      \
          return NAME(ark_mem, &WRAPPER1, &WRAPPER2);                            \
        })

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

  /////////////////////////////////////////////////////////////////////////////
  // ARKODE user-supplied function setters
  /////////////////////////////////////////////////////////////////////////////

  BIND_ARKODE_CALLBACK(ARKodeSetPostprocessStepFn, ARKPostProcessFn,
                       postprocessstepfn, arkode_postprocessstepfn_wrapper);
  BIND_ARKODE_CALLBACK(ARKodeSetPostprocessStageFn, ARKPostProcessFn,
                       postprocessstagefn, arkode_postprocessstagefn_wrapper);
  BIND_ARKODE_CALLBACK(ARKodeSetStagePredictFn, ARKStagePredictFn,
                       stagepredictfn, arkode_stagepredictfn_wrapper);
  BIND_ARKODE_CALLBACK(ARKodeWFtolerances, ARKEwtFn, ewtn, arkode_ewtfn_wrapper);
  BIND_ARKODE_CALLBACK(ARKodeSetNlsRhsFn, ARKRhsFn, nlsfi,
                       arkode_nlsrhsfn_wrapper);
  BIND_ARKODE_CALLBACK(ARKodeSetJacFn, ARKLsJacFn, lsjacfn,
                       arkode_lsjacfn_wrapper);
  BIND_ARKODE_CALLBACK(ARKodeSetMassFn, ARKLsMassFn, lsmassfn,
                       arkode_lsmassfn_wrapper);
  BIND_ARKODE_CALLBACK2(ARKodeSetPreconditioner, ARKLsPrecSetupFn, lsprecsetupfn,
                        arkode_lsprecsetupfn_wrapper, ARKLsPrecSolveFn,
                        lsprecsolvefn, arkode_lsprecsolvefn_wrapper);
  BIND_ARKODE_CALLBACK2(ARKodeSetMassPreconditioner, ARKLsMassPrecSetupFn,
                        lsmassprecsetupfn, arkode_lsmassprecsetupfn_wrapper,
                        ARKLsMassPrecSolveFn, lsmassprecsolvefn,
                        arkode_lsmassprecsolvefn_wrapper);
  BIND_ARKODE_CALLBACK2(ARKodeSetJacTimes, ARKLsJacTimesSetupFn,
                        lsjactimessetupfn, arkode_lsjactimessetupfn_wrapper,
                        ARKLsJacTimesVecFn, lsjactimesvecfn,
                        arkode_lsjactimesvecfn_wrapper);
  BIND_ARKODE_CALLBACK(ARKodeSetJacTimesRhsFn, ARKRhsFn, lsjacrhsfn,
                       arkode_lsjacrhsfn_wrapper);

  // ARKodeSetMassTimes doesnt fit the BIND_ARKODE_CALLBACK macro pattern(s)
  // due to the 4th argument for user data, so we just write it out explicitly.
  m.def("ARKodeSetMassTimes",
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
        });

  BIND_ARKODE_CALLBACK(ARKodeSetLinSysFn, ARKLsLinSysFn, lslinsysfn,
                       arkode_lslinsysfn_wrapper);

  /////////////////////////////////////////////////////////////////////////////
  // Additional functions that litgen cannot generate
  /////////////////////////////////////////////////////////////////////////////

  // These functions have optional arguments which litgen cannot deal with because they are followed by
  // non-optional arguments.
  m.def("ARKodeSetMassLinearSolver", ARKodeSetMassLinearSolver,
        nb::arg("arkode_mem"), nb::arg("LS"), nb::arg("M").none(),
        nb::arg("time_dep"));

  m.def(
    "ARKodeCreateSUNStepper",
    [](void* ark_mem)
    {
      SUNStepper stepper = nullptr;

      int status = ARKodeCreateSUNStepper(ark_mem, &stepper);
      if (status != ARK_SUCCESS)
      {
        throw std::runtime_error(
          "PY-SUNDIALS-OFFICIAL: ARKodeCreateSUNStepper call failed");
      }

      return stepper;
    },
    nb::arg("ark_mem"), nb::rv_policy::reference);

  bind_arkode_arkstep(m);
  bind_arkode_erkstep(m);
  bind_arkode_sprkstep(m);
  bind_arkode_lsrkstep(m);
  bind_arkode_mristep(m);
  bind_arkode_forcingstep(m);
  bind_arkode_splittingstep(m);
}
