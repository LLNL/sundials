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
// #include <sundials/sundials_stepper.hpp>

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

  m.def("ARKodeSetPostprocessStepFn",
        [](void* ark_mem, std::function<std::remove_pointer_t<ARKPostProcessFn>> fn)
        {
          void* user_data = nullptr;

          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
          {
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          }

          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);

          // Set the user-supplied function
          fntable->postprocessstepfn = nb::cast(fn);
          return ARKodeSetPostprocessStepFn(ark_mem,
                                            &arkode_postprocessstepfn_wrapper);
        });

  m.def("ARKodeSetPostprocessStageFn",
        [](void* ark_mem, std::function<std::remove_pointer_t<ARKPostProcessFn>> fn)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->postprocessstagefn = nb::cast(fn);
          return ARKodeSetPostprocessStageFn(ark_mem,
                                             &arkode_postprocessstagefn_wrapper);
        });

  m.def("ARKodeSetStagePredictFn",
        [](void* ark_mem,
           std::function<std::remove_pointer_t<ARKStagePredictFn>> fn)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->stagepredictfn = nb::cast(fn);
          return ARKodeSetStagePredictFn(ark_mem, &arkode_stagepredictfn_wrapper);
        });

  m.def("ARKodeWFtolerances",
        [](void* ark_mem, std::function<std::remove_pointer_t<ARKEwtFn>> fn)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->ewtn = nb::cast(fn);
          return ARKodeWFtolerances(ark_mem, &arkode_ewtfn_wrapper);
        });

  m.def("ARKodeSetNlsRhsFn",
        [](void* ark_mem, std::function<std::remove_pointer_t<ARKRhsFn>> fn)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->nlsfi = nb::cast(fn);
          return ARKodeSetNlsRhsFn(ark_mem, &arkode_nlsrhsfn_wrapper);
        });

  m.def("ARKodeSetJacFn",
        [](void* ark_mem, std::function<std::remove_pointer_t<ARKLsJacFn>> fn)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->lsjacfn = nb::cast(fn);
          return ARKodeSetJacFn(ark_mem, &arkode_lsjacfn_wrapper);
        });

  m.def("ARKodeSetMassFn",
        [](void* ark_mem, std::function<std::remove_pointer_t<ARKLsMassFn>> fn)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->lsmassfn = nb::cast(fn);
          return ARKodeSetMassFn(ark_mem, &arkode_lsmassfn_wrapper);
        });

  m.def("ARKodeSetPreconditioner",
        [](void* ark_mem,
           std::function<std::remove_pointer_t<ARKLsPrecSetupFn>> psetup,
           std::function<std::remove_pointer_t<ARKLsPrecSolveFn>> psolve)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->lsprecsetupfn = nb::cast(psetup);
          fntable->lsprecsolvefn = nb::cast(psolve);
          return ARKodeSetPreconditioner(ark_mem, &arkode_lsprecsetupfn_wrapper,
                                         &arkode_lsprecsolvefn_wrapper);
        });

  m.def("ARKodeSetMassPreconditioner",
        [](void* ark_mem,
           std::function<std::remove_pointer_t<ARKLsMassPrecSetupFn>> psetup,
           std::function<std::remove_pointer_t<ARKLsMassPrecSolveFn>> psolve)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->lsmassprecsetupfn = nb::cast(psetup);
          fntable->lsmassprecsolvefn = nb::cast(psolve);
          return ARKodeSetMassPreconditioner(ark_mem,
                                             &arkode_lsmassprecsetupfn_wrapper,
                                             &arkode_lsmassprecsolvefn_wrapper);
        });

  m.def("ARKodeSetJacTimes",
        [](void* ark_mem,
           std::function<std::remove_pointer_t<ARKLsJacTimesSetupFn>> jtsetup,
           std::function<std::remove_pointer_t<ARKLsJacTimesVecFn>> jtimes)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->lsjactimessetupfn = nb::cast(jtsetup);
          fntable->lsjactimesvecfn   = nb::cast(jtimes);
          return ARKodeSetJacTimes(ark_mem, &arkode_lsjactimessetupfn_wrapper,
                                   &arkode_lsjactimesvecfn_wrapper);
        });

  m.def("ARKodeSetJacTimesRhsFn",
        [](void* ark_mem, std::function<std::remove_pointer_t<ARKRhsFn>> fn)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->lsjacrhsfn = nb::cast(fn);
          return ARKodeSetJacTimesRhsFn(ark_mem, &arkode_lsjacrhsfn_wrapper);
        });

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

  m.def("ARKodeSetLinSysFn",
        [](void* ark_mem, std::function<std::remove_pointer_t<ARKLsLinSysFn>> fn)
        {
          void* user_data = nullptr;
          ARKodeGetUserData(ark_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from ARKODE memory");
          auto fntable = static_cast<arkode_user_supplied_fn_table*>(user_data);
          fntable->lslinsysfn = nb::cast(fn);
          return ARKodeSetLinSysFn(ark_mem, &arkode_lslinsysfn_wrapper);
        });

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
