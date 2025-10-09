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
#include <nanobind/stl/tuple.h>

#include <arkode/arkode.hpp>
#include <arkode/arkode_mristep.hpp>
#include <sundials/sundials_core.hpp>

#include "arkode_mristep_impl.h"

#include "arkode_usersupplied.hpp"

namespace nb = nanobind;

using namespace sundials::experimental;

void bind_arkode_mristep(nb::module_& m)
{
#include "arkode_mristep_generated.hpp"

  nb::class_<MRIStepCouplingView>(m, "MRIStepCouplingView")
    .def_static("Create", &MRIStepCouplingView::Create<MRIStepCouplingView>)
    .def("get", nb::overload_cast<>(&MRIStepCouplingView::get, nb::const_),
         nb::rv_policy::reference);

  // _MRIStepInnerStepper is a opaque/private class forward declared in a public header but
  // defined in a source file elsewhere. As such, we need to declare it here since its
  // not picked up in any header files by the generator.
  nb::class_<_MRIStepInnerStepper>(m, "_MRIStepInnerStepper");

  nb::class_<MRIStepInnerStepperView>(m, "MRIStepInnerStepperView")
    .def_static("Create", &MRIStepInnerStepperView::Create<MRIStepInnerStepper>)
    .def("get", nb::overload_cast<>(&MRIStepInnerStepperView::get, nb::const_),
         nb::rv_policy::reference);

  m.def(
    "MRIStepInnerStepper_Create",
    [](SUNContext sunctx)
    {
      MRIStepInnerStepper stepper = nullptr;

      int status      = MRIStepInnerStepper_Create(sunctx, &stepper);
      auto cb_fns     = mristepinnerstepper_user_supplied_fn_table_alloc();
      stepper->python = static_cast<void*>(cb_fns);

      return std::make_tuple(status, stepper);
    },
    nb::rv_policy::reference);

  m.def(
    "MRIStepInnerStepper_CreateFromSUNStepper",
    [](SUNStepper stepper)
    {
      MRIStepInnerStepper inner_stepper = nullptr;

      int status = MRIStepInnerStepper_CreateFromSUNStepper(stepper,
                                                            &inner_stepper);

      return std::make_tuple(status, stepper);
    },
    nb::rv_policy::reference);

  m.def(
    "ARKodeCreateMRIStepInnerStepper",
    [](void* inner_arkode_mem)
    {
      MRIStepInnerStepper stepper = nullptr;

      int status = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);

      return std::make_tuple(status, stepper);
    },
    nb::rv_policy::reference);

  m.def(
    "MRIStepCreate",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> fse,
       std::function<std::remove_pointer_t<ARKRhsFn>> fsi, sunrealtype t0,
       N_Vector y0, MRIStepInnerStepper stepper, SUNContext sunctx)
    {
      auto fse_wrapper = fse ? mristep_fse_wrapper : nullptr;
      auto fsi_wrapper = fsi ? mristep_fsi_wrapper : nullptr;

      void* ark_mem = MRIStepCreate(fse_wrapper, fsi_wrapper, t0, y0, stepper,
                                    sunctx);
      if (ark_mem == nullptr)
      {
        throw std::runtime_error("Failed to create ARKODE memory");
      }

      // Create the user-supplied function table to store the Python user functions
      auto cb_fns = arkode_user_supplied_fn_table_alloc();

      // Smuggle the user-supplied function table into callback wrappers through the user_data pointer
      int ark_status = ARKodeSetUserData(ark_mem, static_cast<void*>(cb_fns));
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error("Failed to set user data in ARKODE memory");
      }

      // Ensure ARKodeFree will free the user-supplied function table
      ark_status = arkSetOwnUserData(ark_mem, SUNTRUE);
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error(
          "Failed to set user data ownership in ARKODE memory");
      }

      // Finally, set the RHS function
      cb_fns->mristep_fse = nb::cast(fse);
      cb_fns->mristep_fsi = nb::cast(fsi);

      return ark_mem;
    },
    nb::arg("fse").none(), nb::arg("fsi").none(), nb::arg("t0"), nb::arg("y0"),
    nb::arg("inner_stepper"), nb::arg("sunctx"));

  m.def(
    "MRIStepGetCurrentCoupling",
    [](void* ark_mem)
    {
      MRIStepCoupling C = nullptr;

      int status = MRIStepGetCurrentCoupling(ark_mem, &C);

      return std::make_tuple(status, C);
    },
    "WARNING: this function returns a MRIStepCoupling reference, DO NOT WRAP "
    "IT IN A `MRIStepCouplingView`. Doing so will result in a double free "
    "or worse.",
    nb::rv_policy::reference);
}
