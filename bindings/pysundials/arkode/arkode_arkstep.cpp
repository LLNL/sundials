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

#include <arkode/arkode.hpp>
#include <arkode/arkode_arkstep.h>
#include <sundials/sundials_core.hpp>

#include "arkode_usersupplied.hpp"

namespace nb = nanobind;

namespace pysundials {

void bind_arkode_arkstep(nb::module_& m)
{
#include "arkode_arkstep_generated.hpp"

  m.def(
    "ARKStepCreate",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> fe,
       std::function<std::remove_pointer_t<ARKRhsFn>> fi, sunrealtype t0,
       N_Vector y0, SUNContext sunctx)
    {
      auto fe_wrapper = fe ? arkstep_fe_wrapper : nullptr;
      auto fi_wrapper = fi ? arkstep_fi_wrapper : nullptr;

      void* ark_mem = ARKStepCreate(fe_wrapper, fi_wrapper, t0, y0, sunctx);
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

      // Finally, set the RHS functions
      cb_fns->arkstep_fe = nb::cast(fe);
      cb_fns->arkstep_fi = nb::cast(fi);

      return ark_mem;
    }, // .none() must be added to functions that accept nullptr as a valid argument
    nb::arg("fe").none(), nb::arg("fi").none(), nb::arg("t0"), nb::arg("y0"),
    nb::arg("sunctx"));

  m.def(
    "ARKStepCreateAdjointStepper",
    [](void* arkode_mem, std::function<std::remove_pointer_t<SUNAdjRhsFn>> adj_fe,
       std::function<std::remove_pointer_t<SUNAdjRhsFn>> adj_fi, sunrealtype tf,
       N_Vector sf, SUNContext sunctx) -> std::tuple<int, SUNAdjointStepper>
    {
      auto fe_wrapper = adj_fe ? arkstep_adjfe_wrapper : nullptr;
      auto fi_wrapper = adj_fi ? arkstep_adjfi_wrapper : nullptr;

      SUNAdjointStepper adj_stepper = nullptr;
      int ark_status = ARKStepCreateAdjointStepper(arkode_mem, fe_wrapper,
                                                   fi_wrapper, tf, sf, sunctx,
                                                   &adj_stepper);
      if (ark_status != ARK_SUCCESS)
      {
        throw std::runtime_error(
          "Failed to create adjoint stepper in py-sundials memory");
      }

      // Finally, set the RHS functions
      void* user_data = nullptr;
      ark_status      = ARKodeGetUserData(arkode_mem, &user_data);
      if (ark_status != ARK_SUCCESS)
      {
        throw std::runtime_error("Failed to extract ARKODE user data");
      }

      auto cb_fns = static_cast<arkode_user_supplied_fn_table*>(user_data);

      cb_fns->arkstep_adjfe = nb::cast(adj_fe);
      cb_fns->arkstep_adjfi = nb::cast(adj_fi);

      return std::make_tuple(ark_status, adj_stepper);
    },
    nb::arg("arkode_mem"), nb::arg("adj_fe").none(), nb::arg("adj_fi").none(),
    nb::arg("tf"), nb::arg("sf"), nb::arg("sunctx"));

  m.def(
    "ARKStepGetCurrentButcherTables",
    [](void* ark_mem)
    {
      ARKodeButcherTable fe = nullptr;
      ARKodeButcherTable fi = nullptr;

      int status = ARKStepGetCurrentButcherTables(ark_mem, &fi, &fe);

      return std::make_tuple(status, fi, fe);
    },
    "WARNING: this function returns ARKodeButcherTable references, DO NOT WRAP "
    "THEM IN A `ARKodeButcherTableView`. Doing so will result in a double free "
    "or worse.",
    nb::rv_policy::reference);
}

} // namespace pysundials
