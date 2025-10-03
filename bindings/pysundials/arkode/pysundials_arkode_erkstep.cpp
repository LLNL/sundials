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
#include <arkode/arkode_erkstep.h>
#include <sundials/sundials_core.hpp>

#include "pysundials_arkode_usersupplied.hpp"

namespace nb = nanobind;

void bind_arkode_erkstep(nb::module_& m)
{
#include "pysundials_arkode_erkstep_generated.hpp"

  m.def("ERKStepCreate",
        [](std::function<std::remove_pointer_t<ARKRhsFn>> rhs, sunrealtype t0,
           N_Vector y0, SUNContext sunctx)
        {
          void* ark_mem = ERKStepCreate(erkstep_f_wrapper, t0, y0, sunctx);
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
            throw std::runtime_error(
              "Failed to set user data in ARKODE memory");
          }

          // Ensure ARKodeFree will free the user-supplied function table
          ark_status = ARKodeSetOwnUserData(ark_mem, SUNTRUE);
          if (ark_status != ARK_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error(
              "Failed to set user data ownership in ARKODE memory");
          }

          // Finally, set the RHS function
          cb_fns->erkstep_f = nb::cast(rhs);

          return ark_mem;
        });

  m.def(
    "ERKStepCreateAdjointStepper",
    [](void* arkode_mem,
       std::function<std::remove_pointer_t<SUNAdjRhsFn>> adj_f, sunrealtype tf,
       N_Vector sf, SUNContext sunctx) -> std::tuple<int, SUNAdjointStepper>
    {
      auto f_wrapper = adj_f ? erkstep_adjf_wrapper : nullptr;

      SUNAdjointStepper adj_stepper = nullptr;
      int ark_status = ERKStepCreateAdjointStepper(arkode_mem, f_wrapper, tf,
                                                   sf, sunctx, &adj_stepper);
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

      cb_fns->erkstep_adjf = nb::cast(adj_f);

      return std::make_tuple(ark_status, adj_stepper);
    },
    nb::arg("arkode_mem"), nb::arg("adj_f").none(), nb::arg("tf"),
    nb::arg("sf"), nb::arg("sunctx"));

  m.def(
    "ERKStepGetCurrentButcherTable",
    [](void* ark_mem)
    {
      ARKodeButcherTable fe = nullptr;

      int status = ERKStepGetCurrentButcherTable(ark_mem, &fe);

      return std::make_tuple(status, fe);
    },
    "WARNING: this function returns ARKodeButcherTable references, DO NOT WRAP "
    "THEM IN A `ARKodeButcherTableView`. Doing so will result in a double free "
    "or worse.",
    nb::rv_policy::reference);
}
