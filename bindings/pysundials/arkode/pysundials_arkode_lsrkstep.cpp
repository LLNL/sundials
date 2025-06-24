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

#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

#include <arkode/arkode.hpp>
#include <arkode/arkode_lsrkstep.h>
#include <sundials/sundials_core.hpp>

#include "pysundials_arkode_usersupplied.hpp"

namespace nb = nanobind;

void bind_arkode_lsrkstep(nb::module_& m)
{
#include "pysundials_arkode_lsrkstep_generated.hpp"

  m.def(
    "LSRKStepCreateSTS",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> rhs, sunrealtype t0, N_Vector y0, SUNContext sunctx)
    {
      auto rhs_wrapper = rhs ? erkstep_f_wrapper : nullptr;
      void* ark_mem = LSRKStepCreateSTS(rhs_wrapper, t0, y0, sunctx);
      if (ark_mem == nullptr)
      {
        throw std::runtime_error("Failed to create LSRKStep memory");
      }
      auto cb_fns = arkode_user_supplied_fn_table_alloc();
      int ark_status = ARKodeSetUserData(ark_mem, static_cast<void*>(cb_fns));
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error("Failed to set user data in LSRKStep memory");
      }
      ark_status = ARKodeSetOwnUserData(ark_mem, SUNTRUE);
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error("Failed to set user data ownership in LSRKStep memory");
      }
      cb_fns->lsrkstep_rhs = nb::cast(rhs);
      return ark_mem;
    },
    nb::arg("rhs"), nb::arg("t0"), nb::arg("y0"), nb::arg("sunctx"));

  m.def(
    "LSRKStepCreateSSP",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> rhs, sunrealtype t0, N_Vector y0, SUNContext sunctx)
    {
      auto rhs_wrapper = rhs ? erkstep_f_wrapper : nullptr;
      void* ark_mem = LSRKStepCreateSSP(rhs_wrapper, t0, y0, sunctx);
      if (ark_mem == nullptr)
      {
        throw std::runtime_error("Failed to create LSRKStep memory");
      }
      auto cb_fns = arkode_user_supplied_fn_table_alloc();
      int ark_status = ARKodeSetUserData(ark_mem, static_cast<void*>(cb_fns));
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error("Failed to set user data in LSRKStep memory");
      }
      ark_status = ARKodeSetOwnUserData(ark_mem, SUNTRUE);
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error("Failed to set user data ownership in LSRKStep memory");
      }
      cb_fns->lsrkstep_rhs = nb::cast(rhs);
      return ark_mem;
    },
    nb::arg("rhs"), nb::arg("t0"), nb::arg("y0"), nb::arg("sunctx"));
}
