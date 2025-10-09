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
#include <arkode/arkode_sprkstep.h>
#include <sundials/sundials_core.hpp>

#include "arkode_usersupplied.hpp"

namespace nb = nanobind;

namespace pysundials {

void bind_arkode_sprkstep(nb::module_& m)
{
#include "arkode_sprkstep_generated.hpp"

  m.def(
    "SPRKStepCreate",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> f1,
       std::function<std::remove_pointer_t<ARKRhsFn>> f2, sunrealtype t0,
       N_Vector y0, SUNContext sunctx)
    {
      auto f1_wrapper = f1 ? sprkstep_f1_wrapper : nullptr;
      auto f2_wrapper = f2 ? sprkstep_f2_wrapper : nullptr;

      void* ark_mem = SPRKStepCreate(f1_wrapper, f2_wrapper, t0, y0, sunctx);
      if (ark_mem == nullptr)
      {
        throw std::runtime_error("Failed to create SPRKStep memory");
      }

      auto cb_fns    = arkode_user_supplied_fn_table_alloc();
      int ark_status = ARKodeSetUserData(ark_mem, static_cast<void*>(cb_fns));
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error("Failed to set user data in SPRKStep memory");
      }
      ark_status = arkSetOwnUserData(ark_mem, SUNTRUE);
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error(
          "Failed to set user data ownership in SPRKStep memory");
      }
      cb_fns->sprkstep_f1 = nb::cast(f1);
      cb_fns->sprkstep_f2 = nb::cast(f2);
      return ark_mem;
    },
    nb::arg("f1"), nb::arg("f2"), nb::arg("t0"), nb::arg("y0"),
    nb::arg("sunctx"));
}

} // namespace pysundials
