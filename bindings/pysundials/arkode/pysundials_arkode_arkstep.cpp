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
#include <arkode/arkode_arkstep.h>
#include <sundials/sundials_core.hpp>

#include "pysundials_arkode_usersupplied.hpp"

namespace nb = nanobind;

void bind_arkode_arkstep(nb::module_& m)
{
  #include "pysundials_arkode_arkstep_generated.hpp"

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
      // NOTE: we must use malloc since ARKodeFree calls free
      arkstep_user_supplied_fn_table* cb_fns =
        static_cast<arkstep_user_supplied_fn_table*>(
          malloc(sizeof(arkstep_user_supplied_fn_table)));

      // Smuggle the user-supplied function table into callback wrappers through the user_data pointer
      int ark_status = ARKodeSetUserData(ark_mem, static_cast<void*>(cb_fns));
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error("Failed to set user data in ARKODE memory");
      }

      // Ensure ARKodeFree will free the user-supplied function table
      ark_status = ARKodeSetOwnUserData(ark_mem, SUNTRUE);
      if (ark_status != ARK_SUCCESS)
      {
        free(cb_fns);
        throw std::runtime_error(
          "Failed to set user data ownership in ARKODE memory");
      }

      // Finally, set the RHS functions
      cb_fns->fe = nb::cast(fe);
      cb_fns->fi = nb::cast(fi);

      return ark_mem;
    }, // .none() must be added to functions that accept nullptr as a valid argument
    nb::arg("fe").none(), nb::arg("fi").none(), nb::arg("t0"), nb::arg("y0"),
    nb::arg("sunctx"));
}
