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
#include <arkode/arkode_erkstep.h>
#include <sundials/sundials_core.hpp>

#include "pysundials_arkode_usersupplied.hpp"

namespace nb = nanobind;

void bind_arkode_erkstep(nb::module_& m)
{
  m.def("ERKStepCreate",
        [](std::function<ARKRhsStdFn> rhs, sunrealtype t0, N_Vector y0,
           SUNContext sunctx)
        {
          void* ark_mem = ERKStepCreate(erk_rhsfn_wrapper, t0, y0, sunctx);
          if (ark_mem == nullptr)
          {
            throw std::runtime_error("Failed to create ARKODE memory");
          }

          // Create the user-supplied function table to store the Python user functions
          // NOTE: we must use malloc since ARKodeFree calls free
          erkstep_user_supplied_fn_table* cb_fns =
            static_cast<erkstep_user_supplied_fn_table*>(
              malloc(sizeof(erkstep_user_supplied_fn_table)));
         
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
          cb_fns->erkstep_rhsfn = nb::cast(rhs);

          return ark_mem;
        });
}
