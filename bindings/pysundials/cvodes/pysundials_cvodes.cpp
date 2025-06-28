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

#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <sundials/sundials_core.hpp>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes.hpp>
#include <cvodes/cvodes_ls.h>

#include "sundials_adjointcheckpointscheme_impl.h"

#include "pysundials_cvodes_usersupplied.hpp"

namespace nb = nanobind;

using namespace sundials::experimental;

// Forward declarations of functions defined in other translation units

void bind_cvodes(nb::module_& m)
{
#include "pysundials_cvodes_generated.hpp"

  nb::class_<CVodeView>(m, "CVodeView")
    .def_static("Create", &CVodeView::Create<void*>)
    .def("get", nb::overload_cast<>(&CVodeView::get, nb::const_),
         nb::rv_policy::reference);

  m.def("CVodeInit",
        [](void* cv_mem, std::function<std::remove_pointer_t<CVRhsFn>> rhs,
           sunrealtype t0, N_Vector y0)
        {
          int cv_status = CVodeInit(cv_mem, cvode_f_wrapper, t0, y0);

          // Create the user-supplied function table to store the Python user functions
          auto cb_fns = cvode_user_supplied_fn_table_alloc();

          // Smuggle the user-supplied function table into callback wrappers through the user_data pointer
          cv_status = CVodeSetUserData(cv_mem, static_cast<void*>(cb_fns));
          if (cv_status != CV_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error("Failed to set user data in CVODE memory");
          }

          // Ensure CVodeFree will free the user-supplied function table
          cv_status = CVodeSetOwnUserData(cv_mem, SUNTRUE);
          if (cv_status != CV_SUCCESS)
          {
            free(cb_fns);
            throw std::runtime_error(
              "Failed to set user data ownership in CVODE memory");
          }

          // Finally, set the RHS function
          cb_fns->f = nb::cast(rhs);

          return cv_status;
        });

  // We cant use std::function<std::remove_pointer_t<CVQuadSensRhsFn>>
  // because we need to replace the N_Vector* arguments with std::vector<N_Vector>.
  using CVQuadSensRhsStdFn = int(sunrealtype t, N_Vector y,
                                 std::vector<N_Vector> yQS, N_Vector ydot,
                                 std::vector<N_Vector> yQSdot, void* user_data);
  m.def("CVodeQuadSensInit",
        [](void* cv_mem, std::function<CVQuadSensRhsStdFn> fQS,
           std::vector<N_Vector> yQS0)
        {
          void* user_data = nullptr;
          CVodeGetUserData(cv_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvode_user_supplied_fn_table*>(user_data);
          fntable->fQS = nb::cast(fQS);
          return CVodeQuadSensInit(cv_mem, cvode_fqs_wrapper, yQS0.data());
        });

  using CVSensRhsStdFn = int(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                             int iS, std::vector<N_Vector> yS,
                             std::vector<N_Vector> ySdot, void* user_data);
  m.def("CVodeSensInit",
        [](void* cv_mem, int Ns, int ism, std::function<CVSensRhsStdFn> fS,
           std::vector<N_Vector> yS0)
        {
          void* user_data = nullptr;
          CVodeGetUserData(cv_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvode_user_supplied_fn_table*>(user_data);
          fntable->fS  = nb::cast(fS);
          return CVodeSensInit(cv_mem, Ns, ism, cvode_fs_wrapper, yS0.data());
        });

  using CVSensRhs1StdFn = int(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                              int iS, N_Vector yS, N_Vector ySdot,
                              void* user_data, N_Vector tmp1, N_Vector tmp2);
  m.def("CVodeSensInit1",
        [](void* cv_mem, int Ns, std::function<CVSensRhs1StdFn> fS1, int ism,
           std::vector<N_Vector> yS0)
        {
          void* user_data = nullptr;
          CVodeGetUserData(cv_mem, &user_data);
          if (!user_data)
            throw std::runtime_error(
              "Failed to get Python function table from CVODE memory");
          auto fntable = static_cast<cvode_user_supplied_fn_table*>(user_data);
          fntable->fS1 = nb::cast(fS1);
          // void* cvode_mem, int Ns, int ism,
          //  CVSensRhs1Fn fS1, N_Vector* yS0
          return CVodeSensInit1(cv_mem, Ns, ism, cvode_fs1_wrapper, yS0.data());
        });
}
