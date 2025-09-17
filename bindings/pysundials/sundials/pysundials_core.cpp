/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file defines the pysundials.core module. 
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include <sundials/sundials_core.hpp>
#include <sundials/sundials_futils.h>

namespace nb = nanobind;

void bind_nvector(nb::module_& m);
void bind_sunadaptcontroller(nb::module_& m);
void bind_sunadjointcheckpointscheme(nb::module_& m);
void bind_sunadjointstepper(nb::module_& m);
void bind_suncontext(nb::module_& m);
void bind_sundomeigestimator(nb::module_& m);
void bind_sunlinearsolver(nb::module_& m);
void bind_sunlogger(nb::module_& m);
void bind_sunmatrix(nb::module_& m);
void bind_sunmemory(nb::module_& m);
void bind_sunnonlinearsolver(nb::module_& m);
void bind_sunprofiler(nb::module_& m);
void bind_sunstepper(nb::module_& m);

void bind_core(nb::module_& m)
{
#include "pysundials_errors.hpp"
#include "pysundials_types_generated.hpp"

  // handle opening and closing C files
  nb::class_<FILE>(m, "FILE");
  m.def("SUNFileOpen",
        [](const char* filename, const char* modes)
        {
          FILE* tmp = nullptr;
          std::shared_ptr<FILE> fp;
          SUNErrCode status = SUNFileOpen(filename, modes, &tmp);
          if (status) { fp = nullptr; }
          else { fp = std::shared_ptr<FILE>(tmp, std::fclose); }
          return std::make_tuple(status, fp);
        });

  bind_nvector(m);
  bind_sunadaptcontroller(m);
  bind_sunadjointcheckpointscheme(m);
  bind_sunadjointstepper(m);
  bind_suncontext(m);
  bind_sundomeigestimator(m);
  bind_sunlinearsolver(m);
  bind_sunlogger(m);
  bind_sunmatrix(m);
  bind_sunmemory(m);
  bind_sunnonlinearsolver(m);
  bind_sunprofiler(m);
  bind_sunstepper(m);
}