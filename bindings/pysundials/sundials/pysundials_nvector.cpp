/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This file is the entrypoint for the Python binding code for the
 * SUNDIALS N_Vector class. It contains hand-written code for functions
 * that require special treatment, and includes the generated code
 * produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>

#include <sundials/sundials_nvector.hpp>

namespace nb = nanobind;

using namespace sundials::experimental;

void bind_nvector(nb::module_& m)
{
#include "pysundials_nvector_generated.hpp"

  nb::class_<NVectorView>(m, "NVectorView")
    .def_static("Create", &NVectorView::Create<N_Vector>)
    .def("get", nb::overload_cast<>(&NVectorView::get, nb::const_),
         nb::rv_policy::reference);

  // TODO(CJB): allow for nb::numpy to also be other things, in particular nb::dlpack or maybe nb::cupy and nb::pytorch.
  m.def("N_VGetArrayPointer",
        [](N_Vector v)
        {
          auto ptr = N_VGetArrayPointer(v);
          if (!ptr) { throw std::runtime_error("Failed to get array pointer"); }
          auto owner = nb::find(v);
          size_t shape[1]{static_cast<size_t>(N_VGetLength(v))};
          return nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>,
                             nb::c_contig>(ptr, 1, shape, owner);
        });
  m.def("N_VGetDeviceArrayPointer",
        [](N_Vector v)
        {
          auto ptr = N_VGetDeviceArrayPointer(v);
          if (!ptr) { throw std::runtime_error("Failed to get array pointer"); }
          auto owner = nb::find(v);
          size_t shape[1]{static_cast<size_t>(N_VGetLength(v))};
          return nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>,
                             nb::c_contig>(ptr, 1, shape, owner);
        });
  m.def("N_VSetArrayPointer",
        [](nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> arr,
           N_Vector v)
        {
          if (arr.shape(0) != N_VGetLength(v))
          {
            throw std::runtime_error(
              "Array shape does not match vector length");
          }
          N_VSetArrayPointer(arr.data(), v);
        });
}
