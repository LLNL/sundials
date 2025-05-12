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
 * This file is the entrypoint for the Python binding code for the
 * SUNDIALS N_Vector class. It contains hand-written code for functions
 * that require special treatment, and includes the generated code
 * produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials/sundials_nvector.h"
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_nvector.hpp>

namespace nb = nanobind;

void bind_nvector(nb::module_& m)
{
#include "sundials_nvector_generated.cpp"

  nb::class_<sundials::experimental::NVectorView>(m, "NVectorView")
    .def(nb::init<>())
    .def(nb::init<_generic_N_Vector*>())
    .def("get",
         nb::overload_cast<>(&sundials::experimental::NVectorView::get,
                             nb::const_),
         nb::rv_policy::reference);

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
