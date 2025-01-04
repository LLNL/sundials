#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_nvector.hpp>

namespace nb = nanobind;

void bind_nvector(nb::module_& m)
{
  nb::class_<sundials::experimental::NVectorView>(m, "NVectorView")
    .def(nb::init<>())
    .def(nb::init<_generic_N_Vector*>())
    // Option 1: nv.get() must be invoked in Python to convert to N_Vector before calling N_V functions
    .def("get",
         nb::overload_cast<>(&sundials::experimental::NVectorView::get,
                             nb::const_),
         nb::rv_policy::reference)
    // Option 2: nv.GetArrayPointer() must be invoked in Python and we wrap every N_V function as a class method
    .def(
      "GetArrayPointer",
      [](sundials::experimental::NVectorView&& v)
      { return N_VGetArrayPointer(v); },
      nb::rv_policy::reference);

  // I don't think implicit conversion will work unless we make the View classes convertible to the underlying type instead of the pointer type
  // nb::implicitly_convertible<sundials::experimental::NVectorView, _generic_N_Vector*>();

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

  #include "sundials_nvector_generated.cpp"
}
