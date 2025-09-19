#include <nanobind/nanobind.h>
#include <sundials/sundials_core.hpp>
#include <sunmatrix/sunmatrix_band.h>

namespace nb = nanobind;

void bind_sunmatrix_band(nb::module_& m)
{
  m.def("SUNBandMatrix", &SUNBandMatrix, nb::rv_policy::reference);
}
