#include <nanobind/nanobind.h>

#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_band.h>

namespace nb = nanobind;

void bind_sunlinsol_band(nb::module_& m)
{
  m.def("SUNLinSol_Band", &SUNLinSol_Band, nb::rv_policy::reference);
}