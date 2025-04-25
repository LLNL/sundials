#include <nanobind/nanobind.h>

#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_lapackband.h>

namespace nb = nanobind;

void bind_sunlinsol_lapackband(nb::module_& m)
{
  m.def("SUNLinSol_LapackBand", &SUNLinSol_LapackBand, nb::rv_policy::reference);
}