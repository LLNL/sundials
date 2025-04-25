#include <nanobind/nanobind.h>

#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_pcg.h>

namespace nb = nanobind;

void bind_sunlinsol_pcg(nb::module_& m)
{
  m.def("SUNLinSol_PCG", &SUNLinSol_PCG, nb::rv_policy::reference);
}