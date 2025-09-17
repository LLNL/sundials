#include <nanobind/nanobind.h>

#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_spbcgs.h>

namespace nb = nanobind;

void bind_sunlinsol_spbcgs(nb::module_& m)
{
  m.def("SUNLinSol_SPBCGS", &SUNLinSol_SPBCGS, nb::rv_policy::reference);
}