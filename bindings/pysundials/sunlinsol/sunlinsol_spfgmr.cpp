#include <nanobind/nanobind.h>

#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_spfgmr.h>

namespace nb = nanobind;

void bind_sunlinsol_spfgmr(nb::module_& m)
{
  m.def("SUNLinSol_SPFGMR", &SUNLinSol_SPFGMR, nb::rv_policy::reference);
}