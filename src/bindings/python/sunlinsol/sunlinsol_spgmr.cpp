#include <nanobind/nanobind.h>

#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_spgmr.h>

namespace nb = nanobind;

void bind_sunlinsol_spgmr(nb::module_& m)
{
  m.def("SUNLinSol_SPGMR", &SUNLinSol_SPGMR, nb::rv_policy::reference);
}
