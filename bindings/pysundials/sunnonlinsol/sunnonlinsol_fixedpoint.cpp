#include <nanobind/nanobind.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

namespace nb = nanobind;

void bind_sunnonlinsol_fixedpoint(nb::module_& m)
{
  m.def("SUNNonlinSol_FixedPoint", &SUNNonlinSol_FixedPoint,
        nb::rv_policy::reference);
}
