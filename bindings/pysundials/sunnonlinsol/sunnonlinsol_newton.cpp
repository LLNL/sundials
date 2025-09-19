#include <nanobind/nanobind.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

namespace nb = nanobind;

void bind_sunnonlinsol_newton(nb::module_& m)
{
  m.def("SUNNonlinSol_Newton", &SUNNonlinSol_Newton, nb::rv_policy::reference);
}
