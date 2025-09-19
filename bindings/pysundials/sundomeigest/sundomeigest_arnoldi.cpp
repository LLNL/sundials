#include <nanobind/nanobind.h>
#include <sundomeigest/sundomeigest_arnoldi.h>

namespace nb = nanobind;

void bind_sundomeigest_arnoldi(nb::module_& m)
{
  m.def("SUNDomEigEstimator_Arnoldi", &SUNDomEigEstimator_Arnoldi,
        nb::rv_policy::reference);
}
