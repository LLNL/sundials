#include <nanobind/nanobind.h>
#include <sundomeigest/sundomeigest_power.h>

namespace nb = nanobind;

void bind_sundomeigest_power(nb::module_& m)
{
  m.def("SUNDomEigEstimator_Power", &SUNDomEigEstimator_Power,
        nb::rv_policy::reference);
}
