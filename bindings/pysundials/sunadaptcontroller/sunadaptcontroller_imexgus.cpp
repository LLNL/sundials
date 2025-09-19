#include <nanobind/nanobind.h>
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>

namespace nb = nanobind;

void bind_sunadaptcontroller_imexgus(nb::module_& m)
{
  m.def("SUNAdaptController_ImExGus", &SUNAdaptController_ImExGus,
        nb::rv_policy::reference);
}
