#include <nanobind/nanobind.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

namespace nb = nanobind;

void bind_sunadaptcontroller_soderlind(nb::module_& m)
{
  m.def("SUNAdaptController_Soderlind", &SUNAdaptController_Soderlind,
        nb::rv_policy::reference);
}
