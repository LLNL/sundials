#include <nanobind/nanobind.h>
#include <sunadaptcontroller/sunadaptcontroller_mrihtol.h>

namespace nb = nanobind;

void bind_sunadaptcontroller_mrihtol(nb::module_& m)
{
  m.def("SUNAdaptController_MRIHTol", &SUNAdaptController_MRIHTol,
        nb::rv_policy::reference);
}
