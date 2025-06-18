#include <nanobind/nanobind.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_core.h>

namespace nb = nanobind;

void bind_nvector_serial(nb::module_& m)
{
  m.def("N_VNewEmpty_Serial", &N_VNewEmpty_Serial, nb::rv_policy::reference);
  m.def("N_VNew_Serial", &N_VNew_Serial, nb::rv_policy::reference);
  m.def("N_VMake_Serial", &N_VMake_Serial, nb::rv_policy::reference);
}
