#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_core(nb::module_ &m);
void bind_arkode(nb::module_ &m);

void bind_nvector_serial(nb::module_ &m);

void bind_sunlinsol_spgmr(nb::module_ &m);

NB_MODULE(pysundials, m) {

  nb::module_ core_m = m.def_submodule("core", "A submodule of 'pysundials'");
  bind_core(core_m);

  nb::module_ arkode_m = m.def_submodule("arkode", "A submodule of 'pysundials'");
  bind_arkode(arkode_m);

  nb::module_ serial_vector_m = m.def_submodule("core", "A submodule of 'pysundials'");
  bind_nvector_serial(serial_vector_m);
  
  nb::module_ spgmr_m = m.def_submodule("core", "A submodule of 'pysundials'");
  bind_sunlinsol_spgmr(spgmr_m);

}
