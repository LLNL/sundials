#include <nanobind/nanobind.h>

namespace nb = nanobind;

// Forward declarations of all of the binding functions
void bind_core(nb::module_& m);
void bind_arkode(nb::module_& m);
void bind_cvodes(nb::module_& m);
void bind_idas(nb::module_& m);
void bind_nvector_serial(nb::module_& m);
void bind_sunlinsol_spgmr(nb::module_& m);
void bind_sunlinsol_dense(nb::module_& m);
void bind_sunlinsol_band(nb::module_& m);
void bind_sunlinsol_spbcgs(nb::module_& m);
void bind_sunlinsol_spfgmr(nb::module_& m);
void bind_sunlinsol_sptfqmr(nb::module_& m);
void bind_sunlinsol_pcg(nb::module_& m);
void bind_sunlinsol_klu(nb::module_& m);

// void bind_sunlinsol_lapackdense(nb::module_ &m);
// void bind_sunlinsol_lapackband(nb::module_ &m);

NB_MODULE(pysundials, m)
{
  nb::module_ core_m = m.def_submodule("core", "A submodule of 'pysundials'");
  bind_core(core_m);

  nb::module_ arkode_m = m.def_submodule("arkode",
                                         "A submodule of 'pysundials'");
  bind_arkode(arkode_m);

  nb::module_ cvodes_m = m.def_submodule("cvodes",
                                         "A submodule of 'pysundials'");
  bind_cvodes(cvodes_m);

  nb::module_ idas_m = m.def_submodule("idas",
                                         "A submodule of 'pysundials'");
  bind_idas(idas_m);

  {
    // Bind N_Vector implementations
    nb::module_ serial_vector_m =
      m.def_submodule("core", "A submodule of 'pysundials'");
    bind_nvector_serial(serial_vector_m);
  }

  {
    // Bind SUNLinearSolver implementations
    nb::module_ spgmr_m = m.def_submodule("core",
                                          "A submodule of 'pysundials'");
    bind_sunlinsol_spgmr(spgmr_m);

    nb::module_ dense_m = m.def_submodule("core",
                                          "A submodule of 'pysundials'");
    bind_sunlinsol_dense(dense_m);

    nb::module_ band_m = m.def_submodule("core", "A submodule of 'pysundials'");
    bind_sunlinsol_band(band_m);

    nb::module_ spbcgs_m = m.def_submodule("core",
                                           "A submodule of 'pysundials'");
    bind_sunlinsol_spbcgs(spbcgs_m);

    nb::module_ spfgmr_m = m.def_submodule("core",
                                           "A submodule of 'pysundials'");
    bind_sunlinsol_spfgmr(spfgmr_m);

    nb::module_ sptfqmr_m = m.def_submodule("core",
                                            "A submodule of 'pysundials'");
    bind_sunlinsol_sptfqmr(sptfqmr_m);

    nb::module_ pcg_m = m.def_submodule("core", "A submodule of 'pysundials'");
    bind_sunlinsol_pcg(pcg_m);
  }
}
