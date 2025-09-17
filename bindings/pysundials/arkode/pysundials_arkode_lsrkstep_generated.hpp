// #ifndef _LSRKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

m.def("LSRKStepSetSTSMethod", LSRKStepSetSTSMethod, nb::arg("arkode_mem"),
      nb::arg("method"));

m.def("LSRKStepSetSSPMethod", LSRKStepSetSSPMethod, nb::arg("arkode_mem"),
      nb::arg("method"));

m.def("LSRKStepSetSTSMethodByName", LSRKStepSetSTSMethodByName,
      nb::arg("arkode_mem"), nb::arg("emethod"));

m.def("LSRKStepSetSSPMethodByName", LSRKStepSetSSPMethodByName,
      nb::arg("arkode_mem"), nb::arg("emethod"));

m.def("LSRKStepSetDomEigFrequency", LSRKStepSetDomEigFrequency,
      nb::arg("arkode_mem"), nb::arg("nsteps"));

m.def("LSRKStepSetMaxNumStages", LSRKStepSetMaxNumStages, nb::arg("arkode_mem"),
      nb::arg("stage_max_limit"));

m.def("LSRKStepSetDomEigSafetyFactor", LSRKStepSetDomEigSafetyFactor,
      nb::arg("arkode_mem"), nb::arg("dom_eig_safety"));

m.def("LSRKStepSetNumSSPStages", LSRKStepSetNumSSPStages, nb::arg("arkode_mem"),
      nb::arg("num_of_stages"));

m.def(
  "LSRKStepGetNumDomEigUpdates",
  [](void* arkode_mem, long dom_eig_num_evals) -> std::tuple<int, long>
  {
    auto LSRKStepGetNumDomEigUpdates_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, long dom_eig_num_evals) -> std::tuple<int, long>
    {
      long* dom_eig_num_evals_adapt_modifiable = &dom_eig_num_evals;

      int r = LSRKStepGetNumDomEigUpdates(arkode_mem,
                                          dom_eig_num_evals_adapt_modifiable);
      return std::make_tuple(r, dom_eig_num_evals);
    };

    return LSRKStepGetNumDomEigUpdates_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                            dom_eig_num_evals);
  },
  nb::arg("arkode_mem"), nb::arg("dom_eig_num_evals"));

m.def(
  "LSRKStepGetMaxNumStages",
  [](void* arkode_mem, int stage_max) -> std::tuple<int, int>
  {
    auto LSRKStepGetMaxNumStages_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, int stage_max) -> std::tuple<int, int>
    {
      int* stage_max_adapt_modifiable = &stage_max;

      int r = LSRKStepGetMaxNumStages(arkode_mem, stage_max_adapt_modifiable);
      return std::make_tuple(r, stage_max);
    };

    return LSRKStepGetMaxNumStages_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                        stage_max);
  },
  nb::arg("arkode_mem"), nb::arg("stage_max"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
