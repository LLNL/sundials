// #ifndef _LSRKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumARKODE_LSRKMethodType_ =
  nb::enum_<ARKODE_LSRKMethodType_>(m, "ARKODE_LSRKMethodType_",
                                    nb::is_arithmetic(), "")
    .value("ARKODE_LSRK_RKC_2", ARKODE_LSRK_RKC_2, "")
    .value("ARKODE_LSRK_RKL_2", ARKODE_LSRK_RKL_2, "")
    .value("ARKODE_LSRK_SSP_S_2", ARKODE_LSRK_SSP_S_2, "")
    .value("ARKODE_LSRK_SSP_S_3", ARKODE_LSRK_SSP_S_3, "")
    .value("ARKODE_LSRK_SSP_10_4", ARKODE_LSRK_SSP_10_4, "")
    .export_values();

m.def("LSRKStepSetSTSMethod", LSRKStepSetSTSMethod, nb::arg("arkode_mem"),
      nb::arg("method"));

m.def("LSRKStepSetSSPMethod", LSRKStepSetSSPMethod, nb::arg("arkode_mem"),
      nb::arg("method"));

m.def("LSRKStepSetSTSMethodByName", LSRKStepSetSTSMethodByName,
      nb::arg("arkode_mem"), nb::arg("emethod"));

m.def("LSRKStepSetSSPMethodByName", LSRKStepSetSSPMethodByName,
      nb::arg("arkode_mem"), nb::arg("emethod"));

m.def("LSRKStepSetDomEigEstimator", LSRKStepSetDomEigEstimator,
      nb::arg("arkode_mem"), nb::arg("DEE"));

m.def("LSRKStepSetDomEigFrequency", LSRKStepSetDomEigFrequency,
      nb::arg("arkode_mem"), nb::arg("nsteps"));

m.def("LSRKStepSetMaxNumStages", LSRKStepSetMaxNumStages, nb::arg("arkode_mem"),
      nb::arg("stage_max_limit"));

m.def("LSRKStepSetDomEigSafetyFactor", LSRKStepSetDomEigSafetyFactor,
      nb::arg("arkode_mem"), nb::arg("dom_eig_safety"));

m.def("LSRKStepSetNumDomEigEstInitPreprocessIters",
      LSRKStepSetNumDomEigEstInitPreprocessIters, nb::arg("arkode_mem"),
      nb::arg("num_iters"));

m.def("LSRKStepSetNumDomEigEstPreprocessIters",
      LSRKStepSetNumDomEigEstPreprocessIters, nb::arg("arkode_mem"),
      nb::arg("num_iters"));

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

m.def(
  "LSRKStepGetNumDomEigEstRhsEvals",
  [](void* arkode_mem, long nfeDQ) -> std::tuple<int, long>
  {
    auto LSRKStepGetNumDomEigEstRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, long nfeDQ) -> std::tuple<int, long>
    {
      long* nfeDQ_adapt_modifiable = &nfeDQ;

      int r = LSRKStepGetNumDomEigEstRhsEvals(arkode_mem, nfeDQ_adapt_modifiable);
      return std::make_tuple(r, nfeDQ);
    };

    return LSRKStepGetNumDomEigEstRhsEvals_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                                nfeDQ);
  },
  nb::arg("arkode_mem"), nb::arg("nfeDQ"));

m.def(
  "LSRKStepGetNumDomEigEstIters",
  [](void* arkode_mem, long num_iters) -> std::tuple<int, long>
  {
    auto LSRKStepGetNumDomEigEstIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, long num_iters) -> std::tuple<int, long>
    {
      long* num_iters_adapt_modifiable = &num_iters;

      int r = LSRKStepGetNumDomEigEstIters(arkode_mem,
                                           num_iters_adapt_modifiable);
      return std::make_tuple(r, num_iters);
    };

    return LSRKStepGetNumDomEigEstIters_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                             num_iters);
  },
  nb::arg("arkode_mem"), nb::arg("num_iters"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
