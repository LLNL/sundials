// #ifndef _LSRKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumARKODE_LSRKMethodType =
  nb::enum_<ARKODE_LSRKMethodType>(m, "ARKODE_LSRKMethodType",
                                   nb::is_arithmetic(), "")
    .value("ARKODE_LSRK_RKC_2", ARKODE_LSRK_RKC_2, "")
    .value("ARKODE_LSRK_RKL_2", ARKODE_LSRK_RKL_2, "")
    .value("ARKODE_LSRK_SSP_S_2", ARKODE_LSRK_SSP_S_2, "")
    .value("ARKODE_LSRK_SSP_S_3", ARKODE_LSRK_SSP_S_3, "")
    .value("ARKODE_LSRK_SSP_10_4", ARKODE_LSRK_SSP_10_4, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

sundials4py::scoped_def(
  m, "LSRKStepSetSTSMethod",
  [](const ARKodeView& arkode_mem, ARKODE_LSRKMethodType method) -> int
  {
    auto LSRKStepSetSTSMethod_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, ARKODE_LSRKMethodType method) -> int
    {
      auto lambda_result = LSRKStepSetSTSMethod(arkode_mem.get(), method);
      return lambda_result;
    };

    return LSRKStepSetSTSMethod_adapt_solver_memory(arkode_mem, method);
  },
  nb::arg("arkode_mem"), nb::arg("method"));

sundials4py::scoped_def(
  m, "LSRKStepSetSSPMethod",
  [](const ARKodeView& arkode_mem, ARKODE_LSRKMethodType method) -> int
  {
    auto LSRKStepSetSSPMethod_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, ARKODE_LSRKMethodType method) -> int
    {
      auto lambda_result = LSRKStepSetSSPMethod(arkode_mem.get(), method);
      return lambda_result;
    };

    return LSRKStepSetSSPMethod_adapt_solver_memory(arkode_mem, method);
  },
  nb::arg("arkode_mem"), nb::arg("method"));

sundials4py::scoped_def(
  m, "LSRKStepSetSTSMethodByName",
  [](const ARKodeView& arkode_mem, const char* emethod) -> int
  {
    auto LSRKStepSetSTSMethodByName_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, const char* emethod) -> int
    {
      auto lambda_result = LSRKStepSetSTSMethodByName(arkode_mem.get(), emethod);
      return lambda_result;
    };

    return LSRKStepSetSTSMethodByName_adapt_solver_memory(arkode_mem, emethod);
  },
  nb::arg("arkode_mem"), nb::arg("emethod"));

sundials4py::scoped_def(
  m, "LSRKStepSetSSPMethodByName",
  [](const ARKodeView& arkode_mem, const char* emethod) -> int
  {
    auto LSRKStepSetSSPMethodByName_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, const char* emethod) -> int
    {
      auto lambda_result = LSRKStepSetSSPMethodByName(arkode_mem.get(), emethod);
      return lambda_result;
    };

    return LSRKStepSetSSPMethodByName_adapt_solver_memory(arkode_mem, emethod);
  },
  nb::arg("arkode_mem"), nb::arg("emethod"));

sundials4py::scoped_def(
  m, "LSRKStepSetDomEigEstimator",
  [](const ARKodeView& arkode_mem, SUNDomEigEstimator DEE) -> int
  {
    auto LSRKStepSetDomEigEstimator_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, SUNDomEigEstimator DEE) -> int
    {
      auto lambda_result = LSRKStepSetDomEigEstimator(arkode_mem.get(), DEE);
      return lambda_result;
    };

    return LSRKStepSetDomEigEstimator_adapt_solver_memory(arkode_mem, DEE);
  },
  nb::arg("arkode_mem"), nb::arg("DEE"));

sundials4py::scoped_def(
  m, "LSRKStepSetDomEigFrequency",
  [](const ARKodeView& arkode_mem, long nsteps) -> int
  {
    auto LSRKStepSetDomEigFrequency_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, long nsteps) -> int
    {
      auto lambda_result = LSRKStepSetDomEigFrequency(arkode_mem.get(), nsteps);
      return lambda_result;
    };

    return LSRKStepSetDomEigFrequency_adapt_solver_memory(arkode_mem, nsteps);
  },
  nb::arg("arkode_mem"), nb::arg("nsteps"));

sundials4py::scoped_def(
  m, "LSRKStepSetMaxNumStages",
  [](const ARKodeView& arkode_mem, int stage_max_limit) -> int
  {
    auto LSRKStepSetMaxNumStages_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int stage_max_limit) -> int
    {
      auto lambda_result = LSRKStepSetMaxNumStages(arkode_mem.get(),
                                                   stage_max_limit);
      return lambda_result;
    };

    return LSRKStepSetMaxNumStages_adapt_solver_memory(arkode_mem,
                                                       stage_max_limit);
  },
  nb::arg("arkode_mem"), nb::arg("stage_max_limit"));

sundials4py::scoped_def(
  m, "LSRKStepSetDomEigSafetyFactor",
  [](const ARKodeView& arkode_mem, sunrealtype dom_eig_safety) -> int
  {
    auto LSRKStepSetDomEigSafetyFactor_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype dom_eig_safety) -> int
    {
      auto lambda_result = LSRKStepSetDomEigSafetyFactor(arkode_mem.get(),
                                                         dom_eig_safety);
      return lambda_result;
    };

    return LSRKStepSetDomEigSafetyFactor_adapt_solver_memory(arkode_mem,
                                                             dom_eig_safety);
  },
  nb::arg("arkode_mem"), nb::arg("dom_eig_safety"));

sundials4py::scoped_def(
  m, "LSRKStepSetUseAnalyticStabilityRegion",
  [](const ARKodeView& arkode_mem, sunbooleantype analytic_stab_region) -> int
  {
    auto LSRKStepSetUseAnalyticStabilityRegion_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunbooleantype analytic_stab_region) -> int
    {
      auto lambda_result =
        LSRKStepSetUseAnalyticStabilityRegion(arkode_mem.get(),
                                              analytic_stab_region);
      return lambda_result;
    };

    return LSRKStepSetUseAnalyticStabilityRegion_adapt_solver_memory(arkode_mem,
                                                                     analytic_stab_region);
  },
  nb::arg("arkode_mem"), nb::arg("analytic_stab_region"));

sundials4py::scoped_def(
  m, "LSRKStepSetNumDomEigEstInitPreprocessIters",
  [](const ARKodeView& arkode_mem, int num_iters) -> int
  {
    auto LSRKStepSetNumDomEigEstInitPreprocessIters_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int num_iters) -> int
    {
      auto lambda_result =
        LSRKStepSetNumDomEigEstInitPreprocessIters(arkode_mem.get(), num_iters);
      return lambda_result;
    };

    return LSRKStepSetNumDomEigEstInitPreprocessIters_adapt_solver_memory(arkode_mem,
                                                                          num_iters);
  },
  nb::arg("arkode_mem"), nb::arg("num_iters"));

sundials4py::scoped_def(
  m, "LSRKStepSetNumDomEigEstPreprocessIters",
  [](const ARKodeView& arkode_mem, int num_iters) -> int
  {
    auto LSRKStepSetNumDomEigEstPreprocessIters_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int num_iters) -> int
    {
      auto lambda_result =
        LSRKStepSetNumDomEigEstPreprocessIters(arkode_mem.get(), num_iters);
      return lambda_result;
    };

    return LSRKStepSetNumDomEigEstPreprocessIters_adapt_solver_memory(arkode_mem,
                                                                      num_iters);
  },
  nb::arg("arkode_mem"), nb::arg("num_iters"));

sundials4py::scoped_def(
  m, "LSRKStepSetNumSSPStages",
  [](const ARKodeView& arkode_mem, int num_of_stages) -> int
  {
    auto LSRKStepSetNumSSPStages_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int num_of_stages) -> int
    {
      auto lambda_result = LSRKStepSetNumSSPStages(arkode_mem.get(),
                                                   num_of_stages);
      return lambda_result;
    };

    return LSRKStepSetNumSSPStages_adapt_solver_memory(arkode_mem, num_of_stages);
  },
  nb::arg("arkode_mem"), nb::arg("num_of_stages"));

sundials4py::scoped_def(
  m, "LSRKStepGetNumDomEigUpdates",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto LSRKStepGetNumDomEigUpdates_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long dom_eig_num_evals_adapt_modifiable;

      int r = LSRKStepGetNumDomEigUpdates(arkode_mem,
                                          &dom_eig_num_evals_adapt_modifiable);
      return std::make_tuple(r, dom_eig_num_evals_adapt_modifiable);
    };
    auto LSRKStepGetNumDomEigUpdates_adapt_solver_memory =
      [&LSRKStepGetNumDomEigUpdates_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        LSRKStepGetNumDomEigUpdates_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return LSRKStepGetNumDomEigUpdates_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "LSRKStepGetMaxNumStages",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, int>
  {
    auto LSRKStepGetMaxNumStages_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, int>
    {
      int stage_max_adapt_modifiable;

      int r = LSRKStepGetMaxNumStages(arkode_mem, &stage_max_adapt_modifiable);
      return std::make_tuple(r, stage_max_adapt_modifiable);
    };
    auto LSRKStepGetMaxNumStages_adapt_solver_memory =
      [&LSRKStepGetMaxNumStages_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, int>
    {
      auto lambda_result =
        LSRKStepGetMaxNumStages_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return LSRKStepGetMaxNumStages_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "LSRKStepGetNumDomEigEstRhsEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto LSRKStepGetNumDomEigEstRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nfeDQ_adapt_modifiable;

      int r = LSRKStepGetNumDomEigEstRhsEvals(arkode_mem,
                                              &nfeDQ_adapt_modifiable);
      return std::make_tuple(r, nfeDQ_adapt_modifiable);
    };
    auto LSRKStepGetNumDomEigEstRhsEvals_adapt_solver_memory =
      [&LSRKStepGetNumDomEigEstRhsEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        LSRKStepGetNumDomEigEstRhsEvals_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return LSRKStepGetNumDomEigEstRhsEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "LSRKStepGetNumDomEigEstIters",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto LSRKStepGetNumDomEigEstIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long num_iters_adapt_modifiable;

      int r = LSRKStepGetNumDomEigEstIters(arkode_mem,
                                           &num_iters_adapt_modifiable);
      return std::make_tuple(r, num_iters_adapt_modifiable);
    };
    auto LSRKStepGetNumDomEigEstIters_adapt_solver_memory =
      [&LSRKStepGetNumDomEigEstIters_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        LSRKStepGetNumDomEigEstIters_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return LSRKStepGetNumDomEigEstIters_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
