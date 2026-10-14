// #ifndef _ARKODE_SPRKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

sundials4py::scoped_def(
  m, "SPRKStepSetMethod",
  [](const ARKodeView& arkode_mem, ARKodeSPRKTable sprk_storage) -> int
  {
    auto SPRKStepSetMethod_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, ARKodeSPRKTable sprk_storage) -> int
    {
      auto lambda_result = SPRKStepSetMethod(arkode_mem.get(), sprk_storage);
      return lambda_result;
    };

    return SPRKStepSetMethod_adapt_solver_memory(arkode_mem, sprk_storage);
  },
  nb::arg("arkode_mem"), nb::arg("sprk_storage"));

sundials4py::scoped_def(
  m, "SPRKStepSetMethodName",
  [](const ARKodeView& arkode_mem, const char* method) -> int
  {
    auto SPRKStepSetMethodName_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, const char* method) -> int
    {
      auto lambda_result = SPRKStepSetMethodName(arkode_mem.get(), method);
      return lambda_result;
    };

    return SPRKStepSetMethodName_adapt_solver_memory(arkode_mem, method);
  },
  nb::arg("arkode_mem"), nb::arg("method"));

sundials4py::scoped_def(
  m, "SPRKStepGetCurrentMethod",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, ARKodeSPRKTable>
  {
    auto SPRKStepGetCurrentMethod_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, ARKodeSPRKTable>
    {
      ARKodeSPRKTable sprk_storage_adapt_modifiable;

      int r = SPRKStepGetCurrentMethod(arkode_mem,
                                       &sprk_storage_adapt_modifiable);
      return std::make_tuple(r, sprk_storage_adapt_modifiable);
    };
    auto SPRKStepGetCurrentMethod_adapt_solver_memory =
      [&SPRKStepGetCurrentMethod_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, ARKodeSPRKTable>
    {
      auto lambda_result =
        SPRKStepGetCurrentMethod_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return SPRKStepGetCurrentMethod_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
