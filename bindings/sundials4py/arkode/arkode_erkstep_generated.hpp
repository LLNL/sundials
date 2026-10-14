// #ifndef _ERKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

sundials4py::scoped_def(
  m, "ERKStepSetTable",
  [](const ARKodeView& arkode_mem, ARKodeButcherTable B) -> int
  {
    auto ERKStepSetTable_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                                  ARKodeButcherTable B) -> int
    {
      auto lambda_result = ERKStepSetTable(arkode_mem.get(), B);
      return lambda_result;
    };

    return ERKStepSetTable_adapt_solver_memory(arkode_mem, B);
  },
  nb::arg("arkode_mem"), nb::arg("B"));

sundials4py::scoped_def(
  m, "ERKStepSetTableNum",
  [](const ARKodeView& arkode_mem, ARKODE_ERKTableID etable) -> int
  {
    auto ERKStepSetTableNum_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, ARKODE_ERKTableID etable) -> int
    {
      auto lambda_result = ERKStepSetTableNum(arkode_mem.get(), etable);
      return lambda_result;
    };

    return ERKStepSetTableNum_adapt_solver_memory(arkode_mem, etable);
  },
  nb::arg("arkode_mem"), nb::arg("etable"));

sundials4py::scoped_def(
  m, "ERKStepSetTableName",
  [](const ARKodeView& arkode_mem, const char* etable) -> int
  {
    auto ERKStepSetTableName_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, const char* etable) -> int
    {
      auto lambda_result = ERKStepSetTableName(arkode_mem.get(), etable);
      return lambda_result;
    };

    return ERKStepSetTableName_adapt_solver_memory(arkode_mem, etable);
  },
  nb::arg("arkode_mem"), nb::arg("etable"));

sundials4py::scoped_def(
  m, "ERKStepGetCurrentButcherTable",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, ARKodeButcherTable>
  {
    auto ERKStepGetCurrentButcherTable_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, ARKodeButcherTable>
    {
      ARKodeButcherTable B_adapt_modifiable;

      int r = ERKStepGetCurrentButcherTable(arkode_mem, &B_adapt_modifiable);
      return std::make_tuple(r, B_adapt_modifiable);
    };
    auto ERKStepGetCurrentButcherTable_adapt_solver_memory =
      [&ERKStepGetCurrentButcherTable_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, ARKodeButcherTable>
    {
      auto lambda_result =
        ERKStepGetCurrentButcherTable_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ERKStepGetCurrentButcherTable_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"),
  " Optional output functions\n\n nb::rv_policy::reference",
  nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "ERKStepGetTimestepperStats",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long, long, long, long, long>
  {
    auto ERKStepGetTimestepperStats_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long, long, long, long, long>
    {
      long expsteps_adapt_modifiable;
      long accsteps_adapt_modifiable;
      long step_attempts_adapt_modifiable;
      long nfevals_adapt_modifiable;
      long netfails_adapt_modifiable;

      int r = ERKStepGetTimestepperStats(arkode_mem, &expsteps_adapt_modifiable,
                                         &accsteps_adapt_modifiable,
                                         &step_attempts_adapt_modifiable,
                                         &nfevals_adapt_modifiable,
                                         &netfails_adapt_modifiable);
      return std::make_tuple(r, expsteps_adapt_modifiable,
                             accsteps_adapt_modifiable,
                             step_attempts_adapt_modifiable,
                             nfevals_adapt_modifiable, netfails_adapt_modifiable);
    };
    auto ERKStepGetTimestepperStats_adapt_solver_memory =
      [&ERKStepGetTimestepperStats_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long, long, long, long, long>
    {
      auto lambda_result =
        ERKStepGetTimestepperStats_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ERKStepGetTimestepperStats_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"), "Grouped optional output functions");
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
