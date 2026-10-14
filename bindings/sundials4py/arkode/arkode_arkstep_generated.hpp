// #ifndef _ARKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

sundials4py::scoped_def(
  m, "ARKStepSetExplicit",
  [](const ARKodeView& arkode_mem) -> int
  {
    auto ARKStepSetExplicit_adapt_solver_memory =
      [](const ARKodeView& arkode_mem) -> int
    {
      auto lambda_result = ARKStepSetExplicit(arkode_mem.get());
      return lambda_result;
    };

    return ARKStepSetExplicit_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKStepSetImplicit",
  [](const ARKodeView& arkode_mem) -> int
  {
    auto ARKStepSetImplicit_adapt_solver_memory =
      [](const ARKodeView& arkode_mem) -> int
    {
      auto lambda_result = ARKStepSetImplicit(arkode_mem.get());
      return lambda_result;
    };

    return ARKStepSetImplicit_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKStepSetImEx",
  [](const ARKodeView& arkode_mem) -> int
  {
    auto ARKStepSetImEx_adapt_solver_memory = [](const ARKodeView& arkode_mem) -> int
    {
      auto lambda_result = ARKStepSetImEx(arkode_mem.get());
      return lambda_result;
    };

    return ARKStepSetImEx_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKStepSetTables",
  [](const ARKodeView& arkode_mem, int q, int p,
     std::optional<ARKodeButcherTable> Bi = std::nullopt,
     std::optional<ARKodeButcherTable> Be = std::nullopt) -> int
  {
    auto ARKStepSetTables_adapt_optional_arg_with_default_null =
      [](void* arkode_mem, int q, int p,
         std::optional<ARKodeButcherTable> Bi = std::nullopt,
         std::optional<ARKodeButcherTable> Be = std::nullopt) -> int
    {
      ARKodeButcherTable Bi_adapt_default_null = nullptr;
      if (Bi.has_value()) Bi_adapt_default_null = Bi.value();
      ARKodeButcherTable Be_adapt_default_null = nullptr;
      if (Be.has_value()) Be_adapt_default_null = Be.value();

      auto lambda_result = ARKStepSetTables(arkode_mem, q, p,
                                            Bi_adapt_default_null,
                                            Be_adapt_default_null);
      return lambda_result;
    };
    auto ARKStepSetTables_adapt_solver_memory =
      [&ARKStepSetTables_adapt_optional_arg_with_default_null](const ARKodeView& arkode_mem,
                                                               int q, int p,
                                                               std::optional<ARKodeButcherTable>
                                                                 Bi = std::nullopt,
                                                               std::optional<ARKodeButcherTable>
                                                                 Be =
                                                                   std::nullopt) -> int
    {
      auto lambda_result =
        ARKStepSetTables_adapt_optional_arg_with_default_null(arkode_mem.get(),
                                                              q, p, Bi, Be);
      return lambda_result;
    };

    return ARKStepSetTables_adapt_solver_memory(arkode_mem, q, p, Bi, Be);
  },
  nb::arg("arkode_mem"), nb::arg("q"), nb::arg("p"),
  nb::arg("Bi").none() = nb::none(), nb::arg("Be").none() = nb::none());

sundials4py::scoped_def(
  m, "ARKStepSetTableNum",
  [](const ARKodeView& arkode_mem, ARKODE_DIRKTableID itable,
     ARKODE_ERKTableID etable) -> int
  {
    auto ARKStepSetTableNum_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, ARKODE_DIRKTableID itable,
         ARKODE_ERKTableID etable) -> int
    {
      auto lambda_result = ARKStepSetTableNum(arkode_mem.get(), itable, etable);
      return lambda_result;
    };

    return ARKStepSetTableNum_adapt_solver_memory(arkode_mem, itable, etable);
  },
  nb::arg("arkode_mem"), nb::arg("itable"), nb::arg("etable"));

sundials4py::scoped_def(
  m, "ARKStepSetTableName",
  [](const ARKodeView& arkode_mem, const char* itable, const char* etable) -> int
  {
    auto ARKStepSetTableName_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, const char* itable,
         const char* etable) -> int
    {
      auto lambda_result = ARKStepSetTableName(arkode_mem.get(), itable, etable);
      return lambda_result;
    };

    return ARKStepSetTableName_adapt_solver_memory(arkode_mem, itable, etable);
  },
  nb::arg("arkode_mem"), nb::arg("itable"), nb::arg("etable"));

sundials4py::scoped_def(
  m, "ARKStepGetCurrentButcherTables",
  [](const ARKodeView& arkode_mem)
    -> std::tuple<int, ARKodeButcherTable, ARKodeButcherTable>
  {
    auto ARKStepGetCurrentButcherTables_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, ARKodeButcherTable, ARKodeButcherTable>
    {
      ARKodeButcherTable Bi_adapt_modifiable;
      ARKodeButcherTable Be_adapt_modifiable;

      int r = ARKStepGetCurrentButcherTables(arkode_mem, &Bi_adapt_modifiable,
                                             &Be_adapt_modifiable);
      return std::make_tuple(r, Bi_adapt_modifiable, Be_adapt_modifiable);
    };
    auto ARKStepGetCurrentButcherTables_adapt_solver_memory =
      [&ARKStepGetCurrentButcherTables_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem)
      -> std::tuple<int, ARKodeButcherTable, ARKodeButcherTable>
    {
      auto lambda_result =
        ARKStepGetCurrentButcherTables_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKStepGetCurrentButcherTables_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"),
  " Optional output functions\n\n nb::rv_policy::reference",
  nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "ARKStepGetTimestepperStats",
  [](const ARKodeView& arkode_mem)
    -> std::tuple<int, long, long, long, long, long, long, long>
  {
    auto ARKStepGetTimestepperStats_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem)
      -> std::tuple<int, long, long, long, long, long, long, long>
    {
      long expsteps_adapt_modifiable;
      long accsteps_adapt_modifiable;
      long step_attempts_adapt_modifiable;
      long nfe_evals_adapt_modifiable;
      long nfi_evals_adapt_modifiable;
      long nlinsetups_adapt_modifiable;
      long netfails_adapt_modifiable;

      int r = ARKStepGetTimestepperStats(arkode_mem, &expsteps_adapt_modifiable,
                                         &accsteps_adapt_modifiable,
                                         &step_attempts_adapt_modifiable,
                                         &nfe_evals_adapt_modifiable,
                                         &nfi_evals_adapt_modifiable,
                                         &nlinsetups_adapt_modifiable,
                                         &netfails_adapt_modifiable);
      return std::make_tuple(r, expsteps_adapt_modifiable,
                             accsteps_adapt_modifiable,
                             step_attempts_adapt_modifiable,
                             nfe_evals_adapt_modifiable,
                             nfi_evals_adapt_modifiable,
                             nlinsetups_adapt_modifiable,
                             netfails_adapt_modifiable);
    };
    auto ARKStepGetTimestepperStats_adapt_solver_memory =
      [&ARKStepGetTimestepperStats_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem)
      -> std::tuple<int, long, long, long, long, long, long, long>
    {
      auto lambda_result =
        ARKStepGetTimestepperStats_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKStepGetTimestepperStats_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
