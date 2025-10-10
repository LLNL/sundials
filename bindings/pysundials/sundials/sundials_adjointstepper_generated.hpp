// #ifndef _SUNADJOINT_STEPPER_H
//
// #ifdef __cplusplus
//
// #endif
//

m.def("SUNAdjointStepper_ReInit", SUNAdjointStepper_ReInit, nb::arg("adj"),
      nb::arg("t0"), nb::arg("y0"), nb::arg("tf"), nb::arg("sf"));

m.def(
  "SUNAdjointStepper_Evolve",
  [](SUNAdjointStepper adj_stepper, sunrealtype tout,
     N_Vector sens) -> std::tuple<SUNErrCode, sunrealtype>
  {
    auto SUNAdjointStepper_Evolve_adapt_modifiable_immutable_to_return =
      [](SUNAdjointStepper adj_stepper, sunrealtype tout,
         N_Vector sens) -> std::tuple<SUNErrCode, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      SUNErrCode r = SUNAdjointStepper_Evolve(adj_stepper, tout, sens,
                                              &tret_adapt_modifiable);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return SUNAdjointStepper_Evolve_adapt_modifiable_immutable_to_return(adj_stepper,
                                                                         tout,
                                                                         sens);
  },
  nb::arg("adj_stepper"), nb::arg("tout"), nb::arg("sens"));

m.def(
  "SUNAdjointStepper_OneStep",
  [](SUNAdjointStepper adj_stepper, sunrealtype tout,
     N_Vector sens) -> std::tuple<SUNErrCode, sunrealtype>
  {
    auto SUNAdjointStepper_OneStep_adapt_modifiable_immutable_to_return =
      [](SUNAdjointStepper adj_stepper, sunrealtype tout,
         N_Vector sens) -> std::tuple<SUNErrCode, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      SUNErrCode r = SUNAdjointStepper_OneStep(adj_stepper, tout, sens,
                                               &tret_adapt_modifiable);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return SUNAdjointStepper_OneStep_adapt_modifiable_immutable_to_return(adj_stepper,
                                                                          tout,
                                                                          sens);
  },
  nb::arg("adj_stepper"), nb::arg("tout"), nb::arg("sens"));

m.def("SUNAdjointStepper_RecomputeFwd", SUNAdjointStepper_RecomputeFwd,
      nb::arg("adj_stepper"), nb::arg("start_idx"), nb::arg("t0"),
      nb::arg("y0"), nb::arg("tf"));

m.def("SUNAdjointStepper_SetUserData", SUNAdjointStepper_SetUserData,
      nb::arg("param_0"), nb::arg("user_data"));

m.def(
  "SUNAdjointStepper_GetNumSteps",
  [](SUNAdjointStepper adj_stepper) -> std::tuple<SUNErrCode, suncountertype>
  {
    auto SUNAdjointStepper_GetNumSteps_adapt_modifiable_immutable_to_return =
      [](SUNAdjointStepper adj_stepper) -> std::tuple<SUNErrCode, suncountertype>
    {
      suncountertype num_steps_adapt_modifiable;

      SUNErrCode r = SUNAdjointStepper_GetNumSteps(adj_stepper,
                                                   &num_steps_adapt_modifiable);
      return std::make_tuple(r, num_steps_adapt_modifiable);
    };

    return SUNAdjointStepper_GetNumSteps_adapt_modifiable_immutable_to_return(
      adj_stepper);
  },
  nb::arg("adj_stepper"));

m.def(
  "SUNAdjointStepper_GetNumRecompute",
  [](SUNAdjointStepper adj_stepper) -> std::tuple<SUNErrCode, suncountertype>
  {
    auto SUNAdjointStepper_GetNumRecompute_adapt_modifiable_immutable_to_return =
      [](SUNAdjointStepper adj_stepper) -> std::tuple<SUNErrCode, suncountertype>
    {
      suncountertype num_recompute_adapt_modifiable;

      SUNErrCode r =
        SUNAdjointStepper_GetNumRecompute(adj_stepper,
                                          &num_recompute_adapt_modifiable);
      return std::make_tuple(r, num_recompute_adapt_modifiable);
    };

    return SUNAdjointStepper_GetNumRecompute_adapt_modifiable_immutable_to_return(
      adj_stepper);
  },
  nb::arg("adj_stepper"));

m.def("SUNAdjointStepper_PrintAllStats", SUNAdjointStepper_PrintAllStats,
      nb::arg("adj_stepper"), nb::arg("outfile"), nb::arg("fmt"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
