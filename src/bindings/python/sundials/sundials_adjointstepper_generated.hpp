// #ifndef _SUNADJOINT_STEPPER_H
//
// #ifdef __cplusplus
//
// #endif
//

m.def("SUNAdjointStepper_ReInit", SUNAdjointStepper_ReInit, nb::arg("adj"),
      nb::arg("t0"), nb::arg("y0"), nb::arg("tf"), nb::arg("sf"));

m.def("SUNAdjointStepper_Evolve", SUNAdjointStepper_Evolve,
      nb::arg("adj_stepper"), nb::arg("tout"), nb::arg("sens"), nb::arg("tret"));

m.def("SUNAdjointStepper_OneStep", SUNAdjointStepper_OneStep,
      nb::arg("adj_stepper"), nb::arg("tout"), nb::arg("sens"), nb::arg("tret"));

m.def("SUNAdjointStepper_RecomputeFwd", SUNAdjointStepper_RecomputeFwd,
      nb::arg("adj_stepper"), nb::arg("start_idx"), nb::arg("t0"),
      nb::arg("y0"), nb::arg("tf"));

m.def("SUNAdjointStepper_SetUserData", SUNAdjointStepper_SetUserData,
      nb::arg("param_0"), nb::arg("user_data"));

m.def("SUNAdjointStepper_GetNumSteps", SUNAdjointStepper_GetNumSteps,
      nb::arg("adj_stepper"), nb::arg("num_steps"));

m.def("SUNAdjointStepper_GetNumRecompute", SUNAdjointStepper_GetNumRecompute,
      nb::arg("adj_stepper"), nb::arg("num_recompute"));

m.def("SUNAdjointStepper_PrintAllStats", SUNAdjointStepper_PrintAllStats,
      nb::arg("adj_stepper"), nb::arg("outfile"), nb::arg("fmt"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
