// #ifndef _SUNDIALS_STEPPER_H
// 
// #ifdef __cplusplus
// 
// #endif
// 

m.def("SUNStepper_Evolve",
    SUNStepper_Evolve, nb::arg("stepper"), nb::arg("tout"), nb::arg("vret"), nb::arg("tret"));

m.def("SUNStepper_OneStep",
    SUNStepper_OneStep, nb::arg("stepper"), nb::arg("tout"), nb::arg("vret"), nb::arg("tret"));

m.def("SUNStepper_FullRhs",
    SUNStepper_FullRhs, nb::arg("stepper"), nb::arg("t"), nb::arg("v"), nb::arg("f"), nb::arg("mode"));

m.def("SUNStepper_ReInit",
    SUNStepper_ReInit, nb::arg("stepper"), nb::arg("t0"), nb::arg("v0"));

m.def("SUNStepper_Reset",
    SUNStepper_Reset, nb::arg("stepper"), nb::arg("tR"), nb::arg("vR"));

m.def("SUNStepper_ResetCheckpointIndex",
    SUNStepper_ResetCheckpointIndex, nb::arg("stepper"), nb::arg("ckptIdxR"));

m.def("SUNStepper_SetStopTime",
    SUNStepper_SetStopTime, nb::arg("stepper"), nb::arg("tstop"));

m.def("SUNStepper_SetStepDirection",
    SUNStepper_SetStepDirection, nb::arg("stepper"), nb::arg("stepdir"));

m.def("SUNStepper_SetLastFlag",
    SUNStepper_SetLastFlag, nb::arg("stepper"), nb::arg("last_flag"));

m.def("SUNStepper_GetLastFlag",
    SUNStepper_GetLastFlag, nb::arg("stepper"), nb::arg("last_flag"));

m.def("SUNStepper_GetNumSteps",
    SUNStepper_GetNumSteps, nb::arg("stepper"), nb::arg("nst"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif 
