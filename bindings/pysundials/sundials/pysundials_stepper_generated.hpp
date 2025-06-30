// #ifndef _SUNDIALS_STEPPER_H
// 
// #ifdef __cplusplus
// 
// #endif
// 

m.def("SUNStepper_Evolve",
    [](SUNStepper stepper, double tout, N_Vector vret, double tret) -> std::tuple<SUNErrCode, double>
    {
        auto SUNStepper_Evolve_adapt_modifiable_immutable_to_return = [](SUNStepper stepper, double tout, N_Vector vret, double tret) -> std::tuple<SUNErrCode, double>
        {
            double * tret_adapt_modifiable = & tret;

            SUNErrCode r = SUNStepper_Evolve(stepper, tout, vret, tret_adapt_modifiable);
            return std::make_tuple(r, tret);
        };

        return SUNStepper_Evolve_adapt_modifiable_immutable_to_return(stepper, tout, vret, tret);
    },     nb::arg("stepper"), nb::arg("tout"), nb::arg("vret"), nb::arg("tret"));

m.def("SUNStepper_OneStep",
    [](SUNStepper stepper, double tout, N_Vector vret, double tret) -> std::tuple<SUNErrCode, double>
    {
        auto SUNStepper_OneStep_adapt_modifiable_immutable_to_return = [](SUNStepper stepper, double tout, N_Vector vret, double tret) -> std::tuple<SUNErrCode, double>
        {
            double * tret_adapt_modifiable = & tret;

            SUNErrCode r = SUNStepper_OneStep(stepper, tout, vret, tret_adapt_modifiable);
            return std::make_tuple(r, tret);
        };

        return SUNStepper_OneStep_adapt_modifiable_immutable_to_return(stepper, tout, vret, tret);
    },     nb::arg("stepper"), nb::arg("tout"), nb::arg("vret"), nb::arg("tret"));

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
    [](SUNStepper stepper, int last_flag) -> std::tuple<SUNErrCode, int>
    {
        auto SUNStepper_GetLastFlag_adapt_modifiable_immutable_to_return = [](SUNStepper stepper, int last_flag) -> std::tuple<SUNErrCode, int>
        {
            int * last_flag_adapt_modifiable = & last_flag;

            SUNErrCode r = SUNStepper_GetLastFlag(stepper, last_flag_adapt_modifiable);
            return std::make_tuple(r, last_flag);
        };

        return SUNStepper_GetLastFlag_adapt_modifiable_immutable_to_return(stepper, last_flag);
    },     nb::arg("stepper"), nb::arg("last_flag"));

m.def("SUNStepper_GetNumSteps",
    SUNStepper_GetNumSteps, nb::arg("stepper"), nb::arg("nst"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif 
