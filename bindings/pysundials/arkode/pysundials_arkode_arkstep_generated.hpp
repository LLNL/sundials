// #ifndef _ARKSTEP_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("ARKStepReInit",
    ARKStepReInit, nb::arg("arkode_mem"), nb::arg("fe"), nb::arg("fi"), nb::arg("t0"), nb::arg("y0"));

m.def("ARKStepSetExplicit",
    ARKStepSetExplicit, nb::arg("arkode_mem"));

m.def("ARKStepSetImplicit",
    ARKStepSetImplicit, nb::arg("arkode_mem"));

m.def("ARKStepSetImEx",
    ARKStepSetImEx, nb::arg("arkode_mem"));

m.def("ARKStepSetTables",
    ARKStepSetTables, nb::arg("arkode_mem"), nb::arg("q"), nb::arg("p"), nb::arg("Bi"), nb::arg("Be"));

m.def("ARKStepSetTableNum",
    ARKStepSetTableNum, nb::arg("arkode_mem"), nb::arg("itable"), nb::arg("etable"));

m.def("ARKStepSetTableName",
    ARKStepSetTableName, nb::arg("arkode_mem"), nb::arg("itable"), nb::arg("etable"));

m.def("ARKStepGetTimestepperStats",
    ARKStepGetTimestepperStats, nb::arg("arkode_mem"), nb::arg("expsteps"), nb::arg("accsteps"), nb::arg("step_attempts"), nb::arg("nfe_evals"), nb::arg("nfi_evals"), nb::arg("nlinsetups"), nb::arg("netfails"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
