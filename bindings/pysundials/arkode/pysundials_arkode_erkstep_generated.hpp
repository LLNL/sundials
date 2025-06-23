// #ifndef _ERKSTEP_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("ERKStepReInit",
    ERKStepReInit, nb::arg("arkode_mem"), nb::arg("f"), nb::arg("t0"), nb::arg("y0"));

m.def("ERKStepSetTable",
    ERKStepSetTable, nb::arg("arkode_mem"), nb::arg("B"));

m.def("ERKStepSetTableNum",
    ERKStepSetTableNum, nb::arg("arkode_mem"), nb::arg("etable"));

m.def("ERKStepSetTableName",
    ERKStepSetTableName, nb::arg("arkode_mem"), nb::arg("etable"));

m.def("ERKStepGetTimestepperStats",
    ERKStepGetTimestepperStats, nb::arg("arkode_mem"), nb::arg("expsteps"), nb::arg("accsteps"), nb::arg("step_attempts"), nb::arg("nfevals"), nb::arg("netfails"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
