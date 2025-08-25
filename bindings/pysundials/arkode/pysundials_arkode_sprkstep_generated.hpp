// #ifndef _ARKODE_SPRKSTEP_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("SPRKStepSetUseCompensatedSums",
    SPRKStepSetUseCompensatedSums, nb::arg("arkode_mem"), nb::arg("onoff"));

m.def("SPRKStepSetMethod",
    SPRKStepSetMethod, nb::arg("arkode_mem"), nb::arg("sprk_storage"));

m.def("SPRKStepSetMethodName",
    SPRKStepSetMethodName, nb::arg("arkode_mem"), nb::arg("method"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
