// #ifndef _SUNADJOINT_CHECKPOINTSCHEME_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("SUNAdjointCheckpointScheme_NeedsSaving",
    SUNAdjointCheckpointScheme_NeedsSaving, nb::arg("check_scheme"), nb::arg("step_num"), nb::arg("stage_num"), nb::arg("t"), nb::arg("yes_or_no"));

m.def("SUNAdjointCheckpointScheme_EnableDense",
    SUNAdjointCheckpointScheme_EnableDense, nb::arg("check_scheme"), nb::arg("on_or_off"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif 
