// #ifndef _SUNADJOINT_CHECKPOINTSCHEME_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("SUNAdjointCheckpointScheme_NewEmpty",
    SUNAdjointCheckpointScheme_NewEmpty, nb::arg("sunctx"), nb::arg("param_1"));

m.def("SUNAdjointCheckpointScheme_SetNeedsSavingFn",
    SUNAdjointCheckpointScheme_SetNeedsSavingFn, nb::arg("check_scheme"), nb::arg("param_1"));

m.def("SUNAdjointCheckpointScheme_SetInsertVectorFn",
    SUNAdjointCheckpointScheme_SetInsertVectorFn, nb::arg("check_scheme"), nb::arg("param_1"));

m.def("SUNAdjointCheckpointScheme_SetLoadVectorFn",
    SUNAdjointCheckpointScheme_SetLoadVectorFn, nb::arg("check_scheme"), nb::arg("param_1"));

m.def("SUNAdjointCheckpointScheme_SetEnableDenseFn",
    SUNAdjointCheckpointScheme_SetEnableDenseFn, nb::arg("check_scheme"), nb::arg("param_1"));

m.def("SUNAdjointCheckpointScheme_SetContent",
    SUNAdjointCheckpointScheme_SetContent, nb::arg("check_scheme"), nb::arg("content"));

m.def("SUNAdjointCheckpointScheme_NeedsSaving",
    SUNAdjointCheckpointScheme_NeedsSaving, nb::arg("check_scheme"), nb::arg("step_num"), nb::arg("stage_num"), nb::arg("t"), nb::arg("yes_or_no"));

m.def("SUNAdjointCheckpointScheme_InsertVector",
    SUNAdjointCheckpointScheme_InsertVector, nb::arg("check_scheme"), nb::arg("step_num"), nb::arg("stage_num"), nb::arg("t"), nb::arg("state"));

m.def("SUNAdjointCheckpointScheme_LoadVector",
    SUNAdjointCheckpointScheme_LoadVector, nb::arg("check_scheme"), nb::arg("step_num"), nb::arg("stage_num"), nb::arg("peek"), nb::arg("out"), nb::arg("tout"));

m.def("SUNAdjointCheckpointScheme_EnableDense",
    SUNAdjointCheckpointScheme_EnableDense, nb::arg("check_scheme"), nb::arg("on_or_off"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif 
