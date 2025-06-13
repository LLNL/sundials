// #ifndef _SUNDIALS_ADAPTCONTROLLER_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyClass_generic_SUNAdaptController_Ops =
    nb::class_<_generic_SUNAdaptController_Ops>
        (m, "_generic_SUNAdaptController_Ops", "Structure containing function pointers to controller operations")
    .def(nb::init<>()) // implicit default constructor 
    ;


auto pyClass_generic_SUNAdaptController =
    nb::class_<_generic_SUNAdaptController>
        (m, "_generic_SUNAdaptController", " A SUNAdaptController is a structure with an implementation-dependent\n   'content' field, and a pointer to a structure of\n   operations corresponding to that implementation.")
    .def("__init__", [](_generic_SUNAdaptController * self, SUNAdaptController_Ops ops = SUNAdaptController_Ops(), SUNContext sunctx = SUNContext())
    {
        new (self) _generic_SUNAdaptController();  // placement new
        auto r = self;
        r->ops = ops;
        r->sunctx = sunctx;
    },
    nb::arg("ops") = SUNAdaptController_Ops(), nb::arg("sunctx") = SUNContext()
    )
    .def_rw("content", &_generic_SUNAdaptController::content, "")
    .def_rw("ops", &_generic_SUNAdaptController::ops, "")
    .def_rw("sunctx", &_generic_SUNAdaptController::sunctx, "")
    ;


m.def("SUNAdaptController_GetType",
    SUNAdaptController_GetType, nb::arg("C"));

m.def("SUNAdaptController_EstimateStep",
    SUNAdaptController_EstimateStep, nb::arg("C"), nb::arg("h"), nb::arg("p"), nb::arg("dsm"), nb::arg("hnew"));

m.def("SUNAdaptController_EstimateStepTol",
    SUNAdaptController_EstimateStepTol, nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("P"), nb::arg("DSM"), nb::arg("dsm"), nb::arg("Hnew"), nb::arg("tolfacnew"));

m.def("SUNAdaptController_Reset",
    SUNAdaptController_Reset, nb::arg("C"));

m.def("SUNAdaptController_SetDefaults",
    SUNAdaptController_SetDefaults, nb::arg("C"));

m.def("SUNAdaptController_Write",
    SUNAdaptController_Write, nb::arg("C"), nb::arg("fptr"));

m.def("SUNAdaptController_SetErrorBias",
    SUNAdaptController_SetErrorBias, nb::arg("C"), nb::arg("bias"));

m.def("SUNAdaptController_UpdateH",
    SUNAdaptController_UpdateH, nb::arg("C"), nb::arg("h"), nb::arg("dsm"));

m.def("SUNAdaptController_UpdateMRIHTol",
    SUNAdaptController_UpdateMRIHTol, nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("DSM"), nb::arg("dsm"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif 
