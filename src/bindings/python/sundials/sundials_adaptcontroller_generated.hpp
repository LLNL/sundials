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


m.def("SUNAdaptController_NewEmpty",
    SUNAdaptController_NewEmpty, 
    nb::arg("sunctx"), 
    "Function to create an empty SUNAdaptController data structure.");

m.def("SUNAdaptController_GetType",
    SUNAdaptController_GetType, 
    nb::arg("C"), 
    "Function to report the type of a SUNAdaptController object.");

m.def("SUNAdaptController_EstimateStep",
    SUNAdaptController_EstimateStep, 
    nb::arg("C"), nb::arg("h"), nb::arg("p"), nb::arg("dsm"), nb::arg("hnew"), 
    " Main step size controller function.  This is called following\n   a time step with size 'h' and local error factor 'dsm', and the\n   controller should estimate 'hnew' so that the ensuing step\n   will have 'dsm' value JUST BELOW 1.\n\n   Any return value other than SUN_SUCCESS will be treated as\n   an unrecoverable failure.");

m.def("SUNAdaptController_EstimateStepTol",
    SUNAdaptController_EstimateStepTol, 
    nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("P"), nb::arg("DSM"), nb::arg("dsm"), nb::arg("Hnew"), nb::arg("tolfacnew"), 
    " Combined slow step/fast tolerance multirate controller function.\n   This is called following a slow multirate time step with size 'H'\n   and fast/slow relative tolerance ratio 'tolfac', and error factors\n   'DSM' and 'dsm' (slow and fast, resp.).  The controller should\n   estimate slow stepsize 'Hnew' and updated relative tolerance ratio\n   'tolfacnew', so that the ensuing step will have 'DSM' and 'dsm'\n   values JUST BELOW 1 with minimal computational effort.");

m.def("SUNAdaptController_Reset",
    SUNAdaptController_Reset, 
    nb::arg("C"), 
    " Function to reset the controller to its initial state, e.g., if\n   it stores a small number of previous dsm or step size values.");

m.def("SUNAdaptController_SetDefaults",
    SUNAdaptController_SetDefaults, 
    nb::arg("C"), 
    "Function to set the controller parameters to their default values.");

m.def("SUNAdaptController_Write",
    SUNAdaptController_Write, 
    nb::arg("C"), nb::arg("fptr"), 
    " Function to write all controller parameters to the indicated file\n   pointer.");

m.def("SUNAdaptController_SetErrorBias",
    SUNAdaptController_SetErrorBias, 
    nb::arg("C"), nb::arg("bias"), 
    " Function to set an error bias factor to use for scaling the local error\n   'dsm' factors above.");

m.def("SUNAdaptController_UpdateH",
    SUNAdaptController_UpdateH, 
    nb::arg("C"), nb::arg("h"), nb::arg("dsm"), 
    " Function to notify a controller of type SUN_ADAPTCONTROLLER_H that\n   a successful time step was taken with stepsize h and local error factor\n   dsm, indicating that these can be saved for subsequent controller functions.");

m.def("SUNAdaptController_UpdateMRIHTol",
    SUNAdaptController_UpdateMRIHTol, 
    nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("DSM"), nb::arg("dsm"), 
    " Function to notify the controller of a successful multirate time step\n   with size H and fast tolerance factor tolfac, and local error factors\n   DSM and dsm, indicating that the step size, tolerance factor, or local\n   error factors can be saved for subsequent controller functions.");
// #ifdef __cplusplus
// 
// #endif
// 
// #endif 
