// #ifndef _SUNDIALS_CONTEXT_H
//
// #ifdef __cplusplus
// #endif
//

m.def("SUNContext_GetLastError", SUNContext_GetLastError, nb::arg("sunctx"));

m.def("SUNContext_PeekLastError", SUNContext_PeekLastError, nb::arg("sunctx"));

m.def("SUNContext_PushErrHandler", SUNContext_PushErrHandler, nb::arg("sunctx"),
      nb::arg("err_fn"), nb::arg("err_user_data"));

m.def("SUNContext_PopErrHandler", SUNContext_PopErrHandler, nb::arg("sunctx"));

m.def("SUNContext_ClearErrHandlers", SUNContext_ClearErrHandlers,
      nb::arg("sunctx"));

m.def("SUNContext_SetProfiler", SUNContext_SetProfiler, nb::arg("sunctx"),
      nb::arg("profiler"));

m.def("SUNContext_SetLogger", SUNContext_SetLogger, nb::arg("sunctx"),
      nb::arg("logger"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
