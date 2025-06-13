// #ifndef _SUNDIALS_PROFILER_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("SUNProfiler_Begin",
    SUNProfiler_Begin, nb::arg("p"), nb::arg("name"));

m.def("SUNProfiler_End",
    SUNProfiler_End, nb::arg("p"), nb::arg("name"));

m.def("SUNProfiler_GetTimerResolution",
    SUNProfiler_GetTimerResolution, nb::arg("p"), nb::arg("resolution"));

m.def("SUNProfiler_GetElapsedTime",
    SUNProfiler_GetElapsedTime, nb::arg("p"), nb::arg("name"), nb::arg("time"));

m.def("SUNProfiler_Print",
    SUNProfiler_Print, nb::arg("p"), nb::arg("fp"));

m.def("SUNProfiler_Reset",
    SUNProfiler_Reset, nb::arg("p"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif 
