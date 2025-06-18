// #ifndef _SUNDIALS_ERRORS_H
// 


auto pyEnumSUNErrCode_ =
    nb::enum_<SUNErrCode_>(m, "SUNErrCode_", nb::is_arithmetic(), "")
        .value("SUN_ERR_MINIMUM", SUN_ERR_MINIMUM, "")
        .value("SUN_ERR_MAXIMUM", SUN_ERR_MAXIMUM, "")
        .value("SUN_SUCCESS", SUN_SUCCESS, "");
// #ifdef __cplusplus 
// #endif
// 

m.def("SUNLogErrHandlerFn",
    SUNLogErrHandlerFn, nb::arg("line"), nb::arg("func"), nb::arg("file"), nb::arg("msg"), nb::arg("err_code"), nb::arg("err_user_data"), nb::arg("sunctx"));

m.def("SUNAbortErrHandlerFn",
    SUNAbortErrHandlerFn, nb::arg("line"), nb::arg("func"), nb::arg("file"), nb::arg("msg"), nb::arg("err_code"), nb::arg("err_user_data"), nb::arg("sunctx"));

m.def("SUNGetErrMsg",
    SUNGetErrMsg, nb::arg("code"));
// #ifdef __cplusplus 
// #endif
// 
// #endif 
