// #ifndef _SUNDIALS_LOGGER_H
//
// #ifdef __cplusplus
// #endif
//

m.def("SUNLogger_SetErrorFilename", SUNLogger_SetErrorFilename,
      nb::arg("logger"), nb::arg("error_filename"));

m.def("SUNLogger_SetWarningFilename", SUNLogger_SetWarningFilename,
      nb::arg("logger"), nb::arg("warning_filename"));

m.def("SUNLogger_SetDebugFilename", SUNLogger_SetDebugFilename,
      nb::arg("logger"), nb::arg("debug_filename"));

m.def("SUNLogger_SetInfoFilename", SUNLogger_SetInfoFilename, nb::arg("logger"),
      nb::arg("info_filename"));

m.def(
  "SUNLogger_QueueMsg",
  [](SUNLogger logger, SUNLogLevel lvl, const char* scope, const char* label,
     const char* msg_txt) -> SUNErrCode
  {
    auto SUNLogger_QueueMsg_adapt_variadic_format =
      [](SUNLogger logger, SUNLogLevel lvl, const char* scope,
         const char* label, const char* msg_txt) -> SUNErrCode
    {
      auto lambda_result = SUNLogger_QueueMsg(logger, lvl, scope, label, "%s",
                                              msg_txt);
      return lambda_result;
    };

    return SUNLogger_QueueMsg_adapt_variadic_format(logger, lvl, scope, label,
                                                    msg_txt);
  },
  nb::arg("logger"), nb::arg("lvl"), nb::arg("scope"), nb::arg("label"),
  nb::arg("msg_txt"));

m.def("SUNLogger_Flush", SUNLogger_Flush, nb::arg("logger"), nb::arg("lvl"));

m.def("SUNLogger_GetOutputRank", SUNLogger_GetOutputRank, nb::arg("logger"),
      nb::arg("output_rank"));
// #ifdef __cplusplus
// #endif
//
// #endif
