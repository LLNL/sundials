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

m.def(
  "SUNContext_GetProfiler",
  [](SUNContext sunctx, SUNProfiler profiler) -> std::tuple<SUNErrCode, SUNProfiler>
  {
    auto SUNContext_GetProfiler_adapt_modifiable_immutable_to_return =
      [](SUNContext sunctx,
         SUNProfiler profiler) -> std::tuple<SUNErrCode, SUNProfiler>
    {
      SUNProfiler* profiler_adapt_modifiable = &profiler;

      SUNErrCode r = SUNContext_GetProfiler(sunctx, profiler_adapt_modifiable);
      return std::make_tuple(r, profiler);
    };

    return SUNContext_GetProfiler_adapt_modifiable_immutable_to_return(sunctx,
                                                                       profiler);
  },
  nb::arg("sunctx"), nb::arg("profiler"));

m.def("SUNContext_SetProfiler", SUNContext_SetProfiler, nb::arg("sunctx"),
      nb::arg("profiler"));

m.def(
  "SUNContext_GetLogger",
  [](SUNContext sunctx, SUNLogger logger) -> std::tuple<SUNErrCode, SUNLogger>
  {
    auto SUNContext_GetLogger_adapt_modifiable_immutable_to_return =
      [](SUNContext sunctx, SUNLogger logger) -> std::tuple<SUNErrCode, SUNLogger>
    {
      SUNLogger* logger_adapt_modifiable = &logger;

      SUNErrCode r = SUNContext_GetLogger(sunctx, logger_adapt_modifiable);
      return std::make_tuple(r, logger);
    };

    return SUNContext_GetLogger_adapt_modifiable_immutable_to_return(sunctx,
                                                                     logger);
  },
  nb::arg("sunctx"), nb::arg("logger"));

m.def("SUNContext_SetLogger", SUNContext_SetLogger, nb::arg("sunctx"),
      nb::arg("logger"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
