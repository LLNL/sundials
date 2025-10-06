// #ifndef _SUNDIALS_PROFILER_H
//
// #ifdef __cplusplus
// #endif
//

m.def("SUNProfiler_Begin", SUNProfiler_Begin, nb::arg("p"), nb::arg("name"));

m.def("SUNProfiler_End", SUNProfiler_End, nb::arg("p"), nb::arg("name"));

m.def(
  "SUNProfiler_GetTimerResolution",
  [](SUNProfiler p) -> std::tuple<SUNErrCode, double>
  {
    auto SUNProfiler_GetTimerResolution_adapt_modifiable_immutable_to_return =
      [](SUNProfiler p) -> std::tuple<SUNErrCode, double>
    {
      double resolution_adapt_modifiable;

      SUNErrCode r =
        SUNProfiler_GetTimerResolution(p, &resolution_adapt_modifiable);
      return std::make_tuple(r, resolution_adapt_modifiable);
    };

    return SUNProfiler_GetTimerResolution_adapt_modifiable_immutable_to_return(p);
  },
  nb::arg("p"));

m.def(
  "SUNProfiler_GetElapsedTime",
  [](SUNProfiler p, const char* name) -> std::tuple<SUNErrCode, double>
  {
    auto SUNProfiler_GetElapsedTime_adapt_modifiable_immutable_to_return =
      [](SUNProfiler p, const char* name) -> std::tuple<SUNErrCode, double>
    {
      double time_adapt_modifiable;

      SUNErrCode r = SUNProfiler_GetElapsedTime(p, name, &time_adapt_modifiable);
      return std::make_tuple(r, time_adapt_modifiable);
    };

    return SUNProfiler_GetElapsedTime_adapt_modifiable_immutable_to_return(p,
                                                                           name);
  },
  nb::arg("p"), nb::arg("name"));

m.def("SUNProfiler_Print", SUNProfiler_Print, nb::arg("p"), nb::arg("fp"));

m.def("SUNProfiler_Reset", SUNProfiler_Reset, nb::arg("p"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
