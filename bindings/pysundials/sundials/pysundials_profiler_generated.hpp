// #ifndef _SUNDIALS_PROFILER_H
//
// #ifdef __cplusplus
// #endif
//

m.def("SUNProfiler_Begin", SUNProfiler_Begin, nb::arg("p"), nb::arg("name"));

m.def("SUNProfiler_End", SUNProfiler_End, nb::arg("p"), nb::arg("name"));

m.def(
  "SUNProfiler_GetTimerResolution",
  [](SUNProfiler p, double resolution) -> std::tuple<SUNErrCode, double>
  {
    auto SUNProfiler_GetTimerResolution_adapt_modifiable_immutable_to_return =
      [](SUNProfiler p, double resolution) -> std::tuple<SUNErrCode, double>
    {
      double* resolution_adapt_modifiable = &resolution;

      SUNErrCode r =
        SUNProfiler_GetTimerResolution(p, resolution_adapt_modifiable);
      return std::make_tuple(r, resolution);
    };

    return SUNProfiler_GetTimerResolution_adapt_modifiable_immutable_to_return(p,
                                                                               resolution);
  },
  nb::arg("p"), nb::arg("resolution"));

m.def(
  "SUNProfiler_GetElapsedTime",
  [](SUNProfiler p, const char* name, double time) -> std::tuple<SUNErrCode, double>
  {
    auto SUNProfiler_GetElapsedTime_adapt_modifiable_immutable_to_return =
      [](SUNProfiler p, const char* name,
         double time) -> std::tuple<SUNErrCode, double>
    {
      double* time_adapt_modifiable = &time;

      SUNErrCode r = SUNProfiler_GetElapsedTime(p, name, time_adapt_modifiable);
      return std::make_tuple(r, time);
    };

    return SUNProfiler_GetElapsedTime_adapt_modifiable_immutable_to_return(p,
                                                                           name,
                                                                           time);
  },
  nb::arg("p"), nb::arg("name"), nb::arg("time"));

m.def("SUNProfiler_Print", SUNProfiler_Print, nb::arg("p"), nb::arg("fp"));

m.def("SUNProfiler_Reset", SUNProfiler_Reset, nb::arg("p"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
