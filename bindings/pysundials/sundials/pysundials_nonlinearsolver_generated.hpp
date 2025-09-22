// #ifndef _SUNNONLINEARSOLVER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNNonlinearSolver_Type_ =
  nb::enum_<SUNNonlinearSolver_Type_>(m, "SUNNonlinearSolver_Type_",
                                      nb::is_arithmetic(), "")
    .value("SUNNONLINEARSOLVER_ROOTFIND", SUNNONLINEARSOLVER_ROOTFIND, "")
    .value("SUNNONLINEARSOLVER_FIXEDPOINT", SUNNONLINEARSOLVER_FIXEDPOINT, "")
    .export_values();

auto pyClass_generic_SUNNonlinearSolver_Ops =
  nb::class_<_generic_SUNNonlinearSolver_Ops>(m,
                                              "_generic_SUNNonlinearSolver_Ops",
                                              "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_SUNNonlinearSolver =
  nb::class_<_generic_SUNNonlinearSolver>(m, "_generic_SUNNonlinearSolver", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNNonlinSolGetType", SUNNonlinSolGetType, nb::arg("NLS"));

m.def("SUNNonlinSolInitialize", SUNNonlinSolInitialize, nb::arg("NLS"));

m.def("SUNNonlinSolSetMaxIters", SUNNonlinSolSetMaxIters, nb::arg("NLS"),
      nb::arg("maxiters"));

m.def(
  "SUNNonlinSolGetNumIters",
  [](SUNNonlinearSolver NLS, long niters) -> std::tuple<SUNErrCode, long>
  {
    auto SUNNonlinSolGetNumIters_adapt_modifiable_immutable_to_return =
      [](SUNNonlinearSolver NLS, long niters) -> std::tuple<SUNErrCode, long>
    {
      long* niters_adapt_modifiable = &niters;

      SUNErrCode r = SUNNonlinSolGetNumIters(NLS, niters_adapt_modifiable);
      return std::make_tuple(r, niters);
    };

    return SUNNonlinSolGetNumIters_adapt_modifiable_immutable_to_return(NLS,
                                                                        niters);
  },
  nb::arg("NLS"), nb::arg("niters"));

m.def(
  "SUNNonlinSolGetCurIter",
  [](SUNNonlinearSolver NLS, int iter) -> std::tuple<SUNErrCode, int>
  {
    auto SUNNonlinSolGetCurIter_adapt_modifiable_immutable_to_return =
      [](SUNNonlinearSolver NLS, int iter) -> std::tuple<SUNErrCode, int>
    {
      int* iter_adapt_modifiable = &iter;

      SUNErrCode r = SUNNonlinSolGetCurIter(NLS, iter_adapt_modifiable);
      return std::make_tuple(r, iter);
    };

    return SUNNonlinSolGetCurIter_adapt_modifiable_immutable_to_return(NLS, iter);
  },
  nb::arg("NLS"), nb::arg("iter"));

m.def(
  "SUNNonlinSolGetNumConvFails",
  [](SUNNonlinearSolver NLS, long nconvfails) -> std::tuple<SUNErrCode, long>
  {
    auto SUNNonlinSolGetNumConvFails_adapt_modifiable_immutable_to_return =
      [](SUNNonlinearSolver NLS, long nconvfails) -> std::tuple<SUNErrCode, long>
    {
      long* nconvfails_adapt_modifiable = &nconvfails;

      SUNErrCode r = SUNNonlinSolGetNumConvFails(NLS,
                                                 nconvfails_adapt_modifiable);
      return std::make_tuple(r, nconvfails);
    };

    return SUNNonlinSolGetNumConvFails_adapt_modifiable_immutable_to_return(NLS,
                                                                            nconvfails);
  },
  nb::arg("NLS"), nb::arg("nconvfails"));
m.attr("SUN_NLS_CONTINUE")   = +901;
m.attr("SUN_NLS_CONV_RECVR") = +902;
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
