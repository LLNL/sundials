// #ifndef _SUNDOMEIGEST_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClassSUNDomEigEstimator_Ops_ =
  nb::class_<SUNDomEigEstimator_Ops_>(m, "SUNDomEigEstimator_Ops_", "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClassSUNDomEigEstimator_ =
  nb::class_<SUNDomEigEstimator_>(m, "SUNDomEigEstimator_", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNDomEigEstimator_SetMaxIters", SUNDomEigEstimator_SetMaxIters,
      nb::arg("DEE"), nb::arg("max_iters"));

m.def("SUNDomEigEstimator_SetNumPreprocessIters",
      SUNDomEigEstimator_SetNumPreprocessIters, nb::arg("DEE"),
      nb::arg("num_iters"));

m.def("SUNDomEigEstimator_SetRelTol", SUNDomEigEstimator_SetRelTol,
      nb::arg("DEE"), nb::arg("tol"));

m.def("SUNDomEigEstimator_SetInitialGuess", SUNDomEigEstimator_SetInitialGuess,
      nb::arg("DEE"), nb::arg("q"));

m.def("SUNDomEigEstimator_Initialize", SUNDomEigEstimator_Initialize,
      nb::arg("DEE"));

m.def(
  "SUNDomEigEstimator_Estimate",
  [](SUNDomEigEstimator DEE, double lambdaR,
     double lambdaI) -> std::tuple<SUNErrCode, double, double>
  {
    auto SUNDomEigEstimator_Estimate_adapt_modifiable_immutable_to_return =
      [](SUNDomEigEstimator DEE, double lambdaR,
         double lambdaI) -> std::tuple<SUNErrCode, double, double>
    {
      double* lambdaR_adapt_modifiable = &lambdaR;
      double* lambdaI_adapt_modifiable = &lambdaI;

      SUNErrCode r = SUNDomEigEstimator_Estimate(DEE, lambdaR_adapt_modifiable,
                                                 lambdaI_adapt_modifiable);
      return std::make_tuple(r, lambdaR, lambdaI);
    };

    return SUNDomEigEstimator_Estimate_adapt_modifiable_immutable_to_return(DEE,
                                                                            lambdaR,
                                                                            lambdaI);
  },
  nb::arg("DEE"), nb::arg("lambdaR"), nb::arg("lambdaI"));

m.def(
  "SUNDomEigEstimator_GetRes",
  [](SUNDomEigEstimator DEE, double res) -> std::tuple<SUNErrCode, double>
  {
    auto SUNDomEigEstimator_GetRes_adapt_modifiable_immutable_to_return =
      [](SUNDomEigEstimator DEE, double res) -> std::tuple<SUNErrCode, double>
    {
      double* res_adapt_modifiable = &res;

      SUNErrCode r = SUNDomEigEstimator_GetRes(DEE, res_adapt_modifiable);
      return std::make_tuple(r, res);
    };

    return SUNDomEigEstimator_GetRes_adapt_modifiable_immutable_to_return(DEE,
                                                                          res);
  },
  nb::arg("DEE"), nb::arg("res"));

m.def(
  "SUNDomEigEstimator_GetNumIters",
  [](SUNDomEigEstimator DEE, long num_iters) -> std::tuple<SUNErrCode, long>
  {
    auto SUNDomEigEstimator_GetNumIters_adapt_modifiable_immutable_to_return =
      [](SUNDomEigEstimator DEE, long num_iters) -> std::tuple<SUNErrCode, long>
    {
      long* num_iters_adapt_modifiable = &num_iters;

      SUNErrCode r = SUNDomEigEstimator_GetNumIters(DEE,
                                                    num_iters_adapt_modifiable);
      return std::make_tuple(r, num_iters);
    };

    return SUNDomEigEstimator_GetNumIters_adapt_modifiable_immutable_to_return(DEE,
                                                                               num_iters);
  },
  nb::arg("DEE"), nb::arg("num_iters"));

m.def(
  "SUNDomEigEstimator_GetNumATimesCalls",
  [](SUNDomEigEstimator DEE, long num_ATimes) -> std::tuple<SUNErrCode, long>
  {
    auto SUNDomEigEstimator_GetNumATimesCalls_adapt_modifiable_immutable_to_return =
      [](SUNDomEigEstimator DEE, long num_ATimes) -> std::tuple<SUNErrCode, long>
    {
      long* num_ATimes_adapt_modifiable = &num_ATimes;

      SUNErrCode r =
        SUNDomEigEstimator_GetNumATimesCalls(DEE, num_ATimes_adapt_modifiable);
      return std::make_tuple(r, num_ATimes);
    };

    return SUNDomEigEstimator_GetNumATimesCalls_adapt_modifiable_immutable_to_return(DEE,
                                                                                     num_ATimes);
  },
  nb::arg("DEE"), nb::arg("num_ATimes"));

m.def("SUNDomEigEstimator_Write", SUNDomEigEstimator_Write, nb::arg("DEE"),
      nb::arg("outfile"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
