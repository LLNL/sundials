// #ifndef _SUNDIALS_ADAPTCONTROLLER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_generic_SUNAdaptController_Ops =
  nb::class_<_generic_SUNAdaptController_Ops>(m,
                                              "_generic_SUNAdaptController_Ops",
                                              "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_SUNAdaptController =
  nb::class_<_generic_SUNAdaptController>(m, "_generic_SUNAdaptController", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNAdaptController_GetType", SUNAdaptController_GetType, nb::arg("C"));

m.def(
  "SUNAdaptController_EstimateStep",
  [](SUNAdaptController C, double h, int p, double dsm,
     double hnew) -> std::tuple<SUNErrCode, double>
  {
    auto SUNAdaptController_EstimateStep_adapt_modifiable_immutable_to_return =
      [](SUNAdaptController C, double h, int p, double dsm,
         double hnew) -> std::tuple<SUNErrCode, double>
    {
      double* hnew_adapt_modifiable = &hnew;

      SUNErrCode r = SUNAdaptController_EstimateStep(C, h, p, dsm,
                                                     hnew_adapt_modifiable);
      return std::make_tuple(r, hnew);
    };

    return SUNAdaptController_EstimateStep_adapt_modifiable_immutable_to_return(C,
                                                                                h,
                                                                                p,
                                                                                dsm,
                                                                                hnew);
  },
  nb::arg("C"), nb::arg("h"), nb::arg("p"), nb::arg("dsm"), nb::arg("hnew"));

m.def(
  "SUNAdaptController_EstimateStepTol",
  [](SUNAdaptController C, double H, double tolfac, int P, double DSM, double dsm,
     double Hnew, double tolfacnew) -> std::tuple<SUNErrCode, double, double>
  {
    auto SUNAdaptController_EstimateStepTol_adapt_modifiable_immutable_to_return =
      [](SUNAdaptController C, double H, double tolfac, int P, double DSM,
         double dsm, double Hnew,
         double tolfacnew) -> std::tuple<SUNErrCode, double, double>
    {
      double* Hnew_adapt_modifiable      = &Hnew;
      double* tolfacnew_adapt_modifiable = &tolfacnew;

      SUNErrCode r =
        SUNAdaptController_EstimateStepTol(C, H, tolfac, P, DSM, dsm,
                                           Hnew_adapt_modifiable,
                                           tolfacnew_adapt_modifiable);
      return std::make_tuple(r, Hnew, tolfacnew);
    };

    return SUNAdaptController_EstimateStepTol_adapt_modifiable_immutable_to_return(C,
                                                                                   H,
                                                                                   tolfac,
                                                                                   P,
                                                                                   DSM,
                                                                                   dsm,
                                                                                   Hnew,
                                                                                   tolfacnew);
  },
  nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("P"), nb::arg("DSM"),
  nb::arg("dsm"), nb::arg("Hnew"), nb::arg("tolfacnew"));

m.def("SUNAdaptController_Reset", SUNAdaptController_Reset, nb::arg("C"));

m.def("SUNAdaptController_SetDefaults", SUNAdaptController_SetDefaults,
      nb::arg("C"));

m.def("SUNAdaptController_Write", SUNAdaptController_Write, nb::arg("C"),
      nb::arg("fptr"));

m.def("SUNAdaptController_SetErrorBias", SUNAdaptController_SetErrorBias,
      nb::arg("C"), nb::arg("bias"));

m.def("SUNAdaptController_UpdateH", SUNAdaptController_UpdateH, nb::arg("C"),
      nb::arg("h"), nb::arg("dsm"));

m.def("SUNAdaptController_UpdateMRIHTol", SUNAdaptController_UpdateMRIHTol,
      nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("DSM"),
      nb::arg("dsm"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
