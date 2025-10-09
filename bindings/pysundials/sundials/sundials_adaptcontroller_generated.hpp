// #ifndef _SUNDIALS_ADAPTCONTROLLER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNAdaptController_Type =
  nb::enum_<SUNAdaptController_Type>(m, "SUNAdaptController_Type",
                                     nb::is_arithmetic(), "")
    .value("SUN_ADAPTCONTROLLER_NONE", SUN_ADAPTCONTROLLER_NONE, "")
    .value("SUN_ADAPTCONTROLLER_H", SUN_ADAPTCONTROLLER_H, "")
    .value("SUN_ADAPTCONTROLLER_MRI_H_TOL", SUN_ADAPTCONTROLLER_MRI_H_TOL, "")
    .export_values();
// #ifndef SWIG
//
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
  [](SUNAdaptController C, double h, int p,
     double dsm) -> std::tuple<SUNErrCode, double>
  {
    auto SUNAdaptController_EstimateStep_adapt_modifiable_immutable_to_return =
      [](SUNAdaptController C, double h, int p,
         double dsm) -> std::tuple<SUNErrCode, double>
    {
      double hnew_adapt_modifiable;

      SUNErrCode r = SUNAdaptController_EstimateStep(C, h, p, dsm,
                                                     &hnew_adapt_modifiable);
      return std::make_tuple(r, hnew_adapt_modifiable);
    };

    return SUNAdaptController_EstimateStep_adapt_modifiable_immutable_to_return(C,
                                                                                h,
                                                                                p,
                                                                                dsm);
  },
  nb::arg("C"), nb::arg("h"), nb::arg("p"), nb::arg("dsm"));

m.def(
  "SUNAdaptController_EstimateStepTol",
  [](SUNAdaptController C, double H, double tolfac, int P, double DSM,
     double dsm) -> std::tuple<SUNErrCode, double, double>
  {
    auto SUNAdaptController_EstimateStepTol_adapt_modifiable_immutable_to_return =
      [](SUNAdaptController C, double H, double tolfac, int P, double DSM,
         double dsm) -> std::tuple<SUNErrCode, double, double>
    {
      double Hnew_adapt_modifiable;
      double tolfacnew_adapt_modifiable;

      SUNErrCode r =
        SUNAdaptController_EstimateStepTol(C, H, tolfac, P, DSM, dsm,
                                           &Hnew_adapt_modifiable,
                                           &tolfacnew_adapt_modifiable);
      return std::make_tuple(r, Hnew_adapt_modifiable,
                             tolfacnew_adapt_modifiable);
    };

    return SUNAdaptController_EstimateStepTol_adapt_modifiable_immutable_to_return(C,
                                                                                   H,
                                                                                   tolfac,
                                                                                   P,
                                                                                   DSM,
                                                                                   dsm);
  },
  nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("P"), nb::arg("DSM"),
  nb::arg("dsm"));

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
