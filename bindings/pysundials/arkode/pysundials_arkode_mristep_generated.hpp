// #ifndef _MRISTEP_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumMRISTEP_METHOD_TYPE =
  nb::enum_<MRISTEP_METHOD_TYPE>(m, "MRISTEP_METHOD_TYPE", nb::is_arithmetic(), "")
    .value("MRISTEP_EXPLICIT", MRISTEP_EXPLICIT, "")
    .value("MRISTEP_IMPLICIT", MRISTEP_IMPLICIT, "")
    .value("MRISTEP_IMEX", MRISTEP_IMEX, "")
    .value("MRISTEP_MERK", MRISTEP_MERK, "")
    .value("MRISTEP_SR", MRISTEP_SR, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyEnumARKODE_MRITableID =
  nb::enum_<ARKODE_MRITableID>(m, "ARKODE_MRITableID", nb::is_arithmetic(), "")
    .value("ARKODE_MRI_NONE", ARKODE_MRI_NONE, "")
    .value("ARKODE_MIN_MRI_NUM", ARKODE_MIN_MRI_NUM, "")
    .value("ARKODE_MIS_KW3", ARKODE_MIS_KW3, "")
    .value("ARKODE_MRI_GARK_ERK33a", ARKODE_MRI_GARK_ERK33a, "")
    .value("ARKODE_MRI_GARK_ERK45a", ARKODE_MRI_GARK_ERK45a, "")
    .value("ARKODE_MRI_GARK_IRK21a", ARKODE_MRI_GARK_IRK21a, "")
    .value("ARKODE_MRI_GARK_ESDIRK34a", ARKODE_MRI_GARK_ESDIRK34a, "")
    .value("ARKODE_MRI_GARK_ESDIRK46a", ARKODE_MRI_GARK_ESDIRK46a, "")
    .value("ARKODE_IMEX_MRI_GARK3a", ARKODE_IMEX_MRI_GARK3a, "")
    .value("ARKODE_IMEX_MRI_GARK3b", ARKODE_IMEX_MRI_GARK3b, "")
    .value("ARKODE_IMEX_MRI_GARK4", ARKODE_IMEX_MRI_GARK4, "")
    .value("ARKODE_MRI_GARK_FORWARD_EULER", ARKODE_MRI_GARK_FORWARD_EULER, "")
    .value("ARKODE_MRI_GARK_RALSTON2", ARKODE_MRI_GARK_RALSTON2, "")
    .value("ARKODE_MRI_GARK_ERK22a", ARKODE_MRI_GARK_ERK22a, "")
    .value("ARKODE_MRI_GARK_ERK22b", ARKODE_MRI_GARK_ERK22b, "")
    .value("ARKODE_MRI_GARK_RALSTON3", ARKODE_MRI_GARK_RALSTON3, "")
    .value("ARKODE_MRI_GARK_BACKWARD_EULER", ARKODE_MRI_GARK_BACKWARD_EULER, "")
    .value("ARKODE_MRI_GARK_IMPLICIT_MIDPOINT",
           ARKODE_MRI_GARK_IMPLICIT_MIDPOINT, "")
    .value("ARKODE_IMEX_MRI_GARK_EULER", ARKODE_IMEX_MRI_GARK_EULER, "")
    .value("ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL", ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL,
           "")
    .value("ARKODE_IMEX_MRI_GARK_MIDPOINT", ARKODE_IMEX_MRI_GARK_MIDPOINT, "")
    .value("ARKODE_MERK21", ARKODE_MERK21, "")
    .value("ARKODE_MERK32", ARKODE_MERK32, "")
    .value("ARKODE_MERK43", ARKODE_MERK43, "")
    .value("ARKODE_MERK54", ARKODE_MERK54, "")
    .value("ARKODE_IMEX_MRI_SR21", ARKODE_IMEX_MRI_SR21, "")
    .value("ARKODE_IMEX_MRI_SR32", ARKODE_IMEX_MRI_SR32, "")
    .value("ARKODE_IMEX_MRI_SR43", ARKODE_IMEX_MRI_SR43, "")
    .value("ARKODE_MAX_MRI_NUM", ARKODE_MAX_MRI_NUM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClassMRIStepCouplingMem =
  nb::class_<MRIStepCouplingMem>(m, "MRIStepCouplingMem", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("MRIStepCoupling_LoadTable", MRIStepCoupling_LoadTable, nb::arg("method"));

m.def("MRIStepCoupling_LoadTableByName", MRIStepCoupling_LoadTableByName,
      nb::arg("method"));

m.def("MRIStepCoupling_MIStoMRI", MRIStepCoupling_MIStoMRI, nb::arg("B"),
      nb::arg("q"), nb::arg("p"));

m.def("MRIStepCoupling_Copy", MRIStepCoupling_Copy, nb::arg("MRIC"));

m.def("MRIStepCoupling_Write", MRIStepCoupling_Write, nb::arg("MRIC"),
      nb::arg("outfile"));

m.def("MRIStepSetCoupling", MRIStepSetCoupling, nb::arg("arkode_mem"),
      nb::arg("MRIC"));

m.def("MRIStepSetPreInnerFn", MRIStepSetPreInnerFn, nb::arg("arkode_mem"),
      nb::arg("prefn"));

m.def("MRIStepSetPostInnerFn", MRIStepSetPostInnerFn, nb::arg("arkode_mem"),
      nb::arg("postfn"));

m.def(
  "MRIStepGetLastInnerStepFlag",
  [](void* arkode_mem, int flag) -> std::tuple<int, int>
  {
    auto MRIStepGetLastInnerStepFlag_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, int flag) -> std::tuple<int, int>
    {
      int* flag_adapt_modifiable = &flag;

      int r = MRIStepGetLastInnerStepFlag(arkode_mem, flag_adapt_modifiable);
      return std::make_tuple(r, flag);
    };

    return MRIStepGetLastInnerStepFlag_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                            flag);
  },
  nb::arg("arkode_mem"), nb::arg("flag"));

m.def(
  "MRIStepGetNumInnerStepperFails",
  [](void* arkode_mem, long inner_fails) -> std::tuple<int, long>
  {
    auto MRIStepGetNumInnerStepperFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, long inner_fails) -> std::tuple<int, long>
    {
      long* inner_fails_adapt_modifiable = &inner_fails;

      int r = MRIStepGetNumInnerStepperFails(arkode_mem,
                                             inner_fails_adapt_modifiable);
      return std::make_tuple(r, inner_fails);
    };

    return MRIStepGetNumInnerStepperFails_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                               inner_fails);
  },
  nb::arg("arkode_mem"), nb::arg("inner_fails"));

m.def("MRIStepInnerStepper_AddForcing", MRIStepInnerStepper_AddForcing,
      nb::arg("stepper"), nb::arg("t"), nb::arg("f"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
