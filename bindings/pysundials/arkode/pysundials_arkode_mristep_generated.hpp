// #ifndef _MRISTEP_H
//
// #ifdef __cplusplus
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
