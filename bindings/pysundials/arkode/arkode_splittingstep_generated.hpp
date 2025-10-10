// #ifndef ARKODE_SPLITTINGSTEP_H_
//
// #ifdef __cplusplus
// #endif
//

auto pyClassSplittingStepCoefficientsMem =
  nb::class_<SplittingStepCoefficientsMem>(m, "SplittingStepCoefficientsMem", "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyEnumARKODE_SplittingCoefficientsID =
  nb::enum_<ARKODE_SplittingCoefficientsID>(m, "ARKODE_SplittingCoefficientsID",
                                            nb::is_arithmetic(), "")
    .value("ARKODE_SPLITTING_NONE", ARKODE_SPLITTING_NONE, "")
    .value("ARKODE_SPLITTING_LIE_TROTTER_1_1_2",
           ARKODE_SPLITTING_LIE_TROTTER_1_1_2, "")
    .value("ARKODE_MIN_SPLITTING_NUM", ARKODE_MIN_SPLITTING_NUM, "")
    .value("ARKODE_SPLITTING_STRANG_2_2_2", ARKODE_SPLITTING_STRANG_2_2_2, "")
    .value("ARKODE_SPLITTING_BEST_2_2_2", ARKODE_SPLITTING_BEST_2_2_2, "")
    .value("ARKODE_SPLITTING_SUZUKI_3_3_2", ARKODE_SPLITTING_SUZUKI_3_3_2, "")
    .value("ARKODE_SPLITTING_RUTH_3_3_2", ARKODE_SPLITTING_RUTH_3_3_2, "")
    .value("ARKODE_SPLITTING_YOSHIDA_4_4_2", ARKODE_SPLITTING_YOSHIDA_4_4_2, "")
    .value("ARKODE_SPLITTING_YOSHIDA_8_6_2", ARKODE_SPLITTING_YOSHIDA_8_6_2, "")
    .value("ARKODE_MAX_SPLITTING_NUM", ARKODE_MAX_SPLITTING_NUM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

m.def(
  "SplittingStepCoefficients_Create",
  [](int sequential_methods, int stages, int partitions, int order,
     nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> alpha_1d,
     nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> beta_1d)
    -> SplittingStepCoefficients
  {
    auto SplittingStepCoefficients_Create_adapt_arr_ptr_to_std_vector =
      [](int sequential_methods, int stages, int partitions, int order,
         nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> alpha_1d,
         nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> beta_1d)
      -> SplittingStepCoefficients
    {
      sunrealtype* alpha_1d_ptr = reinterpret_cast<sunrealtype*>(alpha_1d.data());
      sunrealtype* beta_1d_ptr = reinterpret_cast<sunrealtype*>(beta_1d.data());

      auto lambda_result =
        SplittingStepCoefficients_Create(sequential_methods, stages, partitions,
                                         order, alpha_1d_ptr, beta_1d_ptr);
      return lambda_result;
    };

    return SplittingStepCoefficients_Create_adapt_arr_ptr_to_std_vector(sequential_methods,
                                                                        stages,
                                                                        partitions,
                                                                        order,
                                                                        alpha_1d,
                                                                        beta_1d);
  },
  nb::arg("sequential_methods"), nb::arg("stages"), nb::arg("partitions"),
  nb::arg("order"), nb::arg("alpha_1d"), nb::arg("beta_1d"),
  nb::rv_policy::reference);

m.def("SplittingStepCoefficients_Copy", SplittingStepCoefficients_Copy,
      nb::arg("coefficients"), nb::rv_policy::reference);

m.def("SplittingStepCoefficients_Write", SplittingStepCoefficients_Write,
      nb::arg("coefficients"), nb::arg("outfile"));

m.def("SplittingStepCoefficients_LoadCoefficients",
      SplittingStepCoefficients_LoadCoefficients, nb::arg("id"),
      nb::rv_policy::reference);

m.def("SplittingStepCoefficients_LoadCoefficientsByName",
      SplittingStepCoefficients_LoadCoefficientsByName, nb::arg("name"),
      nb::rv_policy::reference);

m.def("SplittingStepCoefficients_IDToName", SplittingStepCoefficients_IDToName,
      nb::arg("id"), nb::rv_policy::reference);

m.def("SplittingStepCoefficients_LieTrotter",
      SplittingStepCoefficients_LieTrotter, nb::arg("partitions"),
      nb::rv_policy::reference);

m.def("SplittingStepCoefficients_Strang", SplittingStepCoefficients_Strang,
      nb::arg("partitions"), nb::rv_policy::reference);

m.def("SplittingStepCoefficients_Parallel", SplittingStepCoefficients_Parallel,
      nb::arg("partitions"), nb::rv_policy::reference);

m.def("SplittingStepCoefficients_SymmetricParallel",
      SplittingStepCoefficients_SymmetricParallel, nb::arg("partitions"),
      nb::rv_policy::reference);

m.def("SplittingStepCoefficients_ThirdOrderSuzuki",
      SplittingStepCoefficients_ThirdOrderSuzuki, nb::arg("partitions"),
      nb::rv_policy::reference);

m.def("SplittingStepCoefficients_TripleJump",
      SplittingStepCoefficients_TripleJump, nb::arg("partitions"),
      nb::arg("order"), nb::rv_policy::reference);

m.def("SplittingStepCoefficients_SuzukiFractal",
      SplittingStepCoefficients_SuzukiFractal, nb::arg("partitions"),
      nb::arg("order"), nb::rv_policy::reference);

m.def("SplittingStepSetCoefficients", SplittingStepSetCoefficients,
      nb::arg("arkode_mem"), nb::arg("coefficients"));

m.def(
  "SplittingStepGetNumEvolves",
  [](void* arkode_mem, int partition) -> std::tuple<int, long>
  {
    auto SplittingStepGetNumEvolves_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, int partition) -> std::tuple<int, long>
    {
      long evolves_adapt_modifiable;

      int r = SplittingStepGetNumEvolves(arkode_mem, partition,
                                         &evolves_adapt_modifiable);
      return std::make_tuple(r, evolves_adapt_modifiable);
    };

    return SplittingStepGetNumEvolves_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                           partition);
  },
  nb::arg("arkode_mem"), nb::arg("partition"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
