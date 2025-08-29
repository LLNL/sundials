// #ifndef ARKODE_SPLITTINGSTEP_H_
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyClassSplittingStepCoefficientsMem =
    nb::class_<SplittingStepCoefficientsMem>
        (m, "SplittingStepCoefficientsMem", "")
    .def(nb::init<>()) // implicit default constructor
    ;


auto pyEnumARKODE_SplittingCoefficientsID =
    nb::enum_<ARKODE_SplittingCoefficientsID>(m, "ARKODE_SplittingCoefficientsID", nb::is_arithmetic(), "")
        .value("ARKODE_SPLITTING_NONE", ARKODE_SPLITTING_NONE, "")
        .value("ARKODE_MIN_SPLITTING_NUM", ARKODE_MIN_SPLITTING_NUM, "")
        .value("ARKODE_SPLITTING_LIE_TROTTER_1_1_2", ARKODE_SPLITTING_LIE_TROTTER_1_1_2, "")
        .value("ARKODE_SPLITTING_STRANG_2_2_2", ARKODE_SPLITTING_STRANG_2_2_2, "")
        .value("ARKODE_SPLITTING_BEST_2_2_2", ARKODE_SPLITTING_BEST_2_2_2, "")
        .value("ARKODE_SPLITTING_SUZUKI_3_3_2", ARKODE_SPLITTING_SUZUKI_3_3_2, "")
        .value("ARKODE_SPLITTING_RUTH_3_3_2", ARKODE_SPLITTING_RUTH_3_3_2, "")
        .value("ARKODE_SPLITTING_YOSHIDA_4_4_2", ARKODE_SPLITTING_YOSHIDA_4_4_2, "")
        .value("ARKODE_SPLITTING_YOSHIDA_8_6_2", ARKODE_SPLITTING_YOSHIDA_8_6_2, "")
        .value("ARKODE_MAX_SPLITTING_NUM", ARKODE_MAX_SPLITTING_NUM, "")
    .export_values();


m.def("SplittingStepCoefficients_Create",
    [](int sequential_methods, int stages, int partitions, int order, std::vector<sunrealtype> alpha, std::vector<sunrealtype> beta) -> SplittingStepCoefficients
    {
        auto SplittingStepCoefficients_Create_adapt_arr_ptr_to_std_vector = [](int sequential_methods, int stages, int partitions, int order, std::vector<sunrealtype> alpha, std::vector<sunrealtype> beta) -> SplittingStepCoefficients
        {
            sunrealtype* alpha_ptr = reinterpret_cast<sunrealtype*>( alpha.empty() ? nullptr : alpha.data() );
            sunrealtype* beta_ptr = reinterpret_cast<sunrealtype*>( beta.empty() ? nullptr : beta.data() );

            auto lambda_result = SplittingStepCoefficients_Create(sequential_methods, stages, partitions, order, alpha_ptr, beta_ptr);
            return lambda_result;
        };

        return SplittingStepCoefficients_Create_adapt_arr_ptr_to_std_vector(sequential_methods, stages, partitions, order, alpha, beta);
    },     nb::arg("sequential_methods"), nb::arg("stages"), nb::arg("partitions"), nb::arg("order"), nb::arg("alpha"), nb::arg("beta"));

m.def("SplittingStepCoefficients_Copy",
    SplittingStepCoefficients_Copy, nb::arg("coefficients"));

m.def("SplittingStepCoefficients_LoadCoefficients",
    SplittingStepCoefficients_LoadCoefficients, nb::arg("id"));

m.def("SplittingStepCoefficients_LoadCoefficientsByName",
    SplittingStepCoefficients_LoadCoefficientsByName, nb::arg("name"));

m.def("SplittingStepCoefficients_IDToName",
    SplittingStepCoefficients_IDToName, nb::arg("id"));

m.def("SplittingStepCoefficients_LieTrotter",
    SplittingStepCoefficients_LieTrotter, nb::arg("partitions"));

m.def("SplittingStepCoefficients_Strang",
    SplittingStepCoefficients_Strang, nb::arg("partitions"));

m.def("SplittingStepCoefficients_Parallel",
    SplittingStepCoefficients_Parallel, nb::arg("partitions"));

m.def("SplittingStepCoefficients_SymmetricParallel",
    SplittingStepCoefficients_SymmetricParallel, nb::arg("partitions"));

m.def("SplittingStepCoefficients_ThirdOrderSuzuki",
    SplittingStepCoefficients_ThirdOrderSuzuki, nb::arg("partitions"));

m.def("SplittingStepCoefficients_TripleJump",
    SplittingStepCoefficients_TripleJump, nb::arg("partitions"), nb::arg("order"));

m.def("SplittingStepCoefficients_SuzukiFractal",
    SplittingStepCoefficients_SuzukiFractal, nb::arg("partitions"), nb::arg("order"));

m.def("SplittingStepSetCoefficients",
    SplittingStepSetCoefficients, nb::arg("arkode_mem"), nb::arg("coefficients"));

m.def("SplittingStepGetNumEvolves",
    [](void * arkode_mem, int partition, long evolves) -> std::tuple<int, long>
    {
        auto SplittingStepGetNumEvolves_adapt_modifiable_immutable_to_return = [](void * arkode_mem, int partition, long evolves) -> std::tuple<int, long>
        {
            long * evolves_adapt_modifiable = & evolves;

            int r = SplittingStepGetNumEvolves(arkode_mem, partition, evolves_adapt_modifiable);
            return std::make_tuple(r, evolves);
        };

        return SplittingStepGetNumEvolves_adapt_modifiable_immutable_to_return(arkode_mem, partition, evolves);
    },     nb::arg("arkode_mem"), nb::arg("partition"), nb::arg("evolves"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
