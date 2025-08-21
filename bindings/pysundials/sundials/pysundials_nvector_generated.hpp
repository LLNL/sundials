// #ifndef _NVECTOR_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyEnumN_Vector_ID =
    nb::enum_<N_Vector_ID>(m, "N_Vector_ID", nb::is_arithmetic(), "")
        .value("SUNDIALS_NVEC_SERIAL", SUNDIALS_NVEC_SERIAL, "")
        .value("SUNDIALS_NVEC_PARALLEL", SUNDIALS_NVEC_PARALLEL, "")
        .value("SUNDIALS_NVEC_OPENMP", SUNDIALS_NVEC_OPENMP, "")
        .value("SUNDIALS_NVEC_PTHREADS", SUNDIALS_NVEC_PTHREADS, "")
        .value("SUNDIALS_NVEC_PARHYP", SUNDIALS_NVEC_PARHYP, "")
        .value("SUNDIALS_NVEC_PETSC", SUNDIALS_NVEC_PETSC, "")
        .value("SUNDIALS_NVEC_CUDA", SUNDIALS_NVEC_CUDA, "")
        .value("SUNDIALS_NVEC_HIP", SUNDIALS_NVEC_HIP, "")
        .value("SUNDIALS_NVEC_SYCL", SUNDIALS_NVEC_SYCL, "")
        .value("SUNDIALS_NVEC_RAJA", SUNDIALS_NVEC_RAJA, "")
        .value("SUNDIALS_NVEC_KOKKOS", SUNDIALS_NVEC_KOKKOS, "")
        .value("SUNDIALS_NVEC_OPENMPDEV", SUNDIALS_NVEC_OPENMPDEV, "")
        .value("SUNDIALS_NVEC_TRILINOS", SUNDIALS_NVEC_TRILINOS, "")
        .value("SUNDIALS_NVEC_MANYVECTOR", SUNDIALS_NVEC_MANYVECTOR, "")
        .value("SUNDIALS_NVEC_MPIMANYVECTOR", SUNDIALS_NVEC_MPIMANYVECTOR, "")
        .value("SUNDIALS_NVEC_MPIPLUSX", SUNDIALS_NVEC_MPIPLUSX, "")
        .value("SUNDIALS_NVEC_CUSTOM", SUNDIALS_NVEC_CUSTOM, "")
    .export_values();


auto pyClass_generic_N_Vector_Ops =
    nb::class_<_generic_N_Vector_Ops>
        (m, "_generic_N_Vector_Ops", "")
    .def(nb::init<>()) // implicit default constructor 
    ;


auto pyClass_generic_N_Vector =
    nb::class_<_generic_N_Vector>
        (m, "_generic_N_Vector", "")
    .def("__init__", [](_generic_N_Vector * self, N_Vector_Ops ops = N_Vector_Ops(), SUNContext sunctx = SUNContext())
    {
        new (self) _generic_N_Vector();  // placement new
        auto r_ctor_ = self;
        r_ctor_->ops = ops;
        r_ctor_->sunctx = sunctx;
    },
    nb::arg("ops") = N_Vector_Ops(), nb::arg("sunctx") = SUNContext()
    )
    .def_rw("content", &_generic_N_Vector::content, "")
    .def_rw("ops", &_generic_N_Vector::ops, "")
    .def_rw("sunctx", &_generic_N_Vector::sunctx, "")
    ;


m.def("N_VGetVectorID",
    N_VGetVectorID, nb::arg("w"));

m.def("N_VClone",
    N_VClone, nb::arg("w"));

m.def("N_VCloneEmpty",
    N_VCloneEmpty, nb::arg("w"));

m.def("N_VGetCommunicator",
    N_VGetCommunicator, nb::arg("v"));

m.def("N_VGetLength",
    N_VGetLength, nb::arg("v"));

m.def("N_VGetLocalLength",
    N_VGetLocalLength, nb::arg("v"));

m.def("N_VLinearSum",
    N_VLinearSum, nb::arg("a"), nb::arg("x"), nb::arg("b"), nb::arg("y"), nb::arg("z"));

m.def("N_VConst",
    N_VConst, nb::arg("c"), nb::arg("z"));

m.def("N_VProd",
    N_VProd, nb::arg("x"), nb::arg("y"), nb::arg("z"));

m.def("N_VDiv",
    N_VDiv, nb::arg("x"), nb::arg("y"), nb::arg("z"));

m.def("N_VScale",
    N_VScale, nb::arg("c"), nb::arg("x"), nb::arg("z"));

m.def("N_VAbs",
    N_VAbs, nb::arg("x"), nb::arg("z"));

m.def("N_VInv",
    N_VInv, nb::arg("x"), nb::arg("z"));

m.def("N_VAddConst",
    N_VAddConst, nb::arg("x"), nb::arg("b"), nb::arg("z"));

m.def("N_VDotProd",
    N_VDotProd, nb::arg("x"), nb::arg("y"));

m.def("N_VMaxNorm",
    N_VMaxNorm, nb::arg("x"));

m.def("N_VWrmsNorm",
    N_VWrmsNorm, nb::arg("x"), nb::arg("w"));

m.def("N_VWrmsNormMask",
    N_VWrmsNormMask, nb::arg("x"), nb::arg("w"), nb::arg("id"));

m.def("N_VMin",
    N_VMin, nb::arg("x"));

m.def("N_VWL2Norm",
    N_VWL2Norm, nb::arg("x"), nb::arg("w"));

m.def("N_VL1Norm",
    N_VL1Norm, nb::arg("x"));

m.def("N_VCompare",
    N_VCompare, nb::arg("c"), nb::arg("x"), nb::arg("z"));

m.def("N_VInvTest",
    N_VInvTest, nb::arg("x"), nb::arg("z"));

m.def("N_VConstrMask",
    N_VConstrMask, nb::arg("c"), nb::arg("x"), nb::arg("m"));

m.def("N_VMinQuotient",
    N_VMinQuotient, nb::arg("num"), nb::arg("denom"));

m.def("N_VLinearCombination",
    [](int nvec, std::vector<sunrealtype> c_arr, std::vector<N_Vector> X_arr, N_Vector z) -> SUNErrCode
    {
        auto N_VLinearCombination_adapt_arr_ptr_to_std_vector = [](int nvec, std::vector<sunrealtype> c_arr, std::vector<N_Vector> X_arr, N_Vector z) -> SUNErrCode
        {
            sunrealtype* c_arr_ptr = reinterpret_cast<sunrealtype*>( c_arr.empty() ? nullptr : c_arr.data() );
            N_Vector* X_arr_ptr = reinterpret_cast<N_Vector*>( X_arr.empty() ? nullptr : X_arr.data() );

            auto lambda_result = N_VLinearCombination(nvec, c_arr_ptr, X_arr_ptr, z);
            return lambda_result;
        };

        return N_VLinearCombination_adapt_arr_ptr_to_std_vector(nvec, c_arr, X_arr, z);
    },     nb::arg("nvec"), nb::arg("c_arr"), nb::arg("X_arr"), nb::arg("z"));

m.def("N_VScaleAddMulti",
    [](int nvec, std::vector<sunrealtype> a, N_Vector x, std::vector<N_Vector> Y_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
    {
        auto N_VScaleAddMulti_adapt_arr_ptr_to_std_vector = [](int nvec, std::vector<sunrealtype> a, N_Vector x, std::vector<N_Vector> Y_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
        {
            sunrealtype* a_ptr = reinterpret_cast<sunrealtype*>( a.empty() ? nullptr : a.data() );
            N_Vector* Y_arr_ptr = reinterpret_cast<N_Vector*>( Y_arr.empty() ? nullptr : Y_arr.data() );
            N_Vector* Z_arr_ptr = reinterpret_cast<N_Vector*>( Z_arr.empty() ? nullptr : Z_arr.data() );

            auto lambda_result = N_VScaleAddMulti(nvec, a_ptr, x, Y_arr_ptr, Z_arr_ptr);
            return lambda_result;
        };

        return N_VScaleAddMulti_adapt_arr_ptr_to_std_vector(nvec, a, x, Y_arr, Z_arr);
    },     nb::arg("nvec"), nb::arg("a"), nb::arg("x"), nb::arg("Y_arr"), nb::arg("Z_arr"));

m.def("N_VDotProdMulti",
    [](int nvec, N_Vector x, std::vector<N_Vector> Y_arr, std::vector<sunrealtype> dotprods) -> SUNErrCode
    {
        auto N_VDotProdMulti_adapt_arr_ptr_to_std_vector = [](int nvec, N_Vector x, std::vector<N_Vector> Y_arr, std::vector<sunrealtype> dotprods) -> SUNErrCode
        {
            N_Vector* Y_arr_ptr = reinterpret_cast<N_Vector*>( Y_arr.empty() ? nullptr : Y_arr.data() );
            sunrealtype* dotprods_ptr = reinterpret_cast<sunrealtype*>( dotprods.empty() ? nullptr : dotprods.data() );

            auto lambda_result = N_VDotProdMulti(nvec, x, Y_arr_ptr, dotprods_ptr);
            return lambda_result;
        };

        return N_VDotProdMulti_adapt_arr_ptr_to_std_vector(nvec, x, Y_arr, dotprods);
    },     nb::arg("nvec"), nb::arg("x"), nb::arg("Y_arr"), nb::arg("dotprods"));

m.def("N_VLinearSumVectorArray",
    [](int nvec, double a, std::vector<N_Vector> X_arr, double b, std::vector<N_Vector> Y_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
    {
        auto N_VLinearSumVectorArray_adapt_arr_ptr_to_std_vector = [](int nvec, double a, std::vector<N_Vector> X_arr, double b, std::vector<N_Vector> Y_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
        {
            N_Vector* X_arr_ptr = reinterpret_cast<N_Vector*>( X_arr.empty() ? nullptr : X_arr.data() );
            N_Vector* Y_arr_ptr = reinterpret_cast<N_Vector*>( Y_arr.empty() ? nullptr : Y_arr.data() );
            N_Vector* Z_arr_ptr = reinterpret_cast<N_Vector*>( Z_arr.empty() ? nullptr : Z_arr.data() );

            auto lambda_result = N_VLinearSumVectorArray(nvec, a, X_arr_ptr, b, Y_arr_ptr, Z_arr_ptr);
            return lambda_result;
        };

        return N_VLinearSumVectorArray_adapt_arr_ptr_to_std_vector(nvec, a, X_arr, b, Y_arr, Z_arr);
    },     nb::arg("nvec"), nb::arg("a"), nb::arg("X_arr"), nb::arg("b"), nb::arg("Y_arr"), nb::arg("Z_arr"));

m.def("N_VScaleVectorArray",
    [](int nvec, std::vector<sunrealtype> c, std::vector<N_Vector> X_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
    {
        auto N_VScaleVectorArray_adapt_arr_ptr_to_std_vector = [](int nvec, std::vector<sunrealtype> c, std::vector<N_Vector> X_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
        {
            sunrealtype* c_ptr = reinterpret_cast<sunrealtype*>( c.empty() ? nullptr : c.data() );
            N_Vector* X_arr_ptr = reinterpret_cast<N_Vector*>( X_arr.empty() ? nullptr : X_arr.data() );
            N_Vector* Z_arr_ptr = reinterpret_cast<N_Vector*>( Z_arr.empty() ? nullptr : Z_arr.data() );

            auto lambda_result = N_VScaleVectorArray(nvec, c_ptr, X_arr_ptr, Z_arr_ptr);
            return lambda_result;
        };

        return N_VScaleVectorArray_adapt_arr_ptr_to_std_vector(nvec, c, X_arr, Z_arr);
    },     nb::arg("nvec"), nb::arg("c"), nb::arg("X_arr"), nb::arg("Z_arr"));

m.def("N_VConstVectorArray",
    [](int nvec, double c, std::vector<N_Vector> Z_arr) -> SUNErrCode
    {
        auto N_VConstVectorArray_adapt_arr_ptr_to_std_vector = [](int nvec, double c, std::vector<N_Vector> Z_arr) -> SUNErrCode
        {
            N_Vector* Z_arr_ptr = reinterpret_cast<N_Vector*>( Z_arr.empty() ? nullptr : Z_arr.data() );

            auto lambda_result = N_VConstVectorArray(nvec, c, Z_arr_ptr);
            return lambda_result;
        };

        return N_VConstVectorArray_adapt_arr_ptr_to_std_vector(nvec, c, Z_arr);
    },     nb::arg("nvec"), nb::arg("c"), nb::arg("Z_arr"));

m.def("N_VWrmsNormVectorArray",
    [](int nvec, std::vector<N_Vector> X_arr, std::vector<N_Vector> W_arr, std::vector<sunrealtype> nrm) -> SUNErrCode
    {
        auto N_VWrmsNormVectorArray_adapt_arr_ptr_to_std_vector = [](int nvec, std::vector<N_Vector> X_arr, std::vector<N_Vector> W_arr, std::vector<sunrealtype> nrm) -> SUNErrCode
        {
            N_Vector* X_arr_ptr = reinterpret_cast<N_Vector*>( X_arr.empty() ? nullptr : X_arr.data() );
            N_Vector* W_arr_ptr = reinterpret_cast<N_Vector*>( W_arr.empty() ? nullptr : W_arr.data() );
            sunrealtype* nrm_ptr = reinterpret_cast<sunrealtype*>( nrm.empty() ? nullptr : nrm.data() );

            auto lambda_result = N_VWrmsNormVectorArray(nvec, X_arr_ptr, W_arr_ptr, nrm_ptr);
            return lambda_result;
        };

        return N_VWrmsNormVectorArray_adapt_arr_ptr_to_std_vector(nvec, X_arr, W_arr, nrm);
    },     nb::arg("nvec"), nb::arg("X_arr"), nb::arg("W_arr"), nb::arg("nrm"));

m.def("N_VWrmsNormMaskVectorArray",
    [](int nvec, std::vector<N_Vector> X_arr, std::vector<N_Vector> W_arr, N_Vector id, std::vector<sunrealtype> nrm) -> SUNErrCode
    {
        auto N_VWrmsNormMaskVectorArray_adapt_arr_ptr_to_std_vector = [](int nvec, std::vector<N_Vector> X_arr, std::vector<N_Vector> W_arr, N_Vector id, std::vector<sunrealtype> nrm) -> SUNErrCode
        {
            N_Vector* X_arr_ptr = reinterpret_cast<N_Vector*>( X_arr.empty() ? nullptr : X_arr.data() );
            N_Vector* W_arr_ptr = reinterpret_cast<N_Vector*>( W_arr.empty() ? nullptr : W_arr.data() );
            sunrealtype* nrm_ptr = reinterpret_cast<sunrealtype*>( nrm.empty() ? nullptr : nrm.data() );

            auto lambda_result = N_VWrmsNormMaskVectorArray(nvec, X_arr_ptr, W_arr_ptr, id, nrm_ptr);
            return lambda_result;
        };

        return N_VWrmsNormMaskVectorArray_adapt_arr_ptr_to_std_vector(nvec, X_arr, W_arr, id, nrm);
    },     nb::arg("nvec"), nb::arg("X_arr"), nb::arg("W_arr"), nb::arg("id"), nb::arg("nrm"));

m.def("N_VScaleAddMultiVectorArray",
    [](int nvec, int nsum, std::vector<sunrealtype> a, std::vector<N_Vector> X_arr, std::vector<N_Vector> Y_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
    {
        auto N_VScaleAddMultiVectorArray_adapt_arr_ptr_to_std_vector = [](int nvec, int nsum, std::vector<sunrealtype> a, std::vector<N_Vector> X_arr, std::vector<N_Vector> Y_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
        {
            sunrealtype* a_ptr = reinterpret_cast<sunrealtype*>( a.empty() ? nullptr : a.data() );
            N_Vector* X_arr_ptr = reinterpret_cast<N_Vector*>( X_arr.empty() ? nullptr : X_arr.data() );
            N_Vector** Y_arr_ptr = reinterpret_cast<N_Vector**>( Y_arr.empty() ? nullptr : Y_arr.data() );
            N_Vector** Z_arr_ptr = reinterpret_cast<N_Vector**>( Z_arr.empty() ? nullptr : Z_arr.data() );

            auto lambda_result = N_VScaleAddMultiVectorArray(nvec, nsum, a_ptr, X_arr_ptr, Y_arr_ptr, Z_arr_ptr);
            return lambda_result;
        };

        return N_VScaleAddMultiVectorArray_adapt_arr_ptr_to_std_vector(nvec, nsum, a, X_arr, Y_arr, Z_arr);
    },     nb::arg("nvec"), nb::arg("nsum"), nb::arg("a"), nb::arg("X_arr"), nb::arg("Y_arr"), nb::arg("Z_arr"));

m.def("N_VLinearCombinationVectorArray",
    [](int nvec, int nsum, std::vector<sunrealtype> c_arr, std::vector<N_Vector> X_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
    {
        auto N_VLinearCombinationVectorArray_adapt_arr_ptr_to_std_vector = [](int nvec, int nsum, std::vector<sunrealtype> c_arr, std::vector<N_Vector> X_arr, std::vector<N_Vector> Z_arr) -> SUNErrCode
        {
            sunrealtype* c_arr_ptr = reinterpret_cast<sunrealtype*>( c_arr.empty() ? nullptr : c_arr.data() );
            N_Vector** X_arr_ptr = reinterpret_cast<N_Vector**>( X_arr.empty() ? nullptr : X_arr.data() );
            N_Vector* Z_arr_ptr = reinterpret_cast<N_Vector*>( Z_arr.empty() ? nullptr : Z_arr.data() );

            auto lambda_result = N_VLinearCombinationVectorArray(nvec, nsum, c_arr_ptr, X_arr_ptr, Z_arr_ptr);
            return lambda_result;
        };

        return N_VLinearCombinationVectorArray_adapt_arr_ptr_to_std_vector(nvec, nsum, c_arr, X_arr, Z_arr);
    },     nb::arg("nvec"), nb::arg("nsum"), nb::arg("c_arr"), nb::arg("X_arr"), nb::arg("Z_arr"));

m.def("N_VDotProdLocal",
    N_VDotProdLocal, nb::arg("x"), nb::arg("y"));

m.def("N_VMaxNormLocal",
    N_VMaxNormLocal, nb::arg("x"));

m.def("N_VMinLocal",
    N_VMinLocal, nb::arg("x"));

m.def("N_VL1NormLocal",
    N_VL1NormLocal, nb::arg("x"));

m.def("N_VWSqrSumLocal",
    N_VWSqrSumLocal, nb::arg("x"), nb::arg("w"));

m.def("N_VWSqrSumMaskLocal",
    N_VWSqrSumMaskLocal, nb::arg("x"), nb::arg("w"), nb::arg("id"));

m.def("N_VInvTestLocal",
    N_VInvTestLocal, nb::arg("x"), nb::arg("z"));

m.def("N_VConstrMaskLocal",
    N_VConstrMaskLocal, nb::arg("c"), nb::arg("x"), nb::arg("m"));

m.def("N_VMinQuotientLocal",
    N_VMinQuotientLocal, nb::arg("num"), nb::arg("denom"));

m.def("N_VDotProdMultiLocal",
    [](int nvec, N_Vector x, std::vector<N_Vector> Y_arr, std::vector<sunrealtype> dotprods) -> SUNErrCode
    {
        auto N_VDotProdMultiLocal_adapt_arr_ptr_to_std_vector = [](int nvec, N_Vector x, std::vector<N_Vector> Y_arr, std::vector<sunrealtype> dotprods) -> SUNErrCode
        {
            N_Vector* Y_arr_ptr = reinterpret_cast<N_Vector*>( Y_arr.empty() ? nullptr : Y_arr.data() );
            sunrealtype* dotprods_ptr = reinterpret_cast<sunrealtype*>( dotprods.empty() ? nullptr : dotprods.data() );

            auto lambda_result = N_VDotProdMultiLocal(nvec, x, Y_arr_ptr, dotprods_ptr);
            return lambda_result;
        };

        return N_VDotProdMultiLocal_adapt_arr_ptr_to_std_vector(nvec, x, Y_arr, dotprods);
    },     nb::arg("nvec"), nb::arg("x"), nb::arg("Y_arr"), nb::arg("dotprods"));

m.def("N_VDotProdMultiAllReduce",
    [](int nvec_total, N_Vector x, std::vector<sunrealtype> sum) -> SUNErrCode
    {
        auto N_VDotProdMultiAllReduce_adapt_arr_ptr_to_std_vector = [](int nvec_total, N_Vector x, std::vector<sunrealtype> sum) -> SUNErrCode
        {
            sunrealtype* sum_ptr = reinterpret_cast<sunrealtype*>( sum.empty() ? nullptr : sum.data() );

            auto lambda_result = N_VDotProdMultiAllReduce(nvec_total, x, sum_ptr);
            return lambda_result;
        };

        return N_VDotProdMultiAllReduce_adapt_arr_ptr_to_std_vector(nvec_total, x, sum);
    },     nb::arg("nvec_total"), nb::arg("x"), nb::arg("sum"));

m.def("N_VPrint",
    N_VPrint, nb::arg("v"));

m.def("N_VPrintFile",
    N_VPrintFile, nb::arg("v"), nb::arg("outfile"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
