// #ifndef _SUNDIALS_ITERATIVE_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyEnumSUN_PREC_ID =
    nb::enum_<SUN_PREC_ID>(m, "SUN_PREC_ID", nb::is_arithmetic(), "")
        .value("SUN_PREC_NONE", SUN_PREC_NONE, "")
        .value("SUN_PREC_LEFT", SUN_PREC_LEFT, "")
        .value("SUN_PREC_RIGHT", SUN_PREC_RIGHT, "")
        .value("SUN_PREC_BOTH", SUN_PREC_BOTH, "")
    .export_values();


auto pyEnumSUN_GRAMSCHMIDT_ID =
    nb::enum_<SUN_GRAMSCHMIDT_ID>(m, "SUN_GRAMSCHMIDT_ID", nb::is_arithmetic(), "")
        .value("SUN_MODIFIED_GS", SUN_MODIFIED_GS, "")
        .value("SUN_CLASSICAL_GS", SUN_CLASSICAL_GS, "")
    .export_values();


m.def("SUNModifiedGS",
    [](std::vector<N_Vector> v, sunrealtype2d h, int k, int p, double new_vk_norm) -> std::tuple<SUNErrCode, double>
    {
        auto SUNModifiedGS_adapt_modifiable_immutable_to_return = [](N_Vector * v, sunrealtype2d h, int k, int p, double new_vk_norm) -> std::tuple<SUNErrCode, double>
        {
            double * new_vk_norm_adapt_modifiable = & new_vk_norm;

            SUNErrCode r = SUNModifiedGS(v, h, k, p, new_vk_norm_adapt_modifiable);
            return std::make_tuple(r, new_vk_norm);
        };
        auto SUNModifiedGS_adapt_nvector_ptr_to_vector = [&SUNModifiedGS_adapt_modifiable_immutable_to_return](std::vector<N_Vector> v, sunrealtype2d h, int k, int p, double new_vk_norm) -> std::tuple<SUNErrCode, double>
        {
            N_Vector* v_ptr = v.empty() ? nullptr : v.data();

            auto lambda_result = SUNModifiedGS_adapt_modifiable_immutable_to_return(v_ptr, h, k, p, new_vk_norm);
            return lambda_result;
        };

        return SUNModifiedGS_adapt_nvector_ptr_to_vector(v, h, k, p, new_vk_norm);
    },     nb::arg("v"), nb::arg("h"), nb::arg("k"), nb::arg("p"), nb::arg("new_vk_norm"));

m.def("SUNQRfact",
    SUNQRfact, nb::arg("n"), nb::arg("h"), nb::arg("q"), nb::arg("job"));

m.def("SUNQRsol",
    SUNQRsol, nb::arg("n"), nb::arg("h"), nb::arg("q"), nb::arg("b"));

m.def("SUNQRAdd_MGS",
    [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
    {
        auto SUNQRAdd_MGS_adapt_nvector_ptr_to_vector = [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
        {
            N_Vector* Q_ptr = Q.empty() ? nullptr : Q.data();

            auto lambda_result = SUNQRAdd_MGS(Q_ptr, R, df, m, mMax, QRdata);
            return lambda_result;
        };

        return SUNQRAdd_MGS_adapt_nvector_ptr_to_vector(Q, R, df, m, mMax, QRdata);
    },     nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_ICWY",
    [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
    {
        auto SUNQRAdd_ICWY_adapt_nvector_ptr_to_vector = [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
        {
            N_Vector* Q_ptr = Q.empty() ? nullptr : Q.data();

            auto lambda_result = SUNQRAdd_ICWY(Q_ptr, R, df, m, mMax, QRdata);
            return lambda_result;
        };

        return SUNQRAdd_ICWY_adapt_nvector_ptr_to_vector(Q, R, df, m, mMax, QRdata);
    },     nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_ICWY_SB",
    [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
    {
        auto SUNQRAdd_ICWY_SB_adapt_nvector_ptr_to_vector = [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
        {
            N_Vector* Q_ptr = Q.empty() ? nullptr : Q.data();

            auto lambda_result = SUNQRAdd_ICWY_SB(Q_ptr, R, df, m, mMax, QRdata);
            return lambda_result;
        };

        return SUNQRAdd_ICWY_SB_adapt_nvector_ptr_to_vector(Q, R, df, m, mMax, QRdata);
    },     nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_CGS2",
    [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
    {
        auto SUNQRAdd_CGS2_adapt_nvector_ptr_to_vector = [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
        {
            N_Vector* Q_ptr = Q.empty() ? nullptr : Q.data();

            auto lambda_result = SUNQRAdd_CGS2(Q_ptr, R, df, m, mMax, QRdata);
            return lambda_result;
        };

        return SUNQRAdd_CGS2_adapt_nvector_ptr_to_vector(Q, R, df, m, mMax, QRdata);
    },     nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_DCGS2",
    [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
    {
        auto SUNQRAdd_DCGS2_adapt_nvector_ptr_to_vector = [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
        {
            N_Vector* Q_ptr = Q.empty() ? nullptr : Q.data();

            auto lambda_result = SUNQRAdd_DCGS2(Q_ptr, R, df, m, mMax, QRdata);
            return lambda_result;
        };

        return SUNQRAdd_DCGS2_adapt_nvector_ptr_to_vector(Q, R, df, m, mMax, QRdata);
    },     nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_DCGS2_SB",
    [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
    {
        auto SUNQRAdd_DCGS2_SB_adapt_nvector_ptr_to_vector = [](std::vector<N_Vector> Q, sunrealtype1d R, N_Vector df, int m, int mMax, void * QRdata) -> SUNErrCode
        {
            N_Vector* Q_ptr = Q.empty() ? nullptr : Q.data();

            auto lambda_result = SUNQRAdd_DCGS2_SB(Q_ptr, R, df, m, mMax, QRdata);
            return lambda_result;
        };

        return SUNQRAdd_DCGS2_SB_adapt_nvector_ptr_to_vector(Q, R, df, m, mMax, QRdata);
    },     nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
// #ifndef _SUNLINEARSOLVER_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyClass_generic_SUNLinearSolver_Ops =
    nb::class_<_generic_SUNLinearSolver_Ops>
        (m, "_generic_SUNLinearSolver_Ops", "")
    .def(nb::init<>()) // implicit default constructor 
    ;


auto pyClass_generic_SUNLinearSolver =
    nb::class_<_generic_SUNLinearSolver>
        (m, "_generic_SUNLinearSolver", "")
    .def("__init__", [](_generic_SUNLinearSolver * self, SUNLinearSolver_Ops ops = SUNLinearSolver_Ops(), SUNContext sunctx = SUNContext())
    {
        new (self) _generic_SUNLinearSolver();  // placement new
        auto r_ctor_ = self;
        r_ctor_->ops = ops;
        r_ctor_->sunctx = sunctx;
    },
    nb::arg("ops") = SUNLinearSolver_Ops(), nb::arg("sunctx") = SUNContext()
    )
    .def_rw("content", &_generic_SUNLinearSolver::content, "")
    .def_rw("ops", &_generic_SUNLinearSolver::ops, "")
    .def_rw("sunctx", &_generic_SUNLinearSolver::sunctx, "")
    ;


m.def("SUNLinSolGetType",
    SUNLinSolGetType, nb::arg("S"));

m.def("SUNLinSolGetID",
    SUNLinSolGetID, nb::arg("S"));

m.def("SUNLinSolSetATimes",
    SUNLinSolSetATimes, nb::arg("S"), nb::arg("A_data"), nb::arg("ATimes"));

m.def("SUNLinSolSetPreconditioner",
    SUNLinSolSetPreconditioner, nb::arg("S"), nb::arg("P_data"), nb::arg("Pset"), nb::arg("Psol"));

m.def("SUNLinSolSetScalingVectors",
    SUNLinSolSetScalingVectors, nb::arg("S"), nb::arg("s1"), nb::arg("s2"));

m.def("SUNLinSolSetZeroGuess",
    SUNLinSolSetZeroGuess, nb::arg("S"), nb::arg("onoff"));

m.def("SUNLinSolInitialize",
    SUNLinSolInitialize, nb::arg("S"));

m.def("SUNLinSolSetup",
    SUNLinSolSetup, nb::arg("S"), nb::arg("A"));

m.def("SUNLinSolSolve",
    SUNLinSolSolve, nb::arg("S"), nb::arg("A"), nb::arg("x"), nb::arg("b"), nb::arg("tol"));

m.def("SUNLinSolNumIters",
    SUNLinSolNumIters, nb::arg("S"));

m.def("SUNLinSolResNorm",
    SUNLinSolResNorm, nb::arg("S"));

m.def("SUNLinSolResid",
    SUNLinSolResid, nb::arg("S"));

m.def("SUNLinSolLastFlag",
    SUNLinSolLastFlag, nb::arg("S"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
