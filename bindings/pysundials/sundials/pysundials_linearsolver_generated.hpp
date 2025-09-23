// #ifndef _SUNDIALS_ITERATIVE_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUN_PREC_ID = nb::enum_<SUN_PREC_ID>(m, "SUN_PREC_ID",
                                                nb::is_arithmetic(), "")
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

m.def(
  "SUNModifiedGS",
  [](std::vector<N_Vector> v, std::vector<sunrealtype> h, int k, int p,
     double new_vk_norm) -> std::tuple<SUNErrCode, double>
  {
    auto SUNModifiedGS_adapt_modifiable_immutable_to_return =
      [](N_Vector1d v, sunrealtype2d h, int k, int p,
         double new_vk_norm) -> std::tuple<SUNErrCode, double>
    {
      double* new_vk_norm_adapt_modifiable = &new_vk_norm;

      SUNErrCode r = SUNModifiedGS(v, h, k, p, new_vk_norm_adapt_modifiable);
      return std::make_tuple(r, new_vk_norm);
    };
    auto SUNModifiedGS_adapt_arr_ptr_to_std_vector =
      [&SUNModifiedGS_adapt_modifiable_immutable_to_return](std::vector<N_Vector> v,
                                                            std::vector<sunrealtype> h,
                                                            int k, int p,
                                                            double new_vk_norm)
      -> std::tuple<SUNErrCode, double>
    {
      N_Vector* v_ptr = reinterpret_cast<N_Vector*>(v.empty() ? nullptr
                                                              : v.data());
      sunrealtype** h_ptr =
        reinterpret_cast<sunrealtype**>(h.empty() ? nullptr : h.data());

      auto lambda_result =
        SUNModifiedGS_adapt_modifiable_immutable_to_return(v_ptr, h_ptr, k, p,
                                                           new_vk_norm);
      return lambda_result;
    };

    return SUNModifiedGS_adapt_arr_ptr_to_std_vector(v, h, k, p, new_vk_norm);
  },
  nb::arg("v"), nb::arg("h"), nb::arg("k"), nb::arg("p"), nb::arg("new_vk_norm"));

m.def(
  "SUNClassicalGS",
  [](std::vector<N_Vector> v, std::vector<sunrealtype> h, int k, int p,
     double new_vk_norm, std::vector<sunrealtype> stemp,
     std::vector<N_Vector> vtemp) -> std::tuple<SUNErrCode, double>
  {
    auto SUNClassicalGS_adapt_modifiable_immutable_to_return =
      [](N_Vector1d v, sunrealtype2d h, int k, int p, double new_vk_norm,
         sunrealtype1d stemp, N_Vector1d vtemp) -> std::tuple<SUNErrCode, double>
    {
      double* new_vk_norm_adapt_modifiable = &new_vk_norm;

      SUNErrCode r = SUNClassicalGS(v, h, k, p, new_vk_norm_adapt_modifiable,
                                    stemp, vtemp);
      return std::make_tuple(r, new_vk_norm);
    };
    auto SUNClassicalGS_adapt_arr_ptr_to_std_vector =
      [&SUNClassicalGS_adapt_modifiable_immutable_to_return](std::vector<N_Vector> v,
                                                             std::vector<sunrealtype> h,
                                                             int k, int p,
                                                             double new_vk_norm,
                                                             std::vector<sunrealtype> stemp,
                                                             std::vector<N_Vector> vtemp)
      -> std::tuple<SUNErrCode, double>
    {
      N_Vector* v_ptr = reinterpret_cast<N_Vector*>(v.empty() ? nullptr
                                                              : v.data());
      sunrealtype** h_ptr =
        reinterpret_cast<sunrealtype**>(h.empty() ? nullptr : h.data());
      sunrealtype* stemp_ptr =
        reinterpret_cast<sunrealtype*>(stemp.empty() ? nullptr : stemp.data());
      N_Vector* vtemp_ptr =
        reinterpret_cast<N_Vector*>(vtemp.empty() ? nullptr : vtemp.data());

      auto lambda_result =
        SUNClassicalGS_adapt_modifiable_immutable_to_return(v_ptr, h_ptr, k, p,
                                                            new_vk_norm,
                                                            stemp_ptr, vtemp_ptr);
      return lambda_result;
    };

    return SUNClassicalGS_adapt_arr_ptr_to_std_vector(v, h, k, p, new_vk_norm,
                                                      stemp, vtemp);
  },
  nb::arg("v"), nb::arg("h"), nb::arg("k"), nb::arg("p"),
  nb::arg("new_vk_norm"), nb::arg("stemp"), nb::arg("vtemp"));

m.def(
  "SUNQRfact",
  [](int n, std::vector<sunrealtype> h, std::vector<sunrealtype> q, int job) -> int
  {
    auto SUNQRfact_adapt_arr_ptr_to_std_vector =
      [](int n, std::vector<sunrealtype> h, std::vector<sunrealtype> q,
         int job) -> int
    {
      sunrealtype** h_ptr =
        reinterpret_cast<sunrealtype**>(h.empty() ? nullptr : h.data());
      sunrealtype* q_ptr = reinterpret_cast<sunrealtype*>(q.empty() ? nullptr
                                                                    : q.data());

      auto lambda_result = SUNQRfact(n, h_ptr, q_ptr, job);
      return lambda_result;
    };

    return SUNQRfact_adapt_arr_ptr_to_std_vector(n, h, q, job);
  },
  nb::arg("n"), nb::arg("h"), nb::arg("q"), nb::arg("job"));

m.def(
  "SUNQRsol",
  [](int n, std::vector<sunrealtype> h, std::vector<sunrealtype> q,
     std::vector<sunrealtype> b) -> int
  {
    auto SUNQRsol_adapt_arr_ptr_to_std_vector =
      [](int n, std::vector<sunrealtype> h, std::vector<sunrealtype> q,
         std::vector<sunrealtype> b) -> int
    {
      sunrealtype** h_ptr =
        reinterpret_cast<sunrealtype**>(h.empty() ? nullptr : h.data());
      sunrealtype* q_ptr = reinterpret_cast<sunrealtype*>(q.empty() ? nullptr
                                                                    : q.data());
      sunrealtype* b_ptr = reinterpret_cast<sunrealtype*>(b.empty() ? nullptr
                                                                    : b.data());

      auto lambda_result = SUNQRsol(n, h_ptr, q_ptr, b_ptr);
      return lambda_result;
    };

    return SUNQRsol_adapt_arr_ptr_to_std_vector(n, h, q, b);
  },
  nb::arg("n"), nb::arg("h"), nb::arg("q"), nb::arg("b"));

m.def(
  "SUNQRAdd_MGS",
  [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_MGS_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_ptr    = reinterpret_cast<N_Vector*>(Q.empty() ? nullptr
                                                                 : Q.data());
      sunrealtype* R_ptr = reinterpret_cast<sunrealtype*>(R.empty() ? nullptr
                                                                    : R.data());

      auto lambda_result = SUNQRAdd_MGS(Q_ptr, R_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_MGS_adapt_arr_ptr_to_std_vector(Q, R, df, m, mMax, QRdata);
  },
  nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"),
  nb::arg("QRdata"));

m.def(
  "SUNQRAdd_ICWY",
  [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_ICWY_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_ptr    = reinterpret_cast<N_Vector*>(Q.empty() ? nullptr
                                                                 : Q.data());
      sunrealtype* R_ptr = reinterpret_cast<sunrealtype*>(R.empty() ? nullptr
                                                                    : R.data());

      auto lambda_result = SUNQRAdd_ICWY(Q_ptr, R_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_ICWY_adapt_arr_ptr_to_std_vector(Q, R, df, m, mMax, QRdata);
  },
  nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"),
  nb::arg("QRdata"));

m.def(
  "SUNQRAdd_ICWY_SB",
  [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_ICWY_SB_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_ptr    = reinterpret_cast<N_Vector*>(Q.empty() ? nullptr
                                                                 : Q.data());
      sunrealtype* R_ptr = reinterpret_cast<sunrealtype*>(R.empty() ? nullptr
                                                                    : R.data());

      auto lambda_result = SUNQRAdd_ICWY_SB(Q_ptr, R_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_ICWY_SB_adapt_arr_ptr_to_std_vector(Q, R, df, m, mMax,
                                                        QRdata);
  },
  nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"),
  nb::arg("QRdata"));

m.def(
  "SUNQRAdd_CGS2",
  [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_CGS2_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_ptr    = reinterpret_cast<N_Vector*>(Q.empty() ? nullptr
                                                                 : Q.data());
      sunrealtype* R_ptr = reinterpret_cast<sunrealtype*>(R.empty() ? nullptr
                                                                    : R.data());

      auto lambda_result = SUNQRAdd_CGS2(Q_ptr, R_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_CGS2_adapt_arr_ptr_to_std_vector(Q, R, df, m, mMax, QRdata);
  },
  nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"),
  nb::arg("QRdata"));

m.def(
  "SUNQRAdd_DCGS2",
  [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_DCGS2_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_ptr    = reinterpret_cast<N_Vector*>(Q.empty() ? nullptr
                                                                 : Q.data());
      sunrealtype* R_ptr = reinterpret_cast<sunrealtype*>(R.empty() ? nullptr
                                                                    : R.data());

      auto lambda_result = SUNQRAdd_DCGS2(Q_ptr, R_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_DCGS2_adapt_arr_ptr_to_std_vector(Q, R, df, m, mMax, QRdata);
  },
  nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"),
  nb::arg("QRdata"));

m.def(
  "SUNQRAdd_DCGS2_SB",
  [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_DCGS2_SB_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q, std::vector<sunrealtype> R, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_ptr    = reinterpret_cast<N_Vector*>(Q.empty() ? nullptr
                                                                 : Q.data());
      sunrealtype* R_ptr = reinterpret_cast<sunrealtype*>(R.empty() ? nullptr
                                                                    : R.data());

      auto lambda_result = SUNQRAdd_DCGS2_SB(Q_ptr, R_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_DCGS2_SB_adapt_arr_ptr_to_std_vector(Q, R, df, m, mMax,
                                                         QRdata);
  },
  nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"),
  nb::arg("QRdata"));
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

auto pyEnumSUNLinearSolver_Type =
  nb::enum_<SUNLinearSolver_Type>(m, "SUNLinearSolver_Type",
                                  nb::is_arithmetic(), "")
    .value("SUNLINEARSOLVER_DIRECT", SUNLINEARSOLVER_DIRECT, "")
    .value("SUNLINEARSOLVER_ITERATIVE", SUNLINEARSOLVER_ITERATIVE, "")
    .value("SUNLINEARSOLVER_MATRIX_ITERATIVE", SUNLINEARSOLVER_MATRIX_ITERATIVE,
           "")
    .value("SUNLINEARSOLVER_MATRIX_EMBEDDED", SUNLINEARSOLVER_MATRIX_EMBEDDED, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyEnumSUNLinearSolver_ID =
  nb::enum_<SUNLinearSolver_ID>(m, "SUNLinearSolver_ID", nb::is_arithmetic(), "")
    .value("SUNLINEARSOLVER_BAND", SUNLINEARSOLVER_BAND, "")
    .value("SUNLINEARSOLVER_DENSE", SUNLINEARSOLVER_DENSE, "")
    .value("SUNLINEARSOLVER_KLU", SUNLINEARSOLVER_KLU, "")
    .value("SUNLINEARSOLVER_LAPACKBAND", SUNLINEARSOLVER_LAPACKBAND, "")
    .value("SUNLINEARSOLVER_LAPACKDENSE", SUNLINEARSOLVER_LAPACKDENSE, "")
    .value("SUNLINEARSOLVER_PCG", SUNLINEARSOLVER_PCG, "")
    .value("SUNLINEARSOLVER_SPBCGS", SUNLINEARSOLVER_SPBCGS, "")
    .value("SUNLINEARSOLVER_SPFGMR", SUNLINEARSOLVER_SPFGMR, "")
    .value("SUNLINEARSOLVER_SPGMR", SUNLINEARSOLVER_SPGMR, "")
    .value("SUNLINEARSOLVER_SPTFQMR", SUNLINEARSOLVER_SPTFQMR, "")
    .value("SUNLINEARSOLVER_SUPERLUDIST", SUNLINEARSOLVER_SUPERLUDIST, "")
    .value("SUNLINEARSOLVER_SUPERLUMT", SUNLINEARSOLVER_SUPERLUMT, "")
    .value("SUNLINEARSOLVER_CUSOLVERSP_BATCHQR",
           SUNLINEARSOLVER_CUSOLVERSP_BATCHQR, "")
    .value("SUNLINEARSOLVER_MAGMADENSE", SUNLINEARSOLVER_MAGMADENSE, "")
    .value("SUNLINEARSOLVER_ONEMKLDENSE", SUNLINEARSOLVER_ONEMKLDENSE, "")
    .value("SUNLINEARSOLVER_GINKGO", SUNLINEARSOLVER_GINKGO, "")
    .value("SUNLINEARSOLVER_KOKKOSDENSE", SUNLINEARSOLVER_KOKKOSDENSE, "")
    .value("SUNLINEARSOLVER_CUSTOM", SUNLINEARSOLVER_CUSTOM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClass_generic_SUNLinearSolver_Ops =
  nb::class_<_generic_SUNLinearSolver_Ops>(m, "_generic_SUNLinearSolver_Ops", "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_SUNLinearSolver =
  nb::class_<_generic_SUNLinearSolver>(m, "_generic_SUNLinearSolver", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNLinSolGetType", SUNLinSolGetType, nb::arg("S"));

m.def("SUNLinSolGetID", SUNLinSolGetID, nb::arg("S"));

m.def("SUNLinSolSetATimes", SUNLinSolSetATimes, nb::arg("S"), nb::arg("A_data"),
      nb::arg("ATimes"));

m.def("SUNLinSolSetPreconditioner", SUNLinSolSetPreconditioner, nb::arg("S"),
      nb::arg("P_data"), nb::arg("Pset"), nb::arg("Psol"));

m.def("SUNLinSolSetScalingVectors", SUNLinSolSetScalingVectors, nb::arg("S"),
      nb::arg("s1"), nb::arg("s2"));

m.def("SUNLinSolSetZeroGuess", SUNLinSolSetZeroGuess, nb::arg("S"),
      nb::arg("onoff"));

m.def("SUNLinSolInitialize", SUNLinSolInitialize, nb::arg("S"));

m.def("SUNLinSolSetup", SUNLinSolSetup, nb::arg("S"), nb::arg("A"));

m.def("SUNLinSolSolve", SUNLinSolSolve, nb::arg("S"), nb::arg("A"),
      nb::arg("x"), nb::arg("b"), nb::arg("tol"));

m.def("SUNLinSolNumIters", SUNLinSolNumIters, nb::arg("S"));

m.def("SUNLinSolResNorm", SUNLinSolResNorm, nb::arg("S"));

m.def("SUNLinSolResid", SUNLinSolResid, nb::arg("S"));

m.def("SUNLinSolLastFlag", SUNLinSolLastFlag, nb::arg("S"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
