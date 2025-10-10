// #ifndef _SUNDIALS_ITERATIVE_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNPrecType = nb::enum_<SUNPrecType>(m, "SUNPrecType",
                                                nb::is_arithmetic(), "")
                           .value("SUN_PREC_NONE", SUN_PREC_NONE, "")
                           .value("SUN_PREC_LEFT", SUN_PREC_LEFT, "")
                           .value("SUN_PREC_RIGHT", SUN_PREC_RIGHT, "")
                           .value("SUN_PREC_BOTH", SUN_PREC_BOTH, "")
                           .export_values();

auto pyEnumSUNGramSchmidtType =
  nb::enum_<SUNGramSchmidtType>(m, "SUNGramSchmidtType", nb::is_arithmetic(), "")
    .value("SUN_MODIFIED_GS", SUN_MODIFIED_GS, "")
    .value("SUN_CLASSICAL_GS", SUN_CLASSICAL_GS, "")
    .export_values();

m.def(
  "SUNQRAdd_MGS",
  [](std::vector<N_Vector> Q_1d,
     nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
     N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_MGS_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d,
         nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
         N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr =
        reinterpret_cast<N_Vector*>(Q_1d.empty() ? nullptr : Q_1d.data());
      sunrealtype* R_1d_ptr = reinterpret_cast<sunrealtype*>(R_1d.data());

      auto lambda_result = SUNQRAdd_MGS(Q_1d_ptr, R_1d_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_MGS_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                    QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_ICWY",
  [](std::vector<N_Vector> Q_1d,
     nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
     N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_ICWY_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d,
         nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
         N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr =
        reinterpret_cast<N_Vector*>(Q_1d.empty() ? nullptr : Q_1d.data());
      sunrealtype* R_1d_ptr = reinterpret_cast<sunrealtype*>(R_1d.data());

      auto lambda_result = SUNQRAdd_ICWY(Q_1d_ptr, R_1d_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_ICWY_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                     QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_ICWY_SB",
  [](std::vector<N_Vector> Q_1d,
     nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
     N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_ICWY_SB_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d,
         nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
         N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr =
        reinterpret_cast<N_Vector*>(Q_1d.empty() ? nullptr : Q_1d.data());
      sunrealtype* R_1d_ptr = reinterpret_cast<sunrealtype*>(R_1d.data());

      auto lambda_result = SUNQRAdd_ICWY_SB(Q_1d_ptr, R_1d_ptr, df, m, mMax,
                                            QRdata);
      return lambda_result;
    };

    return SUNQRAdd_ICWY_SB_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                        QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_CGS2",
  [](std::vector<N_Vector> Q_1d,
     nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
     N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_CGS2_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d,
         nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
         N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr =
        reinterpret_cast<N_Vector*>(Q_1d.empty() ? nullptr : Q_1d.data());
      sunrealtype* R_1d_ptr = reinterpret_cast<sunrealtype*>(R_1d.data());

      auto lambda_result = SUNQRAdd_CGS2(Q_1d_ptr, R_1d_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_CGS2_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                     QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_DCGS2",
  [](std::vector<N_Vector> Q_1d,
     nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
     N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_DCGS2_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d,
         nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
         N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr =
        reinterpret_cast<N_Vector*>(Q_1d.empty() ? nullptr : Q_1d.data());
      sunrealtype* R_1d_ptr = reinterpret_cast<sunrealtype*>(R_1d.data());

      auto lambda_result = SUNQRAdd_DCGS2(Q_1d_ptr, R_1d_ptr, df, m, mMax,
                                          QRdata);
      return lambda_result;
    };

    return SUNQRAdd_DCGS2_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                      QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_DCGS2_SB",
  [](std::vector<N_Vector> Q_1d,
     nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
     N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_DCGS2_SB_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d,
         nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> R_1d,
         N_Vector df, int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr =
        reinterpret_cast<N_Vector*>(Q_1d.empty() ? nullptr : Q_1d.data());
      sunrealtype* R_1d_ptr = reinterpret_cast<sunrealtype*>(R_1d.data());

      auto lambda_result = SUNQRAdd_DCGS2_SB(Q_1d_ptr, R_1d_ptr, df, m, mMax,
                                             QRdata);
      return lambda_result;
    };

    return SUNQRAdd_DCGS2_SB_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m,
                                                         mMax, QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));
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
    .value("SUNLINEARSOLVER_GINKGOBATCH", SUNLINEARSOLVER_GINKGOBATCH, "")
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

m.def("SUNLinSolResid", SUNLinSolResid, nb::arg("S"), nb::rv_policy::reference);

m.def("SUNLinSolLastFlag", SUNLinSolLastFlag, nb::arg("S"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
