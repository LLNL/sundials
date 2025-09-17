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
