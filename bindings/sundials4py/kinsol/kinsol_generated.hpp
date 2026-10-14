// #ifndef _KINSOL_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("KIN_SUCCESS")             = 0;
m.attr("KIN_INITIAL_GUESS_OK")    = 1;
m.attr("KIN_STEP_LT_STPTOL")      = 2;
m.attr("KIN_WARNING")             = 99;
m.attr("KIN_MEM_NULL")            = -1;
m.attr("KIN_ILL_INPUT")           = -2;
m.attr("KIN_NO_MALLOC")           = -3;
m.attr("KIN_MEM_FAIL")            = -4;
m.attr("KIN_LINESEARCH_NONCONV")  = -5;
m.attr("KIN_MAXITER_REACHED")     = -6;
m.attr("KIN_MXNEWT_5X_EXCEEDED")  = -7;
m.attr("KIN_LINESEARCH_BCFAIL")   = -8;
m.attr("KIN_LINSOLV_NO_RECOVERY") = -9;
m.attr("KIN_LINIT_FAIL")          = -10;
m.attr("KIN_LSETUP_FAIL")         = -11;
m.attr("KIN_LSOLVE_FAIL")         = -12;
m.attr("KIN_SYSFUNC_FAIL")        = -13;
m.attr("KIN_FIRST_SYSFUNC_ERR")   = -14;
m.attr("KIN_REPTD_SYSFUNC_ERR")   = -15;
m.attr("KIN_VECTOROP_ERR")        = -16;
m.attr("KIN_CONTEXT_ERR")         = -17;
m.attr("KIN_DAMPING_FN_ERR")      = -18;
m.attr("KIN_DEPTH_FN_ERR")        = -19;
m.attr("KIN_ORTH_MGS")            = 0;
m.attr("KIN_ORTH_ICWY")           = 1;
m.attr("KIN_ORTH_CGS2")           = 2;
m.attr("KIN_ORTH_DCGS2")          = 3;
m.attr("KIN_ETACHOICE1")          = 1;
m.attr("KIN_ETACHOICE2")          = 2;
m.attr("KIN_ETACONSTANT")         = 3;
m.attr("KIN_NONE")                = 0;
m.attr("KIN_LINESEARCH")          = 1;
m.attr("KIN_PICARD")              = 2;
m.attr("KIN_FP")                  = 3;

sundials4py::scoped_def(
  m, "KINSetUserData",
  [](const KINView& kinmem, void* user_data) -> int
  {
    auto KINSetUserData_adapt_solver_memory = [](const KINView& kinmem,
                                                 void* user_data) -> int
    {
      auto lambda_result = KINSetUserData(kinmem.get(), user_data);
      return lambda_result;
    };

    return KINSetUserData_adapt_solver_memory(kinmem, user_data);
  },
  nb::arg("kinmem"), nb::arg("user_data"));

sundials4py::scoped_def(
  m, "KINSetDamping",
  [](const KINView& kinmem, sunrealtype beta) -> int
  {
    auto KINSetDamping_adapt_solver_memory = [](const KINView& kinmem,
                                                sunrealtype beta) -> int
    {
      auto lambda_result = KINSetDamping(kinmem.get(), beta);
      return lambda_result;
    };

    return KINSetDamping_adapt_solver_memory(kinmem, beta);
  },
  nb::arg("kinmem"), nb::arg("beta"));

sundials4py::scoped_def(
  m, "KINSetMAA",
  [](const KINView& kinmem, long maa) -> int
  {
    auto KINSetMAA_adapt_solver_memory = [](const KINView& kinmem, long maa) -> int
    {
      auto lambda_result = KINSetMAA(kinmem.get(), maa);
      return lambda_result;
    };

    return KINSetMAA_adapt_solver_memory(kinmem, maa);
  },
  nb::arg("kinmem"), nb::arg("maa"));

sundials4py::scoped_def(
  m, "KINSetOrthAA",
  [](const KINView& kinmem, int orthaa) -> int
  {
    auto KINSetOrthAA_adapt_solver_memory = [](const KINView& kinmem,
                                               int orthaa) -> int
    {
      auto lambda_result = KINSetOrthAA(kinmem.get(), orthaa);
      return lambda_result;
    };

    return KINSetOrthAA_adapt_solver_memory(kinmem, orthaa);
  },
  nb::arg("kinmem"), nb::arg("orthaa"));

sundials4py::scoped_def(
  m, "KINSetDelayAA",
  [](const KINView& kinmem, long delay) -> int
  {
    auto KINSetDelayAA_adapt_solver_memory = [](const KINView& kinmem,
                                                long delay) -> int
    {
      auto lambda_result = KINSetDelayAA(kinmem.get(), delay);
      return lambda_result;
    };

    return KINSetDelayAA_adapt_solver_memory(kinmem, delay);
  },
  nb::arg("kinmem"), nb::arg("delay"));

sundials4py::scoped_def(
  m, "KINSetDampingAA",
  [](const KINView& kinmem, sunrealtype beta) -> int
  {
    auto KINSetDampingAA_adapt_solver_memory = [](const KINView& kinmem,
                                                  sunrealtype beta) -> int
    {
      auto lambda_result = KINSetDampingAA(kinmem.get(), beta);
      return lambda_result;
    };

    return KINSetDampingAA_adapt_solver_memory(kinmem, beta);
  },
  nb::arg("kinmem"), nb::arg("beta"));

sundials4py::scoped_def(
  m, "KINSetReturnNewest",
  [](const KINView& kinmem, sunbooleantype ret_newest) -> int
  {
    auto KINSetReturnNewest_adapt_solver_memory =
      [](const KINView& kinmem, sunbooleantype ret_newest) -> int
    {
      auto lambda_result = KINSetReturnNewest(kinmem.get(), ret_newest);
      return lambda_result;
    };

    return KINSetReturnNewest_adapt_solver_memory(kinmem, ret_newest);
  },
  nb::arg("kinmem"), nb::arg("ret_newest"));

sundials4py::scoped_def(
  m, "KINSetNumMaxIters",
  [](const KINView& kinmem, long mxiter) -> int
  {
    auto KINSetNumMaxIters_adapt_solver_memory = [](const KINView& kinmem,
                                                    long mxiter) -> int
    {
      auto lambda_result = KINSetNumMaxIters(kinmem.get(), mxiter);
      return lambda_result;
    };

    return KINSetNumMaxIters_adapt_solver_memory(kinmem, mxiter);
  },
  nb::arg("kinmem"), nb::arg("mxiter"));

sundials4py::scoped_def(
  m, "KINSetNoInitSetup",
  [](const KINView& kinmem, sunbooleantype noInitSetup) -> int
  {
    auto KINSetNoInitSetup_adapt_solver_memory =
      [](const KINView& kinmem, sunbooleantype noInitSetup) -> int
    {
      auto lambda_result = KINSetNoInitSetup(kinmem.get(), noInitSetup);
      return lambda_result;
    };

    return KINSetNoInitSetup_adapt_solver_memory(kinmem, noInitSetup);
  },
  nb::arg("kinmem"), nb::arg("noInitSetup"));

sundials4py::scoped_def(
  m, "KINSetNoResMon",
  [](const KINView& kinmem, sunbooleantype noNNIResMon) -> int
  {
    auto KINSetNoResMon_adapt_solver_memory =
      [](const KINView& kinmem, sunbooleantype noNNIResMon) -> int
    {
      auto lambda_result = KINSetNoResMon(kinmem.get(), noNNIResMon);
      return lambda_result;
    };

    return KINSetNoResMon_adapt_solver_memory(kinmem, noNNIResMon);
  },
  nb::arg("kinmem"), nb::arg("noNNIResMon"));

sundials4py::scoped_def(
  m, "KINSetMaxSetupCalls",
  [](const KINView& kinmem, long msbset) -> int
  {
    auto KINSetMaxSetupCalls_adapt_solver_memory = [](const KINView& kinmem,
                                                      long msbset) -> int
    {
      auto lambda_result = KINSetMaxSetupCalls(kinmem.get(), msbset);
      return lambda_result;
    };

    return KINSetMaxSetupCalls_adapt_solver_memory(kinmem, msbset);
  },
  nb::arg("kinmem"), nb::arg("msbset"));

sundials4py::scoped_def(
  m, "KINSetMaxSubSetupCalls",
  [](const KINView& kinmem, long msbsetsub) -> int
  {
    auto KINSetMaxSubSetupCalls_adapt_solver_memory = [](const KINView& kinmem,
                                                         long msbsetsub) -> int
    {
      auto lambda_result = KINSetMaxSubSetupCalls(kinmem.get(), msbsetsub);
      return lambda_result;
    };

    return KINSetMaxSubSetupCalls_adapt_solver_memory(kinmem, msbsetsub);
  },
  nb::arg("kinmem"), nb::arg("msbsetsub"));

sundials4py::scoped_def(
  m, "KINSetEtaForm",
  [](const KINView& kinmem, int etachoice) -> int
  {
    auto KINSetEtaForm_adapt_solver_memory = [](const KINView& kinmem,
                                                int etachoice) -> int
    {
      auto lambda_result = KINSetEtaForm(kinmem.get(), etachoice);
      return lambda_result;
    };

    return KINSetEtaForm_adapt_solver_memory(kinmem, etachoice);
  },
  nb::arg("kinmem"), nb::arg("etachoice"));

sundials4py::scoped_def(
  m, "KINSetEtaConstValue",
  [](const KINView& kinmem, sunrealtype eta) -> int
  {
    auto KINSetEtaConstValue_adapt_solver_memory = [](const KINView& kinmem,
                                                      sunrealtype eta) -> int
    {
      auto lambda_result = KINSetEtaConstValue(kinmem.get(), eta);
      return lambda_result;
    };

    return KINSetEtaConstValue_adapt_solver_memory(kinmem, eta);
  },
  nb::arg("kinmem"), nb::arg("eta"));

sundials4py::scoped_def(
  m, "KINSetEtaParams",
  [](const KINView& kinmem, sunrealtype egamma, sunrealtype ealpha) -> int
  {
    auto KINSetEtaParams_adapt_solver_memory =
      [](const KINView& kinmem, sunrealtype egamma, sunrealtype ealpha) -> int
    {
      auto lambda_result = KINSetEtaParams(kinmem.get(), egamma, ealpha);
      return lambda_result;
    };

    return KINSetEtaParams_adapt_solver_memory(kinmem, egamma, ealpha);
  },
  nb::arg("kinmem"), nb::arg("egamma"), nb::arg("ealpha"));

sundials4py::scoped_def(
  m, "KINSetResMonParams",
  [](const KINView& kinmem, sunrealtype omegamin, sunrealtype omegamax) -> int
  {
    auto KINSetResMonParams_adapt_solver_memory =
      [](const KINView& kinmem, sunrealtype omegamin, sunrealtype omegamax) -> int
    {
      auto lambda_result = KINSetResMonParams(kinmem.get(), omegamin, omegamax);
      return lambda_result;
    };

    return KINSetResMonParams_adapt_solver_memory(kinmem, omegamin, omegamax);
  },
  nb::arg("kinmem"), nb::arg("omegamin"), nb::arg("omegamax"));

sundials4py::scoped_def(
  m, "KINSetResMonConstValue",
  [](const KINView& kinmem, sunrealtype omegaconst) -> int
  {
    auto KINSetResMonConstValue_adapt_solver_memory =
      [](const KINView& kinmem, sunrealtype omegaconst) -> int
    {
      auto lambda_result = KINSetResMonConstValue(kinmem.get(), omegaconst);
      return lambda_result;
    };

    return KINSetResMonConstValue_adapt_solver_memory(kinmem, omegaconst);
  },
  nb::arg("kinmem"), nb::arg("omegaconst"));

sundials4py::scoped_def(
  m, "KINSetNoMinEps",
  [](const KINView& kinmem, sunbooleantype noMinEps) -> int
  {
    auto KINSetNoMinEps_adapt_solver_memory = [](const KINView& kinmem,
                                                 sunbooleantype noMinEps) -> int
    {
      auto lambda_result = KINSetNoMinEps(kinmem.get(), noMinEps);
      return lambda_result;
    };

    return KINSetNoMinEps_adapt_solver_memory(kinmem, noMinEps);
  },
  nb::arg("kinmem"), nb::arg("noMinEps"));

sundials4py::scoped_def(
  m, "KINSetMaxNewtonStep",
  [](const KINView& kinmem, sunrealtype mxnewtstep) -> int
  {
    auto KINSetMaxNewtonStep_adapt_solver_memory =
      [](const KINView& kinmem, sunrealtype mxnewtstep) -> int
    {
      auto lambda_result = KINSetMaxNewtonStep(kinmem.get(), mxnewtstep);
      return lambda_result;
    };

    return KINSetMaxNewtonStep_adapt_solver_memory(kinmem, mxnewtstep);
  },
  nb::arg("kinmem"), nb::arg("mxnewtstep"));

sundials4py::scoped_def(
  m, "KINSetMaxBetaFails",
  [](const KINView& kinmem, long mxnbcf) -> int
  {
    auto KINSetMaxBetaFails_adapt_solver_memory = [](const KINView& kinmem,
                                                     long mxnbcf) -> int
    {
      auto lambda_result = KINSetMaxBetaFails(kinmem.get(), mxnbcf);
      return lambda_result;
    };

    return KINSetMaxBetaFails_adapt_solver_memory(kinmem, mxnbcf);
  },
  nb::arg("kinmem"), nb::arg("mxnbcf"));

sundials4py::scoped_def(
  m, "KINSetRelErrFunc",
  [](const KINView& kinmem, sunrealtype relfunc) -> int
  {
    auto KINSetRelErrFunc_adapt_solver_memory = [](const KINView& kinmem,
                                                   sunrealtype relfunc) -> int
    {
      auto lambda_result = KINSetRelErrFunc(kinmem.get(), relfunc);
      return lambda_result;
    };

    return KINSetRelErrFunc_adapt_solver_memory(kinmem, relfunc);
  },
  nb::arg("kinmem"), nb::arg("relfunc"));

sundials4py::scoped_def(
  m, "KINSetFuncNormTol",
  [](const KINView& kinmem, sunrealtype fnormtol) -> int
  {
    auto KINSetFuncNormTol_adapt_solver_memory = [](const KINView& kinmem,
                                                    sunrealtype fnormtol) -> int
    {
      auto lambda_result = KINSetFuncNormTol(kinmem.get(), fnormtol);
      return lambda_result;
    };

    return KINSetFuncNormTol_adapt_solver_memory(kinmem, fnormtol);
  },
  nb::arg("kinmem"), nb::arg("fnormtol"));

sundials4py::scoped_def(
  m, "KINSetScaledStepTol",
  [](const KINView& kinmem, sunrealtype scsteptol) -> int
  {
    auto KINSetScaledStepTol_adapt_solver_memory =
      [](const KINView& kinmem, sunrealtype scsteptol) -> int
    {
      auto lambda_result = KINSetScaledStepTol(kinmem.get(), scsteptol);
      return lambda_result;
    };

    return KINSetScaledStepTol_adapt_solver_memory(kinmem, scsteptol);
  },
  nb::arg("kinmem"), nb::arg("scsteptol"));

sundials4py::scoped_def(
  m, "KINSetConstraints",
  [](const KINView& kinmem, N_Vector constraints) -> int
  {
    auto KINSetConstraints_adapt_solver_memory = [](const KINView& kinmem,
                                                    N_Vector constraints) -> int
    {
      auto lambda_result = KINSetConstraints(kinmem.get(), constraints);
      return lambda_result;
    };

    return KINSetConstraints_adapt_solver_memory(kinmem, constraints);
  },
  nb::arg("kinmem"), nb::arg("constraints"));

sundials4py::scoped_def(
  m, "KINGetNumNonlinSolvIters",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nniters_adapt_modifiable;

      int r = KINGetNumNonlinSolvIters(kinmem, &nniters_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable);
    };
    auto KINGetNumNonlinSolvIters_adapt_solver_memory =
      [&KINGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
          kinmem.get());
      return lambda_result;
    };

    return KINGetNumNonlinSolvIters_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumFuncEvals",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumFuncEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nfevals_adapt_modifiable;

      int r = KINGetNumFuncEvals(kinmem, &nfevals_adapt_modifiable);
      return std::make_tuple(r, nfevals_adapt_modifiable);
    };
    auto KINGetNumFuncEvals_adapt_solver_memory =
      [&KINGetNumFuncEvals_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumFuncEvals_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumFuncEvals_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumBetaCondFails",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumBetaCondFails_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nbcfails_adapt_modifiable;

      int r = KINGetNumBetaCondFails(kinmem, &nbcfails_adapt_modifiable);
      return std::make_tuple(r, nbcfails_adapt_modifiable);
    };
    auto KINGetNumBetaCondFails_adapt_solver_memory =
      [&KINGetNumBetaCondFails_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumBetaCondFails_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumBetaCondFails_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumBacktrackOps",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumBacktrackOps_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nbacktr_adapt_modifiable;

      int r = KINGetNumBacktrackOps(kinmem, &nbacktr_adapt_modifiable);
      return std::make_tuple(r, nbacktr_adapt_modifiable);
    };
    auto KINGetNumBacktrackOps_adapt_solver_memory =
      [&KINGetNumBacktrackOps_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumBacktrackOps_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumBacktrackOps_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetFuncNorm",
  [](const KINView& kinmem) -> std::tuple<int, sunrealtype>
  {
    auto KINGetFuncNorm_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype fnorm_adapt_modifiable;

      int r = KINGetFuncNorm(kinmem, &fnorm_adapt_modifiable);
      return std::make_tuple(r, fnorm_adapt_modifiable);
    };
    auto KINGetFuncNorm_adapt_solver_memory =
      [&KINGetFuncNorm_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        KINGetFuncNorm_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetFuncNorm_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetStepLength",
  [](const KINView& kinmem) -> std::tuple<int, sunrealtype>
  {
    auto KINGetStepLength_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype steplength_adapt_modifiable;

      int r = KINGetStepLength(kinmem, &steplength_adapt_modifiable);
      return std::make_tuple(r, steplength_adapt_modifiable);
    };
    auto KINGetStepLength_adapt_solver_memory =
      [&KINGetStepLength_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        KINGetStepLength_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetStepLength_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINPrintAllStats",
  [](const KINView& kinmem, FILE* outfile, SUNOutputFormat fmt) -> int
  {
    auto KINPrintAllStats_adapt_solver_memory =
      [](const KINView& kinmem, FILE* outfile, SUNOutputFormat fmt) -> int
    {
      auto lambda_result = KINPrintAllStats(kinmem.get(), outfile, fmt);
      return lambda_result;
    };

    return KINPrintAllStats_adapt_solver_memory(kinmem, outfile, fmt);
  },
  nb::arg("kinmem"), nb::arg("outfile"), nb::arg("fmt"));

sundials4py::scoped_def(m, "KINGetReturnFlagName", KINGetReturnFlagName,
                        nb::arg("flag"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _KINLS_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("KINLS_SUCCESS")     = 0;
m.attr("KINLS_MEM_NULL")    = -1;
m.attr("KINLS_LMEM_NULL")   = -2;
m.attr("KINLS_ILL_INPUT")   = -3;
m.attr("KINLS_MEM_FAIL")    = -4;
m.attr("KINLS_PMEM_NULL")   = -5;
m.attr("KINLS_JACFUNC_ERR") = -6;
m.attr("KINLS_SUNMAT_FAIL") = -7;
m.attr("KINLS_SUNLS_FAIL")  = -8;

sundials4py::scoped_def(
  m, "KINSetLinearSolver",
  [](const KINView& kinmem, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto KINSetLinearSolver_adapt_optional_arg_with_default_null =
      [](void* kinmem, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = KINSetLinearSolver(kinmem, LS, A_adapt_default_null);
      return lambda_result;
    };
    auto KINSetLinearSolver_adapt_solver_memory =
      [&KINSetLinearSolver_adapt_optional_arg_with_default_null](const KINView& kinmem,
                                                                 SUNLinearSolver LS,
                                                                 std::optional<SUNMatrix> A =
                                                                   std::nullopt) -> int
    {
      auto lambda_result =
        KINSetLinearSolver_adapt_optional_arg_with_default_null(kinmem.get(),
                                                                LS, A);
      return lambda_result;
    };

    return KINSetLinearSolver_adapt_solver_memory(kinmem, LS, A);
  },
  nb::arg("kinmem"), nb::arg("LS"), nb::arg("A").none() = nb::none());

sundials4py::scoped_def(
  m, "KINGetJac",
  [](const KINView& kinmem) -> std::tuple<int, SUNMatrix>
  {
    auto KINGetJac_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix J_adapt_modifiable;

      int r = KINGetJac(kinmem, &J_adapt_modifiable);
      return std::make_tuple(r, J_adapt_modifiable);
    };
    auto KINGetJac_adapt_solver_memory =
      [&KINGetJac_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, SUNMatrix>
    {
      auto lambda_result =
        KINGetJac_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetJac_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "KINGetJacNumIters",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetJacNumIters_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nni_J_adapt_modifiable;

      int r = KINGetJacNumIters(kinmem, &nni_J_adapt_modifiable);
      return std::make_tuple(r, nni_J_adapt_modifiable);
    };
    auto KINGetJacNumIters_adapt_solver_memory =
      [&KINGetJacNumIters_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetJacNumIters_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetJacNumIters_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumJacEvals",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long njevals_adapt_modifiable;

      int r = KINGetNumJacEvals(kinmem, &njevals_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable);
    };
    auto KINGetNumJacEvals_adapt_solver_memory =
      [&KINGetNumJacEvals_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumJacEvals_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumJacEvals_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumLinFuncEvals",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumLinFuncEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nfevals_adapt_modifiable;

      int r = KINGetNumLinFuncEvals(kinmem, &nfevals_adapt_modifiable);
      return std::make_tuple(r, nfevals_adapt_modifiable);
    };
    auto KINGetNumLinFuncEvals_adapt_solver_memory =
      [&KINGetNumLinFuncEvals_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumLinFuncEvals_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumLinFuncEvals_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumPrecEvals",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long npevals_adapt_modifiable;

      int r = KINGetNumPrecEvals(kinmem, &npevals_adapt_modifiable);
      return std::make_tuple(r, npevals_adapt_modifiable);
    };
    auto KINGetNumPrecEvals_adapt_solver_memory =
      [&KINGetNumPrecEvals_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumPrecEvals_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumPrecEvals_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumPrecSolves",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long npsolves_adapt_modifiable;

      int r = KINGetNumPrecSolves(kinmem, &npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves_adapt_modifiable);
    };
    auto KINGetNumPrecSolves_adapt_solver_memory =
      [&KINGetNumPrecSolves_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumPrecSolves_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumPrecSolves_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumLinIters",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nliters_adapt_modifiable;

      int r = KINGetNumLinIters(kinmem, &nliters_adapt_modifiable);
      return std::make_tuple(r, nliters_adapt_modifiable);
    };
    auto KINGetNumLinIters_adapt_solver_memory =
      [&KINGetNumLinIters_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumLinIters_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumLinIters_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumLinConvFails",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nlcfails_adapt_modifiable;

      int r = KINGetNumLinConvFails(kinmem, &nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails_adapt_modifiable);
    };
    auto KINGetNumLinConvFails_adapt_solver_memory =
      [&KINGetNumLinConvFails_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumLinConvFails_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumLinConvFails_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumJtimesEvals",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long njvevals_adapt_modifiable;

      int r = KINGetNumJtimesEvals(kinmem, &njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals_adapt_modifiable);
    };
    auto KINGetNumJtimesEvals_adapt_solver_memory =
      [&KINGetNumJtimesEvals_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetNumJtimesEvals_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetNumJtimesEvals_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetLastLinFlag",
  [](const KINView& kinmem) -> std::tuple<int, long>
  {
    auto KINGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = KINGetLastLinFlag(kinmem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };
    auto KINGetLastLinFlag_adapt_solver_memory =
      [&KINGetLastLinFlag_adapt_modifiable_immutable_to_return](
        const KINView& kinmem) -> std::tuple<int, long>
    {
      auto lambda_result =
        KINGetLastLinFlag_adapt_modifiable_immutable_to_return(kinmem.get());
      return lambda_result;
    };

    return KINGetLastLinFlag_adapt_solver_memory(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(m, "KINGetLinReturnFlagName", KINGetLinReturnFlagName,
                        nb::arg("flag"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
