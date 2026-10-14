// #ifndef _IDAS_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("IDA_NORMAL")             = 1;
m.attr("IDA_ONE_STEP")           = 2;
m.attr("IDA_YA_YDP_INIT")        = 1;
m.attr("IDA_Y_INIT")             = 2;
m.attr("IDA_SIMULTANEOUS")       = 1;
m.attr("IDA_STAGGERED")          = 2;
m.attr("IDA_CENTERED")           = 1;
m.attr("IDA_FORWARD")            = 2;
m.attr("IDA_HERMITE")            = 1;
m.attr("IDA_POLYNOMIAL")         = 2;
m.attr("IDA_SUCCESS")            = 0;
m.attr("IDA_TSTOP_RETURN")       = 1;
m.attr("IDA_ROOT_RETURN")        = 2;
m.attr("IDA_WARNING")            = 99;
m.attr("IDA_TOO_MUCH_WORK")      = -1;
m.attr("IDA_TOO_MUCH_ACC")       = -2;
m.attr("IDA_ERR_FAIL")           = -3;
m.attr("IDA_CONV_FAIL")          = -4;
m.attr("IDA_LINIT_FAIL")         = -5;
m.attr("IDA_LSETUP_FAIL")        = -6;
m.attr("IDA_LSOLVE_FAIL")        = -7;
m.attr("IDA_RES_FAIL")           = -8;
m.attr("IDA_REP_RES_ERR")        = -9;
m.attr("IDA_RTFUNC_FAIL")        = -10;
m.attr("IDA_CONSTR_FAIL")        = -11;
m.attr("IDA_FIRST_RES_FAIL")     = -12;
m.attr("IDA_LINESEARCH_FAIL")    = -13;
m.attr("IDA_NO_RECOVERY")        = -14;
m.attr("IDA_NLS_INIT_FAIL")      = -15;
m.attr("IDA_NLS_SETUP_FAIL")     = -16;
m.attr("IDA_NLS_FAIL")           = -17;
m.attr("IDA_MEM_NULL")           = -20;
m.attr("IDA_MEM_FAIL")           = -21;
m.attr("IDA_ILL_INPUT")          = -22;
m.attr("IDA_NO_MALLOC")          = -23;
m.attr("IDA_BAD_EWT")            = -24;
m.attr("IDA_BAD_K")              = -25;
m.attr("IDA_BAD_T")              = -26;
m.attr("IDA_BAD_DKY")            = -27;
m.attr("IDA_VECTOROP_ERR")       = -28;
m.attr("IDA_CONTEXT_ERR")        = -29;
m.attr("IDA_NO_QUAD")            = -30;
m.attr("IDA_QRHS_FAIL")          = -31;
m.attr("IDA_FIRST_QRHS_ERR")     = -32;
m.attr("IDA_REP_QRHS_ERR")       = -33;
m.attr("IDA_NO_SENS")            = -40;
m.attr("IDA_SRES_FAIL")          = -41;
m.attr("IDA_REP_SRES_ERR")       = -42;
m.attr("IDA_BAD_IS")             = -43;
m.attr("IDA_NO_QUADSENS")        = -50;
m.attr("IDA_QSRHS_FAIL")         = -51;
m.attr("IDA_FIRST_QSRHS_ERR")    = -52;
m.attr("IDA_REP_QSRHS_ERR")      = -53;
m.attr("IDA_TOO_CLOSE")          = -60;
m.attr("IDA_UNRECOGNIZED_ERROR") = -99;
m.attr("IDA_NO_ADJ")             = -101;
m.attr("IDA_NO_FWD")             = -102;
m.attr("IDA_NO_BCK")             = -103;
m.attr("IDA_BAD_TB0")            = -104;
m.attr("IDA_REIFWD_FAIL")        = -105;
m.attr("IDA_FWD_FAIL")           = -106;
m.attr("IDA_GETY_BADT")          = -107;

sundials4py::scoped_def(
  m, "IDAReInit",
  [](const IDAView& ida_mem, sunrealtype t0, N_Vector yy0, N_Vector yp0) -> int
  {
    auto IDAReInit_adapt_solver_memory = [](const IDAView& ida_mem,
                                            sunrealtype t0, N_Vector yy0,
                                            N_Vector yp0) -> int
    {
      auto lambda_result = IDAReInit(ida_mem.get(), t0, yy0, yp0);
      return lambda_result;
    };

    return IDAReInit_adapt_solver_memory(ida_mem, t0, yy0, yp0);
  },
  nb::arg("ida_mem"), nb::arg("t0"), nb::arg("yy0"), nb::arg("yp0"));

sundials4py::scoped_def(
  m, "IDASStolerances",
  [](const IDAView& ida_mem, sunrealtype reltol, sunrealtype abstol) -> int
  {
    auto IDASStolerances_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype reltol, sunrealtype abstol) -> int
    {
      auto lambda_result = IDASStolerances(ida_mem.get(), reltol, abstol);
      return lambda_result;
    };

    return IDASStolerances_adapt_solver_memory(ida_mem, reltol, abstol);
  },
  nb::arg("ida_mem"), nb::arg("reltol"), nb::arg("abstol"));

sundials4py::scoped_def(
  m, "IDASVtolerances",
  [](const IDAView& ida_mem, sunrealtype reltol, N_Vector abstol) -> int
  {
    auto IDASVtolerances_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype reltol, N_Vector abstol) -> int
    {
      auto lambda_result = IDASVtolerances(ida_mem.get(), reltol, abstol);
      return lambda_result;
    };

    return IDASVtolerances_adapt_solver_memory(ida_mem, reltol, abstol);
  },
  nb::arg("ida_mem"), nb::arg("reltol"), nb::arg("abstol"));

sundials4py::scoped_def(
  m, "IDACalcIC",
  [](const IDAView& ida_mem, int icopt, sunrealtype tout1) -> int
  {
    auto IDACalcIC_adapt_solver_memory = [](const IDAView& ida_mem, int icopt,
                                            sunrealtype tout1) -> int
    {
      auto lambda_result = IDACalcIC(ida_mem.get(), icopt, tout1);
      return lambda_result;
    };

    return IDACalcIC_adapt_solver_memory(ida_mem, icopt, tout1);
  },
  nb::arg("ida_mem"), nb::arg("icopt"), nb::arg("tout1"),
  "Initial condition calculation function");

sundials4py::scoped_def(
  m, "IDASetNonlinConvCoefIC",
  [](const IDAView& ida_mem, sunrealtype epiccon) -> int
  {
    auto IDASetNonlinConvCoefIC_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype epiccon) -> int
    {
      auto lambda_result = IDASetNonlinConvCoefIC(ida_mem.get(), epiccon);
      return lambda_result;
    };

    return IDASetNonlinConvCoefIC_adapt_solver_memory(ida_mem, epiccon);
  },
  nb::arg("ida_mem"), nb::arg("epiccon"));

sundials4py::scoped_def(
  m, "IDASetMaxNumStepsIC",
  [](const IDAView& ida_mem, int maxnh) -> int
  {
    auto IDASetMaxNumStepsIC_adapt_solver_memory = [](const IDAView& ida_mem,
                                                      int maxnh) -> int
    {
      auto lambda_result = IDASetMaxNumStepsIC(ida_mem.get(), maxnh);
      return lambda_result;
    };

    return IDASetMaxNumStepsIC_adapt_solver_memory(ida_mem, maxnh);
  },
  nb::arg("ida_mem"), nb::arg("maxnh"));

sundials4py::scoped_def(
  m, "IDASetMaxNumJacsIC",
  [](const IDAView& ida_mem, int maxnj) -> int
  {
    auto IDASetMaxNumJacsIC_adapt_solver_memory = [](const IDAView& ida_mem,
                                                     int maxnj) -> int
    {
      auto lambda_result = IDASetMaxNumJacsIC(ida_mem.get(), maxnj);
      return lambda_result;
    };

    return IDASetMaxNumJacsIC_adapt_solver_memory(ida_mem, maxnj);
  },
  nb::arg("ida_mem"), nb::arg("maxnj"));

sundials4py::scoped_def(
  m, "IDASetMaxNumItersIC",
  [](const IDAView& ida_mem, int maxnit) -> int
  {
    auto IDASetMaxNumItersIC_adapt_solver_memory = [](const IDAView& ida_mem,
                                                      int maxnit) -> int
    {
      auto lambda_result = IDASetMaxNumItersIC(ida_mem.get(), maxnit);
      return lambda_result;
    };

    return IDASetMaxNumItersIC_adapt_solver_memory(ida_mem, maxnit);
  },
  nb::arg("ida_mem"), nb::arg("maxnit"));

sundials4py::scoped_def(
  m, "IDASetLineSearchOffIC",
  [](const IDAView& ida_mem, sunbooleantype lsoff) -> int
  {
    auto IDASetLineSearchOffIC_adapt_solver_memory =
      [](const IDAView& ida_mem, sunbooleantype lsoff) -> int
    {
      auto lambda_result = IDASetLineSearchOffIC(ida_mem.get(), lsoff);
      return lambda_result;
    };

    return IDASetLineSearchOffIC_adapt_solver_memory(ida_mem, lsoff);
  },
  nb::arg("ida_mem"), nb::arg("lsoff"));

sundials4py::scoped_def(
  m, "IDASetStepToleranceIC",
  [](const IDAView& ida_mem, sunrealtype steptol) -> int
  {
    auto IDASetStepToleranceIC_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype steptol) -> int
    {
      auto lambda_result = IDASetStepToleranceIC(ida_mem.get(), steptol);
      return lambda_result;
    };

    return IDASetStepToleranceIC_adapt_solver_memory(ida_mem, steptol);
  },
  nb::arg("ida_mem"), nb::arg("steptol"));

sundials4py::scoped_def(
  m, "IDASetMaxBacksIC",
  [](const IDAView& ida_mem, int maxbacks) -> int
  {
    auto IDASetMaxBacksIC_adapt_solver_memory = [](const IDAView& ida_mem,
                                                   int maxbacks) -> int
    {
      auto lambda_result = IDASetMaxBacksIC(ida_mem.get(), maxbacks);
      return lambda_result;
    };

    return IDASetMaxBacksIC_adapt_solver_memory(ida_mem, maxbacks);
  },
  nb::arg("ida_mem"), nb::arg("maxbacks"));

sundials4py::scoped_def(m, "IDASetDeltaCjLSetup", IDASetDeltaCjLSetup,
                        nb::arg("ida_max"), nb::arg("dcj"));

sundials4py::scoped_def(
  m, "IDASetMaxOrd",
  [](const IDAView& ida_mem, int maxord) -> int
  {
    auto IDASetMaxOrd_adapt_solver_memory = [](const IDAView& ida_mem,
                                               int maxord) -> int
    {
      auto lambda_result = IDASetMaxOrd(ida_mem.get(), maxord);
      return lambda_result;
    };

    return IDASetMaxOrd_adapt_solver_memory(ida_mem, maxord);
  },
  nb::arg("ida_mem"), nb::arg("maxord"));

sundials4py::scoped_def(
  m, "IDASetMaxNumSteps",
  [](const IDAView& ida_mem, long mxsteps) -> int
  {
    auto IDASetMaxNumSteps_adapt_solver_memory = [](const IDAView& ida_mem,
                                                    long mxsteps) -> int
    {
      auto lambda_result = IDASetMaxNumSteps(ida_mem.get(), mxsteps);
      return lambda_result;
    };

    return IDASetMaxNumSteps_adapt_solver_memory(ida_mem, mxsteps);
  },
  nb::arg("ida_mem"), nb::arg("mxsteps"));

sundials4py::scoped_def(
  m, "IDASetInitStep",
  [](const IDAView& ida_mem, sunrealtype hin) -> int
  {
    auto IDASetInitStep_adapt_solver_memory = [](const IDAView& ida_mem,
                                                 sunrealtype hin) -> int
    {
      auto lambda_result = IDASetInitStep(ida_mem.get(), hin);
      return lambda_result;
    };

    return IDASetInitStep_adapt_solver_memory(ida_mem, hin);
  },
  nb::arg("ida_mem"), nb::arg("hin"));

sundials4py::scoped_def(
  m, "IDASetMaxStep",
  [](const IDAView& ida_mem, sunrealtype hmax) -> int
  {
    auto IDASetMaxStep_adapt_solver_memory = [](const IDAView& ida_mem,
                                                sunrealtype hmax) -> int
    {
      auto lambda_result = IDASetMaxStep(ida_mem.get(), hmax);
      return lambda_result;
    };

    return IDASetMaxStep_adapt_solver_memory(ida_mem, hmax);
  },
  nb::arg("ida_mem"), nb::arg("hmax"));

sundials4py::scoped_def(
  m, "IDASetMinStep",
  [](const IDAView& ida_mem, sunrealtype hmin) -> int
  {
    auto IDASetMinStep_adapt_solver_memory = [](const IDAView& ida_mem,
                                                sunrealtype hmin) -> int
    {
      auto lambda_result = IDASetMinStep(ida_mem.get(), hmin);
      return lambda_result;
    };

    return IDASetMinStep_adapt_solver_memory(ida_mem, hmin);
  },
  nb::arg("ida_mem"), nb::arg("hmin"));

sundials4py::scoped_def(
  m, "IDASetStopTime",
  [](const IDAView& ida_mem, sunrealtype tstop) -> int
  {
    auto IDASetStopTime_adapt_solver_memory = [](const IDAView& ida_mem,
                                                 sunrealtype tstop) -> int
    {
      auto lambda_result = IDASetStopTime(ida_mem.get(), tstop);
      return lambda_result;
    };

    return IDASetStopTime_adapt_solver_memory(ida_mem, tstop);
  },
  nb::arg("ida_mem"), nb::arg("tstop"));

sundials4py::scoped_def(
  m, "IDAClearStopTime",
  [](const IDAView& ida_mem) -> int
  {
    auto IDAClearStopTime_adapt_solver_memory = [](const IDAView& ida_mem) -> int
    {
      auto lambda_result = IDAClearStopTime(ida_mem.get());
      return lambda_result;
    };

    return IDAClearStopTime_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDASetMaxErrTestFails",
  [](const IDAView& ida_mem, int maxnef) -> int
  {
    auto IDASetMaxErrTestFails_adapt_solver_memory = [](const IDAView& ida_mem,
                                                        int maxnef) -> int
    {
      auto lambda_result = IDASetMaxErrTestFails(ida_mem.get(), maxnef);
      return lambda_result;
    };

    return IDASetMaxErrTestFails_adapt_solver_memory(ida_mem, maxnef);
  },
  nb::arg("ida_mem"), nb::arg("maxnef"));

sundials4py::scoped_def(
  m, "IDASetSuppressAlg",
  [](const IDAView& ida_mem, sunbooleantype suppressalg) -> int
  {
    auto IDASetSuppressAlg_adapt_solver_memory =
      [](const IDAView& ida_mem, sunbooleantype suppressalg) -> int
    {
      auto lambda_result = IDASetSuppressAlg(ida_mem.get(), suppressalg);
      return lambda_result;
    };

    return IDASetSuppressAlg_adapt_solver_memory(ida_mem, suppressalg);
  },
  nb::arg("ida_mem"), nb::arg("suppressalg"));

sundials4py::scoped_def(
  m, "IDASetId",
  [](const IDAView& ida_mem, N_Vector id) -> int
  {
    auto IDASetId_adapt_solver_memory = [](const IDAView& ida_mem,
                                           N_Vector id) -> int
    {
      auto lambda_result = IDASetId(ida_mem.get(), id);
      return lambda_result;
    };

    return IDASetId_adapt_solver_memory(ida_mem, id);
  },
  nb::arg("ida_mem"), nb::arg("id"));

sundials4py::scoped_def(
  m, "IDASetConstraints",
  [](const IDAView& ida_mem, N_Vector constraints) -> int
  {
    auto IDASetConstraints_adapt_solver_memory = [](const IDAView& ida_mem,
                                                    N_Vector constraints) -> int
    {
      auto lambda_result = IDASetConstraints(ida_mem.get(), constraints);
      return lambda_result;
    };

    return IDASetConstraints_adapt_solver_memory(ida_mem, constraints);
  },
  nb::arg("ida_mem"), nb::arg("constraints"));

sundials4py::scoped_def(
  m, "IDASetMaxNumConstraintFails",
  [](const IDAView& ida_mem, int max_fails) -> int
  {
    auto IDASetMaxNumConstraintFails_adapt_solver_memory =
      [](const IDAView& ida_mem, int max_fails) -> int
    {
      auto lambda_result = IDASetMaxNumConstraintFails(ida_mem.get(), max_fails);
      return lambda_result;
    };

    return IDASetMaxNumConstraintFails_adapt_solver_memory(ida_mem, max_fails);
  },
  nb::arg("ida_mem"), nb::arg("max_fails"));

sundials4py::scoped_def(
  m, "IDASetEtaFixedStepBounds",
  [](const IDAView& ida_mem, sunrealtype eta_min_fx, sunrealtype eta_max_fx) -> int
  {
    auto IDASetEtaFixedStepBounds_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype eta_min_fx,
         sunrealtype eta_max_fx) -> int
    {
      auto lambda_result = IDASetEtaFixedStepBounds(ida_mem.get(), eta_min_fx,
                                                    eta_max_fx);
      return lambda_result;
    };

    return IDASetEtaFixedStepBounds_adapt_solver_memory(ida_mem, eta_min_fx,
                                                        eta_max_fx);
  },
  nb::arg("ida_mem"), nb::arg("eta_min_fx"), nb::arg("eta_max_fx"));

sundials4py::scoped_def(
  m, "IDASetEtaMin",
  [](const IDAView& ida_mem, sunrealtype eta_min) -> int
  {
    auto IDASetEtaMin_adapt_solver_memory = [](const IDAView& ida_mem,
                                               sunrealtype eta_min) -> int
    {
      auto lambda_result = IDASetEtaMin(ida_mem.get(), eta_min);
      return lambda_result;
    };

    return IDASetEtaMin_adapt_solver_memory(ida_mem, eta_min);
  },
  nb::arg("ida_mem"), nb::arg("eta_min"));

sundials4py::scoped_def(
  m, "IDASetEtaMax",
  [](const IDAView& ida_mem, sunrealtype eta_max) -> int
  {
    auto IDASetEtaMax_adapt_solver_memory = [](const IDAView& ida_mem,
                                               sunrealtype eta_max) -> int
    {
      auto lambda_result = IDASetEtaMax(ida_mem.get(), eta_max);
      return lambda_result;
    };

    return IDASetEtaMax_adapt_solver_memory(ida_mem, eta_max);
  },
  nb::arg("ida_mem"), nb::arg("eta_max"));

sundials4py::scoped_def(
  m, "IDASetEtaLow",
  [](const IDAView& ida_mem, sunrealtype eta_low) -> int
  {
    auto IDASetEtaLow_adapt_solver_memory = [](const IDAView& ida_mem,
                                               sunrealtype eta_low) -> int
    {
      auto lambda_result = IDASetEtaLow(ida_mem.get(), eta_low);
      return lambda_result;
    };

    return IDASetEtaLow_adapt_solver_memory(ida_mem, eta_low);
  },
  nb::arg("ida_mem"), nb::arg("eta_low"));

sundials4py::scoped_def(
  m, "IDASetEtaMinErrFail",
  [](const IDAView& ida_mem, sunrealtype eta_min_ef) -> int
  {
    auto IDASetEtaMinErrFail_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype eta_min_ef) -> int
    {
      auto lambda_result = IDASetEtaMinErrFail(ida_mem.get(), eta_min_ef);
      return lambda_result;
    };

    return IDASetEtaMinErrFail_adapt_solver_memory(ida_mem, eta_min_ef);
  },
  nb::arg("ida_mem"), nb::arg("eta_min_ef"));

sundials4py::scoped_def(
  m, "IDASetEtaConvFail",
  [](const IDAView& ida_mem, sunrealtype eta_cf) -> int
  {
    auto IDASetEtaConvFail_adapt_solver_memory = [](const IDAView& ida_mem,
                                                    sunrealtype eta_cf) -> int
    {
      auto lambda_result = IDASetEtaConvFail(ida_mem.get(), eta_cf);
      return lambda_result;
    };

    return IDASetEtaConvFail_adapt_solver_memory(ida_mem, eta_cf);
  },
  nb::arg("ida_mem"), nb::arg("eta_cf"));

sundials4py::scoped_def(
  m, "IDASetMaxConvFails",
  [](const IDAView& ida_mem, int maxncf) -> int
  {
    auto IDASetMaxConvFails_adapt_solver_memory = [](const IDAView& ida_mem,
                                                     int maxncf) -> int
    {
      auto lambda_result = IDASetMaxConvFails(ida_mem.get(), maxncf);
      return lambda_result;
    };

    return IDASetMaxConvFails_adapt_solver_memory(ida_mem, maxncf);
  },
  nb::arg("ida_mem"), nb::arg("maxncf"));

sundials4py::scoped_def(
  m, "IDASetMaxNonlinIters",
  [](const IDAView& ida_mem, int maxcor) -> int
  {
    auto IDASetMaxNonlinIters_adapt_solver_memory = [](const IDAView& ida_mem,
                                                       int maxcor) -> int
    {
      auto lambda_result = IDASetMaxNonlinIters(ida_mem.get(), maxcor);
      return lambda_result;
    };

    return IDASetMaxNonlinIters_adapt_solver_memory(ida_mem, maxcor);
  },
  nb::arg("ida_mem"), nb::arg("maxcor"));

sundials4py::scoped_def(
  m, "IDASetNonlinConvCoef",
  [](const IDAView& ida_mem, sunrealtype epcon) -> int
  {
    auto IDASetNonlinConvCoef_adapt_solver_memory = [](const IDAView& ida_mem,
                                                       sunrealtype epcon) -> int
    {
      auto lambda_result = IDASetNonlinConvCoef(ida_mem.get(), epcon);
      return lambda_result;
    };

    return IDASetNonlinConvCoef_adapt_solver_memory(ida_mem, epcon);
  },
  nb::arg("ida_mem"), nb::arg("epcon"));

sundials4py::scoped_def(
  m, "IDASetNonlinearSolver",
  [](const IDAView& ida_mem, SUNNonlinearSolver NLS) -> int
  {
    auto IDASetNonlinearSolver_adapt_solver_memory =
      [](const IDAView& ida_mem, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = IDASetNonlinearSolver(ida_mem.get(), NLS);
      return lambda_result;
    };

    return IDASetNonlinearSolver_adapt_solver_memory(ida_mem, NLS);
  },
  nb::arg("ida_mem"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "IDASetRootDirection",
  [](const IDAView& ida_mem, std::vector<int> rootdir_1d) -> int
  {
    auto IDASetRootDirection_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<int> rootdir_1d) -> int
    {
      int* rootdir_1d_ptr = rootdir_1d.empty() ? nullptr : rootdir_1d.data();

      auto lambda_result = IDASetRootDirection(ida_mem, rootdir_1d_ptr);
      return lambda_result;
    };
    auto IDASetRootDirection_adapt_solver_memory =
      [&IDASetRootDirection_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                         std::vector<int> rootdir_1d)
      -> int
    {
      auto lambda_result =
        IDASetRootDirection_adapt_arr_ptr_to_std_vector(ida_mem.get(),
                                                        rootdir_1d);
      return lambda_result;
    };

    return IDASetRootDirection_adapt_solver_memory(ida_mem, rootdir_1d);
  },
  nb::arg("ida_mem"), nb::arg("rootdir_1d"));

sundials4py::scoped_def(
  m, "IDASetNoInactiveRootWarn",
  [](const IDAView& ida_mem) -> int
  {
    auto IDASetNoInactiveRootWarn_adapt_solver_memory =
      [](const IDAView& ida_mem) -> int
    {
      auto lambda_result = IDASetNoInactiveRootWarn(ida_mem.get());
      return lambda_result;
    };

    return IDASetNoInactiveRootWarn_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAComputeY",
  [](const IDAView& ida_mem, N_Vector ycor, N_Vector y) -> int
  {
    auto IDAComputeY_adapt_solver_memory = [](const IDAView& ida_mem,
                                              N_Vector ycor, N_Vector y) -> int
    {
      auto lambda_result = IDAComputeY(ida_mem.get(), ycor, y);
      return lambda_result;
    };

    return IDAComputeY_adapt_solver_memory(ida_mem, ycor, y);
  },
  nb::arg("ida_mem"), nb::arg("ycor"), nb::arg("y"));

sundials4py::scoped_def(
  m, "IDAComputeYp",
  [](const IDAView& ida_mem, N_Vector ycor, N_Vector yp) -> int
  {
    auto IDAComputeYp_adapt_solver_memory = [](const IDAView& ida_mem,
                                               N_Vector ycor, N_Vector yp) -> int
    {
      auto lambda_result = IDAComputeYp(ida_mem.get(), ycor, yp);
      return lambda_result;
    };

    return IDAComputeYp_adapt_solver_memory(ida_mem, ycor, yp);
  },
  nb::arg("ida_mem"), nb::arg("ycor"), nb::arg("yp"));

sundials4py::scoped_def(
  m, "IDAComputeYSens",
  [](const IDAView& ida_mem, std::vector<N_Vector> ycor_1d,
     std::vector<N_Vector> yyS_1d) -> int
  {
    auto IDAComputeYSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> ycor_1d,
         std::vector<N_Vector> yyS_1d) -> int
    {
      N_Vector* ycor_1d_ptr = ycor_1d.empty() ? nullptr : ycor_1d.data();
      N_Vector* yyS_1d_ptr  = yyS_1d.empty() ? nullptr : yyS_1d.data();

      auto lambda_result = IDAComputeYSens(ida_mem, ycor_1d_ptr, yyS_1d_ptr);
      return lambda_result;
    };
    auto IDAComputeYSens_adapt_solver_memory =
      [&IDAComputeYSens_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                     std::vector<N_Vector> ycor_1d,
                                                     std::vector<N_Vector> yyS_1d)
      -> int
    {
      auto lambda_result =
        IDAComputeYSens_adapt_arr_ptr_to_std_vector(ida_mem.get(), ycor_1d,
                                                    yyS_1d);
      return lambda_result;
    };

    return IDAComputeYSens_adapt_solver_memory(ida_mem, ycor_1d, yyS_1d);
  },
  nb::arg("ida_mem"), nb::arg("ycor_1d"), nb::arg("yyS_1d"));

sundials4py::scoped_def(
  m, "IDAComputeYpSens",
  [](const IDAView& ida_mem, std::vector<N_Vector> ycor_1d,
     std::vector<N_Vector> ypS_1d) -> int
  {
    auto IDAComputeYpSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> ycor_1d,
         std::vector<N_Vector> ypS_1d) -> int
    {
      N_Vector* ycor_1d_ptr = ycor_1d.empty() ? nullptr : ycor_1d.data();
      N_Vector* ypS_1d_ptr  = ypS_1d.empty() ? nullptr : ypS_1d.data();

      auto lambda_result = IDAComputeYpSens(ida_mem, ycor_1d_ptr, ypS_1d_ptr);
      return lambda_result;
    };
    auto IDAComputeYpSens_adapt_solver_memory =
      [&IDAComputeYpSens_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                      std::vector<N_Vector> ycor_1d,
                                                      std::vector<N_Vector> ypS_1d)
      -> int
    {
      auto lambda_result =
        IDAComputeYpSens_adapt_arr_ptr_to_std_vector(ida_mem.get(), ycor_1d,
                                                     ypS_1d);
      return lambda_result;
    };

    return IDAComputeYpSens_adapt_solver_memory(ida_mem, ycor_1d, ypS_1d);
  },
  nb::arg("ida_mem"), nb::arg("ycor_1d"), nb::arg("ypS_1d"));

sundials4py::scoped_def(
  m, "IDAGetDky",
  [](const IDAView& ida_mem, sunrealtype t, int k, N_Vector dky) -> int
  {
    auto IDAGetDky_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype t, int k, N_Vector dky) -> int
    {
      auto lambda_result = IDAGetDky(ida_mem.get(), t, k, dky);
      return lambda_result;
    };

    return IDAGetDky_adapt_solver_memory(ida_mem, t, k, dky);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dky"));

sundials4py::scoped_def(
  m, "IDAGetNumSteps",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumSteps_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nsteps_adapt_modifiable;

      int r = IDAGetNumSteps(ida_mem, &nsteps_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable);
    };
    auto IDAGetNumSteps_adapt_solver_memory =
      [&IDAGetNumSteps_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumSteps_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumSteps_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumResEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumResEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nrevals_adapt_modifiable;

      int r = IDAGetNumResEvals(ida_mem, &nrevals_adapt_modifiable);
      return std::make_tuple(r, nrevals_adapt_modifiable);
    };
    auto IDAGetNumResEvals_adapt_solver_memory =
      [&IDAGetNumResEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumResEvals_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumResEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumLinSolvSetups",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nlinsetups_adapt_modifiable;

      int r = IDAGetNumLinSolvSetups(ida_mem, &nlinsetups_adapt_modifiable);
      return std::make_tuple(r, nlinsetups_adapt_modifiable);
    };
    auto IDAGetNumLinSolvSetups_adapt_solver_memory =
      [&IDAGetNumLinSolvSetups_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumLinSolvSetups_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumLinSolvSetups_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumErrTestFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long netfails_adapt_modifiable;

      int r = IDAGetNumErrTestFails(ida_mem, &netfails_adapt_modifiable);
      return std::make_tuple(r, netfails_adapt_modifiable);
    };
    auto IDAGetNumErrTestFails_adapt_solver_memory =
      [&IDAGetNumErrTestFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumErrTestFails_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumErrTestFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumBacktrackOps",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumBacktrackOps_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nbacktr_adapt_modifiable;

      int r = IDAGetNumBacktrackOps(ida_mem, &nbacktr_adapt_modifiable);
      return std::make_tuple(r, nbacktr_adapt_modifiable);
    };
    auto IDAGetNumBacktrackOps_adapt_solver_memory =
      [&IDAGetNumBacktrackOps_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumBacktrackOps_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumBacktrackOps_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetConsistentIC",
  [](const IDAView& ida_mem, N_Vector yy0_mod, N_Vector yp0_mod) -> int
  {
    auto IDAGetConsistentIC_adapt_solver_memory =
      [](const IDAView& ida_mem, N_Vector yy0_mod, N_Vector yp0_mod) -> int
    {
      auto lambda_result = IDAGetConsistentIC(ida_mem.get(), yy0_mod, yp0_mod);
      return lambda_result;
    };

    return IDAGetConsistentIC_adapt_solver_memory(ida_mem, yy0_mod, yp0_mod);
  },
  nb::arg("ida_mem"), nb::arg("yy0_mod"), nb::arg("yp0_mod"));

sundials4py::scoped_def(
  m, "IDAGetLastOrder",
  [](const IDAView& ida_mem) -> std::tuple<int, int>
  {
    auto IDAGetLastOrder_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, int>
    {
      int klast_adapt_modifiable;

      int r = IDAGetLastOrder(ida_mem, &klast_adapt_modifiable);
      return std::make_tuple(r, klast_adapt_modifiable);
    };
    auto IDAGetLastOrder_adapt_solver_memory =
      [&IDAGetLastOrder_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, int>
    {
      auto lambda_result =
        IDAGetLastOrder_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetLastOrder_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetCurrentOrder",
  [](const IDAView& ida_mem) -> std::tuple<int, int>
  {
    auto IDAGetCurrentOrder_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, int>
    {
      int kcur_adapt_modifiable;

      int r = IDAGetCurrentOrder(ida_mem, &kcur_adapt_modifiable);
      return std::make_tuple(r, kcur_adapt_modifiable);
    };
    auto IDAGetCurrentOrder_adapt_solver_memory =
      [&IDAGetCurrentOrder_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, int>
    {
      auto lambda_result =
        IDAGetCurrentOrder_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetCurrentOrder_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetCurrentCj",
  [](const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetCurrentCj_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype cj_adapt_modifiable;

      int r = IDAGetCurrentCj(ida_mem, &cj_adapt_modifiable);
      return std::make_tuple(r, cj_adapt_modifiable);
    };
    auto IDAGetCurrentCj_adapt_solver_memory =
      [&IDAGetCurrentCj_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetCurrentCj_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetCurrentCj_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetCurrentY",
  [](const IDAView& ida_mem) -> std::tuple<int, N_Vector>
  {
    auto IDAGetCurrentY_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, N_Vector>
    {
      N_Vector ycur_adapt_modifiable;

      int r = IDAGetCurrentY(ida_mem, &ycur_adapt_modifiable);
      return std::make_tuple(r, ycur_adapt_modifiable);
    };
    auto IDAGetCurrentY_adapt_solver_memory =
      [&IDAGetCurrentY_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, N_Vector>
    {
      auto lambda_result =
        IDAGetCurrentY_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetCurrentY_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "IDAGetCurrentYp",
  [](const IDAView& ida_mem) -> std::tuple<int, N_Vector>
  {
    auto IDAGetCurrentYp_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, N_Vector>
    {
      N_Vector ypcur_adapt_modifiable;

      int r = IDAGetCurrentYp(ida_mem, &ypcur_adapt_modifiable);
      return std::make_tuple(r, ypcur_adapt_modifiable);
    };
    auto IDAGetCurrentYp_adapt_solver_memory =
      [&IDAGetCurrentYp_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, N_Vector>
    {
      auto lambda_result =
        IDAGetCurrentYp_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetCurrentYp_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "IDAGetActualInitStep",
  [](const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetActualInitStep_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hinused_adapt_modifiable;

      int r = IDAGetActualInitStep(ida_mem, &hinused_adapt_modifiable);
      return std::make_tuple(r, hinused_adapt_modifiable);
    };
    auto IDAGetActualInitStep_adapt_solver_memory =
      [&IDAGetActualInitStep_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetActualInitStep_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetActualInitStep_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetLastStep",
  [](const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetLastStep_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hlast_adapt_modifiable;

      int r = IDAGetLastStep(ida_mem, &hlast_adapt_modifiable);
      return std::make_tuple(r, hlast_adapt_modifiable);
    };
    auto IDAGetLastStep_adapt_solver_memory =
      [&IDAGetLastStep_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetLastStep_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetLastStep_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetCurrentStep",
  [](const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetCurrentStep_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hcur_adapt_modifiable;

      int r = IDAGetCurrentStep(ida_mem, &hcur_adapt_modifiable);
      return std::make_tuple(r, hcur_adapt_modifiable);
    };
    auto IDAGetCurrentStep_adapt_solver_memory =
      [&IDAGetCurrentStep_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetCurrentStep_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetCurrentStep_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetCurrentTime",
  [](const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetCurrentTime_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tcur_adapt_modifiable;

      int r = IDAGetCurrentTime(ida_mem, &tcur_adapt_modifiable);
      return std::make_tuple(r, tcur_adapt_modifiable);
    };
    auto IDAGetCurrentTime_adapt_solver_memory =
      [&IDAGetCurrentTime_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetCurrentTime_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetCurrentTime_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetTolScaleFactor",
  [](const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetTolScaleFactor_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tolsfact_adapt_modifiable;

      int r = IDAGetTolScaleFactor(ida_mem, &tolsfact_adapt_modifiable);
      return std::make_tuple(r, tolsfact_adapt_modifiable);
    };
    auto IDAGetTolScaleFactor_adapt_solver_memory =
      [&IDAGetTolScaleFactor_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetTolScaleFactor_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetTolScaleFactor_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetErrWeights",
  [](const IDAView& ida_mem, N_Vector eweight) -> int
  {
    auto IDAGetErrWeights_adapt_solver_memory = [](const IDAView& ida_mem,
                                                   N_Vector eweight) -> int
    {
      auto lambda_result = IDAGetErrWeights(ida_mem.get(), eweight);
      return lambda_result;
    };

    return IDAGetErrWeights_adapt_solver_memory(ida_mem, eweight);
  },
  nb::arg("ida_mem"), nb::arg("eweight"));

sundials4py::scoped_def(
  m, "IDAGetEstLocalErrors",
  [](const IDAView& ida_mem, N_Vector ele) -> int
  {
    auto IDAGetEstLocalErrors_adapt_solver_memory = [](const IDAView& ida_mem,
                                                       N_Vector ele) -> int
    {
      auto lambda_result = IDAGetEstLocalErrors(ida_mem.get(), ele);
      return lambda_result;
    };

    return IDAGetEstLocalErrors_adapt_solver_memory(ida_mem, ele);
  },
  nb::arg("ida_mem"), nb::arg("ele"));

sundials4py::scoped_def(
  m, "IDAGetNumConstraintFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumConstraintFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long num_fails_out_adapt_modifiable;

      int r = IDAGetNumConstraintFails(ida_mem, &num_fails_out_adapt_modifiable);
      return std::make_tuple(r, num_fails_out_adapt_modifiable);
    };
    auto IDAGetNumConstraintFails_adapt_solver_memory =
      [&IDAGetNumConstraintFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumConstraintFails_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumConstraintFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumConstraintCorrections",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumConstraintCorrections_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long num_corrections_out_adapt_modifiable;

      int r =
        IDAGetNumConstraintCorrections(ida_mem,
                                       &num_corrections_out_adapt_modifiable);
      return std::make_tuple(r, num_corrections_out_adapt_modifiable);
    };
    auto IDAGetNumConstraintCorrections_adapt_solver_memory =
      [&IDAGetNumConstraintCorrections_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumConstraintCorrections_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumConstraintCorrections_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumGEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumGEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long ngevals_adapt_modifiable;

      int r = IDAGetNumGEvals(ida_mem, &ngevals_adapt_modifiable);
      return std::make_tuple(r, ngevals_adapt_modifiable);
    };
    auto IDAGetNumGEvals_adapt_solver_memory =
      [&IDAGetNumGEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumGEvals_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumGEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetRootInfo",
  [](const IDAView& ida_mem, sundials4py::IntArray1d rootsfound_1d) -> int
  {
    auto IDAGetRootInfo_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sundials4py::IntArray1d rootsfound_1d) -> int
    {
      int* rootsfound_1d_ptr = rootsfound_1d.size() == 0 ? nullptr
                                                         : rootsfound_1d.data();

      auto lambda_result = IDAGetRootInfo(ida_mem, rootsfound_1d_ptr);
      return lambda_result;
    };
    auto IDAGetRootInfo_adapt_solver_memory =
      [&IDAGetRootInfo_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                    sundials4py::IntArray1d
                                                      rootsfound_1d) -> int
    {
      auto lambda_result =
        IDAGetRootInfo_adapt_arr_ptr_to_std_vector(ida_mem.get(), rootsfound_1d);
      return lambda_result;
    };

    return IDAGetRootInfo_adapt_solver_memory(ida_mem, rootsfound_1d);
  },
  nb::arg("ida_mem"), nb::arg("rootsfound_1d"));

sundials4py::scoped_def(
  m, "IDAGetIntegratorStats",
  [](const IDAView& ida_mem)
    -> std::tuple<int, long, long, long, long, int, int, sunrealtype,
                  sunrealtype, sunrealtype, sunrealtype>
  {
    auto IDAGetIntegratorStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem)
      -> std::tuple<int, long, long, long, long, int, int, sunrealtype,
                    sunrealtype, sunrealtype, sunrealtype>
    {
      long nsteps_adapt_modifiable;
      long nrevals_adapt_modifiable;
      long nlinsetups_adapt_modifiable;
      long netfails_adapt_modifiable;
      int qlast_adapt_modifiable;
      int qcur_adapt_modifiable;
      sunrealtype hinused_adapt_modifiable;
      sunrealtype hlast_adapt_modifiable;
      sunrealtype hcur_adapt_modifiable;
      sunrealtype tcur_adapt_modifiable;

      int r =
        IDAGetIntegratorStats(ida_mem, &nsteps_adapt_modifiable,
                              &nrevals_adapt_modifiable,
                              &nlinsetups_adapt_modifiable,
                              &netfails_adapt_modifiable,
                              &qlast_adapt_modifiable, &qcur_adapt_modifiable,
                              &hinused_adapt_modifiable, &hlast_adapt_modifiable,
                              &hcur_adapt_modifiable, &tcur_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable, nrevals_adapt_modifiable,
                             nlinsetups_adapt_modifiable,
                             netfails_adapt_modifiable, qlast_adapt_modifiable,
                             qcur_adapt_modifiable, hinused_adapt_modifiable,
                             hlast_adapt_modifiable, hcur_adapt_modifiable,
                             tcur_adapt_modifiable);
    };
    auto IDAGetIntegratorStats_adapt_solver_memory =
      [&IDAGetIntegratorStats_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem)
      -> std::tuple<int, long, long, long, long, int, int, sunrealtype,
                    sunrealtype, sunrealtype, sunrealtype>
    {
      auto lambda_result =
        IDAGetIntegratorStats_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetIntegratorStats_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumNonlinSolvIters",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nniters_adapt_modifiable;

      int r = IDAGetNumNonlinSolvIters(ida_mem, &nniters_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable);
    };
    auto IDAGetNumNonlinSolvIters_adapt_solver_memory =
      [&IDAGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumNonlinSolvIters_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumNonlinSolvConvFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nnfails_adapt_modifiable;

      int r = IDAGetNumNonlinSolvConvFails(ida_mem, &nnfails_adapt_modifiable);
      return std::make_tuple(r, nnfails_adapt_modifiable);
    };
    auto IDAGetNumNonlinSolvConvFails_adapt_solver_memory =
      [&IDAGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumNonlinSolvConvFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNonlinSolvStats",
  [](const IDAView& ida_mem) -> std::tuple<int, long, long>
  {
    auto IDAGetNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long>
    {
      long nniters_adapt_modifiable;
      long nnfails_adapt_modifiable;

      int r = IDAGetNonlinSolvStats(ida_mem, &nniters_adapt_modifiable,
                                    &nnfails_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable,
                             nnfails_adapt_modifiable);
    };
    auto IDAGetNonlinSolvStats_adapt_solver_memory =
      [&IDAGetNonlinSolvStats_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long, long>
    {
      auto lambda_result =
        IDAGetNonlinSolvStats_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNonlinSolvStats_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumStepSolveFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumStepSolveFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nncfails_adapt_modifiable;

      int r = IDAGetNumStepSolveFails(ida_mem, &nncfails_adapt_modifiable);
      return std::make_tuple(r, nncfails_adapt_modifiable);
    };
    auto IDAGetNumStepSolveFails_adapt_solver_memory =
      [&IDAGetNumStepSolveFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumStepSolveFails_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumStepSolveFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAPrintAllStats",
  [](const IDAView& ida_mem, FILE* outfile, SUNOutputFormat fmt) -> int
  {
    auto IDAPrintAllStats_adapt_solver_memory =
      [](const IDAView& ida_mem, FILE* outfile, SUNOutputFormat fmt) -> int
    {
      auto lambda_result = IDAPrintAllStats(ida_mem.get(), outfile, fmt);
      return lambda_result;
    };

    return IDAPrintAllStats_adapt_solver_memory(ida_mem, outfile, fmt);
  },
  nb::arg("ida_mem"), nb::arg("outfile"), nb::arg("fmt"));

sundials4py::scoped_def(m, "IDAGetReturnFlagName", IDAGetReturnFlagName,
                        nb::arg("flag"));

sundials4py::scoped_def(
  m, "IDAQuadReInit",
  [](const IDAView& ida_mem, N_Vector yQ0) -> int
  {
    auto IDAQuadReInit_adapt_solver_memory = [](const IDAView& ida_mem,
                                                N_Vector yQ0) -> int
    {
      auto lambda_result = IDAQuadReInit(ida_mem.get(), yQ0);
      return lambda_result;
    };

    return IDAQuadReInit_adapt_solver_memory(ida_mem, yQ0);
  },
  nb::arg("ida_mem"), nb::arg("yQ0"));

sundials4py::scoped_def(
  m, "IDAQuadSStolerances",
  [](const IDAView& ida_mem, sunrealtype reltolQ, sunrealtype abstolQ) -> int
  {
    auto IDAQuadSStolerances_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype reltolQ, sunrealtype abstolQ) -> int
    {
      auto lambda_result = IDAQuadSStolerances(ida_mem.get(), reltolQ, abstolQ);
      return lambda_result;
    };

    return IDAQuadSStolerances_adapt_solver_memory(ida_mem, reltolQ, abstolQ);
  },
  nb::arg("ida_mem"), nb::arg("reltolQ"), nb::arg("abstolQ"));

sundials4py::scoped_def(
  m, "IDAQuadSVtolerances",
  [](const IDAView& ida_mem, sunrealtype reltolQ, N_Vector abstolQ) -> int
  {
    auto IDAQuadSVtolerances_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype reltolQ, N_Vector abstolQ) -> int
    {
      auto lambda_result = IDAQuadSVtolerances(ida_mem.get(), reltolQ, abstolQ);
      return lambda_result;
    };

    return IDAQuadSVtolerances_adapt_solver_memory(ida_mem, reltolQ, abstolQ);
  },
  nb::arg("ida_mem"), nb::arg("reltolQ"), nb::arg("abstolQ"));

sundials4py::scoped_def(
  m, "IDASetQuadErrCon",
  [](const IDAView& ida_mem, sunbooleantype errconQ) -> int
  {
    auto IDASetQuadErrCon_adapt_solver_memory = [](const IDAView& ida_mem,
                                                   sunbooleantype errconQ) -> int
    {
      auto lambda_result = IDASetQuadErrCon(ida_mem.get(), errconQ);
      return lambda_result;
    };

    return IDASetQuadErrCon_adapt_solver_memory(ida_mem, errconQ);
  },
  nb::arg("ida_mem"), nb::arg("errconQ"));

sundials4py::scoped_def(
  m, "IDAGetQuad",
  [](const IDAView& ida_mem, N_Vector yQout) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetQuad_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, N_Vector yQout) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_adapt_modifiable;

      int r = IDAGetQuad(ida_mem, &t_adapt_modifiable, yQout);
      return std::make_tuple(r, t_adapt_modifiable);
    };
    auto IDAGetQuad_adapt_solver_memory =
      [&IDAGetQuad_adapt_modifiable_immutable_to_return](const IDAView& ida_mem,
                                                         N_Vector yQout)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetQuad_adapt_modifiable_immutable_to_return(ida_mem.get(), yQout);
      return lambda_result;
    };

    return IDAGetQuad_adapt_solver_memory(ida_mem, yQout);
  },
  nb::arg("ida_mem"), nb::arg("yQout"));

sundials4py::scoped_def(
  m, "IDAGetQuadDky",
  [](const IDAView& ida_mem, sunrealtype t, int k, N_Vector dky) -> int
  {
    auto IDAGetQuadDky_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype t, int k, N_Vector dky) -> int
    {
      auto lambda_result = IDAGetQuadDky(ida_mem.get(), t, k, dky);
      return lambda_result;
    };

    return IDAGetQuadDky_adapt_solver_memory(ida_mem, t, k, dky);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dky"));

sundials4py::scoped_def(
  m, "IDAGetQuadNumRhsEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nrhsQevals_adapt_modifiable;

      int r = IDAGetQuadNumRhsEvals(ida_mem, &nrhsQevals_adapt_modifiable);
      return std::make_tuple(r, nrhsQevals_adapt_modifiable);
    };
    auto IDAGetQuadNumRhsEvals_adapt_solver_memory =
      [&IDAGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetQuadNumRhsEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetQuadNumErrTestFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nQetfails_adapt_modifiable;

      int r = IDAGetQuadNumErrTestFails(ida_mem, &nQetfails_adapt_modifiable);
      return std::make_tuple(r, nQetfails_adapt_modifiable);
    };
    auto IDAGetQuadNumErrTestFails_adapt_solver_memory =
      [&IDAGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetQuadNumErrTestFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetQuadErrWeights",
  [](const IDAView& ida_mem, N_Vector eQweight) -> int
  {
    auto IDAGetQuadErrWeights_adapt_solver_memory = [](const IDAView& ida_mem,
                                                       N_Vector eQweight) -> int
    {
      auto lambda_result = IDAGetQuadErrWeights(ida_mem.get(), eQweight);
      return lambda_result;
    };

    return IDAGetQuadErrWeights_adapt_solver_memory(ida_mem, eQweight);
  },
  nb::arg("ida_mem"), nb::arg("eQweight"));

sundials4py::scoped_def(
  m, "IDAGetQuadStats",
  [](const IDAView& ida_mem) -> std::tuple<int, long, long>
  {
    auto IDAGetQuadStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long>
    {
      long nrhsQevals_adapt_modifiable;
      long nQetfails_adapt_modifiable;

      int r = IDAGetQuadStats(ida_mem, &nrhsQevals_adapt_modifiable,
                              &nQetfails_adapt_modifiable);
      return std::make_tuple(r, nrhsQevals_adapt_modifiable,
                             nQetfails_adapt_modifiable);
    };
    auto IDAGetQuadStats_adapt_solver_memory =
      [&IDAGetQuadStats_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long, long>
    {
      auto lambda_result =
        IDAGetQuadStats_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetQuadStats_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDASensReInit",
  [](const IDAView& ida_mem, int ism, std::vector<N_Vector> yS0_1d,
     std::vector<N_Vector> ypS0_1d) -> int
  {
    auto IDASensReInit_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, int ism, std::vector<N_Vector> yS0_1d,
         std::vector<N_Vector> ypS0_1d) -> int
    {
      N_Vector* yS0_1d_ptr  = yS0_1d.empty() ? nullptr : yS0_1d.data();
      N_Vector* ypS0_1d_ptr = ypS0_1d.empty() ? nullptr : ypS0_1d.data();

      auto lambda_result = IDASensReInit(ida_mem, ism, yS0_1d_ptr, ypS0_1d_ptr);
      return lambda_result;
    };
    auto IDASensReInit_adapt_solver_memory =
      [&IDASensReInit_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                   int ism,
                                                   std::vector<N_Vector> yS0_1d,
                                                   std::vector<N_Vector> ypS0_1d) -> int
    {
      auto lambda_result =
        IDASensReInit_adapt_arr_ptr_to_std_vector(ida_mem.get(), ism, yS0_1d,
                                                  ypS0_1d);
      return lambda_result;
    };

    return IDASensReInit_adapt_solver_memory(ida_mem, ism, yS0_1d, ypS0_1d);
  },
  nb::arg("ida_mem"), nb::arg("ism"), nb::arg("yS0_1d"), nb::arg("ypS0_1d"));

sundials4py::scoped_def(
  m, "IDASensSStolerances",
  [](const IDAView& ida_mem, sunrealtype reltolS,
     sundials4py::Array1d abstolS_1d) -> int
  {
    auto IDASensSStolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype reltolS, sundials4py::Array1d abstolS_1d) -> int
    {
      sunrealtype* abstolS_1d_ptr = abstolS_1d.size() == 0 ? nullptr
                                                           : abstolS_1d.data();

      auto lambda_result = IDASensSStolerances(ida_mem, reltolS, abstolS_1d_ptr);
      return lambda_result;
    };
    auto IDASensSStolerances_adapt_solver_memory =
      [&IDASensSStolerances_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                         sunrealtype reltolS,
                                                         sundials4py::Array1d
                                                           abstolS_1d) -> int
    {
      auto lambda_result =
        IDASensSStolerances_adapt_arr_ptr_to_std_vector(ida_mem.get(), reltolS,
                                                        abstolS_1d);
      return lambda_result;
    };

    return IDASensSStolerances_adapt_solver_memory(ida_mem, reltolS, abstolS_1d);
  },
  nb::arg("ida_mem"), nb::arg("reltolS"), nb::arg("abstolS_1d"));

sundials4py::scoped_def(
  m, "IDASensSVtolerances",
  [](const IDAView& ida_mem, sunrealtype reltolS,
     std::vector<N_Vector> abstolS_1d) -> int
  {
    auto IDASensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype reltolS,
         std::vector<N_Vector> abstolS_1d) -> int
    {
      N_Vector* abstolS_1d_ptr = abstolS_1d.empty() ? nullptr : abstolS_1d.data();

      auto lambda_result = IDASensSVtolerances(ida_mem, reltolS, abstolS_1d_ptr);
      return lambda_result;
    };
    auto IDASensSVtolerances_adapt_solver_memory =
      [&IDASensSVtolerances_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                         sunrealtype reltolS,
                                                         std::vector<N_Vector>
                                                           abstolS_1d) -> int
    {
      auto lambda_result =
        IDASensSVtolerances_adapt_arr_ptr_to_std_vector(ida_mem.get(), reltolS,
                                                        abstolS_1d);
      return lambda_result;
    };

    return IDASensSVtolerances_adapt_solver_memory(ida_mem, reltolS, abstolS_1d);
  },
  nb::arg("ida_mem"), nb::arg("reltolS"), nb::arg("abstolS_1d"));

sundials4py::scoped_def(
  m, "IDASensEEtolerances",
  [](const IDAView& ida_mem) -> int
  {
    auto IDASensEEtolerances_adapt_solver_memory = [](const IDAView& ida_mem) -> int
    {
      auto lambda_result = IDASensEEtolerances(ida_mem.get());
      return lambda_result;
    };

    return IDASensEEtolerances_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetSensConsistentIC",
  [](const IDAView& ida_mem, std::vector<N_Vector> yyS0_1d,
     std::vector<N_Vector> ypS0_1d) -> int
  {
    auto IDAGetSensConsistentIC_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> yyS0_1d,
         std::vector<N_Vector> ypS0_1d) -> int
    {
      N_Vector* yyS0_1d_ptr = yyS0_1d.empty() ? nullptr : yyS0_1d.data();
      N_Vector* ypS0_1d_ptr = ypS0_1d.empty() ? nullptr : ypS0_1d.data();

      auto lambda_result = IDAGetSensConsistentIC(ida_mem, yyS0_1d_ptr,
                                                  ypS0_1d_ptr);
      return lambda_result;
    };
    auto IDAGetSensConsistentIC_adapt_solver_memory =
      [&IDAGetSensConsistentIC_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                            std::vector<N_Vector> yyS0_1d,
                                                            std::vector<N_Vector>
                                                              ypS0_1d) -> int
    {
      auto lambda_result =
        IDAGetSensConsistentIC_adapt_arr_ptr_to_std_vector(ida_mem.get(),
                                                           yyS0_1d, ypS0_1d);
      return lambda_result;
    };

    return IDAGetSensConsistentIC_adapt_solver_memory(ida_mem, yyS0_1d, ypS0_1d);
  },
  nb::arg("ida_mem"), nb::arg("yyS0_1d"), nb::arg("ypS0_1d"));

sundials4py::scoped_def(
  m, "IDASetSensDQMethod",
  [](const IDAView& ida_mem, int DQtype, sunrealtype DQrhomax) -> int
  {
    auto IDASetSensDQMethod_adapt_solver_memory =
      [](const IDAView& ida_mem, int DQtype, sunrealtype DQrhomax) -> int
    {
      auto lambda_result = IDASetSensDQMethod(ida_mem.get(), DQtype, DQrhomax);
      return lambda_result;
    };

    return IDASetSensDQMethod_adapt_solver_memory(ida_mem, DQtype, DQrhomax);
  },
  nb::arg("ida_mem"), nb::arg("DQtype"), nb::arg("DQrhomax"));

sundials4py::scoped_def(
  m, "IDASetSensErrCon",
  [](const IDAView& ida_mem, sunbooleantype errconS) -> int
  {
    auto IDASetSensErrCon_adapt_solver_memory = [](const IDAView& ida_mem,
                                                   sunbooleantype errconS) -> int
    {
      auto lambda_result = IDASetSensErrCon(ida_mem.get(), errconS);
      return lambda_result;
    };

    return IDASetSensErrCon_adapt_solver_memory(ida_mem, errconS);
  },
  nb::arg("ida_mem"), nb::arg("errconS"));

sundials4py::scoped_def(
  m, "IDASetSensMaxNonlinIters",
  [](const IDAView& ida_mem, int maxcorS) -> int
  {
    auto IDASetSensMaxNonlinIters_adapt_solver_memory =
      [](const IDAView& ida_mem, int maxcorS) -> int
    {
      auto lambda_result = IDASetSensMaxNonlinIters(ida_mem.get(), maxcorS);
      return lambda_result;
    };

    return IDASetSensMaxNonlinIters_adapt_solver_memory(ida_mem, maxcorS);
  },
  nb::arg("ida_mem"), nb::arg("maxcorS"));

sundials4py::scoped_def(
  m, "IDASetSensParams",
  [](const IDAView& ida_mem, sundials4py::Array1d p_1d,
     sundials4py::Array1d pbar_1d, std::vector<int> plist_1d) -> int
  {
    auto IDASetSensParams_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sundials4py::Array1d p_1d, sundials4py::Array1d pbar_1d,
         std::vector<int> plist_1d) -> int
    {
      sunrealtype* p_1d_ptr    = p_1d.size() == 0 ? nullptr : p_1d.data();
      sunrealtype* pbar_1d_ptr = pbar_1d.size() == 0 ? nullptr : pbar_1d.data();
      int* plist_1d_ptr        = plist_1d.empty() ? nullptr : plist_1d.data();

      auto lambda_result = IDASetSensParams(ida_mem, p_1d_ptr, pbar_1d_ptr,
                                            plist_1d_ptr);
      return lambda_result;
    };
    auto IDASetSensParams_adapt_solver_memory =
      [&IDASetSensParams_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                      sundials4py::Array1d p_1d,
                                                      sundials4py::Array1d pbar_1d,
                                                      std::vector<int> plist_1d) -> int
    {
      auto lambda_result =
        IDASetSensParams_adapt_arr_ptr_to_std_vector(ida_mem.get(), p_1d,
                                                     pbar_1d, plist_1d);
      return lambda_result;
    };

    return IDASetSensParams_adapt_solver_memory(ida_mem, p_1d, pbar_1d, plist_1d);
  },
  nb::arg("ida_mem"), nb::arg("p_1d"), nb::arg("pbar_1d"), nb::arg("plist_1d"));

sundials4py::scoped_def(
  m, "IDASetNonlinearSolverSensSim",
  [](const IDAView& ida_mem, SUNNonlinearSolver NLS) -> int
  {
    auto IDASetNonlinearSolverSensSim_adapt_solver_memory =
      [](const IDAView& ida_mem, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = IDASetNonlinearSolverSensSim(ida_mem.get(), NLS);
      return lambda_result;
    };

    return IDASetNonlinearSolverSensSim_adapt_solver_memory(ida_mem, NLS);
  },
  nb::arg("ida_mem"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "IDASetNonlinearSolverSensStg",
  [](const IDAView& ida_mem, SUNNonlinearSolver NLS) -> int
  {
    auto IDASetNonlinearSolverSensStg_adapt_solver_memory =
      [](const IDAView& ida_mem, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = IDASetNonlinearSolverSensStg(ida_mem.get(), NLS);
      return lambda_result;
    };

    return IDASetNonlinearSolverSensStg_adapt_solver_memory(ida_mem, NLS);
  },
  nb::arg("ida_mem"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "IDASensToggleOff",
  [](const IDAView& ida_mem) -> int
  {
    auto IDASensToggleOff_adapt_solver_memory = [](const IDAView& ida_mem) -> int
    {
      auto lambda_result = IDASensToggleOff(ida_mem.get());
      return lambda_result;
    };

    return IDASensToggleOff_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetSens",
  [](const IDAView& ida_mem,
     std::vector<N_Vector> yySout_1d) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype* tret, std::vector<N_Vector> yySout_1d) -> int
    {
      N_Vector* yySout_1d_ptr = yySout_1d.empty() ? nullptr : yySout_1d.data();

      auto lambda_result = IDAGetSens(ida_mem, tret, yySout_1d_ptr);
      return lambda_result;
    };
    auto IDAGetSens_adapt_modifiable_immutable_to_return =
      [&IDAGetSens_adapt_arr_ptr_to_std_vector](void* ida_mem,
                                                std::vector<N_Vector> yySout_1d)
      -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetSens_adapt_arr_ptr_to_std_vector(ida_mem,
                                                     &tret_adapt_modifiable,
                                                     yySout_1d);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto IDAGetSens_adapt_solver_memory =
      [&IDAGetSens_adapt_modifiable_immutable_to_return](const IDAView& ida_mem,
                                                         std::vector<N_Vector> yySout_1d)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetSens_adapt_modifiable_immutable_to_return(ida_mem.get(), yySout_1d);
      return lambda_result;
    };

    return IDAGetSens_adapt_solver_memory(ida_mem, yySout_1d);
  },
  nb::arg("ida_mem"), nb::arg("yySout_1d"));

sundials4py::scoped_def(
  m, "IDAGetSens1",
  [](const IDAView& ida_mem, int is,
     N_Vector yySret) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetSens1_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int is, N_Vector yySret) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetSens1(ida_mem, &tret_adapt_modifiable, is, yySret);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto IDAGetSens1_adapt_solver_memory =
      [&IDAGetSens1_adapt_modifiable_immutable_to_return](const IDAView& ida_mem,
                                                          int is, N_Vector yySret)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetSens1_adapt_modifiable_immutable_to_return(ida_mem.get(), is,
                                                         yySret);
      return lambda_result;
    };

    return IDAGetSens1_adapt_solver_memory(ida_mem, is, yySret);
  },
  nb::arg("ida_mem"), nb::arg("is_"), nb::arg("yySret"));

sundials4py::scoped_def(
  m, "IDAGetSensDky",
  [](const IDAView& ida_mem, sunrealtype t, int k,
     std::vector<N_Vector> dkyS_1d) -> int
  {
    auto IDAGetSensDky_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype t, int k, std::vector<N_Vector> dkyS_1d) -> int
    {
      N_Vector* dkyS_1d_ptr = dkyS_1d.empty() ? nullptr : dkyS_1d.data();

      auto lambda_result = IDAGetSensDky(ida_mem, t, k, dkyS_1d_ptr);
      return lambda_result;
    };
    auto IDAGetSensDky_adapt_solver_memory =
      [&IDAGetSensDky_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                   sunrealtype t, int k,
                                                   std::vector<N_Vector> dkyS_1d) -> int
    {
      auto lambda_result =
        IDAGetSensDky_adapt_arr_ptr_to_std_vector(ida_mem.get(), t, k, dkyS_1d);
      return lambda_result;
    };

    return IDAGetSensDky_adapt_solver_memory(ida_mem, t, k, dkyS_1d);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyS_1d"));

sundials4py::scoped_def(
  m, "IDAGetSensDky1",
  [](const IDAView& ida_mem, sunrealtype t, int k, int is, N_Vector dkyS) -> int
  {
    auto IDAGetSensDky1_adapt_solver_memory = [](const IDAView& ida_mem,
                                                 sunrealtype t, int k, int is,
                                                 N_Vector dkyS) -> int
    {
      auto lambda_result = IDAGetSensDky1(ida_mem.get(), t, k, is, dkyS);
      return lambda_result;
    };

    return IDAGetSensDky1_adapt_solver_memory(ida_mem, t, k, is, dkyS);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("is_"),
  nb::arg("dkyS"));

sundials4py::scoped_def(
  m, "IDAGetSensNumResEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumResEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nresSevals_adapt_modifiable;

      int r = IDAGetSensNumResEvals(ida_mem, &nresSevals_adapt_modifiable);
      return std::make_tuple(r, nresSevals_adapt_modifiable);
    };
    auto IDAGetSensNumResEvals_adapt_solver_memory =
      [&IDAGetSensNumResEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetSensNumResEvals_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetSensNumResEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumResEvalsSens",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumResEvalsSens_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nresevalsS_adapt_modifiable;

      int r = IDAGetNumResEvalsSens(ida_mem, &nresevalsS_adapt_modifiable);
      return std::make_tuple(r, nresevalsS_adapt_modifiable);
    };
    auto IDAGetNumResEvalsSens_adapt_solver_memory =
      [&IDAGetNumResEvalsSens_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumResEvalsSens_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumResEvalsSens_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetSensNumErrTestFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nSetfails_adapt_modifiable;

      int r = IDAGetSensNumErrTestFails(ida_mem, &nSetfails_adapt_modifiable);
      return std::make_tuple(r, nSetfails_adapt_modifiable);
    };
    auto IDAGetSensNumErrTestFails_adapt_solver_memory =
      [&IDAGetSensNumErrTestFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetSensNumErrTestFails_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetSensNumErrTestFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetSensNumLinSolvSetups",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nlinsetupsS_adapt_modifiable;

      int r = IDAGetSensNumLinSolvSetups(ida_mem, &nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nlinsetupsS_adapt_modifiable);
    };
    auto IDAGetSensNumLinSolvSetups_adapt_solver_memory =
      [&IDAGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetSensNumLinSolvSetups_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetSensErrWeights",
  [](const IDAView& ida_mem, std::vector<N_Vector> eSweight_1d) -> int
  {
    auto IDAGetSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> eSweight_1d) -> int
    {
      N_Vector* eSweight_1d_ptr = eSweight_1d.empty() ? nullptr
                                                      : eSweight_1d.data();

      auto lambda_result = IDAGetSensErrWeights(ida_mem, eSweight_1d_ptr);
      return lambda_result;
    };
    auto IDAGetSensErrWeights_adapt_solver_memory =
      [&IDAGetSensErrWeights_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                          std::vector<N_Vector>
                                                            eSweight_1d) -> int
    {
      auto lambda_result =
        IDAGetSensErrWeights_adapt_arr_ptr_to_std_vector(ida_mem.get(),
                                                         eSweight_1d);
      return lambda_result;
    };

    return IDAGetSensErrWeights_adapt_solver_memory(ida_mem, eSweight_1d);
  },
  nb::arg("ida_mem"), nb::arg("eSweight_1d"));

sundials4py::scoped_def(
  m, "IDAGetSensStats",
  [](const IDAView& ida_mem) -> std::tuple<int, long, long, long, long>
  {
    auto IDAGetSensStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long, long, long>
    {
      long nresSevals_adapt_modifiable;
      long nresevalsS_adapt_modifiable;
      long nSetfails_adapt_modifiable;
      long nlinsetupsS_adapt_modifiable;

      int r = IDAGetSensStats(ida_mem, &nresSevals_adapt_modifiable,
                              &nresevalsS_adapt_modifiable,
                              &nSetfails_adapt_modifiable,
                              &nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nresSevals_adapt_modifiable,
                             nresevalsS_adapt_modifiable,
                             nSetfails_adapt_modifiable,
                             nlinsetupsS_adapt_modifiable);
    };
    auto IDAGetSensStats_adapt_solver_memory =
      [&IDAGetSensStats_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long, long, long, long>
    {
      auto lambda_result =
        IDAGetSensStats_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetSensStats_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetSensNumNonlinSolvIters",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nSniters_adapt_modifiable;

      int r = IDAGetSensNumNonlinSolvIters(ida_mem, &nSniters_adapt_modifiable);
      return std::make_tuple(r, nSniters_adapt_modifiable);
    };
    auto IDAGetSensNumNonlinSolvIters_adapt_solver_memory =
      [&IDAGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetSensNumNonlinSolvIters_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetSensNumNonlinSolvConvFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nSnfails_adapt_modifiable;

      int r = IDAGetSensNumNonlinSolvConvFails(ida_mem,
                                               &nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSnfails_adapt_modifiable);
    };
    auto IDAGetSensNumNonlinSolvConvFails_adapt_solver_memory =
      [&IDAGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetSensNumNonlinSolvConvFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetSensNonlinSolvStats",
  [](const IDAView& ida_mem) -> std::tuple<int, long, long>
  {
    auto IDAGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long>
    {
      long nSniters_adapt_modifiable;
      long nSnfails_adapt_modifiable;

      int r = IDAGetSensNonlinSolvStats(ida_mem, &nSniters_adapt_modifiable,
                                        &nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSniters_adapt_modifiable,
                             nSnfails_adapt_modifiable);
    };
    auto IDAGetSensNonlinSolvStats_adapt_solver_memory =
      [&IDAGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long, long>
    {
      auto lambda_result =
        IDAGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetSensNonlinSolvStats_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumStepSensSolveFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nSncfails_adapt_modifiable;

      int r = IDAGetNumStepSensSolveFails(ida_mem, &nSncfails_adapt_modifiable);
      return std::make_tuple(r, nSncfails_adapt_modifiable);
    };
    auto IDAGetNumStepSensSolveFails_adapt_solver_memory =
      [&IDAGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumStepSensSolveFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAQuadSensReInit",
  [](const IDAView& ida_mem, std::vector<N_Vector> yQS0_1d) -> int
  {
    auto IDAQuadSensReInit_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> yQS0_1d) -> int
    {
      N_Vector* yQS0_1d_ptr = yQS0_1d.empty() ? nullptr : yQS0_1d.data();

      auto lambda_result = IDAQuadSensReInit(ida_mem, yQS0_1d_ptr);
      return lambda_result;
    };
    auto IDAQuadSensReInit_adapt_solver_memory =
      [&IDAQuadSensReInit_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                       std::vector<N_Vector> yQS0_1d)
      -> int
    {
      auto lambda_result =
        IDAQuadSensReInit_adapt_arr_ptr_to_std_vector(ida_mem.get(), yQS0_1d);
      return lambda_result;
    };

    return IDAQuadSensReInit_adapt_solver_memory(ida_mem, yQS0_1d);
  },
  nb::arg("ida_mem"), nb::arg("yQS0_1d"));

sundials4py::scoped_def(
  m, "IDAQuadSensSStolerances",
  [](const IDAView& ida_mem, sunrealtype reltolQS,
     sundials4py::Array1d abstolQS_1d) -> int
  {
    auto IDAQuadSensSStolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype reltolQS,
         sundials4py::Array1d abstolQS_1d) -> int
    {
      sunrealtype* abstolQS_1d_ptr =
        abstolQS_1d.size() == 0 ? nullptr : abstolQS_1d.data();

      auto lambda_result = IDAQuadSensSStolerances(ida_mem, reltolQS,
                                                   abstolQS_1d_ptr);
      return lambda_result;
    };
    auto IDAQuadSensSStolerances_adapt_solver_memory =
      [&IDAQuadSensSStolerances_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                             sunrealtype reltolQS,
                                                             sundials4py::Array1d
                                                               abstolQS_1d) -> int
    {
      auto lambda_result =
        IDAQuadSensSStolerances_adapt_arr_ptr_to_std_vector(ida_mem.get(),
                                                            reltolQS,
                                                            abstolQS_1d);
      return lambda_result;
    };

    return IDAQuadSensSStolerances_adapt_solver_memory(ida_mem, reltolQS,
                                                       abstolQS_1d);
  },
  nb::arg("ida_mem"), nb::arg("reltolQS"), nb::arg("abstolQS_1d"));

sundials4py::scoped_def(
  m, "IDAQuadSensSVtolerances",
  [](const IDAView& ida_mem, sunrealtype reltolQS,
     std::vector<N_Vector> abstolQS_1d) -> int
  {
    auto IDAQuadSensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype reltolQS,
         std::vector<N_Vector> abstolQS_1d) -> int
    {
      N_Vector* abstolQS_1d_ptr = abstolQS_1d.empty() ? nullptr
                                                      : abstolQS_1d.data();

      auto lambda_result = IDAQuadSensSVtolerances(ida_mem, reltolQS,
                                                   abstolQS_1d_ptr);
      return lambda_result;
    };
    auto IDAQuadSensSVtolerances_adapt_solver_memory =
      [&IDAQuadSensSVtolerances_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                             sunrealtype reltolQS,
                                                             std::vector<N_Vector>
                                                               abstolQS_1d) -> int
    {
      auto lambda_result =
        IDAQuadSensSVtolerances_adapt_arr_ptr_to_std_vector(ida_mem.get(),
                                                            reltolQS,
                                                            abstolQS_1d);
      return lambda_result;
    };

    return IDAQuadSensSVtolerances_adapt_solver_memory(ida_mem, reltolQS,
                                                       abstolQS_1d);
  },
  nb::arg("ida_mem"), nb::arg("reltolQS"), nb::arg("abstolQS_1d"));

sundials4py::scoped_def(
  m, "IDAQuadSensEEtolerances",
  [](const IDAView& ida_mem) -> int
  {
    auto IDAQuadSensEEtolerances_adapt_solver_memory =
      [](const IDAView& ida_mem) -> int
    {
      auto lambda_result = IDAQuadSensEEtolerances(ida_mem.get());
      return lambda_result;
    };

    return IDAQuadSensEEtolerances_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDASetQuadSensErrCon",
  [](const IDAView& ida_mem, sunbooleantype errconQS) -> int
  {
    auto IDASetQuadSensErrCon_adapt_solver_memory =
      [](const IDAView& ida_mem, sunbooleantype errconQS) -> int
    {
      auto lambda_result = IDASetQuadSensErrCon(ida_mem.get(), errconQS);
      return lambda_result;
    };

    return IDASetQuadSensErrCon_adapt_solver_memory(ida_mem, errconQS);
  },
  nb::arg("ida_mem"), nb::arg("errconQS"));

sundials4py::scoped_def(
  m, "IDAGetQuadSens",
  [](const IDAView& ida_mem,
     std::vector<N_Vector> yyQSout_1d) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetQuadSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype* tret, std::vector<N_Vector> yyQSout_1d) -> int
    {
      N_Vector* yyQSout_1d_ptr = yyQSout_1d.empty() ? nullptr : yyQSout_1d.data();

      auto lambda_result = IDAGetQuadSens(ida_mem, tret, yyQSout_1d_ptr);
      return lambda_result;
    };
    auto IDAGetQuadSens_adapt_modifiable_immutable_to_return =
      [&IDAGetQuadSens_adapt_arr_ptr_to_std_vector](void* ida_mem,
                                                    std::vector<N_Vector> yyQSout_1d)
      -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetQuadSens_adapt_arr_ptr_to_std_vector(ida_mem,
                                                         &tret_adapt_modifiable,
                                                         yyQSout_1d);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto IDAGetQuadSens_adapt_solver_memory =
      [&IDAGetQuadSens_adapt_modifiable_immutable_to_return](const IDAView& ida_mem,
                                                             std::vector<N_Vector>
                                                               yyQSout_1d)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetQuadSens_adapt_modifiable_immutable_to_return(ida_mem.get(),
                                                            yyQSout_1d);
      return lambda_result;
    };

    return IDAGetQuadSens_adapt_solver_memory(ida_mem, yyQSout_1d);
  },
  nb::arg("ida_mem"), nb::arg("yyQSout_1d"));

sundials4py::scoped_def(
  m, "IDAGetQuadSens1",
  [](const IDAView& ida_mem, int is,
     N_Vector yyQSret) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetQuadSens1_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int is, N_Vector yyQSret) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetQuadSens1(ida_mem, &tret_adapt_modifiable, is, yyQSret);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto IDAGetQuadSens1_adapt_solver_memory =
      [&IDAGetQuadSens1_adapt_modifiable_immutable_to_return](const IDAView& ida_mem,
                                                              int is,
                                                              N_Vector yyQSret)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetQuadSens1_adapt_modifiable_immutable_to_return(ida_mem.get(), is,
                                                             yyQSret);
      return lambda_result;
    };

    return IDAGetQuadSens1_adapt_solver_memory(ida_mem, is, yyQSret);
  },
  nb::arg("ida_mem"), nb::arg("is_"), nb::arg("yyQSret"));

sundials4py::scoped_def(
  m, "IDAGetQuadSensDky",
  [](const IDAView& ida_mem, sunrealtype t, int k,
     std::vector<N_Vector> dkyQS_1d) -> int
  {
    auto IDAGetQuadSensDky_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype t, int k, std::vector<N_Vector> dkyQS_1d) -> int
    {
      N_Vector* dkyQS_1d_ptr = dkyQS_1d.empty() ? nullptr : dkyQS_1d.data();

      auto lambda_result = IDAGetQuadSensDky(ida_mem, t, k, dkyQS_1d_ptr);
      return lambda_result;
    };
    auto IDAGetQuadSensDky_adapt_solver_memory =
      [&IDAGetQuadSensDky_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                       sunrealtype t, int k,
                                                       std::vector<N_Vector>
                                                         dkyQS_1d) -> int
    {
      auto lambda_result =
        IDAGetQuadSensDky_adapt_arr_ptr_to_std_vector(ida_mem.get(), t, k,
                                                      dkyQS_1d);
      return lambda_result;
    };

    return IDAGetQuadSensDky_adapt_solver_memory(ida_mem, t, k, dkyQS_1d);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyQS_1d"));

sundials4py::scoped_def(
  m, "IDAGetQuadSensDky1",
  [](const IDAView& ida_mem, sunrealtype t, int k, int is, N_Vector dkyQS) -> int
  {
    auto IDAGetQuadSensDky1_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype t, int k, int is,
         N_Vector dkyQS) -> int
    {
      auto lambda_result = IDAGetQuadSensDky1(ida_mem.get(), t, k, is, dkyQS);
      return lambda_result;
    };

    return IDAGetQuadSensDky1_adapt_solver_memory(ida_mem, t, k, is, dkyQS);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("is_"),
  nb::arg("dkyQS"));

sundials4py::scoped_def(
  m, "IDAGetQuadSensNumRhsEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nrhsQSevals_adapt_modifiable;

      int r = IDAGetQuadSensNumRhsEvals(ida_mem, &nrhsQSevals_adapt_modifiable);
      return std::make_tuple(r, nrhsQSevals_adapt_modifiable);
    };
    auto IDAGetQuadSensNumRhsEvals_adapt_solver_memory =
      [&IDAGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetQuadSensNumRhsEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetQuadSensNumErrTestFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nQSetfails_adapt_modifiable;

      int r = IDAGetQuadSensNumErrTestFails(ida_mem,
                                            &nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nQSetfails_adapt_modifiable);
    };
    auto IDAGetQuadSensNumErrTestFails_adapt_solver_memory =
      [&IDAGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return(
          ida_mem.get());
      return lambda_result;
    };

    return IDAGetQuadSensNumErrTestFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetQuadSensErrWeights",
  [](const IDAView& ida_mem, std::vector<N_Vector> eQSweight_1d) -> int
  {
    auto IDAGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> eQSweight_1d) -> int
    {
      N_Vector* eQSweight_1d_ptr = eQSweight_1d.empty() ? nullptr
                                                        : eQSweight_1d.data();

      auto lambda_result = IDAGetQuadSensErrWeights(ida_mem, eQSweight_1d_ptr);
      return lambda_result;
    };
    auto IDAGetQuadSensErrWeights_adapt_solver_memory =
      [&IDAGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                              std::vector<N_Vector>
                                                                eQSweight_1d) -> int
    {
      auto lambda_result =
        IDAGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector(ida_mem.get(),
                                                             eQSweight_1d);
      return lambda_result;
    };

    return IDAGetQuadSensErrWeights_adapt_solver_memory(ida_mem, eQSweight_1d);
  },
  nb::arg("ida_mem"), nb::arg("eQSweight_1d"));

sundials4py::scoped_def(
  m, "IDAGetQuadSensStats",
  [](const IDAView& ida_mem) -> std::tuple<int, long, long>
  {
    auto IDAGetQuadSensStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long>
    {
      long nrhsQSevals_adapt_modifiable;
      long nQSetfails_adapt_modifiable;

      int r = IDAGetQuadSensStats(ida_mem, &nrhsQSevals_adapt_modifiable,
                                  &nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nrhsQSevals_adapt_modifiable,
                             nQSetfails_adapt_modifiable);
    };
    auto IDAGetQuadSensStats_adapt_solver_memory =
      [&IDAGetQuadSensStats_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long, long>
    {
      auto lambda_result =
        IDAGetQuadSensStats_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetQuadSensStats_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAAdjInit",
  [](const IDAView& ida_mem, long steps, int interp) -> int
  {
    auto IDAAdjInit_adapt_solver_memory = [](const IDAView& ida_mem, long steps,
                                             int interp) -> int
    {
      auto lambda_result = IDAAdjInit(ida_mem.get(), steps, interp);
      return lambda_result;
    };

    return IDAAdjInit_adapt_solver_memory(ida_mem, steps, interp);
  },
  nb::arg("ida_mem"), nb::arg("steps"), nb::arg("interp"));

sundials4py::scoped_def(
  m, "IDAAdjReInit",
  [](const IDAView& ida_mem) -> int
  {
    auto IDAAdjReInit_adapt_solver_memory = [](const IDAView& ida_mem) -> int
    {
      auto lambda_result = IDAAdjReInit(ida_mem.get());
      return lambda_result;
    };

    return IDAAdjReInit_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDACreateB",
  [](const IDAView& ida_mem) -> std::tuple<int, int>
  {
    auto IDACreateB_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, int>
    {
      int which_adapt_modifiable;

      int r = IDACreateB(ida_mem, &which_adapt_modifiable);
      return std::make_tuple(r, which_adapt_modifiable);
    };
    auto IDACreateB_adapt_solver_memory =
      [&IDACreateB_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, int>
    {
      auto lambda_result =
        IDACreateB_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDACreateB_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAReInitB",
  [](const IDAView& ida_mem, int which, sunrealtype tB0, N_Vector yyB0,
     N_Vector ypB0) -> int
  {
    auto IDAReInitB_adapt_solver_memory = [](const IDAView& ida_mem, int which,
                                             sunrealtype tB0, N_Vector yyB0,
                                             N_Vector ypB0) -> int
    {
      auto lambda_result = IDAReInitB(ida_mem.get(), which, tB0, yyB0, ypB0);
      return lambda_result;
    };

    return IDAReInitB_adapt_solver_memory(ida_mem, which, tB0, yyB0, ypB0);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("tB0"), nb::arg("yyB0"),
  nb::arg("ypB0"));

sundials4py::scoped_def(
  m, "IDASStolerancesB",
  [](const IDAView& ida_mem, int which, sunrealtype relTolB,
     sunrealtype absTolB) -> int
  {
    auto IDASStolerancesB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunrealtype relTolB,
         sunrealtype absTolB) -> int
    {
      auto lambda_result = IDASStolerancesB(ida_mem.get(), which, relTolB,
                                            absTolB);
      return lambda_result;
    };

    return IDASStolerancesB_adapt_solver_memory(ida_mem, which, relTolB, absTolB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("relTolB"), nb::arg("absTolB"));

sundials4py::scoped_def(
  m, "IDASVtolerancesB",
  [](const IDAView& ida_mem, int which, sunrealtype relTolB, N_Vector absTolB) -> int
  {
    auto IDASVtolerancesB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunrealtype relTolB,
         N_Vector absTolB) -> int
    {
      auto lambda_result = IDASVtolerancesB(ida_mem.get(), which, relTolB,
                                            absTolB);
      return lambda_result;
    };

    return IDASVtolerancesB_adapt_solver_memory(ida_mem, which, relTolB, absTolB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("relTolB"), nb::arg("absTolB"));

sundials4py::scoped_def(
  m, "IDAQuadReInitB",
  [](const IDAView& ida_mem, int which, N_Vector yQB0) -> int
  {
    auto IDAQuadReInitB_adapt_solver_memory = [](const IDAView& ida_mem,
                                                 int which, N_Vector yQB0) -> int
    {
      auto lambda_result = IDAQuadReInitB(ida_mem.get(), which, yQB0);
      return lambda_result;
    };

    return IDAQuadReInitB_adapt_solver_memory(ida_mem, which, yQB0);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("yQB0"));

sundials4py::scoped_def(
  m, "IDAQuadSStolerancesB",
  [](const IDAView& ida_mem, int which, sunrealtype reltolQB,
     sunrealtype abstolQB) -> int
  {
    auto IDAQuadSStolerancesB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunrealtype reltolQB,
         sunrealtype abstolQB) -> int
    {
      auto lambda_result = IDAQuadSStolerancesB(ida_mem.get(), which, reltolQB,
                                                abstolQB);
      return lambda_result;
    };

    return IDAQuadSStolerancesB_adapt_solver_memory(ida_mem, which, reltolQB,
                                                    abstolQB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

sundials4py::scoped_def(
  m, "IDAQuadSVtolerancesB",
  [](const IDAView& ida_mem, int which, sunrealtype reltolQB,
     N_Vector abstolQB) -> int
  {
    auto IDAQuadSVtolerancesB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunrealtype reltolQB,
         N_Vector abstolQB) -> int
    {
      auto lambda_result = IDAQuadSVtolerancesB(ida_mem.get(), which, reltolQB,
                                                abstolQB);
      return lambda_result;
    };

    return IDAQuadSVtolerancesB_adapt_solver_memory(ida_mem, which, reltolQB,
                                                    abstolQB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

sundials4py::scoped_def(
  m, "IDACalcICB",
  [](const IDAView& ida_mem, int which, sunrealtype tout1, N_Vector yy0,
     N_Vector yp0) -> int
  {
    auto IDACalcICB_adapt_solver_memory = [](const IDAView& ida_mem, int which,
                                             sunrealtype tout1, N_Vector yy0,
                                             N_Vector yp0) -> int
    {
      auto lambda_result = IDACalcICB(ida_mem.get(), which, tout1, yy0, yp0);
      return lambda_result;
    };

    return IDACalcICB_adapt_solver_memory(ida_mem, which, tout1, yy0, yp0);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("tout1"), nb::arg("yy0"),
  nb::arg("yp0"));

sundials4py::scoped_def(
  m, "IDACalcICBS",
  [](const IDAView& ida_mem, int which, sunrealtype tout1, N_Vector yy0,
     N_Vector yp0, std::vector<N_Vector> yyS0_1d,
     std::vector<N_Vector> ypS0_1d) -> int
  {
    auto IDACalcICBS_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, int which, sunrealtype tout1, N_Vector yy0, N_Vector yp0,
         std::vector<N_Vector> yyS0_1d, std::vector<N_Vector> ypS0_1d) -> int
    {
      N_Vector* yyS0_1d_ptr = yyS0_1d.empty() ? nullptr : yyS0_1d.data();
      N_Vector* ypS0_1d_ptr = ypS0_1d.empty() ? nullptr : ypS0_1d.data();

      auto lambda_result = IDACalcICBS(ida_mem, which, tout1, yy0, yp0,
                                       yyS0_1d_ptr, ypS0_1d_ptr);
      return lambda_result;
    };
    auto IDACalcICBS_adapt_solver_memory =
      [&IDACalcICBS_adapt_arr_ptr_to_std_vector](const IDAView& ida_mem,
                                                 int which, sunrealtype tout1,
                                                 N_Vector yy0, N_Vector yp0,
                                                 std::vector<N_Vector> yyS0_1d,
                                                 std::vector<N_Vector> ypS0_1d) -> int
    {
      auto lambda_result =
        IDACalcICBS_adapt_arr_ptr_to_std_vector(ida_mem.get(), which, tout1,
                                                yy0, yp0, yyS0_1d, ypS0_1d);
      return lambda_result;
    };

    return IDACalcICBS_adapt_solver_memory(ida_mem, which, tout1, yy0, yp0,
                                           yyS0_1d, ypS0_1d);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("tout1"), nb::arg("yy0"),
  nb::arg("yp0"), nb::arg("yyS0_1d"), nb::arg("ypS0_1d"));

sundials4py::scoped_def(
  m, "IDAAdjSetNoSensi",
  [](const IDAView& ida_mem) -> int
  {
    auto IDAAdjSetNoSensi_adapt_solver_memory = [](const IDAView& ida_mem) -> int
    {
      auto lambda_result = IDAAdjSetNoSensi(ida_mem.get());
      return lambda_result;
    };

    return IDAAdjSetNoSensi_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDASetMaxOrdB",
  [](const IDAView& ida_mem, int which, int maxordB) -> int
  {
    auto IDASetMaxOrdB_adapt_solver_memory = [](const IDAView& ida_mem,
                                                int which, int maxordB) -> int
    {
      auto lambda_result = IDASetMaxOrdB(ida_mem.get(), which, maxordB);
      return lambda_result;
    };

    return IDASetMaxOrdB_adapt_solver_memory(ida_mem, which, maxordB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("maxordB"));

sundials4py::scoped_def(
  m, "IDASetMaxNumStepsB",
  [](const IDAView& ida_mem, int which, long mxstepsB) -> int
  {
    auto IDASetMaxNumStepsB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, long mxstepsB) -> int
    {
      auto lambda_result = IDASetMaxNumStepsB(ida_mem.get(), which, mxstepsB);
      return lambda_result;
    };

    return IDASetMaxNumStepsB_adapt_solver_memory(ida_mem, which, mxstepsB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("mxstepsB"));

sundials4py::scoped_def(
  m, "IDASetInitStepB",
  [](const IDAView& ida_mem, int which, sunrealtype hinB) -> int
  {
    auto IDASetInitStepB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunrealtype hinB) -> int
    {
      auto lambda_result = IDASetInitStepB(ida_mem.get(), which, hinB);
      return lambda_result;
    };

    return IDASetInitStepB_adapt_solver_memory(ida_mem, which, hinB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("hinB"));

sundials4py::scoped_def(
  m, "IDASetMaxStepB",
  [](const IDAView& ida_mem, int which, sunrealtype hmaxB) -> int
  {
    auto IDASetMaxStepB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunrealtype hmaxB) -> int
    {
      auto lambda_result = IDASetMaxStepB(ida_mem.get(), which, hmaxB);
      return lambda_result;
    };

    return IDASetMaxStepB_adapt_solver_memory(ida_mem, which, hmaxB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("hmaxB"));

sundials4py::scoped_def(
  m, "IDASetSuppressAlgB",
  [](const IDAView& ida_mem, int which, sunbooleantype suppressalgB) -> int
  {
    auto IDASetSuppressAlgB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunbooleantype suppressalgB) -> int
    {
      auto lambda_result = IDASetSuppressAlgB(ida_mem.get(), which, suppressalgB);
      return lambda_result;
    };

    return IDASetSuppressAlgB_adapt_solver_memory(ida_mem, which, suppressalgB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("suppressalgB"));

sundials4py::scoped_def(
  m, "IDASetIdB",
  [](const IDAView& ida_mem, int which, N_Vector idB) -> int
  {
    auto IDASetIdB_adapt_solver_memory = [](const IDAView& ida_mem, int which,
                                            N_Vector idB) -> int
    {
      auto lambda_result = IDASetIdB(ida_mem.get(), which, idB);
      return lambda_result;
    };

    return IDASetIdB_adapt_solver_memory(ida_mem, which, idB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("idB"));

sundials4py::scoped_def(
  m, "IDASetConstraintsB",
  [](const IDAView& ida_mem, int which, N_Vector constraintsB) -> int
  {
    auto IDASetConstraintsB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, N_Vector constraintsB) -> int
    {
      auto lambda_result = IDASetConstraintsB(ida_mem.get(), which, constraintsB);
      return lambda_result;
    };

    return IDASetConstraintsB_adapt_solver_memory(ida_mem, which, constraintsB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("constraintsB"));

sundials4py::scoped_def(
  m, "IDASetQuadErrConB",
  [](const IDAView& ida_mem, int which, int errconQB) -> int
  {
    auto IDASetQuadErrConB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, int errconQB) -> int
    {
      auto lambda_result = IDASetQuadErrConB(ida_mem.get(), which, errconQB);
      return lambda_result;
    };

    return IDASetQuadErrConB_adapt_solver_memory(ida_mem, which, errconQB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("errconQB"));

sundials4py::scoped_def(
  m, "IDASetNonlinearSolverB",
  [](const IDAView& ida_mem, int which, SUNNonlinearSolver NLS) -> int
  {
    auto IDASetNonlinearSolverB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = IDASetNonlinearSolverB(ida_mem.get(), which, NLS);
      return lambda_result;
    };

    return IDASetNonlinearSolverB_adapt_solver_memory(ida_mem, which, NLS);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "IDAGetB",
  [](const IDAView& ida_mem, int which, N_Vector yy,
     N_Vector yp) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetB_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, N_Vector yy,
         N_Vector yp) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetB(ida_mem, which, &tret_adapt_modifiable, yy, yp);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto IDAGetB_adapt_solver_memory =
      [&IDAGetB_adapt_modifiable_immutable_to_return](const IDAView& ida_mem,
                                                      int which, N_Vector yy,
                                                      N_Vector yp)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetB_adapt_modifiable_immutable_to_return(ida_mem.get(), which, yy,
                                                     yp);
      return lambda_result;
    };

    return IDAGetB_adapt_solver_memory(ida_mem, which, yy, yp);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("yy"), nb::arg("yp"));

sundials4py::scoped_def(
  m, "IDAGetQuadB",
  [](const IDAView& ida_mem, int which, N_Vector qB) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetQuadB_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, N_Vector qB) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetQuadB(ida_mem, which, &tret_adapt_modifiable, qB);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto IDAGetQuadB_adapt_solver_memory =
      [&IDAGetQuadB_adapt_modifiable_immutable_to_return](const IDAView& ida_mem,
                                                          int which, N_Vector qB)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetQuadB_adapt_modifiable_immutable_to_return(ida_mem.get(), which,
                                                         qB);
      return lambda_result;
    };

    return IDAGetQuadB_adapt_solver_memory(ida_mem, which, qB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("qB"));

sundials4py::scoped_def(
  m, "IDAGetAdjIDABmem",
  [](const IDAView& ida_mem, int which) -> void*
  {
    auto IDAGetAdjIDABmem_adapt_solver_memory = [](const IDAView& ida_mem,
                                                   int which) -> void*
    {
      auto lambda_result = IDAGetAdjIDABmem(ida_mem.get(), which);
      return lambda_result;
    };

    return IDAGetAdjIDABmem_adapt_solver_memory(ida_mem, which);
  },
  nb::arg("ida_mem"), nb::arg("which"));

sundials4py::scoped_def(
  m, "IDAGetConsistentICB",
  [](const IDAView& ida_mem, int which, N_Vector yyB0, N_Vector ypB0) -> int
  {
    auto IDAGetConsistentICB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, N_Vector yyB0, N_Vector ypB0) -> int
    {
      auto lambda_result = IDAGetConsistentICB(ida_mem.get(), which, yyB0, ypB0);
      return lambda_result;
    };

    return IDAGetConsistentICB_adapt_solver_memory(ida_mem, which, yyB0, ypB0);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("yyB0"), nb::arg("ypB0"));

sundials4py::scoped_def(
  m, "IDAGetAdjY",
  [](const IDAView& ida_mem, sunrealtype t, N_Vector yy, N_Vector yp) -> int
  {
    auto IDAGetAdjY_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype t, N_Vector yy, N_Vector yp) -> int
    {
      auto lambda_result = IDAGetAdjY(ida_mem.get(), t, yy, yp);
      return lambda_result;
    };

    return IDAGetAdjY_adapt_solver_memory(ida_mem, t, yy, yp);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("yy"), nb::arg("yp"));

sundials4py::scoped_def(
  m, "IDAGetAdjDataPointHermite",
  [](const IDAView& ida_mem, int which, std::optional<N_Vector> yy = std::nullopt,
     std::optional<N_Vector> yd = std::nullopt) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetAdjDataPointHermite_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, N_Vector yy,
         N_Vector yd) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_adapt_modifiable;

      int r = IDAGetAdjDataPointHermite(ida_mem, which, &t_adapt_modifiable, yy,
                                        yd);
      return std::make_tuple(r, t_adapt_modifiable);
    };
    auto IDAGetAdjDataPointHermite_adapt_optional_arg_with_default_null =
      [&IDAGetAdjDataPointHermite_adapt_modifiable_immutable_to_return](void* ida_mem,
                                                                        int which,
                                                                        std::optional<N_Vector>
                                                                          yy =
                                                                            std::nullopt,
                                                                        std::optional<N_Vector>
                                                                          yd =
                                                                            std::nullopt)
      -> std::tuple<int, sunrealtype>
    {
      N_Vector yy_adapt_default_null = nullptr;
      if (yy.has_value()) yy_adapt_default_null = yy.value();
      N_Vector yd_adapt_default_null = nullptr;
      if (yd.has_value()) yd_adapt_default_null = yd.value();

      auto lambda_result =
        IDAGetAdjDataPointHermite_adapt_modifiable_immutable_to_return(ida_mem,
                                                                       which,
                                                                       yy_adapt_default_null,
                                                                       yd_adapt_default_null);
      return lambda_result;
    };
    auto IDAGetAdjDataPointHermite_adapt_solver_memory =
      [&IDAGetAdjDataPointHermite_adapt_optional_arg_with_default_null](const IDAView&
                                                                          ida_mem,
                                                                        int which,
                                                                        std::optional<N_Vector>
                                                                          yy =
                                                                            std::nullopt,
                                                                        std::optional<N_Vector>
                                                                          yd =
                                                                            std::nullopt)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetAdjDataPointHermite_adapt_optional_arg_with_default_null(ida_mem.get(),
                                                                       which,
                                                                       yy, yd);
      return lambda_result;
    };

    return IDAGetAdjDataPointHermite_adapt_solver_memory(ida_mem, which, yy, yd);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("yy").none() = nb::none(),
  nb::arg("yd").none() = nb::none());

sundials4py::scoped_def(
  m, "IDAGetAdjDataPointPolynomial",
  [](const IDAView& ida_mem, int which,
     std::optional<N_Vector> y = std::nullopt) -> std::tuple<int, sunrealtype, int>
  {
    auto IDAGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, N_Vector y) -> std::tuple<int, sunrealtype, int>
    {
      sunrealtype t_adapt_modifiable;
      int order_adapt_modifiable;

      int r = IDAGetAdjDataPointPolynomial(ida_mem, which, &t_adapt_modifiable,
                                           &order_adapt_modifiable, y);
      return std::make_tuple(r, t_adapt_modifiable, order_adapt_modifiable);
    };
    auto IDAGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null =
      [&IDAGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return](void* ida_mem,
                                                                           int which,
                                                                           std::optional<N_Vector>
                                                                             y =
                                                                               std::nullopt)
      -> std::tuple<int, sunrealtype, int>
    {
      N_Vector y_adapt_default_null = nullptr;
      if (y.has_value()) y_adapt_default_null = y.value();

      auto lambda_result =
        IDAGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return(ida_mem,
                                                                          which,
                                                                          y_adapt_default_null);
      return lambda_result;
    };
    auto IDAGetAdjDataPointPolynomial_adapt_solver_memory =
      [&IDAGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null](const IDAView&
                                                                             ida_mem,
                                                                           int which,
                                                                           std::optional<N_Vector>
                                                                             y =
                                                                               std::nullopt)
      -> std::tuple<int, sunrealtype, int>
    {
      auto lambda_result =
        IDAGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null(ida_mem
                                                                            .get(),
                                                                          which,
                                                                          y);
      return lambda_result;
    };

    return IDAGetAdjDataPointPolynomial_adapt_solver_memory(ida_mem, which, y);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("y").none() = nb::none());
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _IDASLS_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("IDALS_SUCCESS")         = 0;
m.attr("IDALS_MEM_NULL")        = -1;
m.attr("IDALS_LMEM_NULL")       = -2;
m.attr("IDALS_ILL_INPUT")       = -3;
m.attr("IDALS_MEM_FAIL")        = -4;
m.attr("IDALS_PMEM_NULL")       = -5;
m.attr("IDALS_JACFUNC_UNRECVR") = -6;
m.attr("IDALS_JACFUNC_RECVR")   = -7;
m.attr("IDALS_SUNMAT_FAIL")     = -8;
m.attr("IDALS_SUNLS_FAIL")      = -9;
m.attr("IDALS_NO_ADJ")          = -101;
m.attr("IDALS_LMEMB_NULL")      = -102;

sundials4py::scoped_def(
  m, "IDASetLinearSolver",
  [](const IDAView& ida_mem, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto IDASetLinearSolver_adapt_optional_arg_with_default_null =
      [](void* ida_mem, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = IDASetLinearSolver(ida_mem, LS, A_adapt_default_null);
      return lambda_result;
    };
    auto IDASetLinearSolver_adapt_solver_memory =
      [&IDASetLinearSolver_adapt_optional_arg_with_default_null](const IDAView& ida_mem,
                                                                 SUNLinearSolver LS,
                                                                 std::optional<SUNMatrix> A =
                                                                   std::nullopt) -> int
    {
      auto lambda_result =
        IDASetLinearSolver_adapt_optional_arg_with_default_null(ida_mem.get(),
                                                                LS, A);
      return lambda_result;
    };

    return IDASetLinearSolver_adapt_solver_memory(ida_mem, LS, A);
  },
  nb::arg("ida_mem"), nb::arg("LS"), nb::arg("A").none() = nb::none());

sundials4py::scoped_def(
  m, "IDASetEpsLin",
  [](const IDAView& ida_mem, sunrealtype eplifac) -> int
  {
    auto IDASetEpsLin_adapt_solver_memory = [](const IDAView& ida_mem,
                                               sunrealtype eplifac) -> int
    {
      auto lambda_result = IDASetEpsLin(ida_mem.get(), eplifac);
      return lambda_result;
    };

    return IDASetEpsLin_adapt_solver_memory(ida_mem, eplifac);
  },
  nb::arg("ida_mem"), nb::arg("eplifac"));

sundials4py::scoped_def(
  m, "IDASetLSNormFactor",
  [](const IDAView& ida_mem, sunrealtype nrmfac) -> int
  {
    auto IDASetLSNormFactor_adapt_solver_memory = [](const IDAView& ida_mem,
                                                     sunrealtype nrmfac) -> int
    {
      auto lambda_result = IDASetLSNormFactor(ida_mem.get(), nrmfac);
      return lambda_result;
    };

    return IDASetLSNormFactor_adapt_solver_memory(ida_mem, nrmfac);
  },
  nb::arg("ida_mem"), nb::arg("nrmfac"));

sundials4py::scoped_def(
  m, "IDASetLinearSolutionScaling",
  [](const IDAView& ida_mem, sunbooleantype onoff) -> int
  {
    auto IDASetLinearSolutionScaling_adapt_solver_memory =
      [](const IDAView& ida_mem, sunbooleantype onoff) -> int
    {
      auto lambda_result = IDASetLinearSolutionScaling(ida_mem.get(), onoff);
      return lambda_result;
    };

    return IDASetLinearSolutionScaling_adapt_solver_memory(ida_mem, onoff);
  },
  nb::arg("ida_mem"), nb::arg("onoff"));

sundials4py::scoped_def(
  m, "IDASetIncrementFactor",
  [](const IDAView& ida_mem, sunrealtype dqincfac) -> int
  {
    auto IDASetIncrementFactor_adapt_solver_memory =
      [](const IDAView& ida_mem, sunrealtype dqincfac) -> int
    {
      auto lambda_result = IDASetIncrementFactor(ida_mem.get(), dqincfac);
      return lambda_result;
    };

    return IDASetIncrementFactor_adapt_solver_memory(ida_mem, dqincfac);
  },
  nb::arg("ida_mem"), nb::arg("dqincfac"));

sundials4py::scoped_def(
  m, "IDAGetJac",
  [](const IDAView& ida_mem) -> std::tuple<int, SUNMatrix>
  {
    auto IDAGetJac_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix J_adapt_modifiable;

      int r = IDAGetJac(ida_mem, &J_adapt_modifiable);
      return std::make_tuple(r, J_adapt_modifiable);
    };
    auto IDAGetJac_adapt_solver_memory =
      [&IDAGetJac_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, SUNMatrix>
    {
      auto lambda_result =
        IDAGetJac_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetJac_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "IDAGetJacCj",
  [](const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetJacCj_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype cj_J_adapt_modifiable;

      int r = IDAGetJacCj(ida_mem, &cj_J_adapt_modifiable);
      return std::make_tuple(r, cj_J_adapt_modifiable);
    };
    auto IDAGetJacCj_adapt_solver_memory =
      [&IDAGetJacCj_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetJacCj_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetJacCj_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetJacTime",
  [](const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetJacTime_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_J_adapt_modifiable;

      int r = IDAGetJacTime(ida_mem, &t_J_adapt_modifiable);
      return std::make_tuple(r, t_J_adapt_modifiable);
    };
    auto IDAGetJacTime_adapt_solver_memory =
      [&IDAGetJacTime_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        IDAGetJacTime_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetJacTime_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetJacNumSteps",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetJacNumSteps_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nst_J_adapt_modifiable;

      int r = IDAGetJacNumSteps(ida_mem, &nst_J_adapt_modifiable);
      return std::make_tuple(r, nst_J_adapt_modifiable);
    };
    auto IDAGetJacNumSteps_adapt_solver_memory =
      [&IDAGetJacNumSteps_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetJacNumSteps_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetJacNumSteps_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumJacEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long njevals_adapt_modifiable;

      int r = IDAGetNumJacEvals(ida_mem, &njevals_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable);
    };
    auto IDAGetNumJacEvals_adapt_solver_memory =
      [&IDAGetNumJacEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumJacEvals_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumJacEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumPrecEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long npevals_adapt_modifiable;

      int r = IDAGetNumPrecEvals(ida_mem, &npevals_adapt_modifiable);
      return std::make_tuple(r, npevals_adapt_modifiable);
    };
    auto IDAGetNumPrecEvals_adapt_solver_memory =
      [&IDAGetNumPrecEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumPrecEvals_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumPrecEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumPrecSolves",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long npsolves_adapt_modifiable;

      int r = IDAGetNumPrecSolves(ida_mem, &npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves_adapt_modifiable);
    };
    auto IDAGetNumPrecSolves_adapt_solver_memory =
      [&IDAGetNumPrecSolves_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumPrecSolves_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumPrecSolves_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumLinIters",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nliters_adapt_modifiable;

      int r = IDAGetNumLinIters(ida_mem, &nliters_adapt_modifiable);
      return std::make_tuple(r, nliters_adapt_modifiable);
    };
    auto IDAGetNumLinIters_adapt_solver_memory =
      [&IDAGetNumLinIters_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumLinIters_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumLinIters_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumLinConvFails",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nlcfails_adapt_modifiable;

      int r = IDAGetNumLinConvFails(ida_mem, &nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails_adapt_modifiable);
    };
    auto IDAGetNumLinConvFails_adapt_solver_memory =
      [&IDAGetNumLinConvFails_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumLinConvFails_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumLinConvFails_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumJTSetupEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumJTSetupEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long njtsetups_adapt_modifiable;

      int r = IDAGetNumJTSetupEvals(ida_mem, &njtsetups_adapt_modifiable);
      return std::make_tuple(r, njtsetups_adapt_modifiable);
    };
    auto IDAGetNumJTSetupEvals_adapt_solver_memory =
      [&IDAGetNumJTSetupEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumJTSetupEvals_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumJTSetupEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumJtimesEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long njvevals_adapt_modifiable;

      int r = IDAGetNumJtimesEvals(ida_mem, &njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals_adapt_modifiable);
    };
    auto IDAGetNumJtimesEvals_adapt_solver_memory =
      [&IDAGetNumJtimesEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumJtimesEvals_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumJtimesEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetNumLinResEvals",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumLinResEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nrevalsLS_adapt_modifiable;

      int r = IDAGetNumLinResEvals(ida_mem, &nrevalsLS_adapt_modifiable);
      return std::make_tuple(r, nrevalsLS_adapt_modifiable);
    };
    auto IDAGetNumLinResEvals_adapt_solver_memory =
      [&IDAGetNumLinResEvals_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetNumLinResEvals_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetNumLinResEvals_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(
  m, "IDAGetLastLinFlag",
  [](const IDAView& ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = IDAGetLastLinFlag(ida_mem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };
    auto IDAGetLastLinFlag_adapt_solver_memory =
      [&IDAGetLastLinFlag_adapt_modifiable_immutable_to_return](
        const IDAView& ida_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        IDAGetLastLinFlag_adapt_modifiable_immutable_to_return(ida_mem.get());
      return lambda_result;
    };

    return IDAGetLastLinFlag_adapt_solver_memory(ida_mem);
  },
  nb::arg("ida_mem"));

sundials4py::scoped_def(m, "IDAGetLinReturnFlagName", IDAGetLinReturnFlagName,
                        nb::arg("flag"));

sundials4py::scoped_def(
  m, "IDASetLinearSolverB",
  [](const IDAView& ida_mem, int which, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto IDASetLinearSolverB_adapt_optional_arg_with_default_null =
      [](void* ida_mem, int which, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = IDASetLinearSolverB(ida_mem, which, LS,
                                               A_adapt_default_null);
      return lambda_result;
    };
    auto IDASetLinearSolverB_adapt_solver_memory =
      [&IDASetLinearSolverB_adapt_optional_arg_with_default_null](const IDAView& ida_mem,
                                                                  int which,
                                                                  SUNLinearSolver LS,
                                                                  std::optional<SUNMatrix> A =
                                                                    std::nullopt) -> int
    {
      auto lambda_result =
        IDASetLinearSolverB_adapt_optional_arg_with_default_null(ida_mem.get(),
                                                                 which, LS, A);
      return lambda_result;
    };

    return IDASetLinearSolverB_adapt_solver_memory(ida_mem, which, LS, A);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("LS"),
  nb::arg("A").none() = nb::none());

sundials4py::scoped_def(
  m, "IDASetEpsLinB",
  [](const IDAView& ida_mem, int which, sunrealtype eplifacB) -> int
  {
    auto IDASetEpsLinB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunrealtype eplifacB) -> int
    {
      auto lambda_result = IDASetEpsLinB(ida_mem.get(), which, eplifacB);
      return lambda_result;
    };

    return IDASetEpsLinB_adapt_solver_memory(ida_mem, which, eplifacB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("eplifacB"));

sundials4py::scoped_def(
  m, "IDASetLSNormFactorB",
  [](const IDAView& ida_mem, int which, sunrealtype nrmfacB) -> int
  {
    auto IDASetLSNormFactorB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunrealtype nrmfacB) -> int
    {
      auto lambda_result = IDASetLSNormFactorB(ida_mem.get(), which, nrmfacB);
      return lambda_result;
    };

    return IDASetLSNormFactorB_adapt_solver_memory(ida_mem, which, nrmfacB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("nrmfacB"));

sundials4py::scoped_def(
  m, "IDASetLinearSolutionScalingB",
  [](const IDAView& ida_mem, int which, sunbooleantype onoffB) -> int
  {
    auto IDASetLinearSolutionScalingB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunbooleantype onoffB) -> int
    {
      auto lambda_result = IDASetLinearSolutionScalingB(ida_mem.get(), which,
                                                        onoffB);
      return lambda_result;
    };

    return IDASetLinearSolutionScalingB_adapt_solver_memory(ida_mem, which,
                                                            onoffB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("onoffB"));

sundials4py::scoped_def(
  m, "IDASetIncrementFactorB",
  [](const IDAView& ida_mem, int which, sunrealtype dqincfacB) -> int
  {
    auto IDASetIncrementFactorB_adapt_solver_memory =
      [](const IDAView& ida_mem, int which, sunrealtype dqincfacB) -> int
    {
      auto lambda_result = IDASetIncrementFactorB(ida_mem.get(), which,
                                                  dqincfacB);
      return lambda_result;
    };

    return IDASetIncrementFactorB_adapt_solver_memory(ida_mem, which, dqincfacB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("dqincfacB"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
