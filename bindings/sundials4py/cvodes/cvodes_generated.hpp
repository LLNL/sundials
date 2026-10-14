// #ifndef _CVODES_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("CV_ADAMS")               = 1;
m.attr("CV_BDF")                 = 2;
m.attr("CV_NORMAL")              = 1;
m.attr("CV_ONE_STEP")            = 2;
m.attr("CV_SIMULTANEOUS")        = 1;
m.attr("CV_STAGGERED")           = 2;
m.attr("CV_STAGGERED1")          = 3;
m.attr("CV_CENTERED")            = 1;
m.attr("CV_FORWARD")             = 2;
m.attr("CV_HERMITE")             = 1;
m.attr("CV_POLYNOMIAL")          = 2;
m.attr("CV_SUCCESS")             = 0;
m.attr("CV_TSTOP_RETURN")        = 1;
m.attr("CV_ROOT_RETURN")         = 2;
m.attr("CV_WARNING")             = 99;
m.attr("CV_TOO_MUCH_WORK")       = -1;
m.attr("CV_TOO_MUCH_ACC")        = -2;
m.attr("CV_ERR_FAILURE")         = -3;
m.attr("CV_CONV_FAILURE")        = -4;
m.attr("CV_LINIT_FAIL")          = -5;
m.attr("CV_LSETUP_FAIL")         = -6;
m.attr("CV_LSOLVE_FAIL")         = -7;
m.attr("CV_RHSFUNC_FAIL")        = -8;
m.attr("CV_FIRST_RHSFUNC_ERR")   = -9;
m.attr("CV_REPTD_RHSFUNC_ERR")   = -10;
m.attr("CV_UNREC_RHSFUNC_ERR")   = -11;
m.attr("CV_RTFUNC_FAIL")         = -12;
m.attr("CV_NLS_INIT_FAIL")       = -13;
m.attr("CV_NLS_SETUP_FAIL")      = -14;
m.attr("CV_CONSTR_FAIL")         = -15;
m.attr("CV_NLS_FAIL")            = -16;
m.attr("CV_MEM_FAIL")            = -20;
m.attr("CV_MEM_NULL")            = -21;
m.attr("CV_ILL_INPUT")           = -22;
m.attr("CV_NO_MALLOC")           = -23;
m.attr("CV_BAD_K")               = -24;
m.attr("CV_BAD_T")               = -25;
m.attr("CV_BAD_DKY")             = -26;
m.attr("CV_TOO_CLOSE")           = -27;
m.attr("CV_VECTOROP_ERR")        = -28;
m.attr("CV_NO_QUAD")             = -30;
m.attr("CV_QRHSFUNC_FAIL")       = -31;
m.attr("CV_FIRST_QRHSFUNC_ERR")  = -32;
m.attr("CV_REPTD_QRHSFUNC_ERR")  = -33;
m.attr("CV_UNREC_QRHSFUNC_ERR")  = -34;
m.attr("CV_NO_SENS")             = -40;
m.attr("CV_SRHSFUNC_FAIL")       = -41;
m.attr("CV_FIRST_SRHSFUNC_ERR")  = -42;
m.attr("CV_REPTD_SRHSFUNC_ERR")  = -43;
m.attr("CV_UNREC_SRHSFUNC_ERR")  = -44;
m.attr("CV_BAD_IS")              = -45;
m.attr("CV_NO_QUADSENS")         = -50;
m.attr("CV_QSRHSFUNC_FAIL")      = -51;
m.attr("CV_FIRST_QSRHSFUNC_ERR") = -52;
m.attr("CV_REPTD_QSRHSFUNC_ERR") = -53;
m.attr("CV_UNREC_QSRHSFUNC_ERR") = -54;
m.attr("CV_CONTEXT_ERR")         = -55;
m.attr("CV_PROJ_MEM_NULL")       = -56;
m.attr("CV_PROJFUNC_FAIL")       = -57;
m.attr("CV_REPTD_PROJFUNC_ERR")  = -58;
m.attr("CV_BAD_TINTERP")         = -59;
m.attr("CV_UNRECOGNIZED_ERR")    = -99;
m.attr("CV_NO_ADJ")              = -101;
m.attr("CV_NO_FWD")              = -102;
m.attr("CV_NO_BCK")              = -103;
m.attr("CV_BAD_TB0")             = -104;
m.attr("CV_REIFWD_FAIL")         = -105;
m.attr("CV_FWD_FAIL")            = -106;
m.attr("CV_GETY_BADT")           = -107;

sundials4py::scoped_def(
  m, "CVodeReInit",
  [](const CVodeView& cvode_mem, sunrealtype t0, N_Vector y0) -> int
  {
    auto CVodeReInit_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                              sunrealtype t0, N_Vector y0) -> int
    {
      auto lambda_result = CVodeReInit(cvode_mem.get(), t0, y0);
      return lambda_result;
    };

    return CVodeReInit_adapt_solver_memory(cvode_mem, t0, y0);
  },
  nb::arg("cvode_mem"), nb::arg("t0"), nb::arg("y0"));

sundials4py::scoped_def(
  m, "CVodeResizeHistory",
  [](const CVodeView& cvode_mem, sundials4py::Array1d t_hist_1d,
     std::vector<N_Vector> y_hist_1d, std::vector<N_Vector> f_hist_1d,
     int num_y_hist, int num_f_hist) -> int
  {
    auto CVodeResizeHistory_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sundials4py::Array1d t_hist_1d,
         std::vector<N_Vector> y_hist_1d, std::vector<N_Vector> f_hist_1d,
         int num_y_hist, int num_f_hist) -> int
    {
      sunrealtype* t_hist_1d_ptr = t_hist_1d.size() == 0 ? nullptr
                                                         : t_hist_1d.data();
      N_Vector* y_hist_1d_ptr = y_hist_1d.empty() ? nullptr : y_hist_1d.data();
      N_Vector* f_hist_1d_ptr = f_hist_1d.empty() ? nullptr : f_hist_1d.data();

      auto lambda_result = CVodeResizeHistory(cvode_mem, t_hist_1d_ptr,
                                              y_hist_1d_ptr, f_hist_1d_ptr,
                                              num_y_hist, num_f_hist);
      return lambda_result;
    };
    auto CVodeResizeHistory_adapt_solver_memory =
      [&CVodeResizeHistory_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                        sundials4py::Array1d t_hist_1d,
                                                        std::vector<N_Vector> y_hist_1d,
                                                        std::vector<N_Vector> f_hist_1d,
                                                        int num_y_hist,
                                                        int num_f_hist) -> int
    {
      auto lambda_result =
        CVodeResizeHistory_adapt_arr_ptr_to_std_vector(cvode_mem.get(), t_hist_1d,
                                                       y_hist_1d, f_hist_1d,
                                                       num_y_hist, num_f_hist);
      return lambda_result;
    };

    return CVodeResizeHistory_adapt_solver_memory(cvode_mem, t_hist_1d,
                                                  y_hist_1d, f_hist_1d,
                                                  num_y_hist, num_f_hist);
  },
  nb::arg("cvode_mem"), nb::arg("t_hist_1d"), nb::arg("y_hist_1d"),
  nb::arg("f_hist_1d"), nb::arg("num_y_hist"), nb::arg("num_f_hist"));

sundials4py::scoped_def(
  m, "CVodeSStolerances",
  [](const CVodeView& cvode_mem, sunrealtype reltol, sunrealtype abstol) -> int
  {
    auto CVodeSStolerances_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                    sunrealtype reltol,
                                                    sunrealtype abstol) -> int
    {
      auto lambda_result = CVodeSStolerances(cvode_mem.get(), reltol, abstol);
      return lambda_result;
    };

    return CVodeSStolerances_adapt_solver_memory(cvode_mem, reltol, abstol);
  },
  nb::arg("cvode_mem"), nb::arg("reltol"), nb::arg("abstol"));

sundials4py::scoped_def(
  m, "CVodeSVtolerances",
  [](const CVodeView& cvode_mem, sunrealtype reltol, N_Vector abstol) -> int
  {
    auto CVodeSVtolerances_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype reltol, N_Vector abstol) -> int
    {
      auto lambda_result = CVodeSVtolerances(cvode_mem.get(), reltol, abstol);
      return lambda_result;
    };

    return CVodeSVtolerances_adapt_solver_memory(cvode_mem, reltol, abstol);
  },
  nb::arg("cvode_mem"), nb::arg("reltol"), nb::arg("abstol"));

sundials4py::scoped_def(
  m, "CVodeSetConstraints",
  [](const CVodeView& cvode_mem, N_Vector constraints) -> int
  {
    auto CVodeSetConstraints_adapt_solver_memory =
      [](const CVodeView& cvode_mem, N_Vector constraints) -> int
    {
      auto lambda_result = CVodeSetConstraints(cvode_mem.get(), constraints);
      return lambda_result;
    };

    return CVodeSetConstraints_adapt_solver_memory(cvode_mem, constraints);
  },
  nb::arg("cvode_mem"), nb::arg("constraints"));

sundials4py::scoped_def(
  m, "CVodeSetMaxNumConstraintFails",
  [](const CVodeView& cvode_mem, int max_fails) -> int
  {
    auto CVodeSetMaxNumConstraintFails_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int max_fails) -> int
    {
      auto lambda_result = CVodeSetMaxNumConstraintFails(cvode_mem.get(),
                                                         max_fails);
      return lambda_result;
    };

    return CVodeSetMaxNumConstraintFails_adapt_solver_memory(cvode_mem,
                                                             max_fails);
  },
  nb::arg("cvode_mem"), nb::arg("max_fails"));

sundials4py::scoped_def(
  m, "CVodeSetDeltaGammaMaxLSetup",
  [](const CVodeView& cvode_mem, sunrealtype dgmax_lsetup) -> int
  {
    auto CVodeSetDeltaGammaMaxLSetup_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype dgmax_lsetup) -> int
    {
      auto lambda_result = CVodeSetDeltaGammaMaxLSetup(cvode_mem.get(),
                                                       dgmax_lsetup);
      return lambda_result;
    };

    return CVodeSetDeltaGammaMaxLSetup_adapt_solver_memory(cvode_mem,
                                                           dgmax_lsetup);
  },
  nb::arg("cvode_mem"), nb::arg("dgmax_lsetup"));

sundials4py::scoped_def(
  m, "CVodeSetInitStep",
  [](const CVodeView& cvode_mem, sunrealtype hin) -> int
  {
    auto CVodeSetInitStep_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                   sunrealtype hin) -> int
    {
      auto lambda_result = CVodeSetInitStep(cvode_mem.get(), hin);
      return lambda_result;
    };

    return CVodeSetInitStep_adapt_solver_memory(cvode_mem, hin);
  },
  nb::arg("cvode_mem"), nb::arg("hin"));

sundials4py::scoped_def(
  m, "CVodeSetLSetupFrequency",
  [](const CVodeView& cvode_mem, long msbp) -> int
  {
    auto CVodeSetLSetupFrequency_adapt_solver_memory =
      [](const CVodeView& cvode_mem, long msbp) -> int
    {
      auto lambda_result = CVodeSetLSetupFrequency(cvode_mem.get(), msbp);
      return lambda_result;
    };

    return CVodeSetLSetupFrequency_adapt_solver_memory(cvode_mem, msbp);
  },
  nb::arg("cvode_mem"), nb::arg("msbp"));

sundials4py::scoped_def(
  m, "CVodeSetMaxConvFails",
  [](const CVodeView& cvode_mem, int maxncf) -> int
  {
    auto CVodeSetMaxConvFails_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int maxncf) -> int
    {
      auto lambda_result = CVodeSetMaxConvFails(cvode_mem.get(), maxncf);
      return lambda_result;
    };

    return CVodeSetMaxConvFails_adapt_solver_memory(cvode_mem, maxncf);
  },
  nb::arg("cvode_mem"), nb::arg("maxncf"));

sundials4py::scoped_def(
  m, "CVodeSetMaxErrTestFails",
  [](const CVodeView& cvode_mem, int maxnef) -> int
  {
    auto CVodeSetMaxErrTestFails_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int maxnef) -> int
    {
      auto lambda_result = CVodeSetMaxErrTestFails(cvode_mem.get(), maxnef);
      return lambda_result;
    };

    return CVodeSetMaxErrTestFails_adapt_solver_memory(cvode_mem, maxnef);
  },
  nb::arg("cvode_mem"), nb::arg("maxnef"));

sundials4py::scoped_def(
  m, "CVodeSetMaxHnilWarns",
  [](const CVodeView& cvode_mem, int mxhnil) -> int
  {
    auto CVodeSetMaxHnilWarns_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int mxhnil) -> int
    {
      auto lambda_result = CVodeSetMaxHnilWarns(cvode_mem.get(), mxhnil);
      return lambda_result;
    };

    return CVodeSetMaxHnilWarns_adapt_solver_memory(cvode_mem, mxhnil);
  },
  nb::arg("cvode_mem"), nb::arg("mxhnil"));

sundials4py::scoped_def(
  m, "CVodeSetMaxNonlinIters",
  [](const CVodeView& cvode_mem, int maxcor) -> int
  {
    auto CVodeSetMaxNonlinIters_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int maxcor) -> int
    {
      auto lambda_result = CVodeSetMaxNonlinIters(cvode_mem.get(), maxcor);
      return lambda_result;
    };

    return CVodeSetMaxNonlinIters_adapt_solver_memory(cvode_mem, maxcor);
  },
  nb::arg("cvode_mem"), nb::arg("maxcor"));

sundials4py::scoped_def(
  m, "CVodeSetMaxNumSteps",
  [](const CVodeView& cvode_mem, long mxsteps) -> int
  {
    auto CVodeSetMaxNumSteps_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                      long mxsteps) -> int
    {
      auto lambda_result = CVodeSetMaxNumSteps(cvode_mem.get(), mxsteps);
      return lambda_result;
    };

    return CVodeSetMaxNumSteps_adapt_solver_memory(cvode_mem, mxsteps);
  },
  nb::arg("cvode_mem"), nb::arg("mxsteps"));

sundials4py::scoped_def(
  m, "CVodeSetMaxOrd",
  [](const CVodeView& cvode_mem, int maxord) -> int
  {
    auto CVodeSetMaxOrd_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                 int maxord) -> int
    {
      auto lambda_result = CVodeSetMaxOrd(cvode_mem.get(), maxord);
      return lambda_result;
    };

    return CVodeSetMaxOrd_adapt_solver_memory(cvode_mem, maxord);
  },
  nb::arg("cvode_mem"), nb::arg("maxord"));

sundials4py::scoped_def(
  m, "CVodeSetMaxStep",
  [](const CVodeView& cvode_mem, sunrealtype hmax) -> int
  {
    auto CVodeSetMaxStep_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                  sunrealtype hmax) -> int
    {
      auto lambda_result = CVodeSetMaxStep(cvode_mem.get(), hmax);
      return lambda_result;
    };

    return CVodeSetMaxStep_adapt_solver_memory(cvode_mem, hmax);
  },
  nb::arg("cvode_mem"), nb::arg("hmax"));

sundials4py::scoped_def(
  m, "CVodeSetMinStep",
  [](const CVodeView& cvode_mem, sunrealtype hmin) -> int
  {
    auto CVodeSetMinStep_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                  sunrealtype hmin) -> int
    {
      auto lambda_result = CVodeSetMinStep(cvode_mem.get(), hmin);
      return lambda_result;
    };

    return CVodeSetMinStep_adapt_solver_memory(cvode_mem, hmin);
  },
  nb::arg("cvode_mem"), nb::arg("hmin"));

sundials4py::scoped_def(
  m, "CVodeSetMonitorFrequency",
  [](const CVodeView& cvode_mem, long nst) -> int
  {
    auto CVodeSetMonitorFrequency_adapt_solver_memory =
      [](const CVodeView& cvode_mem, long nst) -> int
    {
      auto lambda_result = CVodeSetMonitorFrequency(cvode_mem.get(), nst);
      return lambda_result;
    };

    return CVodeSetMonitorFrequency_adapt_solver_memory(cvode_mem, nst);
  },
  nb::arg("cvode_mem"), nb::arg("nst"));

sundials4py::scoped_def(
  m, "CVodeSetNonlinConvCoef",
  [](const CVodeView& cvode_mem, sunrealtype nlscoef) -> int
  {
    auto CVodeSetNonlinConvCoef_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype nlscoef) -> int
    {
      auto lambda_result = CVodeSetNonlinConvCoef(cvode_mem.get(), nlscoef);
      return lambda_result;
    };

    return CVodeSetNonlinConvCoef_adapt_solver_memory(cvode_mem, nlscoef);
  },
  nb::arg("cvode_mem"), nb::arg("nlscoef"));

sundials4py::scoped_def(
  m, "CVodeSetNonlinearSolver",
  [](const CVodeView& cvode_mem, SUNNonlinearSolver NLS) -> int
  {
    auto CVodeSetNonlinearSolver_adapt_solver_memory =
      [](const CVodeView& cvode_mem, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = CVodeSetNonlinearSolver(cvode_mem.get(), NLS);
      return lambda_result;
    };

    return CVodeSetNonlinearSolver_adapt_solver_memory(cvode_mem, NLS);
  },
  nb::arg("cvode_mem"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "CVodeSetStabLimDet",
  [](const CVodeView& cvode_mem, sunbooleantype stldet) -> int
  {
    auto CVodeSetStabLimDet_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunbooleantype stldet) -> int
    {
      auto lambda_result = CVodeSetStabLimDet(cvode_mem.get(), stldet);
      return lambda_result;
    };

    return CVodeSetStabLimDet_adapt_solver_memory(cvode_mem, stldet);
  },
  nb::arg("cvode_mem"), nb::arg("stldet"));

sundials4py::scoped_def(
  m, "CVodeSetStopTime",
  [](const CVodeView& cvode_mem, sunrealtype tstop) -> int
  {
    auto CVodeSetStopTime_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                   sunrealtype tstop) -> int
    {
      auto lambda_result = CVodeSetStopTime(cvode_mem.get(), tstop);
      return lambda_result;
    };

    return CVodeSetStopTime_adapt_solver_memory(cvode_mem, tstop);
  },
  nb::arg("cvode_mem"), nb::arg("tstop"));

sundials4py::scoped_def(
  m, "CVodeSetInterpolateStopTime",
  [](const CVodeView& cvode_mem, sunbooleantype interp) -> int
  {
    auto CVodeSetInterpolateStopTime_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunbooleantype interp) -> int
    {
      auto lambda_result = CVodeSetInterpolateStopTime(cvode_mem.get(), interp);
      return lambda_result;
    };

    return CVodeSetInterpolateStopTime_adapt_solver_memory(cvode_mem, interp);
  },
  nb::arg("cvode_mem"), nb::arg("interp"));

sundials4py::scoped_def(
  m, "CVodeClearStopTime",
  [](const CVodeView& cvode_mem) -> int
  {
    auto CVodeClearStopTime_adapt_solver_memory =
      [](const CVodeView& cvode_mem) -> int
    {
      auto lambda_result = CVodeClearStopTime(cvode_mem.get());
      return lambda_result;
    };

    return CVodeClearStopTime_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeSetEtaFixedStepBounds",
  [](const CVodeView& cvode_mem, sunrealtype eta_min_fx,
     sunrealtype eta_max_fx) -> int
  {
    auto CVodeSetEtaFixedStepBounds_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype eta_min_fx,
         sunrealtype eta_max_fx) -> int
    {
      auto lambda_result = CVodeSetEtaFixedStepBounds(cvode_mem.get(),
                                                      eta_min_fx, eta_max_fx);
      return lambda_result;
    };

    return CVodeSetEtaFixedStepBounds_adapt_solver_memory(cvode_mem, eta_min_fx,
                                                          eta_max_fx);
  },
  nb::arg("cvode_mem"), nb::arg("eta_min_fx"), nb::arg("eta_max_fx"));

sundials4py::scoped_def(
  m, "CVodeSetEtaMaxFirstStep",
  [](const CVodeView& cvode_mem, sunrealtype eta_max_fs) -> int
  {
    auto CVodeSetEtaMaxFirstStep_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype eta_max_fs) -> int
    {
      auto lambda_result = CVodeSetEtaMaxFirstStep(cvode_mem.get(), eta_max_fs);
      return lambda_result;
    };

    return CVodeSetEtaMaxFirstStep_adapt_solver_memory(cvode_mem, eta_max_fs);
  },
  nb::arg("cvode_mem"), nb::arg("eta_max_fs"));

sundials4py::scoped_def(
  m, "CVodeSetEtaMaxEarlyStep",
  [](const CVodeView& cvode_mem, sunrealtype eta_max_es) -> int
  {
    auto CVodeSetEtaMaxEarlyStep_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype eta_max_es) -> int
    {
      auto lambda_result = CVodeSetEtaMaxEarlyStep(cvode_mem.get(), eta_max_es);
      return lambda_result;
    };

    return CVodeSetEtaMaxEarlyStep_adapt_solver_memory(cvode_mem, eta_max_es);
  },
  nb::arg("cvode_mem"), nb::arg("eta_max_es"));

sundials4py::scoped_def(
  m, "CVodeSetNumStepsEtaMaxEarlyStep",
  [](const CVodeView& cvode_mem, long small_nst) -> int
  {
    auto CVodeSetNumStepsEtaMaxEarlyStep_adapt_solver_memory =
      [](const CVodeView& cvode_mem, long small_nst) -> int
    {
      auto lambda_result = CVodeSetNumStepsEtaMaxEarlyStep(cvode_mem.get(),
                                                           small_nst);
      return lambda_result;
    };

    return CVodeSetNumStepsEtaMaxEarlyStep_adapt_solver_memory(cvode_mem,
                                                               small_nst);
  },
  nb::arg("cvode_mem"), nb::arg("small_nst"));

sundials4py::scoped_def(
  m, "CVodeSetEtaMax",
  [](const CVodeView& cvode_mem, sunrealtype eta_max_gs) -> int
  {
    auto CVodeSetEtaMax_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                 sunrealtype eta_max_gs) -> int
    {
      auto lambda_result = CVodeSetEtaMax(cvode_mem.get(), eta_max_gs);
      return lambda_result;
    };

    return CVodeSetEtaMax_adapt_solver_memory(cvode_mem, eta_max_gs);
  },
  nb::arg("cvode_mem"), nb::arg("eta_max_gs"));

sundials4py::scoped_def(
  m, "CVodeSetEtaMin",
  [](const CVodeView& cvode_mem, sunrealtype eta_min) -> int
  {
    auto CVodeSetEtaMin_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                 sunrealtype eta_min) -> int
    {
      auto lambda_result = CVodeSetEtaMin(cvode_mem.get(), eta_min);
      return lambda_result;
    };

    return CVodeSetEtaMin_adapt_solver_memory(cvode_mem, eta_min);
  },
  nb::arg("cvode_mem"), nb::arg("eta_min"));

sundials4py::scoped_def(
  m, "CVodeSetEtaMinErrFail",
  [](const CVodeView& cvode_mem, sunrealtype eta_min_ef) -> int
  {
    auto CVodeSetEtaMinErrFail_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype eta_min_ef) -> int
    {
      auto lambda_result = CVodeSetEtaMinErrFail(cvode_mem.get(), eta_min_ef);
      return lambda_result;
    };

    return CVodeSetEtaMinErrFail_adapt_solver_memory(cvode_mem, eta_min_ef);
  },
  nb::arg("cvode_mem"), nb::arg("eta_min_ef"));

sundials4py::scoped_def(
  m, "CVodeSetEtaMaxErrFail",
  [](const CVodeView& cvode_mem, sunrealtype eta_max_ef) -> int
  {
    auto CVodeSetEtaMaxErrFail_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype eta_max_ef) -> int
    {
      auto lambda_result = CVodeSetEtaMaxErrFail(cvode_mem.get(), eta_max_ef);
      return lambda_result;
    };

    return CVodeSetEtaMaxErrFail_adapt_solver_memory(cvode_mem, eta_max_ef);
  },
  nb::arg("cvode_mem"), nb::arg("eta_max_ef"));

sundials4py::scoped_def(
  m, "CVodeSetNumFailsEtaMaxErrFail",
  [](const CVodeView& cvode_mem, int small_nef) -> int
  {
    auto CVodeSetNumFailsEtaMaxErrFail_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int small_nef) -> int
    {
      auto lambda_result = CVodeSetNumFailsEtaMaxErrFail(cvode_mem.get(),
                                                         small_nef);
      return lambda_result;
    };

    return CVodeSetNumFailsEtaMaxErrFail_adapt_solver_memory(cvode_mem,
                                                             small_nef);
  },
  nb::arg("cvode_mem"), nb::arg("small_nef"));

sundials4py::scoped_def(
  m, "CVodeSetEtaConvFail",
  [](const CVodeView& cvode_mem, sunrealtype eta_cf) -> int
  {
    auto CVodeSetEtaConvFail_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                      sunrealtype eta_cf) -> int
    {
      auto lambda_result = CVodeSetEtaConvFail(cvode_mem.get(), eta_cf);
      return lambda_result;
    };

    return CVodeSetEtaConvFail_adapt_solver_memory(cvode_mem, eta_cf);
  },
  nb::arg("cvode_mem"), nb::arg("eta_cf"));

sundials4py::scoped_def(
  m, "CVodeSetRootDirection",
  [](const CVodeView& cvode_mem, std::vector<int> rootdir_1d) -> int
  {
    auto CVodeSetRootDirection_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<int> rootdir_1d) -> int
    {
      int* rootdir_1d_ptr = rootdir_1d.empty() ? nullptr : rootdir_1d.data();

      auto lambda_result = CVodeSetRootDirection(cvode_mem, rootdir_1d_ptr);
      return lambda_result;
    };
    auto CVodeSetRootDirection_adapt_solver_memory =
      [&CVodeSetRootDirection_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                           std::vector<int>
                                                             rootdir_1d) -> int
    {
      auto lambda_result =
        CVodeSetRootDirection_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                          rootdir_1d);
      return lambda_result;
    };

    return CVodeSetRootDirection_adapt_solver_memory(cvode_mem, rootdir_1d);
  },
  nb::arg("cvode_mem"), nb::arg("rootdir_1d"));

sundials4py::scoped_def(
  m, "CVodeSetNoInactiveRootWarn",
  [](const CVodeView& cvode_mem) -> int
  {
    auto CVodeSetNoInactiveRootWarn_adapt_solver_memory =
      [](const CVodeView& cvode_mem) -> int
    {
      auto lambda_result = CVodeSetNoInactiveRootWarn(cvode_mem.get());
      return lambda_result;
    };

    return CVodeSetNoInactiveRootWarn_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeComputeState",
  [](const CVodeView& cvode_mem, N_Vector ycor, N_Vector y) -> int
  {
    auto CVodeComputeState_adapt_solver_memory =
      [](const CVodeView& cvode_mem, N_Vector ycor, N_Vector y) -> int
    {
      auto lambda_result = CVodeComputeState(cvode_mem.get(), ycor, y);
      return lambda_result;
    };

    return CVodeComputeState_adapt_solver_memory(cvode_mem, ycor, y);
  },
  nb::arg("cvode_mem"), nb::arg("ycor"), nb::arg("y"));

sundials4py::scoped_def(
  m, "CVodeComputeStateSens",
  [](const CVodeView& cvode_mem, std::vector<N_Vector> yScor_1d,
     std::vector<N_Vector> yS_1d) -> int
  {
    auto CVodeComputeStateSens_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<N_Vector> yScor_1d,
         std::vector<N_Vector> yS_1d) -> int
    {
      N_Vector* yScor_1d_ptr = yScor_1d.empty() ? nullptr : yScor_1d.data();
      N_Vector* yS_1d_ptr    = yS_1d.empty() ? nullptr : yS_1d.data();

      auto lambda_result = CVodeComputeStateSens(cvode_mem, yScor_1d_ptr,
                                                 yS_1d_ptr);
      return lambda_result;
    };
    auto CVodeComputeStateSens_adapt_solver_memory =
      [&CVodeComputeStateSens_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                           std::vector<N_Vector> yScor_1d,
                                                           std::vector<N_Vector>
                                                             yS_1d) -> int
    {
      auto lambda_result =
        CVodeComputeStateSens_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                          yScor_1d, yS_1d);
      return lambda_result;
    };

    return CVodeComputeStateSens_adapt_solver_memory(cvode_mem, yScor_1d, yS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("yScor_1d"), nb::arg("yS_1d"));

sundials4py::scoped_def(
  m, "CVodeComputeStateSens1",
  [](const CVodeView& cvode_mem, int idx, N_Vector yScor1, N_Vector yS1) -> int
  {
    auto CVodeComputeStateSens1_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int idx, N_Vector yScor1, N_Vector yS1) -> int
    {
      auto lambda_result = CVodeComputeStateSens1(cvode_mem.get(), idx, yScor1,
                                                  yS1);
      return lambda_result;
    };

    return CVodeComputeStateSens1_adapt_solver_memory(cvode_mem, idx, yScor1,
                                                      yS1);
  },
  nb::arg("cvode_mem"), nb::arg("idx"), nb::arg("yScor1"), nb::arg("yS1"));

sundials4py::scoped_def(
  m, "CVodeGetDky",
  [](const CVodeView& cvode_mem, sunrealtype t, int k, N_Vector dky) -> int
  {
    auto CVodeGetDky_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype t, int k, N_Vector dky) -> int
    {
      auto lambda_result = CVodeGetDky(cvode_mem.get(), t, k, dky);
      return lambda_result;
    };

    return CVodeGetDky_adapt_solver_memory(cvode_mem, t, k, dky);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dky"),
  "Dense output function");

sundials4py::scoped_def(
  m, "CVodeGetNumSteps",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumSteps_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nsteps_adapt_modifiable;

      int r = CVodeGetNumSteps(cvode_mem, &nsteps_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable);
    };
    auto CVodeGetNumSteps_adapt_solver_memory =
      [&CVodeGetNumSteps_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumSteps_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumSteps_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumRhsEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfevals_adapt_modifiable;

      int r = CVodeGetNumRhsEvals(cvode_mem, &nfevals_adapt_modifiable);
      return std::make_tuple(r, nfevals_adapt_modifiable);
    };
    auto CVodeGetNumRhsEvals_adapt_solver_memory =
      [&CVodeGetNumRhsEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumRhsEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumLinSolvSetups",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nlinsetups_adapt_modifiable;

      int r = CVodeGetNumLinSolvSetups(cvode_mem, &nlinsetups_adapt_modifiable);
      return std::make_tuple(r, nlinsetups_adapt_modifiable);
    };
    auto CVodeGetNumLinSolvSetups_adapt_solver_memory =
      [&CVodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumLinSolvSetups_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumErrTestFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long netfails_adapt_modifiable;

      int r = CVodeGetNumErrTestFails(cvode_mem, &netfails_adapt_modifiable);
      return std::make_tuple(r, netfails_adapt_modifiable);
    };
    auto CVodeGetNumErrTestFails_adapt_solver_memory =
      [&CVodeGetNumErrTestFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumErrTestFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumErrTestFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetLastOrder",
  [](const CVodeView& cvode_mem) -> std::tuple<int, int>
  {
    auto CVodeGetLastOrder_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, int>
    {
      int qlast_adapt_modifiable;

      int r = CVodeGetLastOrder(cvode_mem, &qlast_adapt_modifiable);
      return std::make_tuple(r, qlast_adapt_modifiable);
    };
    auto CVodeGetLastOrder_adapt_solver_memory =
      [&CVodeGetLastOrder_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, int>
    {
      auto lambda_result =
        CVodeGetLastOrder_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetLastOrder_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetCurrentOrder",
  [](const CVodeView& cvode_mem) -> std::tuple<int, int>
  {
    auto CVodeGetCurrentOrder_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, int>
    {
      int qcur_adapt_modifiable;

      int r = CVodeGetCurrentOrder(cvode_mem, &qcur_adapt_modifiable);
      return std::make_tuple(r, qcur_adapt_modifiable);
    };
    auto CVodeGetCurrentOrder_adapt_solver_memory =
      [&CVodeGetCurrentOrder_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, int>
    {
      auto lambda_result =
        CVodeGetCurrentOrder_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetCurrentOrder_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetCurrentGamma",
  [](const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetCurrentGamma_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype gamma_adapt_modifiable;

      int r = CVodeGetCurrentGamma(cvode_mem, &gamma_adapt_modifiable);
      return std::make_tuple(r, gamma_adapt_modifiable);
    };
    auto CVodeGetCurrentGamma_adapt_solver_memory =
      [&CVodeGetCurrentGamma_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetCurrentGamma_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetCurrentGamma_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumStabLimOrderReds",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumStabLimOrderReds_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nslred_adapt_modifiable;

      int r = CVodeGetNumStabLimOrderReds(cvode_mem, &nslred_adapt_modifiable);
      return std::make_tuple(r, nslred_adapt_modifiable);
    };
    auto CVodeGetNumStabLimOrderReds_adapt_solver_memory =
      [&CVodeGetNumStabLimOrderReds_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumStabLimOrderReds_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumStabLimOrderReds_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetActualInitStep",
  [](const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetActualInitStep_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hinused_adapt_modifiable;

      int r = CVodeGetActualInitStep(cvode_mem, &hinused_adapt_modifiable);
      return std::make_tuple(r, hinused_adapt_modifiable);
    };
    auto CVodeGetActualInitStep_adapt_solver_memory =
      [&CVodeGetActualInitStep_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetActualInitStep_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetActualInitStep_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetLastStep",
  [](const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetLastStep_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hlast_adapt_modifiable;

      int r = CVodeGetLastStep(cvode_mem, &hlast_adapt_modifiable);
      return std::make_tuple(r, hlast_adapt_modifiable);
    };
    auto CVodeGetLastStep_adapt_solver_memory =
      [&CVodeGetLastStep_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetLastStep_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetLastStep_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetCurrentStep",
  [](const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetCurrentStep_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hcur_adapt_modifiable;

      int r = CVodeGetCurrentStep(cvode_mem, &hcur_adapt_modifiable);
      return std::make_tuple(r, hcur_adapt_modifiable);
    };
    auto CVodeGetCurrentStep_adapt_solver_memory =
      [&CVodeGetCurrentStep_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetCurrentStep_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetCurrentStep_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetCurrentState",
  [](const CVodeView& cvode_mem) -> std::tuple<int, N_Vector>
  {
    auto CVodeGetCurrentState_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, N_Vector>
    {
      N_Vector y_adapt_modifiable;

      int r = CVodeGetCurrentState(cvode_mem, &y_adapt_modifiable);
      return std::make_tuple(r, y_adapt_modifiable);
    };
    auto CVodeGetCurrentState_adapt_solver_memory =
      [&CVodeGetCurrentState_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, N_Vector>
    {
      auto lambda_result =
        CVodeGetCurrentState_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetCurrentState_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "CVodeGetCurrentSensSolveIndex",
  [](const CVodeView& cvode_mem) -> std::tuple<int, int>
  {
    auto CVodeGetCurrentSensSolveIndex_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, int>
    {
      int index_adapt_modifiable;

      int r = CVodeGetCurrentSensSolveIndex(cvode_mem, &index_adapt_modifiable);
      return std::make_tuple(r, index_adapt_modifiable);
    };
    auto CVodeGetCurrentSensSolveIndex_adapt_solver_memory =
      [&CVodeGetCurrentSensSolveIndex_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, int>
    {
      auto lambda_result =
        CVodeGetCurrentSensSolveIndex_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetCurrentSensSolveIndex_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetCurrentTime",
  [](const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetCurrentTime_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tcur_adapt_modifiable;

      int r = CVodeGetCurrentTime(cvode_mem, &tcur_adapt_modifiable);
      return std::make_tuple(r, tcur_adapt_modifiable);
    };
    auto CVodeGetCurrentTime_adapt_solver_memory =
      [&CVodeGetCurrentTime_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetCurrentTime_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetCurrentTime_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetTolScaleFactor",
  [](const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetTolScaleFactor_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tolsfac_adapt_modifiable;

      int r = CVodeGetTolScaleFactor(cvode_mem, &tolsfac_adapt_modifiable);
      return std::make_tuple(r, tolsfac_adapt_modifiable);
    };
    auto CVodeGetTolScaleFactor_adapt_solver_memory =
      [&CVodeGetTolScaleFactor_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetTolScaleFactor_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetTolScaleFactor_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetErrWeights",
  [](const CVodeView& cvode_mem, N_Vector eweight) -> int
  {
    auto CVodeGetErrWeights_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                     N_Vector eweight) -> int
    {
      auto lambda_result = CVodeGetErrWeights(cvode_mem.get(), eweight);
      return lambda_result;
    };

    return CVodeGetErrWeights_adapt_solver_memory(cvode_mem, eweight);
  },
  nb::arg("cvode_mem"), nb::arg("eweight"));

sundials4py::scoped_def(
  m, "CVodeGetEstLocalErrors",
  [](const CVodeView& cvode_mem, N_Vector ele) -> int
  {
    auto CVodeGetEstLocalErrors_adapt_solver_memory =
      [](const CVodeView& cvode_mem, N_Vector ele) -> int
    {
      auto lambda_result = CVodeGetEstLocalErrors(cvode_mem.get(), ele);
      return lambda_result;
    };

    return CVodeGetEstLocalErrors_adapt_solver_memory(cvode_mem, ele);
  },
  nb::arg("cvode_mem"), nb::arg("ele"));

sundials4py::scoped_def(
  m, "CVodeGetNumGEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumGEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long ngevals_adapt_modifiable;

      int r = CVodeGetNumGEvals(cvode_mem, &ngevals_adapt_modifiable);
      return std::make_tuple(r, ngevals_adapt_modifiable);
    };
    auto CVodeGetNumGEvals_adapt_solver_memory =
      [&CVodeGetNumGEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumGEvals_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumGEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetRootInfo",
  [](const CVodeView& cvode_mem, sundials4py::IntArray1d rootsfound_1d) -> int
  {
    auto CVodeGetRootInfo_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sundials4py::IntArray1d rootsfound_1d) -> int
    {
      int* rootsfound_1d_ptr = rootsfound_1d.size() == 0 ? nullptr
                                                         : rootsfound_1d.data();

      auto lambda_result = CVodeGetRootInfo(cvode_mem, rootsfound_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetRootInfo_adapt_solver_memory =
      [&CVodeGetRootInfo_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                      sundials4py::IntArray1d
                                                        rootsfound_1d) -> int
    {
      auto lambda_result =
        CVodeGetRootInfo_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                     rootsfound_1d);
      return lambda_result;
    };

    return CVodeGetRootInfo_adapt_solver_memory(cvode_mem, rootsfound_1d);
  },
  nb::arg("cvode_mem"), nb::arg("rootsfound_1d"));

sundials4py::scoped_def(
  m, "CVodeGetIntegratorStats",
  [](const CVodeView& cvode_mem)
    -> std::tuple<int, long, long, long, long, int, int, sunrealtype,
                  sunrealtype, sunrealtype, sunrealtype>
  {
    auto CVodeGetIntegratorStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem)
      -> std::tuple<int, long, long, long, long, int, int, sunrealtype,
                    sunrealtype, sunrealtype, sunrealtype>
    {
      long nsteps_adapt_modifiable;
      long nfevals_adapt_modifiable;
      long nlinsetups_adapt_modifiable;
      long netfails_adapt_modifiable;
      int qlast_adapt_modifiable;
      int qcur_adapt_modifiable;
      sunrealtype hinused_adapt_modifiable;
      sunrealtype hlast_adapt_modifiable;
      sunrealtype hcur_adapt_modifiable;
      sunrealtype tcur_adapt_modifiable;

      int r =
        CVodeGetIntegratorStats(cvode_mem, &nsteps_adapt_modifiable,
                                &nfevals_adapt_modifiable,
                                &nlinsetups_adapt_modifiable,
                                &netfails_adapt_modifiable,
                                &qlast_adapt_modifiable, &qcur_adapt_modifiable,
                                &hinused_adapt_modifiable,
                                &hlast_adapt_modifiable, &hcur_adapt_modifiable,
                                &tcur_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable, nfevals_adapt_modifiable,
                             nlinsetups_adapt_modifiable,
                             netfails_adapt_modifiable, qlast_adapt_modifiable,
                             qcur_adapt_modifiable, hinused_adapt_modifiable,
                             hlast_adapt_modifiable, hcur_adapt_modifiable,
                             tcur_adapt_modifiable);
    };
    auto CVodeGetIntegratorStats_adapt_solver_memory =
      [&CVodeGetIntegratorStats_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem)
      -> std::tuple<int, long, long, long, long, int, int, sunrealtype,
                    sunrealtype, sunrealtype, sunrealtype>
    {
      auto lambda_result =
        CVodeGetIntegratorStats_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetIntegratorStats_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumNonlinSolvIters",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nniters_adapt_modifiable;

      int r = CVodeGetNumNonlinSolvIters(cvode_mem, &nniters_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable);
    };
    auto CVodeGetNumNonlinSolvIters_adapt_solver_memory =
      [&CVodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumNonlinSolvIters_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumNonlinSolvConvFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nnfails_adapt_modifiable;

      int r = CVodeGetNumNonlinSolvConvFails(cvode_mem,
                                             &nnfails_adapt_modifiable);
      return std::make_tuple(r, nnfails_adapt_modifiable);
    };
    auto CVodeGetNumNonlinSolvConvFails_adapt_solver_memory =
      [&CVodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumNonlinSolvConvFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNonlinSolvStats",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long, long>
  {
    auto CVodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long>
    {
      long nniters_adapt_modifiable;
      long nnfails_adapt_modifiable;

      int r = CVodeGetNonlinSolvStats(cvode_mem, &nniters_adapt_modifiable,
                                      &nnfails_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable,
                             nnfails_adapt_modifiable);
    };
    auto CVodeGetNonlinSolvStats_adapt_solver_memory =
      [&CVodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long, long>
    {
      auto lambda_result =
        CVodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNonlinSolvStats_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumStepSolveFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nncfails_adapt_modifiable;

      int r = CVodeGetNumStepSolveFails(cvode_mem, &nncfails_adapt_modifiable);
      return std::make_tuple(r, nncfails_adapt_modifiable);
    };
    auto CVodeGetNumStepSolveFails_adapt_solver_memory =
      [&CVodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumStepSolveFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumConstraintFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumConstraintFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long num_fails_out_adapt_modifiable;

      int r = CVodeGetNumConstraintFails(cvode_mem,
                                         &num_fails_out_adapt_modifiable);
      return std::make_tuple(r, num_fails_out_adapt_modifiable);
    };
    auto CVodeGetNumConstraintFails_adapt_solver_memory =
      [&CVodeGetNumConstraintFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumConstraintFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumConstraintFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumConstraintCorrections",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumConstraintCorrections_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long num_corrections_out_adapt_modifiable;

      int r =
        CVodeGetNumConstraintCorrections(cvode_mem,
                                         &num_corrections_out_adapt_modifiable);
      return std::make_tuple(r, num_corrections_out_adapt_modifiable);
    };
    auto CVodeGetNumConstraintCorrections_adapt_solver_memory =
      [&CVodeGetNumConstraintCorrections_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumConstraintCorrections_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumConstraintCorrections_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodePrintAllStats",
  [](const CVodeView& cvode_mem, FILE* outfile, SUNOutputFormat fmt) -> int
  {
    auto CVodePrintAllStats_adapt_solver_memory =
      [](const CVodeView& cvode_mem, FILE* outfile, SUNOutputFormat fmt) -> int
    {
      auto lambda_result = CVodePrintAllStats(cvode_mem.get(), outfile, fmt);
      return lambda_result;
    };

    return CVodePrintAllStats_adapt_solver_memory(cvode_mem, outfile, fmt);
  },
  nb::arg("cvode_mem"), nb::arg("outfile"), nb::arg("fmt"));

sundials4py::scoped_def(m, "CVodeGetReturnFlagName", CVodeGetReturnFlagName,
                        nb::arg("flag"));

sundials4py::scoped_def(
  m, "CVodeQuadReInit",
  [](const CVodeView& cvode_mem, N_Vector yQ0) -> int
  {
    auto CVodeQuadReInit_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                  N_Vector yQ0) -> int
    {
      auto lambda_result = CVodeQuadReInit(cvode_mem.get(), yQ0);
      return lambda_result;
    };

    return CVodeQuadReInit_adapt_solver_memory(cvode_mem, yQ0);
  },
  nb::arg("cvode_mem"), nb::arg("yQ0"));

sundials4py::scoped_def(
  m, "CVodeQuadSStolerances",
  [](const CVodeView& cvode_mem, sunrealtype reltolQ, sunrealtype abstolQ) -> int
  {
    auto CVodeQuadSStolerances_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype reltolQ,
         sunrealtype abstolQ) -> int
    {
      auto lambda_result = CVodeQuadSStolerances(cvode_mem.get(), reltolQ,
                                                 abstolQ);
      return lambda_result;
    };

    return CVodeQuadSStolerances_adapt_solver_memory(cvode_mem, reltolQ, abstolQ);
  },
  nb::arg("cvode_mem"), nb::arg("reltolQ"), nb::arg("abstolQ"));

sundials4py::scoped_def(
  m, "CVodeQuadSVtolerances",
  [](const CVodeView& cvode_mem, sunrealtype reltolQ, N_Vector abstolQ) -> int
  {
    auto CVodeQuadSVtolerances_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype reltolQ, N_Vector abstolQ) -> int
    {
      auto lambda_result = CVodeQuadSVtolerances(cvode_mem.get(), reltolQ,
                                                 abstolQ);
      return lambda_result;
    };

    return CVodeQuadSVtolerances_adapt_solver_memory(cvode_mem, reltolQ, abstolQ);
  },
  nb::arg("cvode_mem"), nb::arg("reltolQ"), nb::arg("abstolQ"));

sundials4py::scoped_def(
  m, "CVodeSetQuadErrCon",
  [](const CVodeView& cvode_mem, sunbooleantype errconQ) -> int
  {
    auto CVodeSetQuadErrCon_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunbooleantype errconQ) -> int
    {
      auto lambda_result = CVodeSetQuadErrCon(cvode_mem.get(), errconQ);
      return lambda_result;
    };

    return CVodeSetQuadErrCon_adapt_solver_memory(cvode_mem, errconQ);
  },
  nb::arg("cvode_mem"), nb::arg("errconQ"));

sundials4py::scoped_def(
  m, "CVodeGetQuad",
  [](const CVodeView& cvode_mem, N_Vector yQout) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetQuad_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, N_Vector yQout) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetQuad(cvode_mem, &tret_adapt_modifiable, yQout);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto CVodeGetQuad_adapt_solver_memory =
      [&CVodeGetQuad_adapt_modifiable_immutable_to_return](const CVodeView& cvode_mem,
                                                           N_Vector yQout)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetQuad_adapt_modifiable_immutable_to_return(cvode_mem.get(), yQout);
      return lambda_result;
    };

    return CVodeGetQuad_adapt_solver_memory(cvode_mem, yQout);
  },
  nb::arg("cvode_mem"), nb::arg("yQout"));

sundials4py::scoped_def(
  m, "CVodeGetQuadDky",
  [](const CVodeView& cvode_mem, sunrealtype t, int k, N_Vector dky) -> int
  {
    auto CVodeGetQuadDky_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype t, int k, N_Vector dky) -> int
    {
      auto lambda_result = CVodeGetQuadDky(cvode_mem.get(), t, k, dky);
      return lambda_result;
    };

    return CVodeGetQuadDky_adapt_solver_memory(cvode_mem, t, k, dky);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dky"));

sundials4py::scoped_def(
  m, "CVodeGetQuadNumRhsEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfQevals_adapt_modifiable;

      int r = CVodeGetQuadNumRhsEvals(cvode_mem, &nfQevals_adapt_modifiable);
      return std::make_tuple(r, nfQevals_adapt_modifiable);
    };
    auto CVodeGetQuadNumRhsEvals_adapt_solver_memory =
      [&CVodeGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetQuadNumRhsEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetQuadNumErrTestFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nQetfails_adapt_modifiable;

      int r = CVodeGetQuadNumErrTestFails(cvode_mem, &nQetfails_adapt_modifiable);
      return std::make_tuple(r, nQetfails_adapt_modifiable);
    };
    auto CVodeGetQuadNumErrTestFails_adapt_solver_memory =
      [&CVodeGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetQuadNumErrTestFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetQuadErrWeights",
  [](const CVodeView& cvode_mem, N_Vector eQweight) -> int
  {
    auto CVodeGetQuadErrWeights_adapt_solver_memory =
      [](const CVodeView& cvode_mem, N_Vector eQweight) -> int
    {
      auto lambda_result = CVodeGetQuadErrWeights(cvode_mem.get(), eQweight);
      return lambda_result;
    };

    return CVodeGetQuadErrWeights_adapt_solver_memory(cvode_mem, eQweight);
  },
  nb::arg("cvode_mem"), nb::arg("eQweight"));

sundials4py::scoped_def(
  m, "CVodeGetQuadStats",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long, long>
  {
    auto CVodeGetQuadStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long>
    {
      long nfQevals_adapt_modifiable;
      long nQetfails_adapt_modifiable;

      int r = CVodeGetQuadStats(cvode_mem, &nfQevals_adapt_modifiable,
                                &nQetfails_adapt_modifiable);
      return std::make_tuple(r, nfQevals_adapt_modifiable,
                             nQetfails_adapt_modifiable);
    };
    auto CVodeGetQuadStats_adapt_solver_memory =
      [&CVodeGetQuadStats_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long, long>
    {
      auto lambda_result =
        CVodeGetQuadStats_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetQuadStats_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeSensReInit",
  [](const CVodeView& cvode_mem, int ism, std::vector<N_Vector> yS0_1d) -> int
  {
    auto CVodeSensReInit_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, int ism, std::vector<N_Vector> yS0_1d) -> int
    {
      N_Vector* yS0_1d_ptr = yS0_1d.empty() ? nullptr : yS0_1d.data();

      auto lambda_result = CVodeSensReInit(cvode_mem, ism, yS0_1d_ptr);
      return lambda_result;
    };
    auto CVodeSensReInit_adapt_solver_memory =
      [&CVodeSensReInit_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                     int ism,
                                                     std::vector<N_Vector> yS0_1d)
      -> int
    {
      auto lambda_result =
        CVodeSensReInit_adapt_arr_ptr_to_std_vector(cvode_mem.get(), ism, yS0_1d);
      return lambda_result;
    };

    return CVodeSensReInit_adapt_solver_memory(cvode_mem, ism, yS0_1d);
  },
  nb::arg("cvode_mem"), nb::arg("ism"), nb::arg("yS0_1d"));

sundials4py::scoped_def(
  m, "CVodeSensSStolerances",
  [](const CVodeView& cvode_mem, sunrealtype reltolS,
     sundials4py::Array1d abstolS_1d) -> int
  {
    auto CVodeSensSStolerances_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype reltolS,
         sundials4py::Array1d abstolS_1d) -> int
    {
      sunrealtype* abstolS_1d_ptr = abstolS_1d.size() == 0 ? nullptr
                                                           : abstolS_1d.data();

      auto lambda_result = CVodeSensSStolerances(cvode_mem, reltolS,
                                                 abstolS_1d_ptr);
      return lambda_result;
    };
    auto CVodeSensSStolerances_adapt_solver_memory =
      [&CVodeSensSStolerances_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                           sunrealtype reltolS,
                                                           sundials4py::Array1d
                                                             abstolS_1d) -> int
    {
      auto lambda_result =
        CVodeSensSStolerances_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                          reltolS, abstolS_1d);
      return lambda_result;
    };

    return CVodeSensSStolerances_adapt_solver_memory(cvode_mem, reltolS,
                                                     abstolS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("reltolS"), nb::arg("abstolS_1d"));

sundials4py::scoped_def(
  m, "CVodeSensSVtolerances",
  [](const CVodeView& cvode_mem, sunrealtype reltolS,
     std::vector<N_Vector> abstolS_1d) -> int
  {
    auto CVodeSensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype reltolS,
         std::vector<N_Vector> abstolS_1d) -> int
    {
      N_Vector* abstolS_1d_ptr = abstolS_1d.empty() ? nullptr : abstolS_1d.data();

      auto lambda_result = CVodeSensSVtolerances(cvode_mem, reltolS,
                                                 abstolS_1d_ptr);
      return lambda_result;
    };
    auto CVodeSensSVtolerances_adapt_solver_memory =
      [&CVodeSensSVtolerances_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                           sunrealtype reltolS,
                                                           std::vector<N_Vector>
                                                             abstolS_1d) -> int
    {
      auto lambda_result =
        CVodeSensSVtolerances_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                          reltolS, abstolS_1d);
      return lambda_result;
    };

    return CVodeSensSVtolerances_adapt_solver_memory(cvode_mem, reltolS,
                                                     abstolS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("reltolS"), nb::arg("abstolS_1d"));

sundials4py::scoped_def(
  m, "CVodeSensEEtolerances",
  [](const CVodeView& cvode_mem) -> int
  {
    auto CVodeSensEEtolerances_adapt_solver_memory =
      [](const CVodeView& cvode_mem) -> int
    {
      auto lambda_result = CVodeSensEEtolerances(cvode_mem.get());
      return lambda_result;
    };

    return CVodeSensEEtolerances_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeSetSensDQMethod",
  [](const CVodeView& cvode_mem, int DQtype, sunrealtype DQrhomax) -> int
  {
    auto CVodeSetSensDQMethod_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int DQtype, sunrealtype DQrhomax) -> int
    {
      auto lambda_result = CVodeSetSensDQMethod(cvode_mem.get(), DQtype,
                                                DQrhomax);
      return lambda_result;
    };

    return CVodeSetSensDQMethod_adapt_solver_memory(cvode_mem, DQtype, DQrhomax);
  },
  nb::arg("cvode_mem"), nb::arg("DQtype"), nb::arg("DQrhomax"));

sundials4py::scoped_def(
  m, "CVodeSetSensErrCon",
  [](const CVodeView& cvode_mem, sunbooleantype errconS) -> int
  {
    auto CVodeSetSensErrCon_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunbooleantype errconS) -> int
    {
      auto lambda_result = CVodeSetSensErrCon(cvode_mem.get(), errconS);
      return lambda_result;
    };

    return CVodeSetSensErrCon_adapt_solver_memory(cvode_mem, errconS);
  },
  nb::arg("cvode_mem"), nb::arg("errconS"));

sundials4py::scoped_def(
  m, "CVodeSetSensMaxNonlinIters",
  [](const CVodeView& cvode_mem, int maxcorS) -> int
  {
    auto CVodeSetSensMaxNonlinIters_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int maxcorS) -> int
    {
      auto lambda_result = CVodeSetSensMaxNonlinIters(cvode_mem.get(), maxcorS);
      return lambda_result;
    };

    return CVodeSetSensMaxNonlinIters_adapt_solver_memory(cvode_mem, maxcorS);
  },
  nb::arg("cvode_mem"), nb::arg("maxcorS"));

sundials4py::scoped_def(
  m, "CVodeSetSensParams",
  [](const CVodeView& cvode_mem, sundials4py::Array1d p_1d,
     sundials4py::Array1d pbar_1d, std::vector<int> plist_1d) -> int
  {
    auto CVodeSetSensParams_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sundials4py::Array1d p_1d,
         sundials4py::Array1d pbar_1d, std::vector<int> plist_1d) -> int
    {
      sunrealtype* p_1d_ptr    = p_1d.size() == 0 ? nullptr : p_1d.data();
      sunrealtype* pbar_1d_ptr = pbar_1d.size() == 0 ? nullptr : pbar_1d.data();
      int* plist_1d_ptr        = plist_1d.empty() ? nullptr : plist_1d.data();

      auto lambda_result = CVodeSetSensParams(cvode_mem, p_1d_ptr, pbar_1d_ptr,
                                              plist_1d_ptr);
      return lambda_result;
    };
    auto CVodeSetSensParams_adapt_solver_memory =
      [&CVodeSetSensParams_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                        sundials4py::Array1d p_1d,
                                                        sundials4py::Array1d pbar_1d,
                                                        std::vector<int> plist_1d)
      -> int
    {
      auto lambda_result =
        CVodeSetSensParams_adapt_arr_ptr_to_std_vector(cvode_mem.get(), p_1d,
                                                       pbar_1d, plist_1d);
      return lambda_result;
    };

    return CVodeSetSensParams_adapt_solver_memory(cvode_mem, p_1d, pbar_1d,
                                                  plist_1d);
  },
  nb::arg("cvode_mem"), nb::arg("p_1d"), nb::arg("pbar_1d"), nb::arg("plist_1d"));

sundials4py::scoped_def(
  m, "CVodeSetNonlinearSolverSensSim",
  [](const CVodeView& cvode_mem, SUNNonlinearSolver NLS) -> int
  {
    auto CVodeSetNonlinearSolverSensSim_adapt_solver_memory =
      [](const CVodeView& cvode_mem, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = CVodeSetNonlinearSolverSensSim(cvode_mem.get(), NLS);
      return lambda_result;
    };

    return CVodeSetNonlinearSolverSensSim_adapt_solver_memory(cvode_mem, NLS);
  },
  nb::arg("cvode_mem"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "CVodeSetNonlinearSolverSensStg",
  [](const CVodeView& cvode_mem, SUNNonlinearSolver NLS) -> int
  {
    auto CVodeSetNonlinearSolverSensStg_adapt_solver_memory =
      [](const CVodeView& cvode_mem, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = CVodeSetNonlinearSolverSensStg(cvode_mem.get(), NLS);
      return lambda_result;
    };

    return CVodeSetNonlinearSolverSensStg_adapt_solver_memory(cvode_mem, NLS);
  },
  nb::arg("cvode_mem"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "CVodeSetNonlinearSolverSensStg1",
  [](const CVodeView& cvode_mem, SUNNonlinearSolver NLS) -> int
  {
    auto CVodeSetNonlinearSolverSensStg1_adapt_solver_memory =
      [](const CVodeView& cvode_mem, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = CVodeSetNonlinearSolverSensStg1(cvode_mem.get(), NLS);
      return lambda_result;
    };

    return CVodeSetNonlinearSolverSensStg1_adapt_solver_memory(cvode_mem, NLS);
  },
  nb::arg("cvode_mem"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "CVodeSensToggleOff",
  [](const CVodeView& cvode_mem) -> int
  {
    auto CVodeSensToggleOff_adapt_solver_memory =
      [](const CVodeView& cvode_mem) -> int
    {
      auto lambda_result = CVodeSensToggleOff(cvode_mem.get());
      return lambda_result;
    };

    return CVodeSensToggleOff_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetSens",
  [](const CVodeView& cvode_mem,
     std::vector<N_Vector> ySout_1d) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetSens_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype* tret, std::vector<N_Vector> ySout_1d) -> int
    {
      N_Vector* ySout_1d_ptr = ySout_1d.empty() ? nullptr : ySout_1d.data();

      auto lambda_result = CVodeGetSens(cvode_mem, tret, ySout_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetSens_adapt_modifiable_immutable_to_return =
      [&CVodeGetSens_adapt_arr_ptr_to_std_vector](void* cvode_mem,
                                                  std::vector<N_Vector> ySout_1d)
      -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetSens_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                       &tret_adapt_modifiable,
                                                       ySout_1d);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto CVodeGetSens_adapt_solver_memory =
      [&CVodeGetSens_adapt_modifiable_immutable_to_return](const CVodeView& cvode_mem,
                                                           std::vector<N_Vector> ySout_1d)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetSens_adapt_modifiable_immutable_to_return(cvode_mem.get(),
                                                          ySout_1d);
      return lambda_result;
    };

    return CVodeGetSens_adapt_solver_memory(cvode_mem, ySout_1d);
  },
  nb::arg("cvode_mem"), nb::arg("ySout_1d"));

sundials4py::scoped_def(
  m, "CVodeGetSens1",
  [](const CVodeView& cvode_mem, int is,
     N_Vector ySout) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetSens1_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int is, N_Vector ySout) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetSens1(cvode_mem, &tret_adapt_modifiable, is, ySout);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto CVodeGetSens1_adapt_solver_memory =
      [&CVodeGetSens1_adapt_modifiable_immutable_to_return](const CVodeView& cvode_mem,
                                                            int is,
                                                            N_Vector ySout)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetSens1_adapt_modifiable_immutable_to_return(cvode_mem.get(), is,
                                                           ySout);
      return lambda_result;
    };

    return CVodeGetSens1_adapt_solver_memory(cvode_mem, is, ySout);
  },
  nb::arg("cvode_mem"), nb::arg("is_"), nb::arg("ySout"));

sundials4py::scoped_def(
  m, "CVodeGetSensDky",
  [](const CVodeView& cvode_mem, sunrealtype t, int k,
     std::vector<N_Vector> dkyA_1d) -> int
  {
    auto CVodeGetSensDky_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype t, int k,
         std::vector<N_Vector> dkyA_1d) -> int
    {
      N_Vector* dkyA_1d_ptr = dkyA_1d.empty() ? nullptr : dkyA_1d.data();

      auto lambda_result = CVodeGetSensDky(cvode_mem, t, k, dkyA_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetSensDky_adapt_solver_memory =
      [&CVodeGetSensDky_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                     sunrealtype t, int k,
                                                     std::vector<N_Vector> dkyA_1d)
      -> int
    {
      auto lambda_result =
        CVodeGetSensDky_adapt_arr_ptr_to_std_vector(cvode_mem.get(), t, k,
                                                    dkyA_1d);
      return lambda_result;
    };

    return CVodeGetSensDky_adapt_solver_memory(cvode_mem, t, k, dkyA_1d);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyA_1d"));

sundials4py::scoped_def(
  m, "CVodeGetSensDky1",
  [](const CVodeView& cvode_mem, sunrealtype t, int k, int is, N_Vector dky) -> int
  {
    auto CVodeGetSensDky1_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                   sunrealtype t, int k, int is,
                                                   N_Vector dky) -> int
    {
      auto lambda_result = CVodeGetSensDky1(cvode_mem.get(), t, k, is, dky);
      return lambda_result;
    };

    return CVodeGetSensDky1_adapt_solver_memory(cvode_mem, t, k, is, dky);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("is_"),
  nb::arg("dky"));

sundials4py::scoped_def(
  m, "CVodeGetSensNumRhsEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfSevals_adapt_modifiable;

      int r = CVodeGetSensNumRhsEvals(cvode_mem, &nfSevals_adapt_modifiable);
      return std::make_tuple(r, nfSevals_adapt_modifiable);
    };
    auto CVodeGetSensNumRhsEvals_adapt_solver_memory =
      [&CVodeGetSensNumRhsEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetSensNumRhsEvals_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetSensNumRhsEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumRhsEvalsSens",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumRhsEvalsSens_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfevalsS_adapt_modifiable;

      int r = CVodeGetNumRhsEvalsSens(cvode_mem, &nfevalsS_adapt_modifiable);
      return std::make_tuple(r, nfevalsS_adapt_modifiable);
    };
    auto CVodeGetNumRhsEvalsSens_adapt_solver_memory =
      [&CVodeGetNumRhsEvalsSens_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumRhsEvalsSens_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumRhsEvalsSens_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetSensNumErrTestFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nSetfails_adapt_modifiable;

      int r = CVodeGetSensNumErrTestFails(cvode_mem, &nSetfails_adapt_modifiable);
      return std::make_tuple(r, nSetfails_adapt_modifiable);
    };
    auto CVodeGetSensNumErrTestFails_adapt_solver_memory =
      [&CVodeGetSensNumErrTestFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetSensNumErrTestFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetSensNumErrTestFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetSensNumLinSolvSetups",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nlinsetupsS_adapt_modifiable;

      int r = CVodeGetSensNumLinSolvSetups(cvode_mem,
                                           &nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nlinsetupsS_adapt_modifiable);
    };
    auto CVodeGetSensNumLinSolvSetups_adapt_solver_memory =
      [&CVodeGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetSensNumLinSolvSetups_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetSensErrWeights",
  [](const CVodeView& cvode_mem, std::vector<N_Vector> eSweight_1d) -> int
  {
    auto CVodeGetSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<N_Vector> eSweight_1d) -> int
    {
      N_Vector* eSweight_1d_ptr = eSweight_1d.empty() ? nullptr
                                                      : eSweight_1d.data();

      auto lambda_result = CVodeGetSensErrWeights(cvode_mem, eSweight_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetSensErrWeights_adapt_solver_memory =
      [&CVodeGetSensErrWeights_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                            std::vector<N_Vector>
                                                              eSweight_1d) -> int
    {
      auto lambda_result =
        CVodeGetSensErrWeights_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                           eSweight_1d);
      return lambda_result;
    };

    return CVodeGetSensErrWeights_adapt_solver_memory(cvode_mem, eSweight_1d);
  },
  nb::arg("cvode_mem"), nb::arg("eSweight_1d"));

sundials4py::scoped_def(
  m, "CVodeGetSensStats",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long, long, long, long>
  {
    auto CVodeGetSensStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long, long, long>
    {
      long nfSevals_adapt_modifiable;
      long nfevalsS_adapt_modifiable;
      long nSetfails_adapt_modifiable;
      long nlinsetupsS_adapt_modifiable;

      int r = CVodeGetSensStats(cvode_mem, &nfSevals_adapt_modifiable,
                                &nfevalsS_adapt_modifiable,
                                &nSetfails_adapt_modifiable,
                                &nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nfSevals_adapt_modifiable,
                             nfevalsS_adapt_modifiable,
                             nSetfails_adapt_modifiable,
                             nlinsetupsS_adapt_modifiable);
    };
    auto CVodeGetSensStats_adapt_solver_memory =
      [&CVodeGetSensStats_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long, long, long, long>
    {
      auto lambda_result =
        CVodeGetSensStats_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetSensStats_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetSensNumNonlinSolvIters",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nSniters_adapt_modifiable;

      int r = CVodeGetSensNumNonlinSolvIters(cvode_mem,
                                             &nSniters_adapt_modifiable);
      return std::make_tuple(r, nSniters_adapt_modifiable);
    };
    auto CVodeGetSensNumNonlinSolvIters_adapt_solver_memory =
      [&CVodeGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetSensNumNonlinSolvIters_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetSensNumNonlinSolvConvFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nSnfails_adapt_modifiable;

      int r = CVodeGetSensNumNonlinSolvConvFails(cvode_mem,
                                                 &nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSnfails_adapt_modifiable);
    };
    auto CVodeGetSensNumNonlinSolvConvFails_adapt_solver_memory =
      [&CVodeGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetSensNumNonlinSolvConvFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetSensNonlinSolvStats",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long, long>
  {
    auto CVodeGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long>
    {
      long nSniters_adapt_modifiable;
      long nSnfails_adapt_modifiable;

      int r = CVodeGetSensNonlinSolvStats(cvode_mem, &nSniters_adapt_modifiable,
                                          &nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSniters_adapt_modifiable,
                             nSnfails_adapt_modifiable);
    };
    auto CVodeGetSensNonlinSolvStats_adapt_solver_memory =
      [&CVodeGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long, long>
    {
      auto lambda_result =
        CVodeGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetSensNonlinSolvStats_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumStepSensSolveFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nSncfails_adapt_modifiable;

      int r = CVodeGetNumStepSensSolveFails(cvode_mem,
                                            &nSncfails_adapt_modifiable);
      return std::make_tuple(r, nSncfails_adapt_modifiable);
    };
    auto CVodeGetNumStepSensSolveFails_adapt_solver_memory =
      [&CVodeGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumStepSensSolveFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetStgrSensNumNonlinSolvIters",
  [](const CVodeView& cvode_mem, sundials4py::LongArray1d nSTGR1niters_1d) -> int
  {
    auto CVodeGetStgrSensNumNonlinSolvIters_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sundials4py::LongArray1d nSTGR1niters_1d) -> int
    {
      long* nSTGR1niters_1d_ptr =
        nSTGR1niters_1d.size() == 0 ? nullptr : nSTGR1niters_1d.data();

      auto lambda_result =
        CVodeGetStgrSensNumNonlinSolvIters(cvode_mem, nSTGR1niters_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetStgrSensNumNonlinSolvIters_adapt_solver_memory =
      [&CVodeGetStgrSensNumNonlinSolvIters_adapt_arr_ptr_to_std_vector](const CVodeView&
                                                                          cvode_mem,
                                                                        sundials4py::LongArray1d
                                                                          nSTGR1niters_1d)
      -> int
    {
      auto lambda_result =
        CVodeGetStgrSensNumNonlinSolvIters_adapt_arr_ptr_to_std_vector(cvode_mem
                                                                         .get(),
                                                                       nSTGR1niters_1d);
      return lambda_result;
    };

    return CVodeGetStgrSensNumNonlinSolvIters_adapt_solver_memory(cvode_mem,
                                                                  nSTGR1niters_1d);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1niters_1d"));

sundials4py::scoped_def(
  m, "CVodeGetStgrSensNumNonlinSolvConvFails",
  [](const CVodeView& cvode_mem, sundials4py::LongArray1d nSTGR1nfails_1d) -> int
  {
    auto CVodeGetStgrSensNumNonlinSolvConvFails_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sundials4py::LongArray1d nSTGR1nfails_1d) -> int
    {
      long* nSTGR1nfails_1d_ptr =
        nSTGR1nfails_1d.size() == 0 ? nullptr : nSTGR1nfails_1d.data();

      auto lambda_result =
        CVodeGetStgrSensNumNonlinSolvConvFails(cvode_mem, nSTGR1nfails_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetStgrSensNumNonlinSolvConvFails_adapt_solver_memory =
      [&CVodeGetStgrSensNumNonlinSolvConvFails_adapt_arr_ptr_to_std_vector](const CVodeView&
                                                                              cvode_mem,
                                                                            sundials4py::LongArray1d
                                                                              nSTGR1nfails_1d)
      -> int
    {
      auto lambda_result =
        CVodeGetStgrSensNumNonlinSolvConvFails_adapt_arr_ptr_to_std_vector(cvode_mem
                                                                             .get(),
                                                                           nSTGR1nfails_1d);
      return lambda_result;
    };

    return CVodeGetStgrSensNumNonlinSolvConvFails_adapt_solver_memory(cvode_mem,
                                                                      nSTGR1nfails_1d);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1nfails_1d"));

sundials4py::scoped_def(
  m, "CVodeGetStgrSensNonlinSolvStats",
  [](const CVodeView& cvode_mem, sundials4py::LongArray1d nSTGR1niters_1d,
     sundials4py::LongArray1d nSTGR1nfails_1d) -> int
  {
    auto CVodeGetStgrSensNonlinSolvStats_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sundials4py::LongArray1d nSTGR1niters_1d,
         sundials4py::LongArray1d nSTGR1nfails_1d) -> int
    {
      long* nSTGR1niters_1d_ptr =
        nSTGR1niters_1d.size() == 0 ? nullptr : nSTGR1niters_1d.data();
      long* nSTGR1nfails_1d_ptr =
        nSTGR1nfails_1d.size() == 0 ? nullptr : nSTGR1nfails_1d.data();

      auto lambda_result = CVodeGetStgrSensNonlinSolvStats(cvode_mem,
                                                           nSTGR1niters_1d_ptr,
                                                           nSTGR1nfails_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetStgrSensNonlinSolvStats_adapt_solver_memory =
      [&CVodeGetStgrSensNonlinSolvStats_adapt_arr_ptr_to_std_vector](const CVodeView&
                                                                       cvode_mem,
                                                                     sundials4py::LongArray1d
                                                                       nSTGR1niters_1d,
                                                                     sundials4py::LongArray1d
                                                                       nSTGR1nfails_1d)
      -> int
    {
      auto lambda_result =
        CVodeGetStgrSensNonlinSolvStats_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                                    nSTGR1niters_1d,
                                                                    nSTGR1nfails_1d);
      return lambda_result;
    };

    return CVodeGetStgrSensNonlinSolvStats_adapt_solver_memory(cvode_mem,
                                                               nSTGR1niters_1d,
                                                               nSTGR1nfails_1d);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1niters_1d"), nb::arg("nSTGR1nfails_1d"));

sundials4py::scoped_def(
  m, "CVodeGetNumStepStgrSensSolveFails",
  [](const CVodeView& cvode_mem, sundials4py::LongArray1d nSTGR1ncfails_1d) -> int
  {
    auto CVodeGetNumStepStgrSensSolveFails_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sundials4py::LongArray1d nSTGR1ncfails_1d) -> int
    {
      long* nSTGR1ncfails_1d_ptr =
        nSTGR1ncfails_1d.size() == 0 ? nullptr : nSTGR1ncfails_1d.data();

      auto lambda_result =
        CVodeGetNumStepStgrSensSolveFails(cvode_mem, nSTGR1ncfails_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetNumStepStgrSensSolveFails_adapt_solver_memory =
      [&CVodeGetNumStepStgrSensSolveFails_adapt_arr_ptr_to_std_vector](const CVodeView&
                                                                         cvode_mem,
                                                                       sundials4py::LongArray1d
                                                                         nSTGR1ncfails_1d)
      -> int
    {
      auto lambda_result =
        CVodeGetNumStepStgrSensSolveFails_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                                      nSTGR1ncfails_1d);
      return lambda_result;
    };

    return CVodeGetNumStepStgrSensSolveFails_adapt_solver_memory(cvode_mem,
                                                                 nSTGR1ncfails_1d);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1ncfails_1d"));

sundials4py::scoped_def(
  m, "CVodeQuadSensReInit",
  [](const CVodeView& cvode_mem, std::vector<N_Vector> yQS0_1d) -> int
  {
    auto CVodeQuadSensReInit_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<N_Vector> yQS0_1d) -> int
    {
      N_Vector* yQS0_1d_ptr = yQS0_1d.empty() ? nullptr : yQS0_1d.data();

      auto lambda_result = CVodeQuadSensReInit(cvode_mem, yQS0_1d_ptr);
      return lambda_result;
    };
    auto CVodeQuadSensReInit_adapt_solver_memory =
      [&CVodeQuadSensReInit_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                         std::vector<N_Vector>
                                                           yQS0_1d) -> int
    {
      auto lambda_result =
        CVodeQuadSensReInit_adapt_arr_ptr_to_std_vector(cvode_mem.get(), yQS0_1d);
      return lambda_result;
    };

    return CVodeQuadSensReInit_adapt_solver_memory(cvode_mem, yQS0_1d);
  },
  nb::arg("cvode_mem"), nb::arg("yQS0_1d"));

sundials4py::scoped_def(
  m, "CVodeQuadSensSStolerances",
  [](const CVodeView& cvode_mem, sunrealtype reltolQS,
     sundials4py::Array1d abstolQS_1d) -> int
  {
    auto CVodeQuadSensSStolerances_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype reltolQS,
         sundials4py::Array1d abstolQS_1d) -> int
    {
      sunrealtype* abstolQS_1d_ptr =
        abstolQS_1d.size() == 0 ? nullptr : abstolQS_1d.data();

      auto lambda_result = CVodeQuadSensSStolerances(cvode_mem, reltolQS,
                                                     abstolQS_1d_ptr);
      return lambda_result;
    };
    auto CVodeQuadSensSStolerances_adapt_solver_memory =
      [&CVodeQuadSensSStolerances_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                               sunrealtype reltolQS,
                                                               sundials4py::Array1d
                                                                 abstolQS_1d) -> int
    {
      auto lambda_result =
        CVodeQuadSensSStolerances_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                              reltolQS,
                                                              abstolQS_1d);
      return lambda_result;
    };

    return CVodeQuadSensSStolerances_adapt_solver_memory(cvode_mem, reltolQS,
                                                         abstolQS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("reltolQS"), nb::arg("abstolQS_1d"));

sundials4py::scoped_def(
  m, "CVodeQuadSensSVtolerances",
  [](const CVodeView& cvode_mem, sunrealtype reltolQS,
     std::vector<N_Vector> abstolQS_1d) -> int
  {
    auto CVodeQuadSensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype reltolQS,
         std::vector<N_Vector> abstolQS_1d) -> int
    {
      N_Vector* abstolQS_1d_ptr = abstolQS_1d.empty() ? nullptr
                                                      : abstolQS_1d.data();

      auto lambda_result = CVodeQuadSensSVtolerances(cvode_mem, reltolQS,
                                                     abstolQS_1d_ptr);
      return lambda_result;
    };
    auto CVodeQuadSensSVtolerances_adapt_solver_memory =
      [&CVodeQuadSensSVtolerances_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                               sunrealtype reltolQS,
                                                               std::vector<N_Vector>
                                                                 abstolQS_1d) -> int
    {
      auto lambda_result =
        CVodeQuadSensSVtolerances_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                              reltolQS,
                                                              abstolQS_1d);
      return lambda_result;
    };

    return CVodeQuadSensSVtolerances_adapt_solver_memory(cvode_mem, reltolQS,
                                                         abstolQS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("reltolQS"), nb::arg("abstolQS_1d"));

sundials4py::scoped_def(
  m, "CVodeQuadSensEEtolerances",
  [](const CVodeView& cvode_mem) -> int
  {
    auto CVodeQuadSensEEtolerances_adapt_solver_memory =
      [](const CVodeView& cvode_mem) -> int
    {
      auto lambda_result = CVodeQuadSensEEtolerances(cvode_mem.get());
      return lambda_result;
    };

    return CVodeQuadSensEEtolerances_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeSetQuadSensErrCon",
  [](const CVodeView& cvode_mem, sunbooleantype errconQS) -> int
  {
    auto CVodeSetQuadSensErrCon_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunbooleantype errconQS) -> int
    {
      auto lambda_result = CVodeSetQuadSensErrCon(cvode_mem.get(), errconQS);
      return lambda_result;
    };

    return CVodeSetQuadSensErrCon_adapt_solver_memory(cvode_mem, errconQS);
  },
  nb::arg("cvode_mem"), nb::arg("errconQS"));

sundials4py::scoped_def(
  m, "CVodeGetQuadSens",
  [](const CVodeView& cvode_mem,
     std::vector<N_Vector> yQSout_1d) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetQuadSens_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype* tret, std::vector<N_Vector> yQSout_1d) -> int
    {
      N_Vector* yQSout_1d_ptr = yQSout_1d.empty() ? nullptr : yQSout_1d.data();

      auto lambda_result = CVodeGetQuadSens(cvode_mem, tret, yQSout_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetQuadSens_adapt_modifiable_immutable_to_return =
      [&CVodeGetQuadSens_adapt_arr_ptr_to_std_vector](void* cvode_mem,
                                                      std::vector<N_Vector> yQSout_1d)
      -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetQuadSens_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                           &tret_adapt_modifiable,
                                                           yQSout_1d);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto CVodeGetQuadSens_adapt_solver_memory =
      [&CVodeGetQuadSens_adapt_modifiable_immutable_to_return](const CVodeView& cvode_mem,
                                                               std::vector<N_Vector>
                                                                 yQSout_1d)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetQuadSens_adapt_modifiable_immutable_to_return(cvode_mem.get(),
                                                              yQSout_1d);
      return lambda_result;
    };

    return CVodeGetQuadSens_adapt_solver_memory(cvode_mem, yQSout_1d);
  },
  nb::arg("cvode_mem"), nb::arg("yQSout_1d"));

sundials4py::scoped_def(
  m, "CVodeGetQuadSens1",
  [](const CVodeView& cvode_mem, int is,
     N_Vector yQSout) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetQuadSens1_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int is, N_Vector yQSout) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetQuadSens1(cvode_mem, &tret_adapt_modifiable, is, yQSout);
      return std::make_tuple(r, tret_adapt_modifiable);
    };
    auto CVodeGetQuadSens1_adapt_solver_memory =
      [&CVodeGetQuadSens1_adapt_modifiable_immutable_to_return](const CVodeView& cvode_mem,
                                                                int is,
                                                                N_Vector yQSout)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetQuadSens1_adapt_modifiable_immutable_to_return(cvode_mem.get(),
                                                               is, yQSout);
      return lambda_result;
    };

    return CVodeGetQuadSens1_adapt_solver_memory(cvode_mem, is, yQSout);
  },
  nb::arg("cvode_mem"), nb::arg("is_"), nb::arg("yQSout"));

sundials4py::scoped_def(
  m, "CVodeGetQuadSensDky",
  [](const CVodeView& cvode_mem, sunrealtype t, int k,
     std::vector<N_Vector> dkyQS_all_1d) -> int
  {
    auto CVodeGetQuadSensDky_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype t, int k,
         std::vector<N_Vector> dkyQS_all_1d) -> int
    {
      N_Vector* dkyQS_all_1d_ptr = dkyQS_all_1d.empty() ? nullptr
                                                        : dkyQS_all_1d.data();

      auto lambda_result = CVodeGetQuadSensDky(cvode_mem, t, k, dkyQS_all_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetQuadSensDky_adapt_solver_memory =
      [&CVodeGetQuadSensDky_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                         sunrealtype t, int k,
                                                         std::vector<N_Vector>
                                                           dkyQS_all_1d) -> int
    {
      auto lambda_result =
        CVodeGetQuadSensDky_adapt_arr_ptr_to_std_vector(cvode_mem.get(), t, k,
                                                        dkyQS_all_1d);
      return lambda_result;
    };

    return CVodeGetQuadSensDky_adapt_solver_memory(cvode_mem, t, k, dkyQS_all_1d);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyQS_all_1d"));

sundials4py::scoped_def(
  m, "CVodeGetQuadSensDky1",
  [](const CVodeView& cvode_mem, sunrealtype t, int k, int is, N_Vector dkyQS) -> int
  {
    auto CVodeGetQuadSensDky1_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype t, int k, int is,
         N_Vector dkyQS) -> int
    {
      auto lambda_result = CVodeGetQuadSensDky1(cvode_mem.get(), t, k, is, dkyQS);
      return lambda_result;
    };

    return CVodeGetQuadSensDky1_adapt_solver_memory(cvode_mem, t, k, is, dkyQS);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("is_"),
  nb::arg("dkyQS"));

sundials4py::scoped_def(
  m, "CVodeGetQuadSensNumRhsEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfQSevals_adapt_modifiable;

      int r = CVodeGetQuadSensNumRhsEvals(cvode_mem, &nfQSevals_adapt_modifiable);
      return std::make_tuple(r, nfQSevals_adapt_modifiable);
    };
    auto CVodeGetQuadSensNumRhsEvals_adapt_solver_memory =
      [&CVodeGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetQuadSensNumRhsEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetQuadSensNumErrTestFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nQSetfails_adapt_modifiable;

      int r = CVodeGetQuadSensNumErrTestFails(cvode_mem,
                                              &nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nQSetfails_adapt_modifiable);
    };
    auto CVodeGetQuadSensNumErrTestFails_adapt_solver_memory =
      [&CVodeGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetQuadSensNumErrTestFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetQuadSensErrWeights",
  [](const CVodeView& cvode_mem, std::vector<N_Vector> eQSweight_1d) -> int
  {
    auto CVodeGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<N_Vector> eQSweight_1d) -> int
    {
      N_Vector* eQSweight_1d_ptr = eQSweight_1d.empty() ? nullptr
                                                        : eQSweight_1d.data();

      auto lambda_result = CVodeGetQuadSensErrWeights(cvode_mem,
                                                      eQSweight_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetQuadSensErrWeights_adapt_solver_memory =
      [&CVodeGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector](const CVodeView& cvode_mem,
                                                                std::vector<N_Vector>
                                                                  eQSweight_1d) -> int
    {
      auto lambda_result =
        CVodeGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector(cvode_mem.get(),
                                                               eQSweight_1d);
      return lambda_result;
    };

    return CVodeGetQuadSensErrWeights_adapt_solver_memory(cvode_mem,
                                                          eQSweight_1d);
  },
  nb::arg("cvode_mem"), nb::arg("eQSweight_1d"));

sundials4py::scoped_def(
  m, "CVodeGetQuadSensStats",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long, long>
  {
    auto CVodeGetQuadSensStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long>
    {
      long nfQSevals_adapt_modifiable;
      long nQSetfails_adapt_modifiable;

      int r = CVodeGetQuadSensStats(cvode_mem, &nfQSevals_adapt_modifiable,
                                    &nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nfQSevals_adapt_modifiable,
                             nQSetfails_adapt_modifiable);
    };
    auto CVodeGetQuadSensStats_adapt_solver_memory =
      [&CVodeGetQuadSensStats_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long, long>
    {
      auto lambda_result =
        CVodeGetQuadSensStats_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetQuadSensStats_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeAdjInit",
  [](const CVodeView& cvode_mem, long steps, int interp) -> int
  {
    auto CVodeAdjInit_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                               long steps, int interp) -> int
    {
      auto lambda_result = CVodeAdjInit(cvode_mem.get(), steps, interp);
      return lambda_result;
    };

    return CVodeAdjInit_adapt_solver_memory(cvode_mem, steps, interp);
  },
  nb::arg("cvode_mem"), nb::arg("steps"), nb::arg("interp"));

sundials4py::scoped_def(
  m, "CVodeAdjReInit",
  [](const CVodeView& cvode_mem) -> int
  {
    auto CVodeAdjReInit_adapt_solver_memory = [](const CVodeView& cvode_mem) -> int
    {
      auto lambda_result = CVodeAdjReInit(cvode_mem.get());
      return lambda_result;
    };

    return CVodeAdjReInit_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeCreateB",
  [](const CVodeView& cvode_mem, int lmmB) -> std::tuple<int, int>
  {
    auto CVodeCreateB_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int lmmB) -> std::tuple<int, int>
    {
      int which_adapt_modifiable;

      int r = CVodeCreateB(cvode_mem, lmmB, &which_adapt_modifiable);
      return std::make_tuple(r, which_adapt_modifiable);
    };
    auto CVodeCreateB_adapt_solver_memory =
      [&CVodeCreateB_adapt_modifiable_immutable_to_return](const CVodeView& cvode_mem,
                                                           int lmmB)
      -> std::tuple<int, int>
    {
      auto lambda_result =
        CVodeCreateB_adapt_modifiable_immutable_to_return(cvode_mem.get(), lmmB);
      return lambda_result;
    };

    return CVodeCreateB_adapt_solver_memory(cvode_mem, lmmB);
  },
  nb::arg("cvode_mem"), nb::arg("lmmB"));

sundials4py::scoped_def(
  m, "CVodeReInitB",
  [](const CVodeView& cvode_mem, int which, sunrealtype tB0, N_Vector yB0) -> int
  {
    auto CVodeReInitB_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                               int which, sunrealtype tB0,
                                               N_Vector yB0) -> int
    {
      auto lambda_result = CVodeReInitB(cvode_mem.get(), which, tB0, yB0);
      return lambda_result;
    };

    return CVodeReInitB_adapt_solver_memory(cvode_mem, which, tB0, yB0);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("tB0"), nb::arg("yB0"));

sundials4py::scoped_def(
  m, "CVodeSStolerancesB",
  [](const CVodeView& cvode_mem, int which, sunrealtype reltolB,
     sunrealtype abstolB) -> int
  {
    auto CVodeSStolerancesB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunrealtype reltolB,
         sunrealtype abstolB) -> int
    {
      auto lambda_result = CVodeSStolerancesB(cvode_mem.get(), which, reltolB,
                                              abstolB);
      return lambda_result;
    };

    return CVodeSStolerancesB_adapt_solver_memory(cvode_mem, which, reltolB,
                                                  abstolB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("reltolB"), nb::arg("abstolB"));

sundials4py::scoped_def(
  m, "CVodeSVtolerancesB",
  [](const CVodeView& cvode_mem, int which, sunrealtype reltolB,
     N_Vector abstolB) -> int
  {
    auto CVodeSVtolerancesB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunrealtype reltolB,
         N_Vector abstolB) -> int
    {
      auto lambda_result = CVodeSVtolerancesB(cvode_mem.get(), which, reltolB,
                                              abstolB);
      return lambda_result;
    };

    return CVodeSVtolerancesB_adapt_solver_memory(cvode_mem, which, reltolB,
                                                  abstolB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("reltolB"), nb::arg("abstolB"));

sundials4py::scoped_def(
  m, "CVodeQuadReInitB",
  [](const CVodeView& cvode_mem, int which, N_Vector yQB0) -> int
  {
    auto CVodeQuadReInitB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, N_Vector yQB0) -> int
    {
      auto lambda_result = CVodeQuadReInitB(cvode_mem.get(), which, yQB0);
      return lambda_result;
    };

    return CVodeQuadReInitB_adapt_solver_memory(cvode_mem, which, yQB0);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("yQB0"));

sundials4py::scoped_def(
  m, "CVodeQuadSStolerancesB",
  [](const CVodeView& cvode_mem, int which, sunrealtype reltolQB,
     sunrealtype abstolQB) -> int
  {
    auto CVodeQuadSStolerancesB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunrealtype reltolQB,
         sunrealtype abstolQB) -> int
    {
      auto lambda_result = CVodeQuadSStolerancesB(cvode_mem.get(), which,
                                                  reltolQB, abstolQB);
      return lambda_result;
    };

    return CVodeQuadSStolerancesB_adapt_solver_memory(cvode_mem, which,
                                                      reltolQB, abstolQB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("reltolQB"),
  nb::arg("abstolQB"));

sundials4py::scoped_def(
  m, "CVodeQuadSVtolerancesB",
  [](const CVodeView& cvode_mem, int which, sunrealtype reltolQB,
     N_Vector abstolQB) -> int
  {
    auto CVodeQuadSVtolerancesB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunrealtype reltolQB,
         N_Vector abstolQB) -> int
    {
      auto lambda_result = CVodeQuadSVtolerancesB(cvode_mem.get(), which,
                                                  reltolQB, abstolQB);
      return lambda_result;
    };

    return CVodeQuadSVtolerancesB_adapt_solver_memory(cvode_mem, which,
                                                      reltolQB, abstolQB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("reltolQB"),
  nb::arg("abstolQB"));

sundials4py::scoped_def(
  m, "CVodeSetAdjNoSensi",
  [](const CVodeView& cvode_mem) -> int
  {
    auto CVodeSetAdjNoSensi_adapt_solver_memory =
      [](const CVodeView& cvode_mem) -> int
    {
      auto lambda_result = CVodeSetAdjNoSensi(cvode_mem.get());
      return lambda_result;
    };

    return CVodeSetAdjNoSensi_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeSetMaxOrdB",
  [](const CVodeView& cvode_mem, int which, int maxordB) -> int
  {
    auto CVodeSetMaxOrdB_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                  int which, int maxordB) -> int
    {
      auto lambda_result = CVodeSetMaxOrdB(cvode_mem.get(), which, maxordB);
      return lambda_result;
    };

    return CVodeSetMaxOrdB_adapt_solver_memory(cvode_mem, which, maxordB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("maxordB"));

sundials4py::scoped_def(
  m, "CVodeSetMaxNumStepsB",
  [](const CVodeView& cvode_mem, int which, long mxstepsB) -> int
  {
    auto CVodeSetMaxNumStepsB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, long mxstepsB) -> int
    {
      auto lambda_result = CVodeSetMaxNumStepsB(cvode_mem.get(), which, mxstepsB);
      return lambda_result;
    };

    return CVodeSetMaxNumStepsB_adapt_solver_memory(cvode_mem, which, mxstepsB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("mxstepsB"));

sundials4py::scoped_def(
  m, "CVodeSetStabLimDetB",
  [](const CVodeView& cvode_mem, int which, sunbooleantype stldetB) -> int
  {
    auto CVodeSetStabLimDetB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunbooleantype stldetB) -> int
    {
      auto lambda_result = CVodeSetStabLimDetB(cvode_mem.get(), which, stldetB);
      return lambda_result;
    };

    return CVodeSetStabLimDetB_adapt_solver_memory(cvode_mem, which, stldetB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("stldetB"));

sundials4py::scoped_def(
  m, "CVodeSetInitStepB",
  [](const CVodeView& cvode_mem, int which, sunrealtype hinB) -> int
  {
    auto CVodeSetInitStepB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunrealtype hinB) -> int
    {
      auto lambda_result = CVodeSetInitStepB(cvode_mem.get(), which, hinB);
      return lambda_result;
    };

    return CVodeSetInitStepB_adapt_solver_memory(cvode_mem, which, hinB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("hinB"));

sundials4py::scoped_def(
  m, "CVodeSetMinStepB",
  [](const CVodeView& cvode_mem, int which, sunrealtype hminB) -> int
  {
    auto CVodeSetMinStepB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunrealtype hminB) -> int
    {
      auto lambda_result = CVodeSetMinStepB(cvode_mem.get(), which, hminB);
      return lambda_result;
    };

    return CVodeSetMinStepB_adapt_solver_memory(cvode_mem, which, hminB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("hminB"));

sundials4py::scoped_def(
  m, "CVodeSetMaxStepB",
  [](const CVodeView& cvode_mem, int which, sunrealtype hmaxB) -> int
  {
    auto CVodeSetMaxStepB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunrealtype hmaxB) -> int
    {
      auto lambda_result = CVodeSetMaxStepB(cvode_mem.get(), which, hmaxB);
      return lambda_result;
    };

    return CVodeSetMaxStepB_adapt_solver_memory(cvode_mem, which, hmaxB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("hmaxB"));

sundials4py::scoped_def(
  m, "CVodeSetConstraintsB",
  [](const CVodeView& cvode_mem, int which, N_Vector constraintsB) -> int
  {
    auto CVodeSetConstraintsB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, N_Vector constraintsB) -> int
    {
      auto lambda_result = CVodeSetConstraintsB(cvode_mem.get(), which,
                                                constraintsB);
      return lambda_result;
    };

    return CVodeSetConstraintsB_adapt_solver_memory(cvode_mem, which,
                                                    constraintsB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("constraintsB"));

sundials4py::scoped_def(
  m, "CVodeSetQuadErrConB",
  [](const CVodeView& cvode_mem, int which, sunbooleantype errconQB) -> int
  {
    auto CVodeSetQuadErrConB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunbooleantype errconQB) -> int
    {
      auto lambda_result = CVodeSetQuadErrConB(cvode_mem.get(), which, errconQB);
      return lambda_result;
    };

    return CVodeSetQuadErrConB_adapt_solver_memory(cvode_mem, which, errconQB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("errconQB"));

sundials4py::scoped_def(
  m, "CVodeSetNonlinearSolverB",
  [](const CVodeView& cvode_mem, int which, SUNNonlinearSolver NLS) -> int
  {
    auto CVodeSetNonlinearSolverB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = CVodeSetNonlinearSolverB(cvode_mem.get(), which, NLS);
      return lambda_result;
    };

    return CVodeSetNonlinearSolverB_adapt_solver_memory(cvode_mem, which, NLS);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "CVodeGetB",
  [](const CVodeView& cvode_mem, int which,
     N_Vector yB) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetB_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int which, N_Vector yB) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tBret_adapt_modifiable;

      int r = CVodeGetB(cvode_mem, which, &tBret_adapt_modifiable, yB);
      return std::make_tuple(r, tBret_adapt_modifiable);
    };
    auto CVodeGetB_adapt_solver_memory =
      [&CVodeGetB_adapt_modifiable_immutable_to_return](const CVodeView& cvode_mem,
                                                        int which, N_Vector yB)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetB_adapt_modifiable_immutable_to_return(cvode_mem.get(), which,
                                                       yB);
      return lambda_result;
    };

    return CVodeGetB_adapt_solver_memory(cvode_mem, which, yB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("yB"));

sundials4py::scoped_def(
  m, "CVodeGetQuadB",
  [](const CVodeView& cvode_mem, int which,
     N_Vector qB) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetQuadB_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int which, N_Vector qB) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tBret_adapt_modifiable;

      int r = CVodeGetQuadB(cvode_mem, which, &tBret_adapt_modifiable, qB);
      return std::make_tuple(r, tBret_adapt_modifiable);
    };
    auto CVodeGetQuadB_adapt_solver_memory =
      [&CVodeGetQuadB_adapt_modifiable_immutable_to_return](const CVodeView& cvode_mem,
                                                            int which,
                                                            N_Vector qB)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetQuadB_adapt_modifiable_immutable_to_return(cvode_mem.get(),
                                                           which, qB);
      return lambda_result;
    };

    return CVodeGetQuadB_adapt_solver_memory(cvode_mem, which, qB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("qB"));

sundials4py::scoped_def(
  m, "CVodeGetAdjCVodeBmem",
  [](const CVodeView& cvode_mem, int which) -> void*
  {
    auto CVodeGetAdjCVodeBmem_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which) -> void*
    {
      auto lambda_result = CVodeGetAdjCVodeBmem(cvode_mem.get(), which);
      return lambda_result;
    };

    return CVodeGetAdjCVodeBmem_adapt_solver_memory(cvode_mem, which);
  },
  nb::arg("cvode_mem"), nb::arg("which"));

sundials4py::scoped_def(
  m, "CVodeGetAdjY",
  [](const CVodeView& cvode_mem, sunrealtype t, N_Vector y) -> int
  {
    auto CVodeGetAdjY_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                               sunrealtype t, N_Vector y) -> int
    {
      auto lambda_result = CVodeGetAdjY(cvode_mem.get(), t, y);
      return lambda_result;
    };

    return CVodeGetAdjY_adapt_solver_memory(cvode_mem, t, y);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("y"));

sundials4py::scoped_def(
  m, "CVodeGetAdjDataPointHermite",
  [](const CVodeView& cvode_mem, int which,
     std::optional<N_Vector> y  = std::nullopt,
     std::optional<N_Vector> yd = std::nullopt) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetAdjDataPointHermite_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int which, N_Vector y,
         N_Vector yd) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_adapt_modifiable;

      int r = CVodeGetAdjDataPointHermite(cvode_mem, which, &t_adapt_modifiable,
                                          y, yd);
      return std::make_tuple(r, t_adapt_modifiable);
    };
    auto CVodeGetAdjDataPointHermite_adapt_optional_arg_with_default_null =
      [&CVodeGetAdjDataPointHermite_adapt_modifiable_immutable_to_return](void* cvode_mem,
                                                                          int which,
                                                                          std::optional<N_Vector>
                                                                            y =
                                                                              std::nullopt,
                                                                          std::optional<N_Vector>
                                                                            yd =
                                                                              std::nullopt)
      -> std::tuple<int, sunrealtype>
    {
      N_Vector y_adapt_default_null = nullptr;
      if (y.has_value()) y_adapt_default_null = y.value();
      N_Vector yd_adapt_default_null = nullptr;
      if (yd.has_value()) yd_adapt_default_null = yd.value();

      auto lambda_result =
        CVodeGetAdjDataPointHermite_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                         which,
                                                                         y_adapt_default_null,
                                                                         yd_adapt_default_null);
      return lambda_result;
    };
    auto CVodeGetAdjDataPointHermite_adapt_solver_memory =
      [&CVodeGetAdjDataPointHermite_adapt_optional_arg_with_default_null](const CVodeView&
                                                                            cvode_mem,
                                                                          int which,
                                                                          std::optional<N_Vector>
                                                                            y =
                                                                              std::nullopt,
                                                                          std::optional<N_Vector>
                                                                            yd =
                                                                              std::nullopt)
      -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetAdjDataPointHermite_adapt_optional_arg_with_default_null(cvode_mem
                                                                           .get(),
                                                                         which,
                                                                         y, yd);
      return lambda_result;
    };

    return CVodeGetAdjDataPointHermite_adapt_solver_memory(cvode_mem, which, y,
                                                           yd);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("y").none() = nb::none(),
  nb::arg("yd").none() = nb::none());

sundials4py::scoped_def(
  m, "CVodeGetAdjDataPointPolynomial",
  [](const CVodeView& cvode_mem, int which,
     std::optional<N_Vector> y = std::nullopt) -> std::tuple<int, sunrealtype, int>
  {
    auto CVodeGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int which,
         N_Vector y) -> std::tuple<int, sunrealtype, int>
    {
      sunrealtype t_adapt_modifiable;
      int order_adapt_modifiable;

      int r = CVodeGetAdjDataPointPolynomial(cvode_mem, which,
                                             &t_adapt_modifiable,
                                             &order_adapt_modifiable, y);
      return std::make_tuple(r, t_adapt_modifiable, order_adapt_modifiable);
    };
    auto CVodeGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null =
      [&CVodeGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return](void* cvode_mem,
                                                                             int which,
                                                                             std::optional<N_Vector>
                                                                               y =
                                                                                 std::nullopt)
      -> std::tuple<int, sunrealtype, int>
    {
      N_Vector y_adapt_default_null = nullptr;
      if (y.has_value()) y_adapt_default_null = y.value();

      auto lambda_result =
        CVodeGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                            which,
                                                                            y_adapt_default_null);
      return lambda_result;
    };
    auto CVodeGetAdjDataPointPolynomial_adapt_solver_memory =
      [&CVodeGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null](const CVodeView&
                                                                               cvode_mem,
                                                                             int which,
                                                                             std::optional<N_Vector>
                                                                               y =
                                                                                 std::nullopt)
      -> std::tuple<int, sunrealtype, int>
    {
      auto lambda_result =
        CVodeGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null(cvode_mem
                                                                              .get(),
                                                                            which,
                                                                            y);
      return lambda_result;
    };

    return CVodeGetAdjDataPointPolynomial_adapt_solver_memory(cvode_mem, which,
                                                              y);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("y").none() = nb::none());
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _CVSLS_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("CVLS_SUCCESS")         = 0;
m.attr("CVLS_MEM_NULL")        = -1;
m.attr("CVLS_LMEM_NULL")       = -2;
m.attr("CVLS_ILL_INPUT")       = -3;
m.attr("CVLS_MEM_FAIL")        = -4;
m.attr("CVLS_PMEM_NULL")       = -5;
m.attr("CVLS_JACFUNC_UNRECVR") = -6;
m.attr("CVLS_JACFUNC_RECVR")   = -7;
m.attr("CVLS_SUNMAT_FAIL")     = -8;
m.attr("CVLS_SUNLS_FAIL")      = -9;
m.attr("CVLS_NO_ADJ")          = -101;
m.attr("CVLS_LMEMB_NULL")      = -102;

sundials4py::scoped_def(
  m, "CVodeSetLinearSolver",
  [](const CVodeView& cvode_mem, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto CVodeSetLinearSolver_adapt_optional_arg_with_default_null =
      [](void* cvode_mem, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = CVodeSetLinearSolver(cvode_mem, LS,
                                                A_adapt_default_null);
      return lambda_result;
    };
    auto CVodeSetLinearSolver_adapt_solver_memory =
      [&CVodeSetLinearSolver_adapt_optional_arg_with_default_null](const CVodeView&
                                                                     cvode_mem,
                                                                   SUNLinearSolver LS,
                                                                   std::optional<SUNMatrix>
                                                                     A = std::nullopt)
      -> int
    {
      auto lambda_result =
        CVodeSetLinearSolver_adapt_optional_arg_with_default_null(cvode_mem.get(),
                                                                  LS, A);
      return lambda_result;
    };

    return CVodeSetLinearSolver_adapt_solver_memory(cvode_mem, LS, A);
  },
  nb::arg("cvode_mem"), nb::arg("LS"), nb::arg("A").none() = nb::none());

sundials4py::scoped_def(
  m, "CVodeSetJacEvalFrequency",
  [](const CVodeView& cvode_mem, long msbj) -> int
  {
    auto CVodeSetJacEvalFrequency_adapt_solver_memory =
      [](const CVodeView& cvode_mem, long msbj) -> int
    {
      auto lambda_result = CVodeSetJacEvalFrequency(cvode_mem.get(), msbj);
      return lambda_result;
    };

    return CVodeSetJacEvalFrequency_adapt_solver_memory(cvode_mem, msbj);
  },
  nb::arg("cvode_mem"), nb::arg("msbj"));

sundials4py::scoped_def(
  m, "CVodeSetLinearSolutionScaling",
  [](const CVodeView& cvode_mem, sunbooleantype onoff) -> int
  {
    auto CVodeSetLinearSolutionScaling_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunbooleantype onoff) -> int
    {
      auto lambda_result = CVodeSetLinearSolutionScaling(cvode_mem.get(), onoff);
      return lambda_result;
    };

    return CVodeSetLinearSolutionScaling_adapt_solver_memory(cvode_mem, onoff);
  },
  nb::arg("cvode_mem"), nb::arg("onoff"));

sundials4py::scoped_def(
  m, "CVodeSetDeltaGammaMaxBadJac",
  [](const CVodeView& cvode_mem, sunrealtype dgmax_jbad) -> int
  {
    auto CVodeSetDeltaGammaMaxBadJac_adapt_solver_memory =
      [](const CVodeView& cvode_mem, sunrealtype dgmax_jbad) -> int
    {
      auto lambda_result = CVodeSetDeltaGammaMaxBadJac(cvode_mem.get(),
                                                       dgmax_jbad);
      return lambda_result;
    };

    return CVodeSetDeltaGammaMaxBadJac_adapt_solver_memory(cvode_mem, dgmax_jbad);
  },
  nb::arg("cvode_mem"), nb::arg("dgmax_jbad"));

sundials4py::scoped_def(
  m, "CVodeSetEpsLin",
  [](const CVodeView& cvode_mem, sunrealtype eplifac) -> int
  {
    auto CVodeSetEpsLin_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                 sunrealtype eplifac) -> int
    {
      auto lambda_result = CVodeSetEpsLin(cvode_mem.get(), eplifac);
      return lambda_result;
    };

    return CVodeSetEpsLin_adapt_solver_memory(cvode_mem, eplifac);
  },
  nb::arg("cvode_mem"), nb::arg("eplifac"));

sundials4py::scoped_def(m, "CVodeSetLSNormFactor", CVodeSetLSNormFactor,
                        nb::arg("arkode_mem"), nb::arg("nrmfac"));

sundials4py::scoped_def(
  m, "CVodeGetJac",
  [](const CVodeView& cvode_mem) -> std::tuple<int, SUNMatrix>
  {
    auto CVodeGetJac_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix J_adapt_modifiable;

      int r = CVodeGetJac(cvode_mem, &J_adapt_modifiable);
      return std::make_tuple(r, J_adapt_modifiable);
    };
    auto CVodeGetJac_adapt_solver_memory =
      [&CVodeGetJac_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, SUNMatrix>
    {
      auto lambda_result =
        CVodeGetJac_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetJac_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "CVodeGetJacTime",
  [](const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetJacTime_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_J_adapt_modifiable;

      int r = CVodeGetJacTime(cvode_mem, &t_J_adapt_modifiable);
      return std::make_tuple(r, t_J_adapt_modifiable);
    };
    auto CVodeGetJacTime_adapt_solver_memory =
      [&CVodeGetJacTime_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        CVodeGetJacTime_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetJacTime_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetJacNumSteps",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetJacNumSteps_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nst_J_adapt_modifiable;

      int r = CVodeGetJacNumSteps(cvode_mem, &nst_J_adapt_modifiable);
      return std::make_tuple(r, nst_J_adapt_modifiable);
    };
    auto CVodeGetJacNumSteps_adapt_solver_memory =
      [&CVodeGetJacNumSteps_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetJacNumSteps_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetJacNumSteps_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumJacEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long njevals_adapt_modifiable;

      int r = CVodeGetNumJacEvals(cvode_mem, &njevals_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable);
    };
    auto CVodeGetNumJacEvals_adapt_solver_memory =
      [&CVodeGetNumJacEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumJacEvals_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumJacEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumPrecEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long npevals_adapt_modifiable;

      int r = CVodeGetNumPrecEvals(cvode_mem, &npevals_adapt_modifiable);
      return std::make_tuple(r, npevals_adapt_modifiable);
    };
    auto CVodeGetNumPrecEvals_adapt_solver_memory =
      [&CVodeGetNumPrecEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumPrecEvals_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumPrecEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumPrecSolves",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long npsolves_adapt_modifiable;

      int r = CVodeGetNumPrecSolves(cvode_mem, &npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves_adapt_modifiable);
    };
    auto CVodeGetNumPrecSolves_adapt_solver_memory =
      [&CVodeGetNumPrecSolves_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumPrecSolves_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumPrecSolves_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumLinIters",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nliters_adapt_modifiable;

      int r = CVodeGetNumLinIters(cvode_mem, &nliters_adapt_modifiable);
      return std::make_tuple(r, nliters_adapt_modifiable);
    };
    auto CVodeGetNumLinIters_adapt_solver_memory =
      [&CVodeGetNumLinIters_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumLinIters_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumLinIters_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumLinConvFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nlcfails_adapt_modifiable;

      int r = CVodeGetNumLinConvFails(cvode_mem, &nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails_adapt_modifiable);
    };
    auto CVodeGetNumLinConvFails_adapt_solver_memory =
      [&CVodeGetNumLinConvFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumLinConvFails_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumLinConvFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumJTSetupEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long njtsetups_adapt_modifiable;

      int r = CVodeGetNumJTSetupEvals(cvode_mem, &njtsetups_adapt_modifiable);
      return std::make_tuple(r, njtsetups_adapt_modifiable);
    };
    auto CVodeGetNumJTSetupEvals_adapt_solver_memory =
      [&CVodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumJTSetupEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumJtimesEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long njvevals_adapt_modifiable;

      int r = CVodeGetNumJtimesEvals(cvode_mem, &njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals_adapt_modifiable);
    };
    auto CVodeGetNumJtimesEvals_adapt_solver_memory =
      [&CVodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumJtimesEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumLinRhsEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfevalsLS_adapt_modifiable;

      int r = CVodeGetNumLinRhsEvals(cvode_mem, &nfevalsLS_adapt_modifiable);
      return std::make_tuple(r, nfevalsLS_adapt_modifiable);
    };
    auto CVodeGetNumLinRhsEvals_adapt_solver_memory =
      [&CVodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumLinRhsEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetLinSolveStats",
  [](const CVodeView& cvode_mem)
    -> std::tuple<int, long, long, long, long, long, long, long, long>
  {
    auto CVodeGetLinSolveStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem)
      -> std::tuple<int, long, long, long, long, long, long, long, long>
    {
      long njevals_adapt_modifiable;
      long nfevalsLS_adapt_modifiable;
      long nliters_adapt_modifiable;
      long nlcfails_adapt_modifiable;
      long npevals_adapt_modifiable;
      long npsolves_adapt_modifiable;
      long njtsetups_adapt_modifiable;
      long njtimes_adapt_modifiable;

      int r = CVodeGetLinSolveStats(cvode_mem, &njevals_adapt_modifiable,
                                    &nfevalsLS_adapt_modifiable,
                                    &nliters_adapt_modifiable,
                                    &nlcfails_adapt_modifiable,
                                    &npevals_adapt_modifiable,
                                    &npsolves_adapt_modifiable,
                                    &njtsetups_adapt_modifiable,
                                    &njtimes_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable,
                             nfevalsLS_adapt_modifiable,
                             nliters_adapt_modifiable, nlcfails_adapt_modifiable,
                             npevals_adapt_modifiable, npsolves_adapt_modifiable,
                             njtsetups_adapt_modifiable,
                             njtimes_adapt_modifiable);
    };
    auto CVodeGetLinSolveStats_adapt_solver_memory =
      [&CVodeGetLinSolveStats_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem)
      -> std::tuple<int, long, long, long, long, long, long, long, long>
    {
      auto lambda_result =
        CVodeGetLinSolveStats_adapt_modifiable_immutable_to_return(
          cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetLinSolveStats_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetLastLinFlag",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = CVodeGetLastLinFlag(cvode_mem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };
    auto CVodeGetLastLinFlag_adapt_solver_memory =
      [&CVodeGetLastLinFlag_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetLastLinFlag_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetLastLinFlag_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(m, "CVodeGetLinReturnFlagName",
                        CVodeGetLinReturnFlagName, nb::arg("flag"));

sundials4py::scoped_def(
  m, "CVodeSetLinearSolverB",
  [](const CVodeView& cvode_mem, int which, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto CVodeSetLinearSolverB_adapt_optional_arg_with_default_null =
      [](void* cvode_mem, int which, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = CVodeSetLinearSolverB(cvode_mem, which, LS,
                                                 A_adapt_default_null);
      return lambda_result;
    };
    auto CVodeSetLinearSolverB_adapt_solver_memory =
      [&CVodeSetLinearSolverB_adapt_optional_arg_with_default_null](const CVodeView&
                                                                      cvode_mem,
                                                                    int which,
                                                                    SUNLinearSolver LS,
                                                                    std::optional<SUNMatrix>
                                                                      A = std::nullopt)
      -> int
    {
      auto lambda_result =
        CVodeSetLinearSolverB_adapt_optional_arg_with_default_null(cvode_mem.get(),
                                                                   which, LS, A);
      return lambda_result;
    };

    return CVodeSetLinearSolverB_adapt_solver_memory(cvode_mem, which, LS, A);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("LS"),
  nb::arg("A").none() = nb::none());

sundials4py::scoped_def(
  m, "CVodeSetEpsLinB",
  [](const CVodeView& cvode_mem, int which, sunrealtype eplifacB) -> int
  {
    auto CVodeSetEpsLinB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunrealtype eplifacB) -> int
    {
      auto lambda_result = CVodeSetEpsLinB(cvode_mem.get(), which, eplifacB);
      return lambda_result;
    };

    return CVodeSetEpsLinB_adapt_solver_memory(cvode_mem, which, eplifacB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("eplifacB"));

sundials4py::scoped_def(m, "CVodeSetLSNormFactorB", CVodeSetLSNormFactorB,
                        nb::arg("arkode_mem"), nb::arg("which"),
                        nb::arg("nrmfacB"));

sundials4py::scoped_def(
  m, "CVodeSetLinearSolutionScalingB",
  [](const CVodeView& cvode_mem, int which, sunbooleantype onoffB) -> int
  {
    auto CVodeSetLinearSolutionScalingB_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int which, sunbooleantype onoffB) -> int
    {
      auto lambda_result = CVodeSetLinearSolutionScalingB(cvode_mem.get(),
                                                          which, onoffB);
      return lambda_result;
    };

    return CVodeSetLinearSolutionScalingB_adapt_solver_memory(cvode_mem, which,
                                                              onoffB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("onoffB"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _CVPROJ_H
//
// #ifdef __cplusplus
// #endif
//

sundials4py::scoped_def(
  m, "CVodeSetProjErrEst",
  [](const CVodeView& cvode_mem, sunbooleantype onoff) -> int
  {
    auto CVodeSetProjErrEst_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                     sunbooleantype onoff) -> int
    {
      auto lambda_result = CVodeSetProjErrEst(cvode_mem.get(), onoff);
      return lambda_result;
    };

    return CVodeSetProjErrEst_adapt_solver_memory(cvode_mem, onoff);
  },
  nb::arg("cvode_mem"), nb::arg("onoff"));

sundials4py::scoped_def(
  m, "CVodeSetProjFrequency",
  [](const CVodeView& cvode_mem, long proj_freq) -> int
  {
    auto CVodeSetProjFrequency_adapt_solver_memory =
      [](const CVodeView& cvode_mem, long proj_freq) -> int
    {
      auto lambda_result = CVodeSetProjFrequency(cvode_mem.get(), proj_freq);
      return lambda_result;
    };

    return CVodeSetProjFrequency_adapt_solver_memory(cvode_mem, proj_freq);
  },
  nb::arg("cvode_mem"), nb::arg("proj_freq"));

sundials4py::scoped_def(
  m, "CVodeSetMaxNumProjFails",
  [](const CVodeView& cvode_mem, int max_fails) -> int
  {
    auto CVodeSetMaxNumProjFails_adapt_solver_memory =
      [](const CVodeView& cvode_mem, int max_fails) -> int
    {
      auto lambda_result = CVodeSetMaxNumProjFails(cvode_mem.get(), max_fails);
      return lambda_result;
    };

    return CVodeSetMaxNumProjFails_adapt_solver_memory(cvode_mem, max_fails);
  },
  nb::arg("cvode_mem"), nb::arg("max_fails"));

sundials4py::scoped_def(
  m, "CVodeSetEpsProj",
  [](const CVodeView& cvode_mem, sunrealtype eps) -> int
  {
    auto CVodeSetEpsProj_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                  sunrealtype eps) -> int
    {
      auto lambda_result = CVodeSetEpsProj(cvode_mem.get(), eps);
      return lambda_result;
    };

    return CVodeSetEpsProj_adapt_solver_memory(cvode_mem, eps);
  },
  nb::arg("cvode_mem"), nb::arg("eps"));

sundials4py::scoped_def(
  m, "CVodeSetProjFailEta",
  [](const CVodeView& cvode_mem, sunrealtype eta) -> int
  {
    auto CVodeSetProjFailEta_adapt_solver_memory = [](const CVodeView& cvode_mem,
                                                      sunrealtype eta) -> int
    {
      auto lambda_result = CVodeSetProjFailEta(cvode_mem.get(), eta);
      return lambda_result;
    };

    return CVodeSetProjFailEta_adapt_solver_memory(cvode_mem, eta);
  },
  nb::arg("cvode_mem"), nb::arg("eta"));

sundials4py::scoped_def(
  m, "CVodeGetNumProjEvals",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumProjEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nproj_adapt_modifiable;

      int r = CVodeGetNumProjEvals(cvode_mem, &nproj_adapt_modifiable);
      return std::make_tuple(r, nproj_adapt_modifiable);
    };
    auto CVodeGetNumProjEvals_adapt_solver_memory =
      [&CVodeGetNumProjEvals_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumProjEvals_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumProjEvals_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));

sundials4py::scoped_def(
  m, "CVodeGetNumProjFails",
  [](const CVodeView& cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumProjFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nprf_adapt_modifiable;

      int r = CVodeGetNumProjFails(cvode_mem, &nprf_adapt_modifiable);
      return std::make_tuple(r, nprf_adapt_modifiable);
    };
    auto CVodeGetNumProjFails_adapt_solver_memory =
      [&CVodeGetNumProjFails_adapt_modifiable_immutable_to_return](
        const CVodeView& cvode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        CVodeGetNumProjFails_adapt_modifiable_immutable_to_return(cvode_mem.get());
      return lambda_result;
    };

    return CVodeGetNumProjFails_adapt_solver_memory(cvode_mem);
  },
  nb::arg("cvode_mem"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
