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

m.def("CVodeCreate", CVodeCreate, nb::arg("lmm"), nb::arg("sunctx"));

m.def("CVodeReInit", CVodeReInit, nb::arg("cvode_mem"), nb::arg("t0"),
      nb::arg("y0"));

m.def("CVodeSStolerances", CVodeSStolerances, nb::arg("cvode_mem"),
      nb::arg("reltol"), nb::arg("abstol"));

m.def("CVodeSVtolerances", CVodeSVtolerances, nb::arg("cvode_mem"),
      nb::arg("reltol"), nb::arg("abstol"));

m.def("CVodeWFtolerances", CVodeWFtolerances, nb::arg("cvode_mem"),
      nb::arg("efun"));

m.def("CVodeSetConstraints", CVodeSetConstraints, nb::arg("cvode_mem"),
      nb::arg("constraints"));

m.def("CVodeSetDeltaGammaMaxLSetup", CVodeSetDeltaGammaMaxLSetup,
      nb::arg("cvode_mem"), nb::arg("dgmax_lsetup"));

m.def("CVodeSetInitStep", CVodeSetInitStep, nb::arg("cvode_mem"), nb::arg("hin"));

m.def("CVodeSetLSetupFrequency", CVodeSetLSetupFrequency, nb::arg("cvode_mem"),
      nb::arg("msbp"));

m.def("CVodeSetMaxConvFails", CVodeSetMaxConvFails, nb::arg("cvode_mem"),
      nb::arg("maxncf"));

m.def("CVodeSetMaxErrTestFails", CVodeSetMaxErrTestFails, nb::arg("cvode_mem"),
      nb::arg("maxnef"));

m.def("CVodeSetMaxHnilWarns", CVodeSetMaxHnilWarns, nb::arg("cvode_mem"),
      nb::arg("mxhnil"));

m.def("CVodeSetMaxNonlinIters", CVodeSetMaxNonlinIters, nb::arg("cvode_mem"),
      nb::arg("maxcor"));

m.def("CVodeSetMaxNumSteps", CVodeSetMaxNumSteps, nb::arg("cvode_mem"),
      nb::arg("mxsteps"));

m.def("CVodeSetMaxOrd", CVodeSetMaxOrd, nb::arg("cvode_mem"), nb::arg("maxord"));

m.def("CVodeSetMaxStep", CVodeSetMaxStep, nb::arg("cvode_mem"), nb::arg("hmax"));

m.def("CVodeSetMinStep", CVodeSetMinStep, nb::arg("cvode_mem"), nb::arg("hmin"));

m.def("CVodeSetMonitorFrequency", CVodeSetMonitorFrequency,
      nb::arg("cvode_mem"), nb::arg("nst"));

m.def("CVodeSetNonlinConvCoef", CVodeSetNonlinConvCoef, nb::arg("cvode_mem"),
      nb::arg("nlscoef"));

m.def("CVodeSetNonlinearSolver", CVodeSetNonlinearSolver, nb::arg("cvode_mem"),
      nb::arg("NLS"));

m.def("CVodeSetStabLimDet", CVodeSetStabLimDet, nb::arg("cvode_mem"),
      nb::arg("stldet"));

m.def("CVodeSetStopTime", CVodeSetStopTime, nb::arg("cvode_mem"),
      nb::arg("tstop"));

m.def("CVodeSetInterpolateStopTime", CVodeSetInterpolateStopTime,
      nb::arg("cvode_mem"), nb::arg("interp"));

m.def("CVodeClearStopTime", CVodeClearStopTime, nb::arg("cvode_mem"));

m.def("CVodeSetOwnUserData", CVodeSetOwnUserData, nb::arg("cvode_mem"),
      nb::arg("own_user_data"));

m.def("CVodeSetEtaFixedStepBounds", CVodeSetEtaFixedStepBounds,
      nb::arg("cvode_mem"), nb::arg("eta_min_fx"), nb::arg("eta_max_fx"));

m.def("CVodeSetEtaMaxFirstStep", CVodeSetEtaMaxFirstStep, nb::arg("cvode_mem"),
      nb::arg("eta_max_fs"));

m.def("CVodeSetEtaMaxEarlyStep", CVodeSetEtaMaxEarlyStep, nb::arg("cvode_mem"),
      nb::arg("eta_max_es"));

m.def("CVodeSetNumStepsEtaMaxEarlyStep", CVodeSetNumStepsEtaMaxEarlyStep,
      nb::arg("cvode_mem"), nb::arg("small_nst"));

m.def("CVodeSetEtaMax", CVodeSetEtaMax, nb::arg("cvode_mem"),
      nb::arg("eta_max_gs"));

m.def("CVodeSetEtaMin", CVodeSetEtaMin, nb::arg("cvode_mem"), nb::arg("eta_min"));

m.def("CVodeSetEtaMinErrFail", CVodeSetEtaMinErrFail, nb::arg("cvode_mem"),
      nb::arg("eta_min_ef"));

m.def("CVodeSetEtaMaxErrFail", CVodeSetEtaMaxErrFail, nb::arg("cvode_mem"),
      nb::arg("eta_max_ef"));

m.def("CVodeSetNumFailsEtaMaxErrFail", CVodeSetNumFailsEtaMaxErrFail,
      nb::arg("cvode_mem"), nb::arg("small_nef"));

m.def("CVodeSetEtaConvFail", CVodeSetEtaConvFail, nb::arg("cvode_mem"),
      nb::arg("eta_cf"));

m.def("CVodeRootInit", CVodeRootInit, nb::arg("cvode_mem"), nb::arg("nrtfn"),
      nb::arg("g"));

m.def(
  "CVodeSetRootDirection",
  [](void* cvode_mem, int rootdir) -> std::tuple<int, int>
  {
    auto CVodeSetRootDirection_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int rootdir) -> std::tuple<int, int>
    {
      int* rootdir_adapt_modifiable = &rootdir;

      int r = CVodeSetRootDirection(cvode_mem, rootdir_adapt_modifiable);
      return std::make_tuple(r, rootdir);
    };

    return CVodeSetRootDirection_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                      rootdir);
  },
  nb::arg("cvode_mem"), nb::arg("rootdir"));

m.def("CVodeSetNoInactiveRootWarn", CVodeSetNoInactiveRootWarn,
      nb::arg("cvode_mem"));

m.def("CVode", CVode, nb::arg("cvode_mem"), nb::arg("tout"), nb::arg("yout"),
      nb::arg("tret"), nb::arg("itask"));

m.def("CVodeComputeState", CVodeComputeState, nb::arg("cvode_mem"),
      nb::arg("ycor"), nb::arg("y"));

m.def("CVodeGetDky", CVodeGetDky, nb::arg("cvode_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("dky"));

m.def(
  "CVodeGetNumSteps",
  [](void* cvode_mem, long nsteps) -> std::tuple<int, long>
  {
    auto CVodeGetNumSteps_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nsteps) -> std::tuple<int, long>
    {
      long* nsteps_adapt_modifiable = &nsteps;

      int r = CVodeGetNumSteps(cvode_mem, nsteps_adapt_modifiable);
      return std::make_tuple(r, nsteps);
    };

    return CVodeGetNumSteps_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                 nsteps);
  },
  nb::arg("cvode_mem"), nb::arg("nsteps"));

m.def(
  "CVodeGetNumRhsEvals",
  [](void* cvode_mem, long nfevals) -> std::tuple<int, long>
  {
    auto CVodeGetNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nfevals) -> std::tuple<int, long>
    {
      long* nfevals_adapt_modifiable = &nfevals;

      int r = CVodeGetNumRhsEvals(cvode_mem, nfevals_adapt_modifiable);
      return std::make_tuple(r, nfevals);
    };

    return CVodeGetNumRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                    nfevals);
  },
  nb::arg("cvode_mem"), nb::arg("nfevals"));

m.def(
  "CVodeGetNumLinSolvSetups",
  [](void* cvode_mem, long nlinsetups) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nlinsetups) -> std::tuple<int, long>
    {
      long* nlinsetups_adapt_modifiable = &nlinsetups;

      int r = CVodeGetNumLinSolvSetups(cvode_mem, nlinsetups_adapt_modifiable);
      return std::make_tuple(r, nlinsetups);
    };

    return CVodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                         nlinsetups);
  },
  nb::arg("cvode_mem"), nb::arg("nlinsetups"));

m.def(
  "CVodeGetNumErrTestFails",
  [](void* cvode_mem, long netfails) -> std::tuple<int, long>
  {
    auto CVodeGetNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long netfails) -> std::tuple<int, long>
    {
      long* netfails_adapt_modifiable = &netfails;

      int r = CVodeGetNumErrTestFails(cvode_mem, netfails_adapt_modifiable);
      return std::make_tuple(r, netfails);
    };

    return CVodeGetNumErrTestFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                        netfails);
  },
  nb::arg("cvode_mem"), nb::arg("netfails"));

m.def(
  "CVodeGetLastOrder",
  [](void* cvode_mem, int qlast) -> std::tuple<int, int>
  {
    auto CVodeGetLastOrder_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int qlast) -> std::tuple<int, int>
    {
      int* qlast_adapt_modifiable = &qlast;

      int r = CVodeGetLastOrder(cvode_mem, qlast_adapt_modifiable);
      return std::make_tuple(r, qlast);
    };

    return CVodeGetLastOrder_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                  qlast);
  },
  nb::arg("cvode_mem"), nb::arg("qlast"));

m.def(
  "CVodeGetCurrentOrder",
  [](void* cvode_mem, int qcur) -> std::tuple<int, int>
  {
    auto CVodeGetCurrentOrder_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int qcur) -> std::tuple<int, int>
    {
      int* qcur_adapt_modifiable = &qcur;

      int r = CVodeGetCurrentOrder(cvode_mem, qcur_adapt_modifiable);
      return std::make_tuple(r, qcur);
    };

    return CVodeGetCurrentOrder_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                     qcur);
  },
  nb::arg("cvode_mem"), nb::arg("qcur"));

m.def("CVodeGetCurrentGamma", CVodeGetCurrentGamma, nb::arg("cvode_mem"),
      nb::arg("gamma"));

m.def(
  "CVodeGetNumStabLimOrderReds",
  [](void* cvode_mem, long nslred) -> std::tuple<int, long>
  {
    auto CVodeGetNumStabLimOrderReds_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nslred) -> std::tuple<int, long>
    {
      long* nslred_adapt_modifiable = &nslred;

      int r = CVodeGetNumStabLimOrderReds(cvode_mem, nslred_adapt_modifiable);
      return std::make_tuple(r, nslred);
    };

    return CVodeGetNumStabLimOrderReds_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                            nslred);
  },
  nb::arg("cvode_mem"), nb::arg("nslred"));

m.def("CVodeGetActualInitStep", CVodeGetActualInitStep, nb::arg("cvode_mem"),
      nb::arg("hinused"));

m.def("CVodeGetLastStep", CVodeGetLastStep, nb::arg("cvode_mem"),
      nb::arg("hlast"));

m.def("CVodeGetCurrentStep", CVodeGetCurrentStep, nb::arg("cvode_mem"),
      nb::arg("hcur"));

m.def(
  "CVodeGetCurrentState",
  [](void* cvode_mem, std::vector<N_Vector> y) -> int
  {
    auto CVodeGetCurrentState_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, std::vector<N_Vector> y) -> int
    {
      N_Vector* y_ptr = y.empty() ? nullptr : y.data();

      auto lambda_result = CVodeGetCurrentState(cvode_mem, y_ptr);
      return lambda_result;
    };

    return CVodeGetCurrentState_adapt_nvector_ptr_to_vector(cvode_mem, y);
  },
  nb::arg("cvode_mem"), nb::arg("y"));

m.def(
  "CVodeGetCurrentSensSolveIndex",
  [](void* cvode_mem, int index) -> std::tuple<int, int>
  {
    auto CVodeGetCurrentSensSolveIndex_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int index) -> std::tuple<int, int>
    {
      int* index_adapt_modifiable = &index;

      int r = CVodeGetCurrentSensSolveIndex(cvode_mem, index_adapt_modifiable);
      return std::make_tuple(r, index);
    };

    return CVodeGetCurrentSensSolveIndex_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                              index);
  },
  nb::arg("cvode_mem"), nb::arg("index"));

m.def("CVodeGetCurrentTime", CVodeGetCurrentTime, nb::arg("cvode_mem"),
      nb::arg("tcur"));

m.def("CVodeGetTolScaleFactor", CVodeGetTolScaleFactor, nb::arg("cvode_mem"),
      nb::arg("tolsfac"));

m.def("CVodeGetErrWeights", CVodeGetErrWeights, nb::arg("cvode_mem"),
      nb::arg("eweight"));

m.def("CVodeGetEstLocalErrors", CVodeGetEstLocalErrors, nb::arg("cvode_mem"),
      nb::arg("ele"));

m.def(
  "CVodeGetNumGEvals",
  [](void* cvode_mem, long ngevals) -> std::tuple<int, long>
  {
    auto CVodeGetNumGEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long ngevals) -> std::tuple<int, long>
    {
      long* ngevals_adapt_modifiable = &ngevals;

      int r = CVodeGetNumGEvals(cvode_mem, ngevals_adapt_modifiable);
      return std::make_tuple(r, ngevals);
    };

    return CVodeGetNumGEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                  ngevals);
  },
  nb::arg("cvode_mem"), nb::arg("ngevals"));

m.def(
  "CVodeGetRootInfo",
  [](void* cvode_mem, int rootsfound) -> std::tuple<int, int>
  {
    auto CVodeGetRootInfo_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int rootsfound) -> std::tuple<int, int>
    {
      int* rootsfound_adapt_modifiable = &rootsfound;

      int r = CVodeGetRootInfo(cvode_mem, rootsfound_adapt_modifiable);
      return std::make_tuple(r, rootsfound);
    };

    return CVodeGetRootInfo_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                 rootsfound);
  },
  nb::arg("cvode_mem"), nb::arg("rootsfound"));

m.def(
  "CVodeGetIntegratorStats",
  [](void* cvode_mem, long nsteps, long nfevals, long nlinsetups, long netfails,
     int qlast, int qcur, sunrealtype* hinused, sunrealtype* hlast,
     sunrealtype* hcur,
     sunrealtype* tcur) -> std::tuple<int, long, long, long, long, int, int>
  {
    auto CVodeGetIntegratorStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nsteps, long nfevals, long nlinsetups,
         long netfails, int qlast, int qcur, sunrealtype* hinused,
         sunrealtype* hlast, sunrealtype* hcur,
         sunrealtype* tcur) -> std::tuple<int, long, long, long, long, int, int>
    {
      long* nsteps_adapt_modifiable     = &nsteps;
      long* nfevals_adapt_modifiable    = &nfevals;
      long* nlinsetups_adapt_modifiable = &nlinsetups;
      long* netfails_adapt_modifiable   = &netfails;
      int* qlast_adapt_modifiable       = &qlast;
      int* qcur_adapt_modifiable        = &qcur;

      int r = CVodeGetIntegratorStats(cvode_mem, nsteps_adapt_modifiable,
                                      nfevals_adapt_modifiable,
                                      nlinsetups_adapt_modifiable,
                                      netfails_adapt_modifiable,
                                      qlast_adapt_modifiable,
                                      qcur_adapt_modifiable, hinused, hlast,
                                      hcur, tcur);
      return std::make_tuple(r, nsteps, nfevals, nlinsetups, netfails, qlast,
                             qcur);
    };

    return CVodeGetIntegratorStats_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                        nsteps,
                                                                        nfevals,
                                                                        nlinsetups,
                                                                        netfails,
                                                                        qlast,
                                                                        qcur,
                                                                        hinused,
                                                                        hlast,
                                                                        hcur,
                                                                        tcur);
  },
  nb::arg("cvode_mem"), nb::arg("nsteps"), nb::arg("nfevals"),
  nb::arg("nlinsetups"), nb::arg("netfails"), nb::arg("qlast"), nb::arg("qcur"),
  nb::arg("hinused"), nb::arg("hlast"), nb::arg("hcur"), nb::arg("tcur"));

m.def(
  "CVodeGetNumNonlinSolvIters",
  [](void* cvode_mem, long nniters) -> std::tuple<int, long>
  {
    auto CVodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nniters) -> std::tuple<int, long>
    {
      long* nniters_adapt_modifiable = &nniters;

      int r = CVodeGetNumNonlinSolvIters(cvode_mem, nniters_adapt_modifiable);
      return std::make_tuple(r, nniters);
    };

    return CVodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                           nniters);
  },
  nb::arg("cvode_mem"), nb::arg("nniters"));

m.def(
  "CVodeGetNumNonlinSolvConvFails",
  [](void* cvode_mem, long nnfails) -> std::tuple<int, long>
  {
    auto CVodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nnfails) -> std::tuple<int, long>
    {
      long* nnfails_adapt_modifiable = &nnfails;

      int r = CVodeGetNumNonlinSolvConvFails(cvode_mem, nnfails_adapt_modifiable);
      return std::make_tuple(r, nnfails);
    };

    return CVodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                               nnfails);
  },
  nb::arg("cvode_mem"), nb::arg("nnfails"));

m.def(
  "CVodeGetNonlinSolvStats",
  [](void* cvode_mem, long nniters, long nnfails) -> std::tuple<int, long, long>
  {
    auto CVodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nniters,
         long nnfails) -> std::tuple<int, long, long>
    {
      long* nniters_adapt_modifiable = &nniters;
      long* nnfails_adapt_modifiable = &nnfails;

      int r = CVodeGetNonlinSolvStats(cvode_mem, nniters_adapt_modifiable,
                                      nnfails_adapt_modifiable);
      return std::make_tuple(r, nniters, nnfails);
    };

    return CVodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                        nniters,
                                                                        nnfails);
  },
  nb::arg("cvode_mem"), nb::arg("nniters"), nb::arg("nnfails"));

m.def(
  "CVodeGetNumStepSolveFails",
  [](void* cvode_mem, long nncfails) -> std::tuple<int, long>
  {
    auto CVodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nncfails) -> std::tuple<int, long>
    {
      long* nncfails_adapt_modifiable = &nncfails;

      int r = CVodeGetNumStepSolveFails(cvode_mem, nncfails_adapt_modifiable);
      return std::make_tuple(r, nncfails);
    };

    return CVodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                          nncfails);
  },
  nb::arg("cvode_mem"), nb::arg("nncfails"));

m.def("CVodePrintAllStats", CVodePrintAllStats, nb::arg("cvode_mem"),
      nb::arg("outfile"), nb::arg("fmt"));

m.def("CVodeGetReturnFlagName", CVodeGetReturnFlagName, nb::arg("flag"));

m.def("CVodeQuadReInit", CVodeQuadReInit, nb::arg("cvode_mem"), nb::arg("yQ0"));

m.def("CVodeQuadSStolerances", CVodeQuadSStolerances, nb::arg("cvode_mem"),
      nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("CVodeQuadSVtolerances", CVodeQuadSVtolerances, nb::arg("cvode_mem"),
      nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("CVodeSetQuadErrCon", CVodeSetQuadErrCon, nb::arg("cvode_mem"),
      nb::arg("errconQ"));

m.def("CVodeGetQuad", CVodeGetQuad, nb::arg("cvode_mem"), nb::arg("tret"),
      nb::arg("yQout"));

m.def("CVodeGetQuadDky", CVodeGetQuadDky, nb::arg("cvode_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("dky"));

m.def(
  "CVodeGetQuadNumRhsEvals",
  [](void* cvode_mem, long nfQevals) -> std::tuple<int, long>
  {
    auto CVodeGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nfQevals) -> std::tuple<int, long>
    {
      long* nfQevals_adapt_modifiable = &nfQevals;

      int r = CVodeGetQuadNumRhsEvals(cvode_mem, nfQevals_adapt_modifiable);
      return std::make_tuple(r, nfQevals);
    };

    return CVodeGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                        nfQevals);
  },
  nb::arg("cvode_mem"), nb::arg("nfQevals"));

m.def(
  "CVodeGetQuadNumErrTestFails",
  [](void* cvode_mem, long nQetfails) -> std::tuple<int, long>
  {
    auto CVodeGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nQetfails) -> std::tuple<int, long>
    {
      long* nQetfails_adapt_modifiable = &nQetfails;

      int r = CVodeGetQuadNumErrTestFails(cvode_mem, nQetfails_adapt_modifiable);
      return std::make_tuple(r, nQetfails);
    };

    return CVodeGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                            nQetfails);
  },
  nb::arg("cvode_mem"), nb::arg("nQetfails"));

m.def("CVodeGetQuadErrWeights", CVodeGetQuadErrWeights, nb::arg("cvode_mem"),
      nb::arg("eQweight"));

m.def(
  "CVodeGetQuadStats",
  [](void* cvode_mem, long nfQevals, long nQetfails) -> std::tuple<int, long, long>
  {
    auto CVodeGetQuadStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nfQevals,
         long nQetfails) -> std::tuple<int, long, long>
    {
      long* nfQevals_adapt_modifiable  = &nfQevals;
      long* nQetfails_adapt_modifiable = &nQetfails;

      int r = CVodeGetQuadStats(cvode_mem, nfQevals_adapt_modifiable,
                                nQetfails_adapt_modifiable);
      return std::make_tuple(r, nfQevals, nQetfails);
    };

    return CVodeGetQuadStats_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                  nfQevals,
                                                                  nQetfails);
  },
  nb::arg("cvode_mem"), nb::arg("nfQevals"), nb::arg("nQetfails"));

m.def(
  "CVodeSensReInit",
  [](void* cvode_mem, int ism, std::vector<N_Vector> yS0) -> int
  {
    auto CVodeSensReInit_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, int ism, std::vector<N_Vector> yS0) -> int
    {
      N_Vector* yS0_ptr = yS0.empty() ? nullptr : yS0.data();

      auto lambda_result = CVodeSensReInit(cvode_mem, ism, yS0_ptr);
      return lambda_result;
    };

    return CVodeSensReInit_adapt_nvector_ptr_to_vector(cvode_mem, ism, yS0);
  },
  nb::arg("cvode_mem"), nb::arg("ism"), nb::arg("yS0"));

m.def("CVodeSensSStolerances", CVodeSensSStolerances, nb::arg("cvode_mem"),
      nb::arg("reltolS"), nb::arg("abstolS"));

m.def(
  "CVodeSensSVtolerances",
  [](void* cvode_mem, sunrealtype reltolS, std::vector<N_Vector> abstolS) -> int
  {
    auto CVodeSensSVtolerances_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, sunrealtype reltolS, std::vector<N_Vector> abstolS) -> int
    {
      N_Vector* abstolS_ptr = abstolS.empty() ? nullptr : abstolS.data();

      auto lambda_result = CVodeSensSVtolerances(cvode_mem, reltolS, abstolS_ptr);
      return lambda_result;
    };

    return CVodeSensSVtolerances_adapt_nvector_ptr_to_vector(cvode_mem, reltolS,
                                                             abstolS);
  },
  nb::arg("cvode_mem"), nb::arg("reltolS"), nb::arg("abstolS"));

m.def("CVodeSensEEtolerances", CVodeSensEEtolerances, nb::arg("cvode_mem"));

m.def("CVodeSetSensDQMethod", CVodeSetSensDQMethod, nb::arg("cvode_mem"),
      nb::arg("DQtype"), nb::arg("DQrhomax"));

m.def("CVodeSetSensErrCon", CVodeSetSensErrCon, nb::arg("cvode_mem"),
      nb::arg("errconS"));

m.def("CVodeSetSensMaxNonlinIters", CVodeSetSensMaxNonlinIters,
      nb::arg("cvode_mem"), nb::arg("maxcorS"));

m.def(
  "CVodeSetSensParams",
  [](void* cvode_mem, sunrealtype* p, sunrealtype* pbar,
     int plist) -> std::tuple<int, int>
  {
    auto CVodeSetSensParams_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, sunrealtype* p, sunrealtype* pbar,
         int plist) -> std::tuple<int, int>
    {
      int* plist_adapt_modifiable = &plist;

      int r = CVodeSetSensParams(cvode_mem, p, pbar, plist_adapt_modifiable);
      return std::make_tuple(r, plist);
    };

    return CVodeSetSensParams_adapt_modifiable_immutable_to_return(cvode_mem, p,
                                                                   pbar, plist);
  },
  nb::arg("cvode_mem"), nb::arg("p"), nb::arg("pbar"), nb::arg("plist"));

m.def("CVodeSetNonlinearSolverSensSim", CVodeSetNonlinearSolverSensSim,
      nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSetNonlinearSolverSensStg", CVodeSetNonlinearSolverSensStg,
      nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSetNonlinearSolverSensStg1", CVodeSetNonlinearSolverSensStg1,
      nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSensToggleOff", CVodeSensToggleOff, nb::arg("cvode_mem"));

m.def(
  "CVodeGetSens",
  [](void* cvode_mem, sunrealtype* tret, std::vector<N_Vector> ySout) -> int
  {
    auto CVodeGetSens_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, sunrealtype* tret, std::vector<N_Vector> ySout) -> int
    {
      N_Vector* ySout_ptr = ySout.empty() ? nullptr : ySout.data();

      auto lambda_result = CVodeGetSens(cvode_mem, tret, ySout_ptr);
      return lambda_result;
    };

    return CVodeGetSens_adapt_nvector_ptr_to_vector(cvode_mem, tret, ySout);
  },
  nb::arg("cvode_mem"), nb::arg("tret"), nb::arg("ySout"));

m.def("CVodeGetSens1", CVodeGetSens1, nb::arg("cvode_mem"), nb::arg("tret"),
      nb::arg("is_"), nb::arg("ySout"));

m.def(
  "CVodeGetSensDky",
  [](void* cvode_mem, sunrealtype t, int k, std::vector<N_Vector> dkyA) -> int
  {
    auto CVodeGetSensDky_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, sunrealtype t, int k, std::vector<N_Vector> dkyA) -> int
    {
      N_Vector* dkyA_ptr = dkyA.empty() ? nullptr : dkyA.data();

      auto lambda_result = CVodeGetSensDky(cvode_mem, t, k, dkyA_ptr);
      return lambda_result;
    };

    return CVodeGetSensDky_adapt_nvector_ptr_to_vector(cvode_mem, t, k, dkyA);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyA"));

m.def("CVodeGetSensDky1", CVodeGetSensDky1, nb::arg("cvode_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("is_"), nb::arg("dky"));

m.def(
  "CVodeGetSensNumRhsEvals",
  [](void* cvode_mem, long nfSevals) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nfSevals) -> std::tuple<int, long>
    {
      long* nfSevals_adapt_modifiable = &nfSevals;

      int r = CVodeGetSensNumRhsEvals(cvode_mem, nfSevals_adapt_modifiable);
      return std::make_tuple(r, nfSevals);
    };

    return CVodeGetSensNumRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                        nfSevals);
  },
  nb::arg("cvode_mem"), nb::arg("nfSevals"));

m.def(
  "CVodeGetNumRhsEvalsSens",
  [](void* cvode_mem, long nfevalsS) -> std::tuple<int, long>
  {
    auto CVodeGetNumRhsEvalsSens_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nfevalsS) -> std::tuple<int, long>
    {
      long* nfevalsS_adapt_modifiable = &nfevalsS;

      int r = CVodeGetNumRhsEvalsSens(cvode_mem, nfevalsS_adapt_modifiable);
      return std::make_tuple(r, nfevalsS);
    };

    return CVodeGetNumRhsEvalsSens_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                        nfevalsS);
  },
  nb::arg("cvode_mem"), nb::arg("nfevalsS"));

m.def(
  "CVodeGetSensNumErrTestFails",
  [](void* cvode_mem, long nSetfails) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nSetfails) -> std::tuple<int, long>
    {
      long* nSetfails_adapt_modifiable = &nSetfails;

      int r = CVodeGetSensNumErrTestFails(cvode_mem, nSetfails_adapt_modifiable);
      return std::make_tuple(r, nSetfails);
    };

    return CVodeGetSensNumErrTestFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                            nSetfails);
  },
  nb::arg("cvode_mem"), nb::arg("nSetfails"));

m.def(
  "CVodeGetSensNumLinSolvSetups",
  [](void* cvode_mem, long nlinsetupsS) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nlinsetupsS) -> std::tuple<int, long>
    {
      long* nlinsetupsS_adapt_modifiable = &nlinsetupsS;

      int r = CVodeGetSensNumLinSolvSetups(cvode_mem,
                                           nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nlinsetupsS);
    };

    return CVodeGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                             nlinsetupsS);
  },
  nb::arg("cvode_mem"), nb::arg("nlinsetupsS"));

m.def(
  "CVodeGetSensErrWeights",
  [](void* cvode_mem, std::vector<N_Vector> eSweight) -> int
  {
    auto CVodeGetSensErrWeights_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, std::vector<N_Vector> eSweight) -> int
    {
      N_Vector* eSweight_ptr = eSweight.empty() ? nullptr : eSweight.data();

      auto lambda_result = CVodeGetSensErrWeights(cvode_mem, eSweight_ptr);
      return lambda_result;
    };

    return CVodeGetSensErrWeights_adapt_nvector_ptr_to_vector(cvode_mem,
                                                              eSweight);
  },
  nb::arg("cvode_mem"), nb::arg("eSweight"));

m.def(
  "CVodeGetSensStats",
  [](void* cvode_mem, long nfSevals, long nfevalsS, long nSetfails,
     long nlinsetupsS) -> std::tuple<int, long, long, long, long>
  {
    auto CVodeGetSensStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nfSevals, long nfevalsS, long nSetfails,
         long nlinsetupsS) -> std::tuple<int, long, long, long, long>
    {
      long* nfSevals_adapt_modifiable    = &nfSevals;
      long* nfevalsS_adapt_modifiable    = &nfevalsS;
      long* nSetfails_adapt_modifiable   = &nSetfails;
      long* nlinsetupsS_adapt_modifiable = &nlinsetupsS;

      int r = CVodeGetSensStats(cvode_mem, nfSevals_adapt_modifiable,
                                nfevalsS_adapt_modifiable,
                                nSetfails_adapt_modifiable,
                                nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nfSevals, nfevalsS, nSetfails, nlinsetupsS);
    };

    return CVodeGetSensStats_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                  nfSevals,
                                                                  nfevalsS,
                                                                  nSetfails,
                                                                  nlinsetupsS);
  },
  nb::arg("cvode_mem"), nb::arg("nfSevals"), nb::arg("nfevalsS"),
  nb::arg("nSetfails"), nb::arg("nlinsetupsS"));

m.def(
  "CVodeGetSensNumNonlinSolvIters",
  [](void* cvode_mem, long nSniters) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nSniters) -> std::tuple<int, long>
    {
      long* nSniters_adapt_modifiable = &nSniters;

      int r = CVodeGetSensNumNonlinSolvIters(cvode_mem,
                                             nSniters_adapt_modifiable);
      return std::make_tuple(r, nSniters);
    };

    return CVodeGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                               nSniters);
  },
  nb::arg("cvode_mem"), nb::arg("nSniters"));

m.def(
  "CVodeGetSensNumNonlinSolvConvFails",
  [](void* cvode_mem, long nSnfails) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nSnfails) -> std::tuple<int, long>
    {
      long* nSnfails_adapt_modifiable = &nSnfails;

      int r = CVodeGetSensNumNonlinSolvConvFails(cvode_mem,
                                                 nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSnfails);
    };

    return CVodeGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                                   nSnfails);
  },
  nb::arg("cvode_mem"), nb::arg("nSnfails"));

m.def(
  "CVodeGetSensNonlinSolvStats",
  [](void* cvode_mem, long nSniters, long nSnfails) -> std::tuple<int, long, long>
  {
    auto CVodeGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nSniters,
         long nSnfails) -> std::tuple<int, long, long>
    {
      long* nSniters_adapt_modifiable = &nSniters;
      long* nSnfails_adapt_modifiable = &nSnfails;

      int r = CVodeGetSensNonlinSolvStats(cvode_mem, nSniters_adapt_modifiable,
                                          nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSniters, nSnfails);
    };

    return CVodeGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                            nSniters,
                                                                            nSnfails);
  },
  nb::arg("cvode_mem"), nb::arg("nSniters"), nb::arg("nSnfails"));

m.def(
  "CVodeGetNumStepSensSolveFails",
  [](void* cvode_mem, long nSncfails) -> std::tuple<int, long>
  {
    auto CVodeGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nSncfails) -> std::tuple<int, long>
    {
      long* nSncfails_adapt_modifiable = &nSncfails;

      int r = CVodeGetNumStepSensSolveFails(cvode_mem,
                                            nSncfails_adapt_modifiable);
      return std::make_tuple(r, nSncfails);
    };

    return CVodeGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                              nSncfails);
  },
  nb::arg("cvode_mem"), nb::arg("nSncfails"));

m.def(
  "CVodeGetStgrSensNumNonlinSolvIters",
  [](void* cvode_mem, long nSTGR1niters) -> std::tuple<int, long>
  {
    auto CVodeGetStgrSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nSTGR1niters) -> std::tuple<int, long>
    {
      long* nSTGR1niters_adapt_modifiable = &nSTGR1niters;

      int r = CVodeGetStgrSensNumNonlinSolvIters(cvode_mem,
                                                 nSTGR1niters_adapt_modifiable);
      return std::make_tuple(r, nSTGR1niters);
    };

    return CVodeGetStgrSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                                   nSTGR1niters);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1niters"));

m.def(
  "CVodeGetStgrSensNumNonlinSolvConvFails",
  [](void* cvode_mem, long nSTGR1nfails) -> std::tuple<int, long>
  {
    auto CVodeGetStgrSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nSTGR1nfails) -> std::tuple<int, long>
    {
      long* nSTGR1nfails_adapt_modifiable = &nSTGR1nfails;

      int r =
        CVodeGetStgrSensNumNonlinSolvConvFails(cvode_mem,
                                               nSTGR1nfails_adapt_modifiable);
      return std::make_tuple(r, nSTGR1nfails);
    };

    return CVodeGetStgrSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                                       nSTGR1nfails);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1nfails"));

m.def(
  "CVodeGetStgrSensNonlinSolvStats",
  [](void* cvode_mem, long nSTGR1niters,
     long nSTGR1nfails) -> std::tuple<int, long, long>
  {
    auto CVodeGetStgrSensNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nSTGR1niters,
         long nSTGR1nfails) -> std::tuple<int, long, long>
    {
      long* nSTGR1niters_adapt_modifiable = &nSTGR1niters;
      long* nSTGR1nfails_adapt_modifiable = &nSTGR1nfails;

      int r = CVodeGetStgrSensNonlinSolvStats(cvode_mem,
                                              nSTGR1niters_adapt_modifiable,
                                              nSTGR1nfails_adapt_modifiable);
      return std::make_tuple(r, nSTGR1niters, nSTGR1nfails);
    };

    return CVodeGetStgrSensNonlinSolvStats_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                                nSTGR1niters,
                                                                                nSTGR1nfails);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1niters"), nb::arg("nSTGR1nfails"));

m.def(
  "CVodeGetNumStepStgrSensSolveFails",
  [](void* cvode_mem, long nSTGR1ncfails) -> std::tuple<int, long>
  {
    auto CVodeGetNumStepStgrSensSolveFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nSTGR1ncfails) -> std::tuple<int, long>
    {
      long* nSTGR1ncfails_adapt_modifiable = &nSTGR1ncfails;

      int r = CVodeGetNumStepStgrSensSolveFails(cvode_mem,
                                                nSTGR1ncfails_adapt_modifiable);
      return std::make_tuple(r, nSTGR1ncfails);
    };

    return CVodeGetNumStepStgrSensSolveFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                                  nSTGR1ncfails);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1ncfails"));

m.def(
  "CVodeQuadSensReInit",
  [](void* cvode_mem, std::vector<N_Vector> yQS0) -> int
  {
    auto CVodeQuadSensReInit_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, std::vector<N_Vector> yQS0) -> int
    {
      N_Vector* yQS0_ptr = yQS0.empty() ? nullptr : yQS0.data();

      auto lambda_result = CVodeQuadSensReInit(cvode_mem, yQS0_ptr);
      return lambda_result;
    };

    return CVodeQuadSensReInit_adapt_nvector_ptr_to_vector(cvode_mem, yQS0);
  },
  nb::arg("cvode_mem"), nb::arg("yQS0"));

m.def("CVodeQuadSensSStolerances", CVodeQuadSensSStolerances,
      nb::arg("cvode_mem"), nb::arg("reltolQS"), nb::arg("abstolQS"));

m.def(
  "CVodeQuadSensSVtolerances",
  [](void* cvode_mem, sunrealtype reltolQS, std::vector<N_Vector> abstolQS) -> int
  {
    auto CVodeQuadSensSVtolerances_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, sunrealtype reltolQS,
         std::vector<N_Vector> abstolQS) -> int
    {
      N_Vector* abstolQS_ptr = abstolQS.empty() ? nullptr : abstolQS.data();

      auto lambda_result = CVodeQuadSensSVtolerances(cvode_mem, reltolQS,
                                                     abstolQS_ptr);
      return lambda_result;
    };

    return CVodeQuadSensSVtolerances_adapt_nvector_ptr_to_vector(cvode_mem,
                                                                 reltolQS,
                                                                 abstolQS);
  },
  nb::arg("cvode_mem"), nb::arg("reltolQS"), nb::arg("abstolQS"));

m.def("CVodeQuadSensEEtolerances", CVodeQuadSensEEtolerances,
      nb::arg("cvode_mem"));

m.def("CVodeSetQuadSensErrCon", CVodeSetQuadSensErrCon, nb::arg("cvode_mem"),
      nb::arg("errconQS"));

m.def(
  "CVodeGetQuadSens",
  [](void* cvode_mem, sunrealtype* tret, std::vector<N_Vector> yQSout) -> int
  {
    auto CVodeGetQuadSens_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, sunrealtype* tret, std::vector<N_Vector> yQSout) -> int
    {
      N_Vector* yQSout_ptr = yQSout.empty() ? nullptr : yQSout.data();

      auto lambda_result = CVodeGetQuadSens(cvode_mem, tret, yQSout_ptr);
      return lambda_result;
    };

    return CVodeGetQuadSens_adapt_nvector_ptr_to_vector(cvode_mem, tret, yQSout);
  },
  nb::arg("cvode_mem"), nb::arg("tret"), nb::arg("yQSout"));

m.def("CVodeGetQuadSens1", CVodeGetQuadSens1, nb::arg("cvode_mem"),
      nb::arg("tret"), nb::arg("is_"), nb::arg("yQSout"));

m.def(
  "CVodeGetQuadSensDky",
  [](void* cvode_mem, sunrealtype t, int k, std::vector<N_Vector> dkyQS_all) -> int
  {
    auto CVodeGetQuadSensDky_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, sunrealtype t, int k,
         std::vector<N_Vector> dkyQS_all) -> int
    {
      N_Vector* dkyQS_all_ptr = dkyQS_all.empty() ? nullptr : dkyQS_all.data();

      auto lambda_result = CVodeGetQuadSensDky(cvode_mem, t, k, dkyQS_all_ptr);
      return lambda_result;
    };

    return CVodeGetQuadSensDky_adapt_nvector_ptr_to_vector(cvode_mem, t, k,
                                                           dkyQS_all);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyQS_all"));

m.def("CVodeGetQuadSensDky1", CVodeGetQuadSensDky1, nb::arg("cvode_mem"),
      nb::arg("t"), nb::arg("k"), nb::arg("is_"), nb::arg("dkyQS"));

m.def(
  "CVodeGetQuadSensNumRhsEvals",
  [](void* cvode_mem, long nfQSevals) -> std::tuple<int, long>
  {
    auto CVodeGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nfQSevals) -> std::tuple<int, long>
    {
      long* nfQSevals_adapt_modifiable = &nfQSevals;

      int r = CVodeGetQuadSensNumRhsEvals(cvode_mem, nfQSevals_adapt_modifiable);
      return std::make_tuple(r, nfQSevals);
    };

    return CVodeGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                            nfQSevals);
  },
  nb::arg("cvode_mem"), nb::arg("nfQSevals"));

m.def(
  "CVodeGetQuadSensNumErrTestFails",
  [](void* cvode_mem, long nQSetfails) -> std::tuple<int, long>
  {
    auto CVodeGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nQSetfails) -> std::tuple<int, long>
    {
      long* nQSetfails_adapt_modifiable = &nQSetfails;

      int r = CVodeGetQuadSensNumErrTestFails(cvode_mem,
                                              nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nQSetfails);
    };

    return CVodeGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                                nQSetfails);
  },
  nb::arg("cvode_mem"), nb::arg("nQSetfails"));

m.def(
  "CVodeGetQuadSensErrWeights",
  [](void* cvode_mem, std::vector<N_Vector> eQSweight) -> int
  {
    auto CVodeGetQuadSensErrWeights_adapt_nvector_ptr_to_vector =
      [](void* cvode_mem, std::vector<N_Vector> eQSweight) -> int
    {
      N_Vector* eQSweight_ptr = eQSweight.empty() ? nullptr : eQSweight.data();

      auto lambda_result = CVodeGetQuadSensErrWeights(cvode_mem, eQSweight_ptr);
      return lambda_result;
    };

    return CVodeGetQuadSensErrWeights_adapt_nvector_ptr_to_vector(cvode_mem,
                                                                  eQSweight);
  },
  nb::arg("cvode_mem"), nb::arg("eQSweight"));

m.def(
  "CVodeGetQuadSensStats",
  [](void* cvode_mem, long nfQSevals,
     long nQSetfails) -> std::tuple<int, long, long>
  {
    auto CVodeGetQuadSensStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nfQSevals,
         long nQSetfails) -> std::tuple<int, long, long>
    {
      long* nfQSevals_adapt_modifiable  = &nfQSevals;
      long* nQSetfails_adapt_modifiable = &nQSetfails;

      int r = CVodeGetQuadSensStats(cvode_mem, nfQSevals_adapt_modifiable,
                                    nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nfQSevals, nQSetfails);
    };

    return CVodeGetQuadSensStats_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                      nfQSevals,
                                                                      nQSetfails);
  },
  nb::arg("cvode_mem"), nb::arg("nfQSevals"), nb::arg("nQSetfails"));

m.def("CVodeAdjInit", CVodeAdjInit, nb::arg("cvode_mem"), nb::arg("steps"),
      nb::arg("interp"));

m.def("CVodeAdjReInit", CVodeAdjReInit, nb::arg("cvode_mem"));

m.def(
  "CVodeCreateB",
  [](void* cvode_mem, int lmmB, int which) -> std::tuple<int, int>
  {
    auto CVodeCreateB_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int lmmB, int which) -> std::tuple<int, int>
    {
      int* which_adapt_modifiable = &which;

      int r = CVodeCreateB(cvode_mem, lmmB, which_adapt_modifiable);
      return std::make_tuple(r, which);
    };

    return CVodeCreateB_adapt_modifiable_immutable_to_return(cvode_mem, lmmB,
                                                             which);
  },
  nb::arg("cvode_mem"), nb::arg("lmmB"), nb::arg("which"));

m.def("CVodeReInitB", CVodeReInitB, nb::arg("cvode_mem"), nb::arg("which"),
      nb::arg("tB0"), nb::arg("yB0"));

m.def("CVodeSStolerancesB", CVodeSStolerancesB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("reltolB"), nb::arg("abstolB"));

m.def("CVodeSVtolerancesB", CVodeSVtolerancesB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("reltolB"), nb::arg("abstolB"));

m.def("CVodeQuadReInitB", CVodeQuadReInitB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("yQB0"));

m.def("CVodeQuadSStolerancesB", CVodeQuadSStolerancesB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def("CVodeQuadSVtolerancesB", CVodeQuadSVtolerancesB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def(
  "CVodeF",
  [](void* cvode_mem, sunrealtype tout, N_Vector yout, sunrealtype* tret,
     int itask, int ncheckPtr) -> std::tuple<int, int>
  {
    auto CVodeF_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, sunrealtype tout, N_Vector yout, sunrealtype* tret,
         int itask, int ncheckPtr) -> std::tuple<int, int>
    {
      int* ncheckPtr_adapt_modifiable = &ncheckPtr;

      int r = CVodeF(cvode_mem, tout, yout, tret, itask,
                     ncheckPtr_adapt_modifiable);
      return std::make_tuple(r, ncheckPtr);
    };

    return CVodeF_adapt_modifiable_immutable_to_return(cvode_mem, tout, yout,
                                                       tret, itask, ncheckPtr);
  },
  nb::arg("cvode_mem"), nb::arg("tout"), nb::arg("yout"), nb::arg("tret"),
  nb::arg("itask"), nb::arg("ncheckPtr"));

m.def("CVodeB", CVodeB, nb::arg("cvode_mem"), nb::arg("tBout"),
      nb::arg("itaskB"));

m.def("CVodeSetAdjNoSensi", CVodeSetAdjNoSensi, nb::arg("cvode_mem"));

m.def("CVodeSetOwnUserDataB", CVodeSetOwnUserDataB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("own_user_data"));

m.def("CVodeSetMaxOrdB", CVodeSetMaxOrdB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("maxordB"));

m.def("CVodeSetMaxNumStepsB", CVodeSetMaxNumStepsB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("mxstepsB"));

m.def("CVodeSetStabLimDetB", CVodeSetStabLimDetB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("stldetB"));

m.def("CVodeSetInitStepB", CVodeSetInitStepB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("hinB"));

m.def("CVodeSetMinStepB", CVodeSetMinStepB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("hminB"));

m.def("CVodeSetMaxStepB", CVodeSetMaxStepB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("hmaxB"));

m.def("CVodeSetConstraintsB", CVodeSetConstraintsB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("constraintsB"));

m.def("CVodeSetQuadErrConB", CVodeSetQuadErrConB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("errconQB"));

m.def("CVodeSetNonlinearSolverB", CVodeSetNonlinearSolverB,
      nb::arg("cvode_mem"), nb::arg("which"), nb::arg("NLS"));

m.def("CVodeGetB", CVodeGetB, nb::arg("cvode_mem"), nb::arg("which"),
      nb::arg("tBret"), nb::arg("yB"));

m.def("CVodeGetQuadB", CVodeGetQuadB, nb::arg("cvode_mem"), nb::arg("which"),
      nb::arg("tBret"), nb::arg("qB"));

m.def("CVodeGetAdjCVodeBmem", CVodeGetAdjCVodeBmem, nb::arg("cvode_mem"),
      nb::arg("which"));

m.def("CVodeGetAdjY", CVodeGetAdjY, nb::arg("cvode_mem"), nb::arg("t"),
      nb::arg("y"));

m.def("CVodeGetAdjCheckPointsInfo", CVodeGetAdjCheckPointsInfo,
      nb::arg("cvode_mem"), nb::arg("ckpnt"));

m.def("CVodeGetAdjDataPointHermite", CVodeGetAdjDataPointHermite,
      nb::arg("cvode_mem"), nb::arg("which"), nb::arg("t"), nb::arg("y"),
      nb::arg("yd"));

m.def(
  "CVodeGetAdjDataPointPolynomial",
  [](void* cvode_mem, int which, sunrealtype* t, int order,
     N_Vector y) -> std::tuple<int, int>
  {
    auto CVodeGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int which, sunrealtype* t, int order,
         N_Vector y) -> std::tuple<int, int>
    {
      int* order_adapt_modifiable = &order;

      int r = CVodeGetAdjDataPointPolynomial(cvode_mem, which, t,
                                             order_adapt_modifiable, y);
      return std::make_tuple(r, order);
    };

    return CVodeGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                               which,
                                                                               t,
                                                                               order,
                                                                               y);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("t"), nb::arg("order"),
  nb::arg("y"));
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

m.def(
  "CVodeSetLinearSolver",
  [](void* cvode_mem, SUNLinearSolver LS,
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

    return CVodeSetLinearSolver_adapt_optional_arg_with_default_null(cvode_mem,
                                                                     LS, A);
  },
  nb::arg("cvode_mem"), nb::arg("LS"), nb::arg("A") = nb::none());

m.def("CVodeSetJacEvalFrequency", CVodeSetJacEvalFrequency,
      nb::arg("cvode_mem"), nb::arg("msbj"));

m.def("CVodeSetLinearSolutionScaling", CVodeSetLinearSolutionScaling,
      nb::arg("cvode_mem"), nb::arg("onoff"));

m.def("CVodeSetDeltaGammaMaxBadJac", CVodeSetDeltaGammaMaxBadJac,
      nb::arg("cvode_mem"), nb::arg("dgmax_jbad"));

m.def("CVodeSetEpsLin", CVodeSetEpsLin, nb::arg("cvode_mem"), nb::arg("eplifac"));

m.def("CVodeSetLSNormFactor", CVodeSetLSNormFactor, nb::arg("arkode_mem"),
      nb::arg("nrmfac"));

m.def("CVodeGetJacTime", CVodeGetJacTime, nb::arg("cvode_mem"), nb::arg("t_J"));

m.def(
  "CVodeGetJacNumSteps",
  [](void* cvode_mem, long nst_J) -> std::tuple<int, long>
  {
    auto CVodeGetJacNumSteps_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nst_J) -> std::tuple<int, long>
    {
      long* nst_J_adapt_modifiable = &nst_J;

      int r = CVodeGetJacNumSteps(cvode_mem, nst_J_adapt_modifiable);
      return std::make_tuple(r, nst_J);
    };

    return CVodeGetJacNumSteps_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                    nst_J);
  },
  nb::arg("cvode_mem"), nb::arg("nst_J"));

m.def(
  "CVodeGetNumJacEvals",
  [](void* cvode_mem, long njevals) -> std::tuple<int, long>
  {
    auto CVodeGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long njevals) -> std::tuple<int, long>
    {
      long* njevals_adapt_modifiable = &njevals;

      int r = CVodeGetNumJacEvals(cvode_mem, njevals_adapt_modifiable);
      return std::make_tuple(r, njevals);
    };

    return CVodeGetNumJacEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                    njevals);
  },
  nb::arg("cvode_mem"), nb::arg("njevals"));

m.def(
  "CVodeGetNumPrecEvals",
  [](void* cvode_mem, long npevals) -> std::tuple<int, long>
  {
    auto CVodeGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long npevals) -> std::tuple<int, long>
    {
      long* npevals_adapt_modifiable = &npevals;

      int r = CVodeGetNumPrecEvals(cvode_mem, npevals_adapt_modifiable);
      return std::make_tuple(r, npevals);
    };

    return CVodeGetNumPrecEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                     npevals);
  },
  nb::arg("cvode_mem"), nb::arg("npevals"));

m.def(
  "CVodeGetNumPrecSolves",
  [](void* cvode_mem, long npsolves) -> std::tuple<int, long>
  {
    auto CVodeGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long npsolves) -> std::tuple<int, long>
    {
      long* npsolves_adapt_modifiable = &npsolves;

      int r = CVodeGetNumPrecSolves(cvode_mem, npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves);
    };

    return CVodeGetNumPrecSolves_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                      npsolves);
  },
  nb::arg("cvode_mem"), nb::arg("npsolves"));

m.def(
  "CVodeGetNumLinIters",
  [](void* cvode_mem, long nliters) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nliters) -> std::tuple<int, long>
    {
      long* nliters_adapt_modifiable = &nliters;

      int r = CVodeGetNumLinIters(cvode_mem, nliters_adapt_modifiable);
      return std::make_tuple(r, nliters);
    };

    return CVodeGetNumLinIters_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                    nliters);
  },
  nb::arg("cvode_mem"), nb::arg("nliters"));

m.def(
  "CVodeGetNumLinConvFails",
  [](void* cvode_mem, long nlcfails) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nlcfails) -> std::tuple<int, long>
    {
      long* nlcfails_adapt_modifiable = &nlcfails;

      int r = CVodeGetNumLinConvFails(cvode_mem, nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails);
    };

    return CVodeGetNumLinConvFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                        nlcfails);
  },
  nb::arg("cvode_mem"), nb::arg("nlcfails"));

m.def(
  "CVodeGetNumJTSetupEvals",
  [](void* cvode_mem, long njtsetups) -> std::tuple<int, long>
  {
    auto CVodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long njtsetups) -> std::tuple<int, long>
    {
      long* njtsetups_adapt_modifiable = &njtsetups;

      int r = CVodeGetNumJTSetupEvals(cvode_mem, njtsetups_adapt_modifiable);
      return std::make_tuple(r, njtsetups);
    };

    return CVodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                        njtsetups);
  },
  nb::arg("cvode_mem"), nb::arg("njtsetups"));

m.def(
  "CVodeGetNumJtimesEvals",
  [](void* cvode_mem, long njvevals) -> std::tuple<int, long>
  {
    auto CVodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long njvevals) -> std::tuple<int, long>
    {
      long* njvevals_adapt_modifiable = &njvevals;

      int r = CVodeGetNumJtimesEvals(cvode_mem, njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals);
    };

    return CVodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                       njvevals);
  },
  nb::arg("cvode_mem"), nb::arg("njvevals"));

m.def(
  "CVodeGetNumLinRhsEvals",
  [](void* cvode_mem, long nfevalsLS) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nfevalsLS) -> std::tuple<int, long>
    {
      long* nfevalsLS_adapt_modifiable = &nfevalsLS;

      int r = CVodeGetNumLinRhsEvals(cvode_mem, nfevalsLS_adapt_modifiable);
      return std::make_tuple(r, nfevalsLS);
    };

    return CVodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                       nfevalsLS);
  },
  nb::arg("cvode_mem"), nb::arg("nfevalsLS"));

m.def(
  "CVodeGetLinSolveStats",
  [](void* cvode_mem, long njevals, long nfevalsLS, long nliters, long nlcfails,
     long npevals, long npsolves, long njtsetups, long njtimes)
    -> std::tuple<int, long, long, long, long, long, long, long, long>
  {
    auto CVodeGetLinSolveStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long njevals, long nfevalsLS, long nliters,
         long nlcfails, long npevals, long npsolves, long njtsetups, long njtimes)
      -> std::tuple<int, long, long, long, long, long, long, long, long>
    {
      long* njevals_adapt_modifiable   = &njevals;
      long* nfevalsLS_adapt_modifiable = &nfevalsLS;
      long* nliters_adapt_modifiable   = &nliters;
      long* nlcfails_adapt_modifiable  = &nlcfails;
      long* npevals_adapt_modifiable   = &npevals;
      long* npsolves_adapt_modifiable  = &npsolves;
      long* njtsetups_adapt_modifiable = &njtsetups;
      long* njtimes_adapt_modifiable   = &njtimes;

      int r = CVodeGetLinSolveStats(cvode_mem, njevals_adapt_modifiable,
                                    nfevalsLS_adapt_modifiable,
                                    nliters_adapt_modifiable,
                                    nlcfails_adapt_modifiable,
                                    npevals_adapt_modifiable,
                                    npsolves_adapt_modifiable,
                                    njtsetups_adapt_modifiable,
                                    njtimes_adapt_modifiable);
      return std::make_tuple(r, njevals, nfevalsLS, nliters, nlcfails, npevals,
                             npsolves, njtsetups, njtimes);
    };

    return CVodeGetLinSolveStats_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                      njevals,
                                                                      nfevalsLS,
                                                                      nliters,
                                                                      nlcfails,
                                                                      npevals,
                                                                      npsolves,
                                                                      njtsetups,
                                                                      njtimes);
  },
  nb::arg("cvode_mem"), nb::arg("njevals"), nb::arg("nfevalsLS"),
  nb::arg("nliters"), nb::arg("nlcfails"), nb::arg("npevals"),
  nb::arg("npsolves"), nb::arg("njtsetups"), nb::arg("njtimes"));

m.def(
  "CVodeGetLastLinFlag",
  [](void* cvode_mem, long flag) -> std::tuple<int, long>
  {
    auto CVodeGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long flag) -> std::tuple<int, long>
    {
      long* flag_adapt_modifiable = &flag;

      int r = CVodeGetLastLinFlag(cvode_mem, flag_adapt_modifiable);
      return std::make_tuple(r, flag);
    };

    return CVodeGetLastLinFlag_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                    flag);
  },
  nb::arg("cvode_mem"), nb::arg("flag"));

m.def("CVodeGetLinReturnFlagName", CVodeGetLinReturnFlagName, nb::arg("flag"));

m.def(
  "CVodeSetLinearSolverB",
  [](void* cvode_mem, int which, SUNLinearSolver LS,
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

    return CVodeSetLinearSolverB_adapt_optional_arg_with_default_null(cvode_mem,
                                                                      which, LS,
                                                                      A);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("LS"),
  nb::arg("A") = nb::none());

m.def("CVodeSetEpsLinB", CVodeSetEpsLinB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("eplifacB"));

m.def("CVodeSetLSNormFactorB", CVodeSetLSNormFactorB, nb::arg("arkode_mem"),
      nb::arg("which"), nb::arg("nrmfacB"));

m.def("CVodeSetLinearSolutionScalingB", CVodeSetLinearSolutionScalingB,
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

m.def("CVodeSetProjErrEst", CVodeSetProjErrEst, nb::arg("cvode_mem"),
      nb::arg("onoff"));

m.def("CVodeSetProjFrequency", CVodeSetProjFrequency, nb::arg("cvode_mem"),
      nb::arg("proj_freq"));

m.def("CVodeSetMaxNumProjFails", CVodeSetMaxNumProjFails, nb::arg("cvode_mem"),
      nb::arg("max_fails"));

m.def("CVodeSetEpsProj", CVodeSetEpsProj, nb::arg("cvode_mem"), nb::arg("eps"));

m.def("CVodeSetProjFailEta", CVodeSetProjFailEta, nb::arg("cvode_mem"),
      nb::arg("eta"));

m.def(
  "CVodeGetNumProjEvals",
  [](void* cvode_mem, long nproj) -> std::tuple<int, long>
  {
    auto CVodeGetNumProjEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nproj) -> std::tuple<int, long>
    {
      long* nproj_adapt_modifiable = &nproj;

      int r = CVodeGetNumProjEvals(cvode_mem, nproj_adapt_modifiable);
      return std::make_tuple(r, nproj);
    };

    return CVodeGetNumProjEvals_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                     nproj);
  },
  nb::arg("cvode_mem"), nb::arg("nproj"));

m.def(
  "CVodeGetNumProjFails",
  [](void* cvode_mem, long nprf) -> std::tuple<int, long>
  {
    auto CVodeGetNumProjFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, long nprf) -> std::tuple<int, long>
    {
      long* nprf_adapt_modifiable = &nprf;

      int r = CVodeGetNumProjFails(cvode_mem, nprf_adapt_modifiable);
      return std::make_tuple(r, nprf);
    };

    return CVodeGetNumProjFails_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                     nprf);
  },
  nb::arg("cvode_mem"), nb::arg("nprf"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
