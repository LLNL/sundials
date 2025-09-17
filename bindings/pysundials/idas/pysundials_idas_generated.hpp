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
m.attr("IDA_UNRECOGNIZED_ERROR") = -99;
m.attr("IDA_NO_ADJ")             = -101;
m.attr("IDA_NO_FWD")             = -102;
m.attr("IDA_NO_BCK")             = -103;
m.attr("IDA_BAD_TB0")            = -104;
m.attr("IDA_REIFWD_FAIL")        = -105;
m.attr("IDA_FWD_FAIL")           = -106;
m.attr("IDA_GETY_BADT")          = -107;

m.def("IDACreate", IDACreate, nb::arg("sunctx"));

m.def("IDAReInit", IDAReInit, nb::arg("ida_mem"), nb::arg("t0"), nb::arg("yy0"),
      nb::arg("yp0"));

m.def("IDASStolerances", IDASStolerances, nb::arg("ida_mem"), nb::arg("reltol"),
      nb::arg("abstol"));

m.def("IDASVtolerances", IDASVtolerances, nb::arg("ida_mem"), nb::arg("reltol"),
      nb::arg("abstol"));

m.def("IDAWFtolerances", IDAWFtolerances, nb::arg("ida_mem"), nb::arg("efun"));

m.def("IDACalcIC", IDACalcIC, nb::arg("ida_mem"), nb::arg("icopt"),
      nb::arg("tout1"));

m.def("IDASetNonlinConvCoefIC", IDASetNonlinConvCoefIC, nb::arg("ida_mem"),
      nb::arg("epiccon"));

m.def("IDASetMaxNumStepsIC", IDASetMaxNumStepsIC, nb::arg("ida_mem"),
      nb::arg("maxnh"));

m.def("IDASetMaxNumJacsIC", IDASetMaxNumJacsIC, nb::arg("ida_mem"),
      nb::arg("maxnj"));

m.def("IDASetMaxNumItersIC", IDASetMaxNumItersIC, nb::arg("ida_mem"),
      nb::arg("maxnit"));

m.def("IDASetLineSearchOffIC", IDASetLineSearchOffIC, nb::arg("ida_mem"),
      nb::arg("lsoff"));

m.def("IDASetStepToleranceIC", IDASetStepToleranceIC, nb::arg("ida_mem"),
      nb::arg("steptol"));

m.def("IDASetMaxBacksIC", IDASetMaxBacksIC, nb::arg("ida_mem"),
      nb::arg("maxbacks"));

m.def("IDASetDeltaCjLSetup", IDASetDeltaCjLSetup, nb::arg("ida_max"),
      nb::arg("dcj"));

m.def("IDASetOwnUserData", IDASetOwnUserData, nb::arg("ida_mem"),
      nb::arg("own_user_data"));

m.def("IDASetMaxOrd", IDASetMaxOrd, nb::arg("ida_mem"), nb::arg("maxord"));

m.def("IDASetMaxNumSteps", IDASetMaxNumSteps, nb::arg("ida_mem"),
      nb::arg("mxsteps"));

m.def("IDASetInitStep", IDASetInitStep, nb::arg("ida_mem"), nb::arg("hin"));

m.def("IDASetMaxStep", IDASetMaxStep, nb::arg("ida_mem"), nb::arg("hmax"));

m.def("IDASetMinStep", IDASetMinStep, nb::arg("ida_mem"), nb::arg("hmin"));

m.def("IDASetStopTime", IDASetStopTime, nb::arg("ida_mem"), nb::arg("tstop"));

m.def("IDAClearStopTime", IDAClearStopTime, nb::arg("ida_mem"));

m.def("IDASetMaxErrTestFails", IDASetMaxErrTestFails, nb::arg("ida_mem"),
      nb::arg("maxnef"));

m.def("IDASetSuppressAlg", IDASetSuppressAlg, nb::arg("ida_mem"),
      nb::arg("suppressalg"));

m.def("IDASetId", IDASetId, nb::arg("ida_mem"), nb::arg("id"));

m.def("IDASetConstraints", IDASetConstraints, nb::arg("ida_mem"),
      nb::arg("constraints"));

m.def("IDASetEtaFixedStepBounds", IDASetEtaFixedStepBounds, nb::arg("ida_mem"),
      nb::arg("eta_min_fx"), nb::arg("eta_max_fx"));

m.def("IDASetEtaMin", IDASetEtaMin, nb::arg("ida_mem"), nb::arg("eta_min"));

m.def("IDASetEtaMax", IDASetEtaMax, nb::arg("ida_mem"), nb::arg("eta_max"));

m.def("IDASetEtaLow", IDASetEtaLow, nb::arg("ida_mem"), nb::arg("eta_low"));

m.def("IDASetEtaMinErrFail", IDASetEtaMinErrFail, nb::arg("ida_mem"),
      nb::arg("eta_min_ef"));

m.def("IDASetEtaConvFail", IDASetEtaConvFail, nb::arg("ida_mem"),
      nb::arg("eta_cf"));

m.def("IDASetMaxConvFails", IDASetMaxConvFails, nb::arg("ida_mem"),
      nb::arg("maxncf"));

m.def("IDASetMaxNonlinIters", IDASetMaxNonlinIters, nb::arg("ida_mem"),
      nb::arg("maxcor"));

m.def("IDASetNonlinConvCoef", IDASetNonlinConvCoef, nb::arg("ida_mem"),
      nb::arg("epcon"));

m.def("IDASetNonlinearSolver", IDASetNonlinearSolver, nb::arg("ida_mem"),
      nb::arg("NLS"));

m.def("IDARootInit", IDARootInit, nb::arg("ida_mem"), nb::arg("nrtfn"),
      nb::arg("g"));

m.def(
  "IDASetRootDirection",
  [](void* ida_mem, int rootdir) -> std::tuple<int, int>
  {
    auto IDASetRootDirection_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int rootdir) -> std::tuple<int, int>
    {
      int* rootdir_adapt_modifiable = &rootdir;

      int r = IDASetRootDirection(ida_mem, rootdir_adapt_modifiable);
      return std::make_tuple(r, rootdir);
    };

    return IDASetRootDirection_adapt_modifiable_immutable_to_return(ida_mem,
                                                                    rootdir);
  },
  nb::arg("ida_mem"), nb::arg("rootdir"));

m.def("IDASetNoInactiveRootWarn", IDASetNoInactiveRootWarn, nb::arg("ida_mem"));

m.def(
  "IDASolve",
  [](void* ida_mem, double tout, double tret, N_Vector yret, N_Vector ypret,
     int itask) -> std::tuple<int, double>
  {
    auto IDASolve_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double tout, double tret, N_Vector yret, N_Vector ypret,
         int itask) -> std::tuple<int, double>
    {
      double* tret_adapt_modifiable = &tret;

      int r = IDASolve(ida_mem, tout, tret_adapt_modifiable, yret, ypret, itask);
      return std::make_tuple(r, tret);
    };

    return IDASolve_adapt_modifiable_immutable_to_return(ida_mem, tout, tret,
                                                         yret, ypret, itask);
  },
  nb::arg("ida_mem"), nb::arg("tout"), nb::arg("tret"), nb::arg("yret"),
  nb::arg("ypret"), nb::arg("itask"));

m.def("IDAComputeY", IDAComputeY, nb::arg("ida_mem"), nb::arg("ycor"),
      nb::arg("y"));

m.def("IDAComputeYp", IDAComputeYp, nb::arg("ida_mem"), nb::arg("ycor"),
      nb::arg("yp"));

m.def(
  "IDAComputeYSens",
  [](void* ida_mem, std::vector<N_Vector> ycor, std::vector<N_Vector> yyS) -> int
  {
    auto IDAComputeYSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> ycor,
         std::vector<N_Vector> yyS) -> int
    {
      N_Vector* ycor_ptr =
        reinterpret_cast<N_Vector*>(ycor.empty() ? nullptr : ycor.data());
      N_Vector* yyS_ptr = reinterpret_cast<N_Vector*>(yyS.empty() ? nullptr
                                                                  : yyS.data());

      auto lambda_result = IDAComputeYSens(ida_mem, ycor_ptr, yyS_ptr);
      return lambda_result;
    };

    return IDAComputeYSens_adapt_arr_ptr_to_std_vector(ida_mem, ycor, yyS);
  },
  nb::arg("ida_mem"), nb::arg("ycor"), nb::arg("yyS"));

m.def(
  "IDAComputeYpSens",
  [](void* ida_mem, std::vector<N_Vector> ycor, std::vector<N_Vector> ypS) -> int
  {
    auto IDAComputeYpSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> ycor,
         std::vector<N_Vector> ypS) -> int
    {
      N_Vector* ycor_ptr =
        reinterpret_cast<N_Vector*>(ycor.empty() ? nullptr : ycor.data());
      N_Vector* ypS_ptr = reinterpret_cast<N_Vector*>(ypS.empty() ? nullptr
                                                                  : ypS.data());

      auto lambda_result = IDAComputeYpSens(ida_mem, ycor_ptr, ypS_ptr);
      return lambda_result;
    };

    return IDAComputeYpSens_adapt_arr_ptr_to_std_vector(ida_mem, ycor, ypS);
  },
  nb::arg("ida_mem"), nb::arg("ycor"), nb::arg("ypS"));

m.def("IDAGetDky", IDAGetDky, nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"),
      nb::arg("dky"));

m.def(
  "IDAGetNumSteps",
  [](void* ida_mem, long nsteps) -> std::tuple<int, long>
  {
    auto IDAGetNumSteps_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nsteps) -> std::tuple<int, long>
    {
      long* nsteps_adapt_modifiable = &nsteps;

      int r = IDAGetNumSteps(ida_mem, nsteps_adapt_modifiable);
      return std::make_tuple(r, nsteps);
    };

    return IDAGetNumSteps_adapt_modifiable_immutable_to_return(ida_mem, nsteps);
  },
  nb::arg("ida_mem"), nb::arg("nsteps"));

m.def(
  "IDAGetNumResEvals",
  [](void* ida_mem, long nrevals) -> std::tuple<int, long>
  {
    auto IDAGetNumResEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nrevals) -> std::tuple<int, long>
    {
      long* nrevals_adapt_modifiable = &nrevals;

      int r = IDAGetNumResEvals(ida_mem, nrevals_adapt_modifiable);
      return std::make_tuple(r, nrevals);
    };

    return IDAGetNumResEvals_adapt_modifiable_immutable_to_return(ida_mem,
                                                                  nrevals);
  },
  nb::arg("ida_mem"), nb::arg("nrevals"));

m.def(
  "IDAGetNumLinSolvSetups",
  [](void* ida_mem, long nlinsetups) -> std::tuple<int, long>
  {
    auto IDAGetNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nlinsetups) -> std::tuple<int, long>
    {
      long* nlinsetups_adapt_modifiable = &nlinsetups;

      int r = IDAGetNumLinSolvSetups(ida_mem, nlinsetups_adapt_modifiable);
      return std::make_tuple(r, nlinsetups);
    };

    return IDAGetNumLinSolvSetups_adapt_modifiable_immutable_to_return(ida_mem,
                                                                       nlinsetups);
  },
  nb::arg("ida_mem"), nb::arg("nlinsetups"));

m.def(
  "IDAGetNumErrTestFails",
  [](void* ida_mem, long netfails) -> std::tuple<int, long>
  {
    auto IDAGetNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long netfails) -> std::tuple<int, long>
    {
      long* netfails_adapt_modifiable = &netfails;

      int r = IDAGetNumErrTestFails(ida_mem, netfails_adapt_modifiable);
      return std::make_tuple(r, netfails);
    };

    return IDAGetNumErrTestFails_adapt_modifiable_immutable_to_return(ida_mem,
                                                                      netfails);
  },
  nb::arg("ida_mem"), nb::arg("netfails"));

m.def(
  "IDAGetNumBacktrackOps",
  [](void* ida_mem, long nbacktr) -> std::tuple<int, long>
  {
    auto IDAGetNumBacktrackOps_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nbacktr) -> std::tuple<int, long>
    {
      long* nbacktr_adapt_modifiable = &nbacktr;

      int r = IDAGetNumBacktrackOps(ida_mem, nbacktr_adapt_modifiable);
      return std::make_tuple(r, nbacktr);
    };

    return IDAGetNumBacktrackOps_adapt_modifiable_immutable_to_return(ida_mem,
                                                                      nbacktr);
  },
  nb::arg("ida_mem"), nb::arg("nbacktr"));

m.def("IDAGetConsistentIC", IDAGetConsistentIC, nb::arg("ida_mem"),
      nb::arg("yy0_mod"), nb::arg("yp0_mod"));

m.def(
  "IDAGetLastOrder",
  [](void* ida_mem, int klast) -> std::tuple<int, int>
  {
    auto IDAGetLastOrder_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int klast) -> std::tuple<int, int>
    {
      int* klast_adapt_modifiable = &klast;

      int r = IDAGetLastOrder(ida_mem, klast_adapt_modifiable);
      return std::make_tuple(r, klast);
    };

    return IDAGetLastOrder_adapt_modifiable_immutable_to_return(ida_mem, klast);
  },
  nb::arg("ida_mem"), nb::arg("klast"));

m.def(
  "IDAGetCurrentOrder",
  [](void* ida_mem, int kcur) -> std::tuple<int, int>
  {
    auto IDAGetCurrentOrder_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int kcur) -> std::tuple<int, int>
    {
      int* kcur_adapt_modifiable = &kcur;

      int r = IDAGetCurrentOrder(ida_mem, kcur_adapt_modifiable);
      return std::make_tuple(r, kcur);
    };

    return IDAGetCurrentOrder_adapt_modifiable_immutable_to_return(ida_mem, kcur);
  },
  nb::arg("ida_mem"), nb::arg("kcur"));

m.def(
  "IDAGetCurrentCj",
  [](void* ida_mem, double cj) -> std::tuple<int, double>
  {
    auto IDAGetCurrentCj_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double cj) -> std::tuple<int, double>
    {
      double* cj_adapt_modifiable = &cj;

      int r = IDAGetCurrentCj(ida_mem, cj_adapt_modifiable);
      return std::make_tuple(r, cj);
    };

    return IDAGetCurrentCj_adapt_modifiable_immutable_to_return(ida_mem, cj);
  },
  nb::arg("ida_mem"), nb::arg("cj"));

m.def(
  "IDAGetActualInitStep",
  [](void* ida_mem, double hinused) -> std::tuple<int, double>
  {
    auto IDAGetActualInitStep_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double hinused) -> std::tuple<int, double>
    {
      double* hinused_adapt_modifiable = &hinused;

      int r = IDAGetActualInitStep(ida_mem, hinused_adapt_modifiable);
      return std::make_tuple(r, hinused);
    };

    return IDAGetActualInitStep_adapt_modifiable_immutable_to_return(ida_mem,
                                                                     hinused);
  },
  nb::arg("ida_mem"), nb::arg("hinused"));

m.def(
  "IDAGetLastStep",
  [](void* ida_mem, double hlast) -> std::tuple<int, double>
  {
    auto IDAGetLastStep_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double hlast) -> std::tuple<int, double>
    {
      double* hlast_adapt_modifiable = &hlast;

      int r = IDAGetLastStep(ida_mem, hlast_adapt_modifiable);
      return std::make_tuple(r, hlast);
    };

    return IDAGetLastStep_adapt_modifiable_immutable_to_return(ida_mem, hlast);
  },
  nb::arg("ida_mem"), nb::arg("hlast"));

m.def(
  "IDAGetCurrentStep",
  [](void* ida_mem, double hcur) -> std::tuple<int, double>
  {
    auto IDAGetCurrentStep_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double hcur) -> std::tuple<int, double>
    {
      double* hcur_adapt_modifiable = &hcur;

      int r = IDAGetCurrentStep(ida_mem, hcur_adapt_modifiable);
      return std::make_tuple(r, hcur);
    };

    return IDAGetCurrentStep_adapt_modifiable_immutable_to_return(ida_mem, hcur);
  },
  nb::arg("ida_mem"), nb::arg("hcur"));

m.def(
  "IDAGetCurrentTime",
  [](void* ida_mem, double tcur) -> std::tuple<int, double>
  {
    auto IDAGetCurrentTime_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double tcur) -> std::tuple<int, double>
    {
      double* tcur_adapt_modifiable = &tcur;

      int r = IDAGetCurrentTime(ida_mem, tcur_adapt_modifiable);
      return std::make_tuple(r, tcur);
    };

    return IDAGetCurrentTime_adapt_modifiable_immutable_to_return(ida_mem, tcur);
  },
  nb::arg("ida_mem"), nb::arg("tcur"));

m.def(
  "IDAGetTolScaleFactor",
  [](void* ida_mem, double tolsfact) -> std::tuple<int, double>
  {
    auto IDAGetTolScaleFactor_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double tolsfact) -> std::tuple<int, double>
    {
      double* tolsfact_adapt_modifiable = &tolsfact;

      int r = IDAGetTolScaleFactor(ida_mem, tolsfact_adapt_modifiable);
      return std::make_tuple(r, tolsfact);
    };

    return IDAGetTolScaleFactor_adapt_modifiable_immutable_to_return(ida_mem,
                                                                     tolsfact);
  },
  nb::arg("ida_mem"), nb::arg("tolsfact"));

m.def("IDAGetErrWeights", IDAGetErrWeights, nb::arg("ida_mem"),
      nb::arg("eweight"));

m.def("IDAGetEstLocalErrors", IDAGetEstLocalErrors, nb::arg("ida_mem"),
      nb::arg("ele"));

m.def(
  "IDAGetNumGEvals",
  [](void* ida_mem, long ngevals) -> std::tuple<int, long>
  {
    auto IDAGetNumGEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long ngevals) -> std::tuple<int, long>
    {
      long* ngevals_adapt_modifiable = &ngevals;

      int r = IDAGetNumGEvals(ida_mem, ngevals_adapt_modifiable);
      return std::make_tuple(r, ngevals);
    };

    return IDAGetNumGEvals_adapt_modifiable_immutable_to_return(ida_mem, ngevals);
  },
  nb::arg("ida_mem"), nb::arg("ngevals"));

m.def(
  "IDAGetRootInfo",
  [](void* ida_mem, int rootsfound) -> std::tuple<int, int>
  {
    auto IDAGetRootInfo_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int rootsfound) -> std::tuple<int, int>
    {
      int* rootsfound_adapt_modifiable = &rootsfound;

      int r = IDAGetRootInfo(ida_mem, rootsfound_adapt_modifiable);
      return std::make_tuple(r, rootsfound);
    };

    return IDAGetRootInfo_adapt_modifiable_immutable_to_return(ida_mem,
                                                               rootsfound);
  },
  nb::arg("ida_mem"), nb::arg("rootsfound"));

m.def(
  "IDAGetIntegratorStats",
  [](void* ida_mem, long nsteps, long nrevals, long nlinsetups, long netfails,
     int qlast, int qcur, double hinused, double hlast, double hcur, double tcur)
    -> std::tuple<int, long, long, long, long, int, int, double, double, double, double>
  {
    auto IDAGetIntegratorStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nsteps, long nrevals, long nlinsetups, long netfails,
         int qlast, int qcur, double hinused, double hlast, double hcur,
         double tcur) -> std::tuple<int, long, long, long, long, int, int,
                                    double, double, double, double>
    {
      long* nsteps_adapt_modifiable     = &nsteps;
      long* nrevals_adapt_modifiable    = &nrevals;
      long* nlinsetups_adapt_modifiable = &nlinsetups;
      long* netfails_adapt_modifiable   = &netfails;
      int* qlast_adapt_modifiable       = &qlast;
      int* qcur_adapt_modifiable        = &qcur;
      double* hinused_adapt_modifiable  = &hinused;
      double* hlast_adapt_modifiable    = &hlast;
      double* hcur_adapt_modifiable     = &hcur;
      double* tcur_adapt_modifiable     = &tcur;

      int r =
        IDAGetIntegratorStats(ida_mem, nsteps_adapt_modifiable,
                              nrevals_adapt_modifiable,
                              nlinsetups_adapt_modifiable,
                              netfails_adapt_modifiable, qlast_adapt_modifiable,
                              qcur_adapt_modifiable, hinused_adapt_modifiable,
                              hlast_adapt_modifiable, hcur_adapt_modifiable,
                              tcur_adapt_modifiable);
      return std::make_tuple(r, nsteps, nrevals, nlinsetups, netfails, qlast,
                             qcur, hinused, hlast, hcur, tcur);
    };

    return IDAGetIntegratorStats_adapt_modifiable_immutable_to_return(ida_mem,
                                                                      nsteps,
                                                                      nrevals,
                                                                      nlinsetups,
                                                                      netfails,
                                                                      qlast, qcur,
                                                                      hinused,
                                                                      hlast,
                                                                      hcur, tcur);
  },
  nb::arg("ida_mem"), nb::arg("nsteps"), nb::arg("nrevals"),
  nb::arg("nlinsetups"), nb::arg("netfails"), nb::arg("qlast"), nb::arg("qcur"),
  nb::arg("hinused"), nb::arg("hlast"), nb::arg("hcur"), nb::arg("tcur"));

m.def(
  "IDAGetNumNonlinSolvIters",
  [](void* ida_mem, long nniters) -> std::tuple<int, long>
  {
    auto IDAGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nniters) -> std::tuple<int, long>
    {
      long* nniters_adapt_modifiable = &nniters;

      int r = IDAGetNumNonlinSolvIters(ida_mem, nniters_adapt_modifiable);
      return std::make_tuple(r, nniters);
    };

    return IDAGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(ida_mem,
                                                                         nniters);
  },
  nb::arg("ida_mem"), nb::arg("nniters"));

m.def(
  "IDAGetNumNonlinSolvConvFails",
  [](void* ida_mem, long nnfails) -> std::tuple<int, long>
  {
    auto IDAGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nnfails) -> std::tuple<int, long>
    {
      long* nnfails_adapt_modifiable = &nnfails;

      int r = IDAGetNumNonlinSolvConvFails(ida_mem, nnfails_adapt_modifiable);
      return std::make_tuple(r, nnfails);
    };

    return IDAGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(ida_mem,
                                                                             nnfails);
  },
  nb::arg("ida_mem"), nb::arg("nnfails"));

m.def(
  "IDAGetNonlinSolvStats",
  [](void* ida_mem, long nniters, long nnfails) -> std::tuple<int, long, long>
  {
    auto IDAGetNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nniters, long nnfails) -> std::tuple<int, long, long>
    {
      long* nniters_adapt_modifiable = &nniters;
      long* nnfails_adapt_modifiable = &nnfails;

      int r = IDAGetNonlinSolvStats(ida_mem, nniters_adapt_modifiable,
                                    nnfails_adapt_modifiable);
      return std::make_tuple(r, nniters, nnfails);
    };

    return IDAGetNonlinSolvStats_adapt_modifiable_immutable_to_return(ida_mem,
                                                                      nniters,
                                                                      nnfails);
  },
  nb::arg("ida_mem"), nb::arg("nniters"), nb::arg("nnfails"));

m.def(
  "IDAGetNumStepSolveFails",
  [](void* ida_mem, long nncfails) -> std::tuple<int, long>
  {
    auto IDAGetNumStepSolveFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nncfails) -> std::tuple<int, long>
    {
      long* nncfails_adapt_modifiable = &nncfails;

      int r = IDAGetNumStepSolveFails(ida_mem, nncfails_adapt_modifiable);
      return std::make_tuple(r, nncfails);
    };

    return IDAGetNumStepSolveFails_adapt_modifiable_immutable_to_return(ida_mem,
                                                                        nncfails);
  },
  nb::arg("ida_mem"), nb::arg("nncfails"));

m.def("IDAPrintAllStats", IDAPrintAllStats, nb::arg("ida_mem"),
      nb::arg("outfile"), nb::arg("fmt"));

m.def("IDAGetReturnFlagName", IDAGetReturnFlagName, nb::arg("flag"));

m.def("IDAQuadReInit", IDAQuadReInit, nb::arg("ida_mem"), nb::arg("yQ0"));

m.def("IDAQuadSStolerances", IDAQuadSStolerances, nb::arg("ida_mem"),
      nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("IDAQuadSVtolerances", IDAQuadSVtolerances, nb::arg("ida_mem"),
      nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("IDASetQuadErrCon", IDASetQuadErrCon, nb::arg("ida_mem"),
      nb::arg("errconQ"));

m.def(
  "IDAGetQuad",
  [](void* ida_mem, double t, N_Vector yQout) -> std::tuple<int, double>
  {
    auto IDAGetQuad_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double t, N_Vector yQout) -> std::tuple<int, double>
    {
      double* t_adapt_modifiable = &t;

      int r = IDAGetQuad(ida_mem, t_adapt_modifiable, yQout);
      return std::make_tuple(r, t);
    };

    return IDAGetQuad_adapt_modifiable_immutable_to_return(ida_mem, t, yQout);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("yQout"));

m.def("IDAGetQuadDky", IDAGetQuadDky, nb::arg("ida_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("dky"));

m.def(
  "IDAGetQuadNumRhsEvals",
  [](void* ida_mem, long nrhsQevals) -> std::tuple<int, long>
  {
    auto IDAGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nrhsQevals) -> std::tuple<int, long>
    {
      long* nrhsQevals_adapt_modifiable = &nrhsQevals;

      int r = IDAGetQuadNumRhsEvals(ida_mem, nrhsQevals_adapt_modifiable);
      return std::make_tuple(r, nrhsQevals);
    };

    return IDAGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return(ida_mem,
                                                                      nrhsQevals);
  },
  nb::arg("ida_mem"), nb::arg("nrhsQevals"));

m.def(
  "IDAGetQuadNumErrTestFails",
  [](void* ida_mem, long nQetfails) -> std::tuple<int, long>
  {
    auto IDAGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nQetfails) -> std::tuple<int, long>
    {
      long* nQetfails_adapt_modifiable = &nQetfails;

      int r = IDAGetQuadNumErrTestFails(ida_mem, nQetfails_adapt_modifiable);
      return std::make_tuple(r, nQetfails);
    };

    return IDAGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return(ida_mem,
                                                                          nQetfails);
  },
  nb::arg("ida_mem"), nb::arg("nQetfails"));

m.def("IDAGetQuadErrWeights", IDAGetQuadErrWeights, nb::arg("ida_mem"),
      nb::arg("eQweight"));

m.def(
  "IDAGetQuadStats",
  [](void* ida_mem, long nrhsQevals, long nQetfails) -> std::tuple<int, long, long>
  {
    auto IDAGetQuadStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nrhsQevals,
         long nQetfails) -> std::tuple<int, long, long>
    {
      long* nrhsQevals_adapt_modifiable = &nrhsQevals;
      long* nQetfails_adapt_modifiable  = &nQetfails;

      int r = IDAGetQuadStats(ida_mem, nrhsQevals_adapt_modifiable,
                              nQetfails_adapt_modifiable);
      return std::make_tuple(r, nrhsQevals, nQetfails);
    };

    return IDAGetQuadStats_adapt_modifiable_immutable_to_return(ida_mem,
                                                                nrhsQevals,
                                                                nQetfails);
  },
  nb::arg("ida_mem"), nb::arg("nrhsQevals"), nb::arg("nQetfails"));

m.def(
  "IDASensReInit",
  [](void* ida_mem, int ism, std::vector<N_Vector> yS0,
     std::vector<N_Vector> ypS0) -> int
  {
    auto IDASensReInit_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, int ism, std::vector<N_Vector> yS0,
         std::vector<N_Vector> ypS0) -> int
    {
      N_Vector* yS0_ptr = reinterpret_cast<N_Vector*>(yS0.empty() ? nullptr
                                                                  : yS0.data());
      N_Vector* ypS0_ptr =
        reinterpret_cast<N_Vector*>(ypS0.empty() ? nullptr : ypS0.data());

      auto lambda_result = IDASensReInit(ida_mem, ism, yS0_ptr, ypS0_ptr);
      return lambda_result;
    };

    return IDASensReInit_adapt_arr_ptr_to_std_vector(ida_mem, ism, yS0, ypS0);
  },
  nb::arg("ida_mem"), nb::arg("ism"), nb::arg("yS0"), nb::arg("ypS0"));

m.def(
  "IDASensSStolerances",
  [](void* ida_mem, double reltolS, double abstolS) -> std::tuple<int, double>
  {
    auto IDASensSStolerances_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double reltolS, double abstolS) -> std::tuple<int, double>
    {
      double* abstolS_adapt_modifiable = &abstolS;

      int r = IDASensSStolerances(ida_mem, reltolS, abstolS_adapt_modifiable);
      return std::make_tuple(r, abstolS);
    };

    return IDASensSStolerances_adapt_modifiable_immutable_to_return(ida_mem,
                                                                    reltolS,
                                                                    abstolS);
  },
  nb::arg("ida_mem"), nb::arg("reltolS"), nb::arg("abstolS"));

m.def(
  "IDASensSVtolerances",
  [](void* ida_mem, double reltolS, std::vector<N_Vector> abstolS) -> int
  {
    auto IDASensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, double reltolS, std::vector<N_Vector> abstolS) -> int
    {
      N_Vector* abstolS_ptr =
        reinterpret_cast<N_Vector*>(abstolS.empty() ? nullptr : abstolS.data());

      auto lambda_result = IDASensSVtolerances(ida_mem, reltolS, abstolS_ptr);
      return lambda_result;
    };

    return IDASensSVtolerances_adapt_arr_ptr_to_std_vector(ida_mem, reltolS,
                                                           abstolS);
  },
  nb::arg("ida_mem"), nb::arg("reltolS"), nb::arg("abstolS"));

m.def("IDASensEEtolerances", IDASensEEtolerances, nb::arg("ida_mem"));

m.def(
  "IDAGetSensConsistentIC",
  [](void* ida_mem, std::vector<N_Vector> yyS0, std::vector<N_Vector> ypS0) -> int
  {
    auto IDAGetSensConsistentIC_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> yyS0,
         std::vector<N_Vector> ypS0) -> int
    {
      N_Vector* yyS0_ptr =
        reinterpret_cast<N_Vector*>(yyS0.empty() ? nullptr : yyS0.data());
      N_Vector* ypS0_ptr =
        reinterpret_cast<N_Vector*>(ypS0.empty() ? nullptr : ypS0.data());

      auto lambda_result = IDAGetSensConsistentIC(ida_mem, yyS0_ptr, ypS0_ptr);
      return lambda_result;
    };

    return IDAGetSensConsistentIC_adapt_arr_ptr_to_std_vector(ida_mem, yyS0,
                                                              ypS0);
  },
  nb::arg("ida_mem"), nb::arg("yyS0"), nb::arg("ypS0"));

m.def("IDASetSensDQMethod", IDASetSensDQMethod, nb::arg("ida_mem"),
      nb::arg("DQtype"), nb::arg("DQrhomax"));

m.def("IDASetSensErrCon", IDASetSensErrCon, nb::arg("ida_mem"),
      nb::arg("errconS"));

m.def("IDASetSensMaxNonlinIters", IDASetSensMaxNonlinIters, nb::arg("ida_mem"),
      nb::arg("maxcorS"));

m.def(
  "IDASetSensParams",
  [](void* ida_mem, std::vector<sunrealtype> p, std::vector<sunrealtype> pbar,
     std::vector<int> plist) -> int
  {
    auto IDASetSensParams_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<sunrealtype> p,
         std::vector<sunrealtype> pbar, std::vector<int> plist) -> int
    {
      sunrealtype* p_ptr = reinterpret_cast<sunrealtype*>(p.empty() ? nullptr
                                                                    : p.data());
      sunrealtype* pbar_ptr =
        reinterpret_cast<sunrealtype*>(pbar.empty() ? nullptr : pbar.data());
      int* plist_ptr = reinterpret_cast<int*>(plist.empty() ? nullptr
                                                            : plist.data());

      auto lambda_result = IDASetSensParams(ida_mem, p_ptr, pbar_ptr, plist_ptr);
      return lambda_result;
    };

    return IDASetSensParams_adapt_arr_ptr_to_std_vector(ida_mem, p, pbar, plist);
  },
  nb::arg("ida_mem"), nb::arg("p"), nb::arg("pbar"), nb::arg("plist"));

m.def("IDASetNonlinearSolverSensSim", IDASetNonlinearSolverSensSim,
      nb::arg("ida_mem"), nb::arg("NLS"));

m.def("IDASetNonlinearSolverSensStg", IDASetNonlinearSolverSensStg,
      nb::arg("ida_mem"), nb::arg("NLS"));

m.def("IDASensToggleOff", IDASensToggleOff, nb::arg("ida_mem"));

m.def(
  "IDAGetSens",
  [](void* ida_mem, double tret,
     std::vector<N_Vector> yySout) -> std::tuple<int, double>
  {
    auto IDAGetSens_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double tret, N_Vector1d yySout) -> std::tuple<int, double>
    {
      double* tret_adapt_modifiable = &tret;

      int r = IDAGetSens(ida_mem, tret_adapt_modifiable, yySout);
      return std::make_tuple(r, tret);
    };
    auto IDAGetSens_adapt_arr_ptr_to_std_vector =
      [&IDAGetSens_adapt_modifiable_immutable_to_return](void* ida_mem,
                                                         double tret,
                                                         std::vector<N_Vector> yySout)
      -> std::tuple<int, double>
    {
      N_Vector* yySout_ptr =
        reinterpret_cast<N_Vector*>(yySout.empty() ? nullptr : yySout.data());

      auto lambda_result =
        IDAGetSens_adapt_modifiable_immutable_to_return(ida_mem, tret,
                                                        yySout_ptr);
      return lambda_result;
    };

    return IDAGetSens_adapt_arr_ptr_to_std_vector(ida_mem, tret, yySout);
  },
  nb::arg("ida_mem"), nb::arg("tret"), nb::arg("yySout"));

m.def(
  "IDAGetSens1",
  [](void* ida_mem, double tret, int is, N_Vector yySret) -> std::tuple<int, double>
  {
    auto IDAGetSens1_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double tret, int is,
         N_Vector yySret) -> std::tuple<int, double>
    {
      double* tret_adapt_modifiable = &tret;

      int r = IDAGetSens1(ida_mem, tret_adapt_modifiable, is, yySret);
      return std::make_tuple(r, tret);
    };

    return IDAGetSens1_adapt_modifiable_immutable_to_return(ida_mem, tret, is,
                                                            yySret);
  },
  nb::arg("ida_mem"), nb::arg("tret"), nb::arg("is_"), nb::arg("yySret"));

m.def(
  "IDAGetSensDky",
  [](void* ida_mem, double t, int k, std::vector<N_Vector> dkyS) -> int
  {
    auto IDAGetSensDky_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, double t, int k, std::vector<N_Vector> dkyS) -> int
    {
      N_Vector* dkyS_ptr =
        reinterpret_cast<N_Vector*>(dkyS.empty() ? nullptr : dkyS.data());

      auto lambda_result = IDAGetSensDky(ida_mem, t, k, dkyS_ptr);
      return lambda_result;
    };

    return IDAGetSensDky_adapt_arr_ptr_to_std_vector(ida_mem, t, k, dkyS);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyS"));

m.def("IDAGetSensDky1", IDAGetSensDky1, nb::arg("ida_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("is_"), nb::arg("dkyS"));

m.def(
  "IDAGetSensNumResEvals",
  [](void* ida_mem, long nresSevals) -> std::tuple<int, long>
  {
    auto IDAGetSensNumResEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nresSevals) -> std::tuple<int, long>
    {
      long* nresSevals_adapt_modifiable = &nresSevals;

      int r = IDAGetSensNumResEvals(ida_mem, nresSevals_adapt_modifiable);
      return std::make_tuple(r, nresSevals);
    };

    return IDAGetSensNumResEvals_adapt_modifiable_immutable_to_return(ida_mem,
                                                                      nresSevals);
  },
  nb::arg("ida_mem"), nb::arg("nresSevals"));

m.def(
  "IDAGetNumResEvalsSens",
  [](void* ida_mem, long nresevalsS) -> std::tuple<int, long>
  {
    auto IDAGetNumResEvalsSens_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nresevalsS) -> std::tuple<int, long>
    {
      long* nresevalsS_adapt_modifiable = &nresevalsS;

      int r = IDAGetNumResEvalsSens(ida_mem, nresevalsS_adapt_modifiable);
      return std::make_tuple(r, nresevalsS);
    };

    return IDAGetNumResEvalsSens_adapt_modifiable_immutable_to_return(ida_mem,
                                                                      nresevalsS);
  },
  nb::arg("ida_mem"), nb::arg("nresevalsS"));

m.def(
  "IDAGetSensNumErrTestFails",
  [](void* ida_mem, long nSetfails) -> std::tuple<int, long>
  {
    auto IDAGetSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nSetfails) -> std::tuple<int, long>
    {
      long* nSetfails_adapt_modifiable = &nSetfails;

      int r = IDAGetSensNumErrTestFails(ida_mem, nSetfails_adapt_modifiable);
      return std::make_tuple(r, nSetfails);
    };

    return IDAGetSensNumErrTestFails_adapt_modifiable_immutable_to_return(ida_mem,
                                                                          nSetfails);
  },
  nb::arg("ida_mem"), nb::arg("nSetfails"));

m.def(
  "IDAGetSensNumLinSolvSetups",
  [](void* ida_mem, long nlinsetupsS) -> std::tuple<int, long>
  {
    auto IDAGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nlinsetupsS) -> std::tuple<int, long>
    {
      long* nlinsetupsS_adapt_modifiable = &nlinsetupsS;

      int r = IDAGetSensNumLinSolvSetups(ida_mem, nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nlinsetupsS);
    };

    return IDAGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return(ida_mem,
                                                                           nlinsetupsS);
  },
  nb::arg("ida_mem"), nb::arg("nlinsetupsS"));

m.def(
  "IDAGetSensErrWeights",
  [](void* ida_mem, std::vector<N_Vector> eSweight) -> int
  {
    auto IDAGetSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> eSweight) -> int
    {
      N_Vector* eSweight_ptr = reinterpret_cast<N_Vector*>(
        eSweight.empty() ? nullptr : eSweight.data());

      auto lambda_result = IDAGetSensErrWeights(ida_mem, eSweight_ptr);
      return lambda_result;
    };

    return IDAGetSensErrWeights_adapt_arr_ptr_to_std_vector(ida_mem, eSweight);
  },
  nb::arg("ida_mem"), nb::arg("eSweight"));

m.def(
  "IDAGetSensStats",
  [](void* ida_mem, long nresSevals, long nresevalsS, long nSetfails,
     long nlinsetupsS) -> std::tuple<int, long, long, long, long>
  {
    auto IDAGetSensStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nresSevals, long nresevalsS, long nSetfails,
         long nlinsetupsS) -> std::tuple<int, long, long, long, long>
    {
      long* nresSevals_adapt_modifiable  = &nresSevals;
      long* nresevalsS_adapt_modifiable  = &nresevalsS;
      long* nSetfails_adapt_modifiable   = &nSetfails;
      long* nlinsetupsS_adapt_modifiable = &nlinsetupsS;

      int r = IDAGetSensStats(ida_mem, nresSevals_adapt_modifiable,
                              nresevalsS_adapt_modifiable,
                              nSetfails_adapt_modifiable,
                              nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nresSevals, nresevalsS, nSetfails, nlinsetupsS);
    };

    return IDAGetSensStats_adapt_modifiable_immutable_to_return(ida_mem,
                                                                nresSevals,
                                                                nresevalsS,
                                                                nSetfails,
                                                                nlinsetupsS);
  },
  nb::arg("ida_mem"), nb::arg("nresSevals"), nb::arg("nresevalsS"),
  nb::arg("nSetfails"), nb::arg("nlinsetupsS"));

m.def(
  "IDAGetSensNumNonlinSolvIters",
  [](void* ida_mem, long nSniters) -> std::tuple<int, long>
  {
    auto IDAGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nSniters) -> std::tuple<int, long>
    {
      long* nSniters_adapt_modifiable = &nSniters;

      int r = IDAGetSensNumNonlinSolvIters(ida_mem, nSniters_adapt_modifiable);
      return std::make_tuple(r, nSniters);
    };

    return IDAGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return(ida_mem,
                                                                             nSniters);
  },
  nb::arg("ida_mem"), nb::arg("nSniters"));

m.def(
  "IDAGetSensNumNonlinSolvConvFails",
  [](void* ida_mem, long nSnfails) -> std::tuple<int, long>
  {
    auto IDAGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nSnfails) -> std::tuple<int, long>
    {
      long* nSnfails_adapt_modifiable = &nSnfails;

      int r = IDAGetSensNumNonlinSolvConvFails(ida_mem,
                                               nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSnfails);
    };

    return IDAGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(ida_mem,
                                                                                 nSnfails);
  },
  nb::arg("ida_mem"), nb::arg("nSnfails"));

m.def(
  "IDAGetSensNonlinSolvStats",
  [](void* ida_mem, long nSniters, long nSnfails) -> std::tuple<int, long, long>
  {
    auto IDAGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nSniters,
         long nSnfails) -> std::tuple<int, long, long>
    {
      long* nSniters_adapt_modifiable = &nSniters;
      long* nSnfails_adapt_modifiable = &nSnfails;

      int r = IDAGetSensNonlinSolvStats(ida_mem, nSniters_adapt_modifiable,
                                        nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSniters, nSnfails);
    };

    return IDAGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return(ida_mem,
                                                                          nSniters,
                                                                          nSnfails);
  },
  nb::arg("ida_mem"), nb::arg("nSniters"), nb::arg("nSnfails"));

m.def(
  "IDAGetNumStepSensSolveFails",
  [](void* ida_mem, long nSncfails) -> std::tuple<int, long>
  {
    auto IDAGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nSncfails) -> std::tuple<int, long>
    {
      long* nSncfails_adapt_modifiable = &nSncfails;

      int r = IDAGetNumStepSensSolveFails(ida_mem, nSncfails_adapt_modifiable);
      return std::make_tuple(r, nSncfails);
    };

    return IDAGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return(ida_mem,
                                                                            nSncfails);
  },
  nb::arg("ida_mem"), nb::arg("nSncfails"));

m.def(
  "IDAQuadSensReInit",
  [](void* ida_mem, std::vector<N_Vector> yQS0) -> int
  {
    auto IDAQuadSensReInit_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> yQS0) -> int
    {
      N_Vector* yQS0_ptr =
        reinterpret_cast<N_Vector*>(yQS0.empty() ? nullptr : yQS0.data());

      auto lambda_result = IDAQuadSensReInit(ida_mem, yQS0_ptr);
      return lambda_result;
    };

    return IDAQuadSensReInit_adapt_arr_ptr_to_std_vector(ida_mem, yQS0);
  },
  nb::arg("ida_mem"), nb::arg("yQS0"));

m.def(
  "IDAQuadSensSStolerances",
  [](void* ida_mem, double reltolQS, double abstolQS) -> std::tuple<int, double>
  {
    auto IDAQuadSensSStolerances_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double reltolQS,
         double abstolQS) -> std::tuple<int, double>
    {
      double* abstolQS_adapt_modifiable = &abstolQS;

      int r = IDAQuadSensSStolerances(ida_mem, reltolQS,
                                      abstolQS_adapt_modifiable);
      return std::make_tuple(r, abstolQS);
    };

    return IDAQuadSensSStolerances_adapt_modifiable_immutable_to_return(ida_mem,
                                                                        reltolQS,
                                                                        abstolQS);
  },
  nb::arg("ida_mem"), nb::arg("reltolQS"), nb::arg("abstolQS"));

m.def(
  "IDAQuadSensSVtolerances",
  [](void* ida_mem, double reltolQS, std::vector<N_Vector> abstolQS) -> int
  {
    auto IDAQuadSensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, double reltolQS, std::vector<N_Vector> abstolQS) -> int
    {
      N_Vector* abstolQS_ptr = reinterpret_cast<N_Vector*>(
        abstolQS.empty() ? nullptr : abstolQS.data());

      auto lambda_result = IDAQuadSensSVtolerances(ida_mem, reltolQS,
                                                   abstolQS_ptr);
      return lambda_result;
    };

    return IDAQuadSensSVtolerances_adapt_arr_ptr_to_std_vector(ida_mem, reltolQS,
                                                               abstolQS);
  },
  nb::arg("ida_mem"), nb::arg("reltolQS"), nb::arg("abstolQS"));

m.def("IDAQuadSensEEtolerances", IDAQuadSensEEtolerances, nb::arg("ida_mem"));

m.def("IDASetQuadSensErrCon", IDASetQuadSensErrCon, nb::arg("ida_mem"),
      nb::arg("errconQS"));

m.def(
  "IDAGetQuadSens",
  [](void* ida_mem, double tret,
     std::vector<N_Vector> yyQSout) -> std::tuple<int, double>
  {
    auto IDAGetQuadSens_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double tret, N_Vector1d yyQSout) -> std::tuple<int, double>
    {
      double* tret_adapt_modifiable = &tret;

      int r = IDAGetQuadSens(ida_mem, tret_adapt_modifiable, yyQSout);
      return std::make_tuple(r, tret);
    };
    auto IDAGetQuadSens_adapt_arr_ptr_to_std_vector =
      [&IDAGetQuadSens_adapt_modifiable_immutable_to_return](void* ida_mem,
                                                             double tret,
                                                             std::vector<N_Vector> yyQSout)
      -> std::tuple<int, double>
    {
      N_Vector* yyQSout_ptr =
        reinterpret_cast<N_Vector*>(yyQSout.empty() ? nullptr : yyQSout.data());

      auto lambda_result =
        IDAGetQuadSens_adapt_modifiable_immutable_to_return(ida_mem, tret,
                                                            yyQSout_ptr);
      return lambda_result;
    };

    return IDAGetQuadSens_adapt_arr_ptr_to_std_vector(ida_mem, tret, yyQSout);
  },
  nb::arg("ida_mem"), nb::arg("tret"), nb::arg("yyQSout"));

m.def(
  "IDAGetQuadSens1",
  [](void* ida_mem, double tret, int is,
     N_Vector yyQSret) -> std::tuple<int, double>
  {
    auto IDAGetQuadSens1_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double tret, int is,
         N_Vector yyQSret) -> std::tuple<int, double>
    {
      double* tret_adapt_modifiable = &tret;

      int r = IDAGetQuadSens1(ida_mem, tret_adapt_modifiable, is, yyQSret);
      return std::make_tuple(r, tret);
    };

    return IDAGetQuadSens1_adapt_modifiable_immutable_to_return(ida_mem, tret,
                                                                is, yyQSret);
  },
  nb::arg("ida_mem"), nb::arg("tret"), nb::arg("is_"), nb::arg("yyQSret"));

m.def(
  "IDAGetQuadSensDky",
  [](void* ida_mem, double t, int k, std::vector<N_Vector> dkyQS) -> int
  {
    auto IDAGetQuadSensDky_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, double t, int k, std::vector<N_Vector> dkyQS) -> int
    {
      N_Vector* dkyQS_ptr =
        reinterpret_cast<N_Vector*>(dkyQS.empty() ? nullptr : dkyQS.data());

      auto lambda_result = IDAGetQuadSensDky(ida_mem, t, k, dkyQS_ptr);
      return lambda_result;
    };

    return IDAGetQuadSensDky_adapt_arr_ptr_to_std_vector(ida_mem, t, k, dkyQS);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyQS"));

m.def("IDAGetQuadSensDky1", IDAGetQuadSensDky1, nb::arg("ida_mem"),
      nb::arg("t"), nb::arg("k"), nb::arg("is_"), nb::arg("dkyQS"));

m.def(
  "IDAGetQuadSensNumRhsEvals",
  [](void* ida_mem, long nrhsQSevals) -> std::tuple<int, long>
  {
    auto IDAGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nrhsQSevals) -> std::tuple<int, long>
    {
      long* nrhsQSevals_adapt_modifiable = &nrhsQSevals;

      int r = IDAGetQuadSensNumRhsEvals(ida_mem, nrhsQSevals_adapt_modifiable);
      return std::make_tuple(r, nrhsQSevals);
    };

    return IDAGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return(ida_mem,
                                                                          nrhsQSevals);
  },
  nb::arg("ida_mem"), nb::arg("nrhsQSevals"));

m.def(
  "IDAGetQuadSensNumErrTestFails",
  [](void* ida_mem, long nQSetfails) -> std::tuple<int, long>
  {
    auto IDAGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nQSetfails) -> std::tuple<int, long>
    {
      long* nQSetfails_adapt_modifiable = &nQSetfails;

      int r = IDAGetQuadSensNumErrTestFails(ida_mem, nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nQSetfails);
    };

    return IDAGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return(ida_mem,
                                                                              nQSetfails);
  },
  nb::arg("ida_mem"), nb::arg("nQSetfails"));

m.def(
  "IDAGetQuadSensErrWeights",
  [](void* ida_mem, std::vector<N_Vector> eQSweight) -> int
  {
    auto IDAGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> eQSweight) -> int
    {
      N_Vector* eQSweight_ptr = reinterpret_cast<N_Vector*>(
        eQSweight.empty() ? nullptr : eQSweight.data());

      auto lambda_result = IDAGetQuadSensErrWeights(ida_mem, eQSweight_ptr);
      return lambda_result;
    };

    return IDAGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector(ida_mem,
                                                                eQSweight);
  },
  nb::arg("ida_mem"), nb::arg("eQSweight"));

m.def(
  "IDAGetQuadSensStats",
  [](void* ida_mem, long nrhsQSevals,
     long nQSetfails) -> std::tuple<int, long, long>
  {
    auto IDAGetQuadSensStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nrhsQSevals,
         long nQSetfails) -> std::tuple<int, long, long>
    {
      long* nrhsQSevals_adapt_modifiable = &nrhsQSevals;
      long* nQSetfails_adapt_modifiable  = &nQSetfails;

      int r = IDAGetQuadSensStats(ida_mem, nrhsQSevals_adapt_modifiable,
                                  nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nrhsQSevals, nQSetfails);
    };

    return IDAGetQuadSensStats_adapt_modifiable_immutable_to_return(ida_mem,
                                                                    nrhsQSevals,
                                                                    nQSetfails);
  },
  nb::arg("ida_mem"), nb::arg("nrhsQSevals"), nb::arg("nQSetfails"));

m.def("IDAAdjInit", IDAAdjInit, nb::arg("ida_mem"), nb::arg("steps"),
      nb::arg("interp"));

m.def("IDAAdjReInit", IDAAdjReInit, nb::arg("ida_mem"));

m.def(
  "IDACreateB",
  [](void* ida_mem, int which) -> std::tuple<int, int>
  {
    auto IDACreateB_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which) -> std::tuple<int, int>
    {
      int* which_adapt_modifiable = &which;

      int r = IDACreateB(ida_mem, which_adapt_modifiable);
      return std::make_tuple(r, which);
    };

    return IDACreateB_adapt_modifiable_immutable_to_return(ida_mem, which);
  },
  nb::arg("ida_mem"), nb::arg("which"));

m.def("IDAReInitB", IDAReInitB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("tB0"), nb::arg("yyB0"), nb::arg("ypB0"));

m.def("IDASStolerancesB", IDASStolerancesB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("relTolB"), nb::arg("absTolB"));

m.def("IDASVtolerancesB", IDASVtolerancesB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("relTolB"), nb::arg("absTolB"));

m.def("IDAQuadReInitB", IDAQuadReInitB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("yQB0"));

m.def("IDAQuadSStolerancesB", IDAQuadSStolerancesB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def("IDAQuadSVtolerancesB", IDAQuadSVtolerancesB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def("IDACalcICB", IDACalcICB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("tout1"), nb::arg("yy0"), nb::arg("yp0"));

m.def(
  "IDACalcICBS",
  [](void* ida_mem, int which, double tout1, N_Vector yy0, N_Vector yp0,
     std::vector<N_Vector> yyS0, std::vector<N_Vector> ypS0) -> int
  {
    auto IDACalcICBS_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, int which, double tout1, N_Vector yy0, N_Vector yp0,
         std::vector<N_Vector> yyS0, std::vector<N_Vector> ypS0) -> int
    {
      N_Vector* yyS0_ptr =
        reinterpret_cast<N_Vector*>(yyS0.empty() ? nullptr : yyS0.data());
      N_Vector* ypS0_ptr =
        reinterpret_cast<N_Vector*>(ypS0.empty() ? nullptr : ypS0.data());

      auto lambda_result = IDACalcICBS(ida_mem, which, tout1, yy0, yp0,
                                       yyS0_ptr, ypS0_ptr);
      return lambda_result;
    };

    return IDACalcICBS_adapt_arr_ptr_to_std_vector(ida_mem, which, tout1, yy0,
                                                   yp0, yyS0, ypS0);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("tout1"), nb::arg("yy0"),
  nb::arg("yp0"), nb::arg("yyS0"), nb::arg("ypS0"));

m.def(
  "IDASolveF",
  [](void* ida_mem, double tout, double tret, N_Vector yret, N_Vector ypret,
     int itask, int ncheckPtr) -> std::tuple<int, double, int>
  {
    auto IDASolveF_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double tout, double tret, N_Vector yret, N_Vector ypret,
         int itask, int ncheckPtr) -> std::tuple<int, double, int>
    {
      double* tret_adapt_modifiable   = &tret;
      int* ncheckPtr_adapt_modifiable = &ncheckPtr;

      int r = IDASolveF(ida_mem, tout, tret_adapt_modifiable, yret, ypret,
                        itask, ncheckPtr_adapt_modifiable);
      return std::make_tuple(r, tret, ncheckPtr);
    };

    return IDASolveF_adapt_modifiable_immutable_to_return(ida_mem, tout, tret,
                                                          yret, ypret, itask,
                                                          ncheckPtr);
  },
  nb::arg("ida_mem"), nb::arg("tout"), nb::arg("tret"), nb::arg("yret"),
  nb::arg("ypret"), nb::arg("itask"), nb::arg("ncheckPtr"));

m.def("IDASolveB", IDASolveB, nb::arg("ida_mem"), nb::arg("tBout"),
      nb::arg("itaskB"));

m.def("IDAAdjSetNoSensi", IDAAdjSetNoSensi, nb::arg("ida_mem"));

m.def("IDASetMaxOrdB", IDASetMaxOrdB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("maxordB"));

m.def("IDASetMaxNumStepsB", IDASetMaxNumStepsB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("mxstepsB"));

m.def("IDASetInitStepB", IDASetInitStepB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("hinB"));

m.def("IDASetMaxStepB", IDASetMaxStepB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("hmaxB"));

m.def("IDASetSuppressAlgB", IDASetSuppressAlgB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("suppressalgB"));

m.def("IDASetIdB", IDASetIdB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("idB"));

m.def("IDASetConstraintsB", IDASetConstraintsB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("constraintsB"));

m.def("IDASetQuadErrConB", IDASetQuadErrConB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("errconQB"));

m.def("IDASetNonlinearSolverB", IDASetNonlinearSolverB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("NLS"));

m.def(
  "IDAGetB",
  [](void* ida_mem, int which, double tret, N_Vector yy,
     N_Vector yp) -> std::tuple<int, double>
  {
    auto IDAGetB_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, double tret, N_Vector yy,
         N_Vector yp) -> std::tuple<int, double>
    {
      double* tret_adapt_modifiable = &tret;

      int r = IDAGetB(ida_mem, which, tret_adapt_modifiable, yy, yp);
      return std::make_tuple(r, tret);
    };

    return IDAGetB_adapt_modifiable_immutable_to_return(ida_mem, which, tret,
                                                        yy, yp);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("tret"), nb::arg("yy"),
  nb::arg("yp"));

m.def(
  "IDAGetQuadB",
  [](void* ida_mem, int which, double tret, N_Vector qB) -> std::tuple<int, double>
  {
    auto IDAGetQuadB_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, double tret,
         N_Vector qB) -> std::tuple<int, double>
    {
      double* tret_adapt_modifiable = &tret;

      int r = IDAGetQuadB(ida_mem, which, tret_adapt_modifiable, qB);
      return std::make_tuple(r, tret);
    };

    return IDAGetQuadB_adapt_modifiable_immutable_to_return(ida_mem, which,
                                                            tret, qB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("tret"), nb::arg("qB"));

m.def("IDAGetAdjIDABmem", IDAGetAdjIDABmem, nb::arg("ida_mem"), nb::arg("which"));

m.def("IDAGetConsistentICB", IDAGetConsistentICB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("yyB0"), nb::arg("ypB0"));

m.def("IDAGetAdjY", IDAGetAdjY, nb::arg("ida_mem"), nb::arg("t"), nb::arg("yy"),
      nb::arg("yp"));

m.def("IDAGetAdjCheckPointsInfo", IDAGetAdjCheckPointsInfo, nb::arg("ida_mem"),
      nb::arg("ckpnt"));

m.def(
  "IDAGetAdjDataPointHermite",
  [](void* ida_mem, int which, double t, N_Vector yy,
     N_Vector yd) -> std::tuple<int, double>
  {
    auto IDAGetAdjDataPointHermite_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, double t, N_Vector yy,
         N_Vector yd) -> std::tuple<int, double>
    {
      double* t_adapt_modifiable = &t;

      int r = IDAGetAdjDataPointHermite(ida_mem, which, t_adapt_modifiable, yy,
                                        yd);
      return std::make_tuple(r, t);
    };

    return IDAGetAdjDataPointHermite_adapt_modifiable_immutable_to_return(ida_mem,
                                                                          which,
                                                                          t, yy,
                                                                          yd);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("t"), nb::arg("yy"),
  nb::arg("yd"));

m.def(
  "IDAGetAdjDataPointPolynomial",
  [](void* ida_mem, int which, double t, int order,
     N_Vector y) -> std::tuple<int, double, int>
  {
    auto IDAGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, double t, int order,
         N_Vector y) -> std::tuple<int, double, int>
    {
      double* t_adapt_modifiable  = &t;
      int* order_adapt_modifiable = &order;

      int r = IDAGetAdjDataPointPolynomial(ida_mem, which, t_adapt_modifiable,
                                           order_adapt_modifiable, y);
      return std::make_tuple(r, t, order);
    };

    return IDAGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return(ida_mem,
                                                                             which,
                                                                             t,
                                                                             order,
                                                                             y);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("t"), nb::arg("order"),
  nb::arg("y"));
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

m.def(
  "IDASetLinearSolver",
  [](void* ida_mem, SUNLinearSolver LS,
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

    return IDASetLinearSolver_adapt_optional_arg_with_default_null(ida_mem, LS,
                                                                   A);
  },
  nb::arg("ida_mem"), nb::arg("LS"), nb::arg("A") = nb::none());

m.def("IDASetEpsLin", IDASetEpsLin, nb::arg("ida_mem"), nb::arg("eplifac"));

m.def("IDASetLSNormFactor", IDASetLSNormFactor, nb::arg("ida_mem"),
      nb::arg("nrmfac"));

m.def("IDASetLinearSolutionScaling", IDASetLinearSolutionScaling,
      nb::arg("ida_mem"), nb::arg("onoff"));

m.def("IDASetIncrementFactor", IDASetIncrementFactor, nb::arg("ida_mem"),
      nb::arg("dqincfac"));

m.def(
  "IDAGetJacCj",
  [](void* ida_mem, double cj_J) -> std::tuple<int, double>
  {
    auto IDAGetJacCj_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double cj_J) -> std::tuple<int, double>
    {
      double* cj_J_adapt_modifiable = &cj_J;

      int r = IDAGetJacCj(ida_mem, cj_J_adapt_modifiable);
      return std::make_tuple(r, cj_J);
    };

    return IDAGetJacCj_adapt_modifiable_immutable_to_return(ida_mem, cj_J);
  },
  nb::arg("ida_mem"), nb::arg("cj_J"));

m.def(
  "IDAGetJacTime",
  [](void* ida_mem, double t_J) -> std::tuple<int, double>
  {
    auto IDAGetJacTime_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, double t_J) -> std::tuple<int, double>
    {
      double* t_J_adapt_modifiable = &t_J;

      int r = IDAGetJacTime(ida_mem, t_J_adapt_modifiable);
      return std::make_tuple(r, t_J);
    };

    return IDAGetJacTime_adapt_modifiable_immutable_to_return(ida_mem, t_J);
  },
  nb::arg("ida_mem"), nb::arg("t_J"));

m.def(
  "IDAGetJacNumSteps",
  [](void* ida_mem, long nst_J) -> std::tuple<int, long>
  {
    auto IDAGetJacNumSteps_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nst_J) -> std::tuple<int, long>
    {
      long* nst_J_adapt_modifiable = &nst_J;

      int r = IDAGetJacNumSteps(ida_mem, nst_J_adapt_modifiable);
      return std::make_tuple(r, nst_J);
    };

    return IDAGetJacNumSteps_adapt_modifiable_immutable_to_return(ida_mem, nst_J);
  },
  nb::arg("ida_mem"), nb::arg("nst_J"));

m.def(
  "IDAGetNumJacEvals",
  [](void* ida_mem, long njevals) -> std::tuple<int, long>
  {
    auto IDAGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long njevals) -> std::tuple<int, long>
    {
      long* njevals_adapt_modifiable = &njevals;

      int r = IDAGetNumJacEvals(ida_mem, njevals_adapt_modifiable);
      return std::make_tuple(r, njevals);
    };

    return IDAGetNumJacEvals_adapt_modifiable_immutable_to_return(ida_mem,
                                                                  njevals);
  },
  nb::arg("ida_mem"), nb::arg("njevals"));

m.def(
  "IDAGetNumPrecEvals",
  [](void* ida_mem, long npevals) -> std::tuple<int, long>
  {
    auto IDAGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long npevals) -> std::tuple<int, long>
    {
      long* npevals_adapt_modifiable = &npevals;

      int r = IDAGetNumPrecEvals(ida_mem, npevals_adapt_modifiable);
      return std::make_tuple(r, npevals);
    };

    return IDAGetNumPrecEvals_adapt_modifiable_immutable_to_return(ida_mem,
                                                                   npevals);
  },
  nb::arg("ida_mem"), nb::arg("npevals"));

m.def(
  "IDAGetNumPrecSolves",
  [](void* ida_mem, long npsolves) -> std::tuple<int, long>
  {
    auto IDAGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long npsolves) -> std::tuple<int, long>
    {
      long* npsolves_adapt_modifiable = &npsolves;

      int r = IDAGetNumPrecSolves(ida_mem, npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves);
    };

    return IDAGetNumPrecSolves_adapt_modifiable_immutable_to_return(ida_mem,
                                                                    npsolves);
  },
  nb::arg("ida_mem"), nb::arg("npsolves"));

m.def(
  "IDAGetNumLinIters",
  [](void* ida_mem, long nliters) -> std::tuple<int, long>
  {
    auto IDAGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nliters) -> std::tuple<int, long>
    {
      long* nliters_adapt_modifiable = &nliters;

      int r = IDAGetNumLinIters(ida_mem, nliters_adapt_modifiable);
      return std::make_tuple(r, nliters);
    };

    return IDAGetNumLinIters_adapt_modifiable_immutable_to_return(ida_mem,
                                                                  nliters);
  },
  nb::arg("ida_mem"), nb::arg("nliters"));

m.def(
  "IDAGetNumLinConvFails",
  [](void* ida_mem, long nlcfails) -> std::tuple<int, long>
  {
    auto IDAGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nlcfails) -> std::tuple<int, long>
    {
      long* nlcfails_adapt_modifiable = &nlcfails;

      int r = IDAGetNumLinConvFails(ida_mem, nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails);
    };

    return IDAGetNumLinConvFails_adapt_modifiable_immutable_to_return(ida_mem,
                                                                      nlcfails);
  },
  nb::arg("ida_mem"), nb::arg("nlcfails"));

m.def(
  "IDAGetNumJTSetupEvals",
  [](void* ida_mem, long njtsetups) -> std::tuple<int, long>
  {
    auto IDAGetNumJTSetupEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long njtsetups) -> std::tuple<int, long>
    {
      long* njtsetups_adapt_modifiable = &njtsetups;

      int r = IDAGetNumJTSetupEvals(ida_mem, njtsetups_adapt_modifiable);
      return std::make_tuple(r, njtsetups);
    };

    return IDAGetNumJTSetupEvals_adapt_modifiable_immutable_to_return(ida_mem,
                                                                      njtsetups);
  },
  nb::arg("ida_mem"), nb::arg("njtsetups"));

m.def(
  "IDAGetNumJtimesEvals",
  [](void* ida_mem, long njvevals) -> std::tuple<int, long>
  {
    auto IDAGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long njvevals) -> std::tuple<int, long>
    {
      long* njvevals_adapt_modifiable = &njvevals;

      int r = IDAGetNumJtimesEvals(ida_mem, njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals);
    };

    return IDAGetNumJtimesEvals_adapt_modifiable_immutable_to_return(ida_mem,
                                                                     njvevals);
  },
  nb::arg("ida_mem"), nb::arg("njvevals"));

m.def(
  "IDAGetNumLinResEvals",
  [](void* ida_mem, long nrevalsLS) -> std::tuple<int, long>
  {
    auto IDAGetNumLinResEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long nrevalsLS) -> std::tuple<int, long>
    {
      long* nrevalsLS_adapt_modifiable = &nrevalsLS;

      int r = IDAGetNumLinResEvals(ida_mem, nrevalsLS_adapt_modifiable);
      return std::make_tuple(r, nrevalsLS);
    };

    return IDAGetNumLinResEvals_adapt_modifiable_immutable_to_return(ida_mem,
                                                                     nrevalsLS);
  },
  nb::arg("ida_mem"), nb::arg("nrevalsLS"));

m.def(
  "IDAGetLastLinFlag",
  [](void* ida_mem, long flag) -> std::tuple<int, long>
  {
    auto IDAGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, long flag) -> std::tuple<int, long>
    {
      long* flag_adapt_modifiable = &flag;

      int r = IDAGetLastLinFlag(ida_mem, flag_adapt_modifiable);
      return std::make_tuple(r, flag);
    };

    return IDAGetLastLinFlag_adapt_modifiable_immutable_to_return(ida_mem, flag);
  },
  nb::arg("ida_mem"), nb::arg("flag"));

m.def("IDAGetLinReturnFlagName", IDAGetLinReturnFlagName, nb::arg("flag"));

m.def(
  "IDASetLinearSolverB",
  [](void* ida_mem, int which, SUNLinearSolver LS,
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

    return IDASetLinearSolverB_adapt_optional_arg_with_default_null(ida_mem,
                                                                    which, LS, A);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("LS"), nb::arg("A") = nb::none());

m.def("IDASetEpsLinB", IDASetEpsLinB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("eplifacB"));

m.def("IDASetLSNormFactorB", IDASetLSNormFactorB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("nrmfacB"));

m.def("IDASetLinearSolutionScalingB", IDASetLinearSolutionScalingB,
      nb::arg("ida_mem"), nb::arg("which"), nb::arg("onoffB"));

m.def("IDASetIncrementFactorB", IDASetIncrementFactorB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("dqincfacB"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
