// #ifndef _CVODES_H
// 
// #ifdef __cplusplus 
// #endif
// 
m.attr("CV_ADAMS") = 1;
m.attr("CV_BDF") = 2;
m.attr("CV_NORMAL") = 1;
m.attr("CV_ONE_STEP") = 2;
m.attr("CV_SIMULTANEOUS") = 1;
m.attr("CV_STAGGERED") = 2;
m.attr("CV_STAGGERED1") = 3;
m.attr("CV_CENTERED") = 1;
m.attr("CV_FORWARD") = 2;
m.attr("CV_HERMITE") = 1;
m.attr("CV_POLYNOMIAL") = 2;
m.attr("CV_SUCCESS") = 0;
m.attr("CV_TSTOP_RETURN") = 1;
m.attr("CV_ROOT_RETURN") = 2;
m.attr("CV_WARNING") = 99;
m.attr("CV_TOO_MUCH_WORK") = -1;
m.attr("CV_TOO_MUCH_ACC") = -2;
m.attr("CV_ERR_FAILURE") = -3;
m.attr("CV_CONV_FAILURE") = -4;
m.attr("CV_LINIT_FAIL") = -5;
m.attr("CV_LSETUP_FAIL") = -6;
m.attr("CV_LSOLVE_FAIL") = -7;
m.attr("CV_RHSFUNC_FAIL") = -8;
m.attr("CV_FIRST_RHSFUNC_ERR") = -9;
m.attr("CV_REPTD_RHSFUNC_ERR") = -10;
m.attr("CV_UNREC_RHSFUNC_ERR") = -11;
m.attr("CV_RTFUNC_FAIL") = -12;
m.attr("CV_NLS_INIT_FAIL") = -13;
m.attr("CV_NLS_SETUP_FAIL") = -14;
m.attr("CV_CONSTR_FAIL") = -15;
m.attr("CV_NLS_FAIL") = -16;
m.attr("CV_MEM_FAIL") = -20;
m.attr("CV_MEM_NULL") = -21;
m.attr("CV_ILL_INPUT") = -22;
m.attr("CV_NO_MALLOC") = -23;
m.attr("CV_BAD_K") = -24;
m.attr("CV_BAD_T") = -25;
m.attr("CV_BAD_DKY") = -26;
m.attr("CV_TOO_CLOSE") = -27;
m.attr("CV_VECTOROP_ERR") = -28;
m.attr("CV_NO_QUAD") = -30;
m.attr("CV_QRHSFUNC_FAIL") = -31;
m.attr("CV_FIRST_QRHSFUNC_ERR") = -32;
m.attr("CV_REPTD_QRHSFUNC_ERR") = -33;
m.attr("CV_UNREC_QRHSFUNC_ERR") = -34;
m.attr("CV_NO_SENS") = -40;
m.attr("CV_SRHSFUNC_FAIL") = -41;
m.attr("CV_FIRST_SRHSFUNC_ERR") = -42;
m.attr("CV_REPTD_SRHSFUNC_ERR") = -43;
m.attr("CV_UNREC_SRHSFUNC_ERR") = -44;
m.attr("CV_BAD_IS") = -45;
m.attr("CV_NO_QUADSENS") = -50;
m.attr("CV_QSRHSFUNC_FAIL") = -51;
m.attr("CV_FIRST_QSRHSFUNC_ERR") = -52;
m.attr("CV_REPTD_QSRHSFUNC_ERR") = -53;
m.attr("CV_UNREC_QSRHSFUNC_ERR") = -54;
m.attr("CV_CONTEXT_ERR") = -55;
m.attr("CV_PROJ_MEM_NULL") = -56;
m.attr("CV_PROJFUNC_FAIL") = -57;
m.attr("CV_REPTD_PROJFUNC_ERR") = -58;
m.attr("CV_BAD_TINTERP") = -59;
m.attr("CV_UNRECOGNIZED_ERR") = -99;
m.attr("CV_NO_ADJ") = -101;
m.attr("CV_NO_FWD") = -102;
m.attr("CV_NO_BCK") = -103;
m.attr("CV_BAD_TB0") = -104;
m.attr("CV_REIFWD_FAIL") = -105;
m.attr("CV_FWD_FAIL") = -106;
m.attr("CV_GETY_BADT") = -107;

m.def("CVodeInit",
    CVodeInit, nb::arg("cvode_mem"), nb::arg("f"), nb::arg("t0"), nb::arg("y0"));

m.def("CVodeReInit",
    CVodeReInit, nb::arg("cvode_mem"), nb::arg("t0"), nb::arg("y0"));

m.def("CVodeSStolerances",
    CVodeSStolerances, nb::arg("cvode_mem"), nb::arg("reltol"), nb::arg("abstol"));

m.def("CVodeSVtolerances",
    CVodeSVtolerances, nb::arg("cvode_mem"), nb::arg("reltol"), nb::arg("abstol"));

m.def("CVodeWFtolerances",
    CVodeWFtolerances, nb::arg("cvode_mem"), nb::arg("efun"));

m.def("CVodeSetConstraints",
    CVodeSetConstraints, nb::arg("cvode_mem"), nb::arg("constraints"));

m.def("CVodeSetDeltaGammaMaxLSetup",
    CVodeSetDeltaGammaMaxLSetup, nb::arg("cvode_mem"), nb::arg("dgmax_lsetup"));

m.def("CVodeSetInitStep",
    CVodeSetInitStep, nb::arg("cvode_mem"), nb::arg("hin"));

m.def("CVodeSetLSetupFrequency",
    CVodeSetLSetupFrequency, nb::arg("cvode_mem"), nb::arg("msbp"));

m.def("CVodeSetMaxConvFails",
    CVodeSetMaxConvFails, nb::arg("cvode_mem"), nb::arg("maxncf"));

m.def("CVodeSetMaxErrTestFails",
    CVodeSetMaxErrTestFails, nb::arg("cvode_mem"), nb::arg("maxnef"));

m.def("CVodeSetMaxHnilWarns",
    CVodeSetMaxHnilWarns, nb::arg("cvode_mem"), nb::arg("mxhnil"));

m.def("CVodeSetMaxNonlinIters",
    CVodeSetMaxNonlinIters, nb::arg("cvode_mem"), nb::arg("maxcor"));

m.def("CVodeSetMaxNumSteps",
    CVodeSetMaxNumSteps, nb::arg("cvode_mem"), nb::arg("mxsteps"));

m.def("CVodeSetMaxOrd",
    CVodeSetMaxOrd, nb::arg("cvode_mem"), nb::arg("maxord"));

m.def("CVodeSetMaxStep",
    CVodeSetMaxStep, nb::arg("cvode_mem"), nb::arg("hmax"));

m.def("CVodeSetMinStep",
    CVodeSetMinStep, nb::arg("cvode_mem"), nb::arg("hmin"));

m.def("CVodeSetMonitorFrequency",
    CVodeSetMonitorFrequency, nb::arg("cvode_mem"), nb::arg("nst"));

m.def("CVodeSetNonlinConvCoef",
    CVodeSetNonlinConvCoef, nb::arg("cvode_mem"), nb::arg("nlscoef"));

m.def("CVodeSetNonlinearSolver",
    CVodeSetNonlinearSolver, nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSetStabLimDet",
    CVodeSetStabLimDet, nb::arg("cvode_mem"), nb::arg("stldet"));

m.def("CVodeSetStopTime",
    CVodeSetStopTime, nb::arg("cvode_mem"), nb::arg("tstop"));

m.def("CVodeSetInterpolateStopTime",
    CVodeSetInterpolateStopTime, nb::arg("cvode_mem"), nb::arg("interp"));

m.def("CVodeClearStopTime",
    CVodeClearStopTime, nb::arg("cvode_mem"));

m.def("CVodeSetOwnUserData",
    CVodeSetOwnUserData, nb::arg("cvode_mem"), nb::arg("own_user_data"));

m.def("CVodeSetEtaFixedStepBounds",
    CVodeSetEtaFixedStepBounds, nb::arg("cvode_mem"), nb::arg("eta_min_fx"), nb::arg("eta_max_fx"));

m.def("CVodeSetEtaMaxFirstStep",
    CVodeSetEtaMaxFirstStep, nb::arg("cvode_mem"), nb::arg("eta_max_fs"));

m.def("CVodeSetEtaMaxEarlyStep",
    CVodeSetEtaMaxEarlyStep, nb::arg("cvode_mem"), nb::arg("eta_max_es"));

m.def("CVodeSetNumStepsEtaMaxEarlyStep",
    CVodeSetNumStepsEtaMaxEarlyStep, nb::arg("cvode_mem"), nb::arg("small_nst"));

m.def("CVodeSetEtaMax",
    CVodeSetEtaMax, nb::arg("cvode_mem"), nb::arg("eta_max_gs"));

m.def("CVodeSetEtaMin",
    CVodeSetEtaMin, nb::arg("cvode_mem"), nb::arg("eta_min"));

m.def("CVodeSetEtaMinErrFail",
    CVodeSetEtaMinErrFail, nb::arg("cvode_mem"), nb::arg("eta_min_ef"));

m.def("CVodeSetEtaMaxErrFail",
    CVodeSetEtaMaxErrFail, nb::arg("cvode_mem"), nb::arg("eta_max_ef"));

m.def("CVodeSetNumFailsEtaMaxErrFail",
    CVodeSetNumFailsEtaMaxErrFail, nb::arg("cvode_mem"), nb::arg("small_nef"));

m.def("CVodeSetEtaConvFail",
    CVodeSetEtaConvFail, nb::arg("cvode_mem"), nb::arg("eta_cf"));

m.def("CVodeRootInit",
    CVodeRootInit, nb::arg("cvode_mem"), nb::arg("nrtfn"), nb::arg("g"));

m.def("CVodeSetRootDirection",
    CVodeSetRootDirection, nb::arg("cvode_mem"), nb::arg("rootdir"));

m.def("CVodeSetNoInactiveRootWarn",
    CVodeSetNoInactiveRootWarn, nb::arg("cvode_mem"));

m.def("CVode",
    CVode, nb::arg("cvode_mem"), nb::arg("tout"), nb::arg("yout"), nb::arg("tret"), nb::arg("itask"));

m.def("CVodeComputeState",
    CVodeComputeState, nb::arg("cvode_mem"), nb::arg("ycor"), nb::arg("y"));

m.def("CVodeGetDky",
    CVodeGetDky, nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dky"));

m.def("CVodeGetNumSteps",
    CVodeGetNumSteps, nb::arg("cvode_mem"), nb::arg("nsteps"));

m.def("CVodeGetNumRhsEvals",
    CVodeGetNumRhsEvals, nb::arg("cvode_mem"), nb::arg("nfevals"));

m.def("CVodeGetNumLinSolvSetups",
    CVodeGetNumLinSolvSetups, nb::arg("cvode_mem"), nb::arg("nlinsetups"));

m.def("CVodeGetNumErrTestFails",
    CVodeGetNumErrTestFails, nb::arg("cvode_mem"), nb::arg("netfails"));

m.def("CVodeGetLastOrder",
    CVodeGetLastOrder, nb::arg("cvode_mem"), nb::arg("qlast"));

m.def("CVodeGetCurrentOrder",
    CVodeGetCurrentOrder, nb::arg("cvode_mem"), nb::arg("qcur"));

m.def("CVodeGetCurrentGamma",
    CVodeGetCurrentGamma, nb::arg("cvode_mem"), nb::arg("gamma"));

m.def("CVodeGetNumStabLimOrderReds",
    CVodeGetNumStabLimOrderReds, nb::arg("cvode_mem"), nb::arg("nslred"));

m.def("CVodeGetActualInitStep",
    CVodeGetActualInitStep, nb::arg("cvode_mem"), nb::arg("hinused"));

m.def("CVodeGetLastStep",
    CVodeGetLastStep, nb::arg("cvode_mem"), nb::arg("hlast"));

m.def("CVodeGetCurrentStep",
    CVodeGetCurrentStep, nb::arg("cvode_mem"), nb::arg("hcur"));

m.def("CVodeGetCurrentSensSolveIndex",
    CVodeGetCurrentSensSolveIndex, nb::arg("cvode_mem"), nb::arg("index"));

m.def("CVodeGetCurrentTime",
    CVodeGetCurrentTime, nb::arg("cvode_mem"), nb::arg("tcur"));

m.def("CVodeGetTolScaleFactor",
    CVodeGetTolScaleFactor, nb::arg("cvode_mem"), nb::arg("tolsfac"));

m.def("CVodeGetErrWeights",
    CVodeGetErrWeights, nb::arg("cvode_mem"), nb::arg("eweight"));

m.def("CVodeGetEstLocalErrors",
    CVodeGetEstLocalErrors, nb::arg("cvode_mem"), nb::arg("ele"));

m.def("CVodeGetNumGEvals",
    CVodeGetNumGEvals, nb::arg("cvode_mem"), nb::arg("ngevals"));

m.def("CVodeGetRootInfo",
    CVodeGetRootInfo, nb::arg("cvode_mem"), nb::arg("rootsfound"));

m.def("CVodeGetIntegratorStats",
    CVodeGetIntegratorStats, nb::arg("cvode_mem"), nb::arg("nsteps"), nb::arg("nfevals"), nb::arg("nlinsetups"), nb::arg("netfails"), nb::arg("qlast"), nb::arg("qcur"), nb::arg("hinused"), nb::arg("hlast"), nb::arg("hcur"), nb::arg("tcur"));

m.def("CVodeGetNumNonlinSolvIters",
    CVodeGetNumNonlinSolvIters, nb::arg("cvode_mem"), nb::arg("nniters"));

m.def("CVodeGetNumNonlinSolvConvFails",
    CVodeGetNumNonlinSolvConvFails, nb::arg("cvode_mem"), nb::arg("nnfails"));

m.def("CVodeGetNonlinSolvStats",
    CVodeGetNonlinSolvStats, nb::arg("cvode_mem"), nb::arg("nniters"), nb::arg("nnfails"));

m.def("CVodeGetNumStepSolveFails",
    CVodeGetNumStepSolveFails, nb::arg("cvode_mem"), nb::arg("nncfails"));

m.def("CVodePrintAllStats",
    CVodePrintAllStats, nb::arg("cvode_mem"), nb::arg("outfile"), nb::arg("fmt"));

m.def("CVodeGetReturnFlagName",
    CVodeGetReturnFlagName, nb::arg("flag"));

m.def("CVodeQuadInit",
    CVodeQuadInit, nb::arg("cvode_mem"), nb::arg("fQ"), nb::arg("yQ0"));

m.def("CVodeQuadReInit",
    CVodeQuadReInit, nb::arg("cvode_mem"), nb::arg("yQ0"));

m.def("CVodeQuadSStolerances",
    CVodeQuadSStolerances, nb::arg("cvode_mem"), nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("CVodeQuadSVtolerances",
    CVodeQuadSVtolerances, nb::arg("cvode_mem"), nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("CVodeSetQuadErrCon",
    CVodeSetQuadErrCon, nb::arg("cvode_mem"), nb::arg("errconQ"));

m.def("CVodeGetQuad",
    CVodeGetQuad, nb::arg("cvode_mem"), nb::arg("tret"), nb::arg("yQout"));

m.def("CVodeGetQuadDky",
    CVodeGetQuadDky, nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dky"));

m.def("CVodeGetQuadNumRhsEvals",
    CVodeGetQuadNumRhsEvals, nb::arg("cvode_mem"), nb::arg("nfQevals"));

m.def("CVodeGetQuadNumErrTestFails",
    CVodeGetQuadNumErrTestFails, nb::arg("cvode_mem"), nb::arg("nQetfails"));

m.def("CVodeGetQuadErrWeights",
    CVodeGetQuadErrWeights, nb::arg("cvode_mem"), nb::arg("eQweight"));

m.def("CVodeGetQuadStats",
    CVodeGetQuadStats, nb::arg("cvode_mem"), nb::arg("nfQevals"), nb::arg("nQetfails"));

m.def("CVodeSensSStolerances",
    CVodeSensSStolerances, nb::arg("cvode_mem"), nb::arg("reltolS"), nb::arg("abstolS"));

m.def("CVodeSensEEtolerances",
    CVodeSensEEtolerances, nb::arg("cvode_mem"));

m.def("CVodeSetSensDQMethod",
    CVodeSetSensDQMethod, nb::arg("cvode_mem"), nb::arg("DQtype"), nb::arg("DQrhomax"));

m.def("CVodeSetSensErrCon",
    CVodeSetSensErrCon, nb::arg("cvode_mem"), nb::arg("errconS"));

m.def("CVodeSetSensMaxNonlinIters",
    CVodeSetSensMaxNonlinIters, nb::arg("cvode_mem"), nb::arg("maxcorS"));

m.def("CVodeSetSensParams",
    CVodeSetSensParams, nb::arg("cvode_mem"), nb::arg("p"), nb::arg("pbar"), nb::arg("plist"));

m.def("CVodeSetNonlinearSolverSensSim",
    CVodeSetNonlinearSolverSensSim, nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSetNonlinearSolverSensStg",
    CVodeSetNonlinearSolverSensStg, nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSetNonlinearSolverSensStg1",
    CVodeSetNonlinearSolverSensStg1, nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSensToggleOff",
    CVodeSensToggleOff, nb::arg("cvode_mem"));

m.def("CVodeGetNumRhsEvalsSens",
    CVodeGetNumRhsEvalsSens, nb::arg("cvode_mem"), nb::arg("nfevalsS"));

m.def("CVodeGetNumStepSensSolveFails",
    CVodeGetNumStepSensSolveFails, nb::arg("cvode_mem"), nb::arg("nSncfails"));

m.def("CVodeGetStgrSensNumNonlinSolvIters",
    CVodeGetStgrSensNumNonlinSolvIters, nb::arg("cvode_mem"), nb::arg("nSTGR1niters"));

m.def("CVodeGetStgrSensNumNonlinSolvConvFails",
    CVodeGetStgrSensNumNonlinSolvConvFails, nb::arg("cvode_mem"), nb::arg("nSTGR1nfails"));

m.def("CVodeGetStgrSensNonlinSolvStats",
    CVodeGetStgrSensNonlinSolvStats, nb::arg("cvode_mem"), nb::arg("nSTGR1niters"), nb::arg("nSTGR1nfails"));

m.def("CVodeGetNumStepStgrSensSolveFails",
    CVodeGetNumStepStgrSensSolveFails, nb::arg("cvode_mem"), nb::arg("nSTGR1ncfails"));

m.def("CVodeQuadSensSStolerances",
    CVodeQuadSensSStolerances, nb::arg("cvode_mem"), nb::arg("reltolQS"), nb::arg("abstolQS"));

m.def("CVodeQuadSensEEtolerances",
    CVodeQuadSensEEtolerances, nb::arg("cvode_mem"));

m.def("CVodeSetQuadSensErrCon",
    CVodeSetQuadSensErrCon, nb::arg("cvode_mem"), nb::arg("errconQS"));

m.def("CVodeAdjInit",
    CVodeAdjInit, nb::arg("cvode_mem"), nb::arg("steps"), nb::arg("interp"));

m.def("CVodeAdjReInit",
    CVodeAdjReInit, nb::arg("cvode_mem"));

m.def("CVodeInitB",
    CVodeInitB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("fB"), nb::arg("tB0"), nb::arg("yB0"));

m.def("CVodeInitBS",
    CVodeInitBS, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("fBs"), nb::arg("tB0"), nb::arg("yB0"));

m.def("CVodeReInitB",
    CVodeReInitB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("tB0"), nb::arg("yB0"));

m.def("CVodeSStolerancesB",
    CVodeSStolerancesB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("reltolB"), nb::arg("abstolB"));

m.def("CVodeSVtolerancesB",
    CVodeSVtolerancesB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("reltolB"), nb::arg("abstolB"));

m.def("CVodeQuadInitB",
    CVodeQuadInitB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("fQB"), nb::arg("yQB0"));

m.def("CVodeQuadInitBS",
    CVodeQuadInitBS, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("fQBs"), nb::arg("yQB0"));

m.def("CVodeQuadReInitB",
    CVodeQuadReInitB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("yQB0"));

m.def("CVodeQuadSStolerancesB",
    CVodeQuadSStolerancesB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def("CVodeQuadSVtolerancesB",
    CVodeQuadSVtolerancesB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def("CVodeF",
    CVodeF, nb::arg("cvode_mem"), nb::arg("tout"), nb::arg("yout"), nb::arg("tret"), nb::arg("itask"), nb::arg("ncheckPtr"));

m.def("CVodeB",
    CVodeB, nb::arg("cvode_mem"), nb::arg("tBout"), nb::arg("itaskB"));

m.def("CVodeSetAdjNoSensi",
    CVodeSetAdjNoSensi, nb::arg("cvode_mem"));

m.def("CVodeSetMaxOrdB",
    CVodeSetMaxOrdB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("maxordB"));

m.def("CVodeSetMaxNumStepsB",
    CVodeSetMaxNumStepsB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("mxstepsB"));

m.def("CVodeSetStabLimDetB",
    CVodeSetStabLimDetB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("stldetB"));

m.def("CVodeSetInitStepB",
    CVodeSetInitStepB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("hinB"));

m.def("CVodeSetMinStepB",
    CVodeSetMinStepB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("hminB"));

m.def("CVodeSetMaxStepB",
    CVodeSetMaxStepB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("hmaxB"));

m.def("CVodeSetConstraintsB",
    CVodeSetConstraintsB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("constraintsB"));

m.def("CVodeSetQuadErrConB",
    CVodeSetQuadErrConB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("errconQB"));

m.def("CVodeSetNonlinearSolverB",
    CVodeSetNonlinearSolverB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("NLS"));

m.def("CVodeGetB",
    CVodeGetB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("tBret"), nb::arg("yB"));

m.def("CVodeGetQuadB",
    CVodeGetQuadB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("tBret"), nb::arg("qB"));

m.def("CVodeGetAdjCVodeBmem",
    CVodeGetAdjCVodeBmem, nb::arg("cvode_mem"), nb::arg("which"));

m.def("CVodeGetAdjY",
    CVodeGetAdjY, nb::arg("cvode_mem"), nb::arg("t"), nb::arg("y"));

m.def("CVodeGetAdjCheckPointsInfo",
    CVodeGetAdjCheckPointsInfo, nb::arg("cvode_mem"), nb::arg("ckpnt"));

m.def("CVodeGetAdjDataPointHermite",
    CVodeGetAdjDataPointHermite, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("t"), nb::arg("y"), nb::arg("yd"));

m.def("CVodeGetAdjDataPointPolynomial",
    CVodeGetAdjDataPointPolynomial, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("t"), nb::arg("order"), nb::arg("y"));
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

m.def("CVodeSetLinearSolver",
    CVodeSetLinearSolver, nb::arg("cvode_mem"), nb::arg("LS"), nb::arg("A"));

m.def("CVodeSetJacEvalFrequency",
    CVodeSetJacEvalFrequency, nb::arg("cvode_mem"), nb::arg("msbj"));

m.def("CVodeSetLinearSolutionScaling",
    CVodeSetLinearSolutionScaling, nb::arg("cvode_mem"), nb::arg("onoff"));

m.def("CVodeSetDeltaGammaMaxBadJac",
    CVodeSetDeltaGammaMaxBadJac, nb::arg("cvode_mem"), nb::arg("dgmax_jbad"));

m.def("CVodeSetEpsLin",
    CVodeSetEpsLin, nb::arg("cvode_mem"), nb::arg("eplifac"));

m.def("CVodeSetLSNormFactor",
    CVodeSetLSNormFactor, nb::arg("arkode_mem"), nb::arg("nrmfac"));

m.def("CVodeGetNumJacEvals",
    CVodeGetNumJacEvals, nb::arg("cvode_mem"), nb::arg("njevals"));

m.def("CVodeGetNumPrecEvals",
    CVodeGetNumPrecEvals, nb::arg("cvode_mem"), nb::arg("npevals"));

m.def("CVodeGetNumPrecSolves",
    CVodeGetNumPrecSolves, nb::arg("cvode_mem"), nb::arg("npsolves"));

m.def("CVodeGetNumLinIters",
    CVodeGetNumLinIters, nb::arg("cvode_mem"), nb::arg("nliters"));

m.def("CVodeGetNumLinConvFails",
    CVodeGetNumLinConvFails, nb::arg("cvode_mem"), nb::arg("nlcfails"));

m.def("CVodeGetNumJTSetupEvals",
    CVodeGetNumJTSetupEvals, nb::arg("cvode_mem"), nb::arg("njtsetups"));

m.def("CVodeGetNumJtimesEvals",
    CVodeGetNumJtimesEvals, nb::arg("cvode_mem"), nb::arg("njvevals"));

m.def("CVodeGetNumLinRhsEvals",
    CVodeGetNumLinRhsEvals, nb::arg("cvode_mem"), nb::arg("nfevalsLS"));

m.def("CVodeGetLinSolveStats",
    CVodeGetLinSolveStats, nb::arg("cvode_mem"), nb::arg("njevals"), nb::arg("nfevalsLS"), nb::arg("nliters"), nb::arg("nlcfails"), nb::arg("npevals"), nb::arg("npsolves"), nb::arg("njtsetups"), nb::arg("njtimes"));

m.def("CVodeGetLastLinFlag",
    CVodeGetLastLinFlag, nb::arg("cvode_mem"), nb::arg("flag"));

m.def("CVodeGetLinReturnFlagName",
    CVodeGetLinReturnFlagName, nb::arg("flag"));

m.def("CVodeSetLinearSolverB",
    CVodeSetLinearSolverB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("LS"), nb::arg("A"));

m.def("CVodeSetEpsLinB",
    CVodeSetEpsLinB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("eplifacB"));

m.def("CVodeSetLSNormFactorB",
    CVodeSetLSNormFactorB, nb::arg("arkode_mem"), nb::arg("which"), nb::arg("nrmfacB"));

m.def("CVodeSetLinearSolutionScalingB",
    CVodeSetLinearSolutionScalingB, nb::arg("cvode_mem"), nb::arg("which"), nb::arg("onoffB"));
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

m.def("CVodeSetProjErrEst",
    CVodeSetProjErrEst, nb::arg("cvode_mem"), nb::arg("onoff"));

m.def("CVodeSetProjFrequency",
    CVodeSetProjFrequency, nb::arg("cvode_mem"), nb::arg("proj_freq"));

m.def("CVodeSetMaxNumProjFails",
    CVodeSetMaxNumProjFails, nb::arg("cvode_mem"), nb::arg("max_fails"));

m.def("CVodeSetEpsProj",
    CVodeSetEpsProj, nb::arg("cvode_mem"), nb::arg("eps"));

m.def("CVodeSetProjFailEta",
    CVodeSetProjFailEta, nb::arg("cvode_mem"), nb::arg("eta"));

m.def("CVodeGetNumProjEvals",
    CVodeGetNumProjEvals, nb::arg("cvode_mem"), nb::arg("nproj"));

m.def("CVodeGetNumProjFails",
    CVodeGetNumProjFails, nb::arg("cvode_mem"), nb::arg("nprf"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
