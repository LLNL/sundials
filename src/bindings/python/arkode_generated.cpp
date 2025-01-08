// #ifndef _ARKODE_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("ARKodeResize",
    ARKodeResize, nb::arg("arkode_mem"), nb::arg("ynew"), nb::arg("hscale"), nb::arg("t0"), nb::arg("resize"), nb::arg("resize_data"));

m.def("ARKodeReset",
    ARKodeReset, nb::arg("arkode_mem"), nb::arg("tR"), nb::arg("yR"));

m.def("ARKodeCreateMRIStepInnerStepper",
    ARKodeCreateMRIStepInnerStepper, 
    nb::arg("arkode_mem"), nb::arg("stepper"), 
    "Utility to wrap ARKODE as an MRIStepInnerStepper");

m.def("ARKodeSStolerances",
    ARKodeSStolerances, nb::arg("arkode_mem"), nb::arg("reltol"), nb::arg("abstol"));

m.def("ARKodeSVtolerances",
    ARKodeSVtolerances, nb::arg("arkode_mem"), nb::arg("reltol"), nb::arg("abstol"));

m.def("ARKodeWFtolerances",
    ARKodeWFtolerances, nb::arg("arkode_mem"), nb::arg("efun"));

m.def("ARKodeResStolerance",
    ARKodeResStolerance, nb::arg("arkode_mem"), nb::arg("rabstol"));

m.def("ARKodeResVtolerance",
    ARKodeResVtolerance, nb::arg("arkode_mem"), nb::arg("rabstol"));

m.def("ARKodeResFtolerance",
    ARKodeResFtolerance, nb::arg("arkode_mem"), nb::arg("rfun"));

m.def("ARKodeRootInit",
    ARKodeRootInit, nb::arg("arkode_mem"), nb::arg("nrtfn"), nb::arg("g"));

m.def("ARKodeSetRootDirection",
    ARKodeSetRootDirection, nb::arg("arkode_mem"), nb::arg("rootdir"));

m.def("ARKodeSetNoInactiveRootWarn",
    ARKodeSetNoInactiveRootWarn, nb::arg("arkode_mem"));

m.def("ARKodeSetDefaults",
    ARKodeSetDefaults, nb::arg("arkode_mem"));

m.def("ARKodeSetOrder",
    ARKodeSetOrder, nb::arg("arkode_mem"), nb::arg("maxord"));

m.def("ARKodeSetInterpolantType",
    ARKodeSetInterpolantType, nb::arg("arkode_mem"), nb::arg("itype"));

m.def("ARKodeSetInterpolantDegree",
    ARKodeSetInterpolantDegree, nb::arg("arkode_mem"), nb::arg("degree"));

m.def("ARKodeSetMaxNumSteps",
    ARKodeSetMaxNumSteps, nb::arg("arkode_mem"), nb::arg("mxsteps"));

m.def("ARKodeSetInterpolateStopTime",
    ARKodeSetInterpolateStopTime, nb::arg("arkode_mem"), nb::arg("interp"));

m.def("ARKodeSetStopTime",
    ARKodeSetStopTime, nb::arg("arkode_mem"), nb::arg("tstop"));

m.def("ARKodeClearStopTime",
    ARKodeClearStopTime, nb::arg("arkode_mem"));

m.def("ARKodeSetFixedStep",
    ARKodeSetFixedStep, nb::arg("arkode_mem"), nb::arg("hfixed"));

m.def("ARKodeSetStepDirection",
    ARKodeSetStepDirection, nb::arg("arkode_mem"), nb::arg("stepdir"));

m.def("ARKodeSetUserData",
    ARKodeSetUserData, nb::arg("arkode_mem"), nb::arg("user_data"));

m.def("ARKodeSetPostprocessStepFn",
    ARKodeSetPostprocessStepFn, nb::arg("arkode_mem"), nb::arg("ProcessStep"));

m.def("ARKodeSetPostprocessStageFn",
    ARKodeSetPostprocessStageFn, nb::arg("arkode_mem"), nb::arg("ProcessStage"));

m.def("ARKodeSetNonlinearSolver",
    ARKodeSetNonlinearSolver, nb::arg("arkode_mem"), nb::arg("NLS"));

m.def("ARKodeSetLinear",
    ARKodeSetLinear, nb::arg("arkode_mem"), nb::arg("timedepend"));

m.def("ARKodeSetNonlinear",
    ARKodeSetNonlinear, nb::arg("arkode_mem"));

m.def("ARKodeSetAutonomous",
    ARKodeSetAutonomous, nb::arg("arkode_mem"), nb::arg("autonomous"));

m.def("ARKodeSetNlsRhsFn",
    ARKodeSetNlsRhsFn, nb::arg("arkode_mem"), nb::arg("nls_fi"));

m.def("ARKodeSetDeduceImplicitRhs",
    ARKodeSetDeduceImplicitRhs, nb::arg("arkode_mem"), nb::arg("deduce"));

m.def("ARKodeSetNonlinCRDown",
    ARKodeSetNonlinCRDown, nb::arg("arkode_mem"), nb::arg("crdown"));

m.def("ARKodeSetNonlinRDiv",
    ARKodeSetNonlinRDiv, nb::arg("arkode_mem"), nb::arg("rdiv"));

m.def("ARKodeSetDeltaGammaMax",
    ARKodeSetDeltaGammaMax, nb::arg("arkode_mem"), nb::arg("dgmax"));

m.def("ARKodeSetLSetupFrequency",
    ARKodeSetLSetupFrequency, nb::arg("arkode_mem"), nb::arg("msbp"));

m.def("ARKodeSetPredictorMethod",
    ARKodeSetPredictorMethod, nb::arg("arkode_mem"), nb::arg("method"));

m.def("ARKodeSetMaxNonlinIters",
    ARKodeSetMaxNonlinIters, nb::arg("arkode_mem"), nb::arg("maxcor"));

m.def("ARKodeSetMaxConvFails",
    ARKodeSetMaxConvFails, nb::arg("arkode_mem"), nb::arg("maxncf"));

m.def("ARKodeSetNonlinConvCoef",
    ARKodeSetNonlinConvCoef, nb::arg("arkode_mem"), nb::arg("nlscoef"));

m.def("ARKodeSetStagePredictFn",
    ARKodeSetStagePredictFn, nb::arg("arkode_mem"), nb::arg("PredictStage"));

m.def("ARKodeSetAdaptController",
    ARKodeSetAdaptController, nb::arg("arkode_mem"), nb::arg("C"));

m.def("ARKodeSetAdaptivityAdjustment",
    ARKodeSetAdaptivityAdjustment, nb::arg("arkode_mem"), nb::arg("adjust"));

m.def("ARKodeSetCFLFraction",
    ARKodeSetCFLFraction, nb::arg("arkode_mem"), nb::arg("cfl_frac"));

m.def("ARKodeSetErrorBias",
    ARKodeSetErrorBias, nb::arg("arkode_mem"), nb::arg("bias"));

m.def("ARKodeSetSafetyFactor",
    ARKodeSetSafetyFactor, nb::arg("arkode_mem"), nb::arg("safety"));

m.def("ARKodeSetMaxGrowth",
    ARKodeSetMaxGrowth, nb::arg("arkode_mem"), nb::arg("mx_growth"));

m.def("ARKodeSetMinReduction",
    ARKodeSetMinReduction, nb::arg("arkode_mem"), nb::arg("eta_min"));

m.def("ARKodeSetFixedStepBounds",
    ARKodeSetFixedStepBounds, nb::arg("arkode_mem"), nb::arg("lb"), nb::arg("ub"));

m.def("ARKodeSetMaxFirstGrowth",
    ARKodeSetMaxFirstGrowth, nb::arg("arkode_mem"), nb::arg("etamx1"));

m.def("ARKodeSetMaxEFailGrowth",
    ARKodeSetMaxEFailGrowth, nb::arg("arkode_mem"), nb::arg("etamxf"));

m.def("ARKodeSetSmallNumEFails",
    ARKodeSetSmallNumEFails, nb::arg("arkode_mem"), nb::arg("small_nef"));

m.def("ARKodeSetMaxCFailGrowth",
    ARKodeSetMaxCFailGrowth, nb::arg("arkode_mem"), nb::arg("etacf"));

m.def("ARKodeSetStabilityFn",
    ARKodeSetStabilityFn, nb::arg("arkode_mem"), nb::arg("EStab"), nb::arg("estab_data"));

m.def("ARKodeSetMaxErrTestFails",
    ARKodeSetMaxErrTestFails, nb::arg("arkode_mem"), nb::arg("maxnef"));

m.def("ARKodeSetConstraints",
    ARKodeSetConstraints, nb::arg("arkode_mem"), nb::arg("constraints"));

m.def("ARKodeSetMaxHnilWarns",
    ARKodeSetMaxHnilWarns, nb::arg("arkode_mem"), nb::arg("mxhnil"));

m.def("ARKodeSetInitStep",
    ARKodeSetInitStep, nb::arg("arkode_mem"), nb::arg("hin"));

m.def("ARKodeSetMinStep",
    ARKodeSetMinStep, nb::arg("arkode_mem"), nb::arg("hmin"));

m.def("ARKodeSetMaxStep",
    ARKodeSetMaxStep, nb::arg("arkode_mem"), nb::arg("hmax"));

m.def("ARKodeSetMaxNumConstrFails",
    ARKodeSetMaxNumConstrFails, nb::arg("arkode_mem"), nb::arg("maxfails"));

m.def("ARKodeSetAccumulatedErrorType",
    ARKodeSetAccumulatedErrorType, nb::arg("arkode_mem"), nb::arg("accum_type"));

m.def("ARKodeResetAccumulatedError",
    ARKodeResetAccumulatedError, nb::arg("arkode_mem"));

m.def("ARKodeEvolve",
    ARKodeEvolve, 
    nb::arg("arkode_mem"), nb::arg("tout"), nb::arg("yout"), nb::arg("tret"), nb::arg("itask"), 
    "Integrate the ODE over an interval in t");

m.def("ARKodeGetDky",
    ARKodeGetDky, 
    nb::arg("arkode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dky"), 
    "Computes the kth derivative of the y function at time t");

m.def("ARKodeComputeState",
    ARKodeComputeState, 
    nb::arg("arkode_mem"), nb::arg("zcor"), nb::arg("z"), 
    "Utility function to update/compute y based on zcor");

m.def("ARKodeGetNumRhsEvals",
    ARKodeGetNumRhsEvals, nb::arg("arkode_mem"), nb::arg("partition_index"), nb::arg("num_rhs_evals"));

m.def("ARKodeGetNumStepAttempts",
    ARKodeGetNumStepAttempts, nb::arg("arkode_mem"), nb::arg("step_attempts"));

m.def("ARKodeGetWorkSpace",
    ARKodeGetWorkSpace, nb::arg("arkode_mem"), nb::arg("lenrw"), nb::arg("leniw"));

m.def("ARKodeGetNumSteps",
    ARKodeGetNumSteps, nb::arg("arkode_mem"), nb::arg("nsteps"));

m.def("ARKodeGetLastStep",
    ARKodeGetLastStep, nb::arg("arkode_mem"), nb::arg("hlast"));

m.def("ARKodeGetCurrentStep",
    ARKodeGetCurrentStep, nb::arg("arkode_mem"), nb::arg("hcur"));

m.def("ARKodeGetStepDirection",
    ARKodeGetStepDirection, nb::arg("arkode_mem"), nb::arg("stepdir"));

m.def("ARKodeGetErrWeights",
    ARKodeGetErrWeights, nb::arg("arkode_mem"), nb::arg("eweight"));

m.def("ARKodeGetNumGEvals",
    ARKodeGetNumGEvals, nb::arg("arkode_mem"), nb::arg("ngevals"));

m.def("ARKodeGetRootInfo",
    ARKodeGetRootInfo, nb::arg("arkode_mem"), nb::arg("rootsfound"));

m.def("ARKodePrintAllStats",
    ARKodePrintAllStats, nb::arg("arkode_mem"), nb::arg("outfile"), nb::arg("fmt"));

m.def("ARKodeGetReturnFlagName",
    ARKodeGetReturnFlagName, nb::arg("flag"));

m.def("ARKodeWriteParameters",
    ARKodeWriteParameters, nb::arg("arkode_mem"), nb::arg("fp"));

m.def("ARKodeGetNumExpSteps",
    ARKodeGetNumExpSteps, nb::arg("arkode_mem"), nb::arg("expsteps"));

m.def("ARKodeGetNumAccSteps",
    ARKodeGetNumAccSteps, nb::arg("arkode_mem"), nb::arg("accsteps"));

m.def("ARKodeGetNumErrTestFails",
    ARKodeGetNumErrTestFails, nb::arg("arkode_mem"), nb::arg("netfails"));

m.def("ARKodeGetEstLocalErrors",
    ARKodeGetEstLocalErrors, nb::arg("arkode_mem"), nb::arg("ele"));

m.def("ARKodeGetActualInitStep",
    ARKodeGetActualInitStep, nb::arg("arkode_mem"), nb::arg("hinused"));

m.def("ARKodeGetTolScaleFactor",
    ARKodeGetTolScaleFactor, nb::arg("arkode_mem"), nb::arg("tolsfac"));

m.def("ARKodeGetNumConstrFails",
    ARKodeGetNumConstrFails, nb::arg("arkode_mem"), nb::arg("nconstrfails"));

m.def("ARKodeGetStepStats",
    ARKodeGetStepStats, nb::arg("arkode_mem"), nb::arg("nsteps"), nb::arg("hinused"), nb::arg("hlast"), nb::arg("hcur"), nb::arg("tcur"));

m.def("ARKodeGetAccumulatedError",
    ARKodeGetAccumulatedError, nb::arg("arkode_mem"), nb::arg("accum_error"));

m.def("ARKodeGetNumLinSolvSetups",
    ARKodeGetNumLinSolvSetups, nb::arg("arkode_mem"), nb::arg("nlinsetups"));

m.def("ARKodeGetCurrentTime",
    ARKodeGetCurrentTime, nb::arg("arkode_mem"), nb::arg("tcur"));

m.def("ARKodeGetCurrentState",
    ARKodeGetCurrentState, nb::arg("arkode_mem"), nb::arg("state"));

m.def("ARKodeGetCurrentGamma",
    ARKodeGetCurrentGamma, nb::arg("arkode_mem"), nb::arg("gamma"));

m.def("ARKodeGetNumNonlinSolvIters",
    ARKodeGetNumNonlinSolvIters, nb::arg("arkode_mem"), nb::arg("nniters"));

m.def("ARKodeGetNumNonlinSolvConvFails",
    ARKodeGetNumNonlinSolvConvFails, nb::arg("arkode_mem"), nb::arg("nnfails"));

m.def("ARKodeGetNonlinSolvStats",
    ARKodeGetNonlinSolvStats, nb::arg("arkode_mem"), nb::arg("nniters"), nb::arg("nnfails"));

m.def("ARKodeGetNumStepSolveFails",
    ARKodeGetNumStepSolveFails, nb::arg("arkode_mem"), nb::arg("nncfails"));

m.def("ARKodeGetJac",
    ARKodeGetJac, nb::arg("arkode_mem"), nb::arg("J"));

m.def("ARKodeGetJacTime",
    ARKodeGetJacTime, nb::arg("arkode_mem"), nb::arg("t_J"));

m.def("ARKodeGetJacNumSteps",
    ARKodeGetJacNumSteps, nb::arg("arkode_mem"), nb::arg("nst_J"));

m.def("ARKodeGetLinWorkSpace",
    ARKodeGetLinWorkSpace, nb::arg("arkode_mem"), nb::arg("lenrwLS"), nb::arg("leniwLS"));

m.def("ARKodeGetNumJacEvals",
    ARKodeGetNumJacEvals, nb::arg("arkode_mem"), nb::arg("njevals"));

m.def("ARKodeGetNumPrecEvals",
    ARKodeGetNumPrecEvals, nb::arg("arkode_mem"), nb::arg("npevals"));

m.def("ARKodeGetNumPrecSolves",
    ARKodeGetNumPrecSolves, nb::arg("arkode_mem"), nb::arg("npsolves"));

m.def("ARKodeGetNumLinIters",
    ARKodeGetNumLinIters, nb::arg("arkode_mem"), nb::arg("nliters"));

m.def("ARKodeGetNumLinConvFails",
    ARKodeGetNumLinConvFails, nb::arg("arkode_mem"), nb::arg("nlcfails"));

m.def("ARKodeGetNumJTSetupEvals",
    ARKodeGetNumJTSetupEvals, nb::arg("arkode_mem"), nb::arg("njtsetups"));

m.def("ARKodeGetNumJtimesEvals",
    ARKodeGetNumJtimesEvals, nb::arg("arkode_mem"), nb::arg("njvevals"));

m.def("ARKodeGetNumLinRhsEvals",
    ARKodeGetNumLinRhsEvals, nb::arg("arkode_mem"), nb::arg("nfevalsLS"));

m.def("ARKodeGetLastLinFlag",
    ARKodeGetLastLinFlag, nb::arg("arkode_mem"), nb::arg("flag"));

m.def("ARKodeGetLinReturnFlagName",
    ARKodeGetLinReturnFlagName, nb::arg("flag"));

m.def("ARKodeGetCurrentMassMatrix",
    ARKodeGetCurrentMassMatrix, nb::arg("arkode_mem"), nb::arg("M"));

m.def("ARKodeGetResWeights",
    ARKodeGetResWeights, nb::arg("arkode_mem"), nb::arg("rweight"));

m.def("ARKodeGetMassWorkSpace",
    ARKodeGetMassWorkSpace, nb::arg("arkode_mem"), nb::arg("lenrwMLS"), nb::arg("leniwMLS"));

m.def("ARKodeGetNumMassSetups",
    ARKodeGetNumMassSetups, nb::arg("arkode_mem"), nb::arg("nmsetups"));

m.def("ARKodeGetNumMassMultSetups",
    ARKodeGetNumMassMultSetups, nb::arg("arkode_mem"), nb::arg("nmvsetups"));

m.def("ARKodeGetNumMassMult",
    ARKodeGetNumMassMult, nb::arg("arkode_mem"), nb::arg("nmvevals"));

m.def("ARKodeGetNumMassSolves",
    ARKodeGetNumMassSolves, nb::arg("arkode_mem"), nb::arg("nmsolves"));

m.def("ARKodeGetNumMassPrecEvals",
    ARKodeGetNumMassPrecEvals, nb::arg("arkode_mem"), nb::arg("nmpevals"));

m.def("ARKodeGetNumMassPrecSolves",
    ARKodeGetNumMassPrecSolves, nb::arg("arkode_mem"), nb::arg("nmpsolves"));

m.def("ARKodeGetNumMassIters",
    ARKodeGetNumMassIters, nb::arg("arkode_mem"), nb::arg("nmiters"));

m.def("ARKodeGetNumMassConvFails",
    ARKodeGetNumMassConvFails, nb::arg("arkode_mem"), nb::arg("nmcfails"));

m.def("ARKodeGetNumMTSetups",
    ARKodeGetNumMTSetups, nb::arg("arkode_mem"), nb::arg("nmtsetups"));

m.def("ARKodeGetLastMassFlag",
    ARKodeGetLastMassFlag, nb::arg("arkode_mem"), nb::arg("flag"));

m.def("ARKodePrintMem",
    ARKodePrintMem, 
    nb::arg("arkode_mem"), nb::arg("outfile"), 
    "Output the ARKODE memory structure (useful when debugging)");

m.def("ARKodeSetRelaxFn",
    ARKodeSetRelaxFn, nb::arg("arkode_mem"), nb::arg("rfn"), nb::arg("rjac"));

m.def("ARKodeSetRelaxEtaFail",
    ARKodeSetRelaxEtaFail, nb::arg("arkode_mem"), nb::arg("eta_rf"));

m.def("ARKodeSetRelaxLowerBound",
    ARKodeSetRelaxLowerBound, nb::arg("arkode_mem"), nb::arg("lower"));

m.def("ARKodeSetRelaxMaxFails",
    ARKodeSetRelaxMaxFails, nb::arg("arkode_mem"), nb::arg("max_fails"));

m.def("ARKodeSetRelaxMaxIters",
    ARKodeSetRelaxMaxIters, nb::arg("arkode_mem"), nb::arg("max_iters"));

m.def("ARKodeSetRelaxSolver",
    ARKodeSetRelaxSolver, nb::arg("arkode_mem"), nb::arg("solver"));

m.def("ARKodeSetRelaxResTol",
    ARKodeSetRelaxResTol, nb::arg("arkode_mem"), nb::arg("res_tol"));

m.def("ARKodeSetRelaxTol",
    ARKodeSetRelaxTol, nb::arg("arkode_mem"), nb::arg("rel_tol"), nb::arg("abs_tol"));

m.def("ARKodeSetRelaxUpperBound",
    ARKodeSetRelaxUpperBound, nb::arg("arkode_mem"), nb::arg("upper"));

m.def("ARKodeGetNumRelaxFnEvals",
    ARKodeGetNumRelaxFnEvals, nb::arg("arkode_mem"), nb::arg("r_evals"));

m.def("ARKodeGetNumRelaxJacEvals",
    ARKodeGetNumRelaxJacEvals, nb::arg("arkode_mem"), nb::arg("J_evals"));

m.def("ARKodeGetNumRelaxFails",
    ARKodeGetNumRelaxFails, nb::arg("arkode_mem"), nb::arg("relax_fails"));

m.def("ARKodeGetNumRelaxBoundFails",
    ARKodeGetNumRelaxBoundFails, nb::arg("arkode_mem"), nb::arg("fails"));

m.def("ARKodeGetNumRelaxSolveFails",
    ARKodeGetNumRelaxSolveFails, nb::arg("arkode_mem"), nb::arg("fails"));

m.def("ARKodeGetNumRelaxSolveIters",
    ARKodeGetNumRelaxSolveIters, nb::arg("arkode_mem"), nb::arg("iters"));

m.def("ARKodeCreateSUNStepper",
    ARKodeCreateSUNStepper, 
    nb::arg("arkode_mem"), nb::arg("stepper"), 
    "SUNStepper functions");
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
