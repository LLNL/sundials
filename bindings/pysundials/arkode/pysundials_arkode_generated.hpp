// #ifndef _ARKODE_H
// 
// #ifdef __cplusplus 
// #endif
// 
m.attr("ARK_NORMAL") = 1;
m.attr("ARK_ONE_STEP") = 2;
m.attr("ARK_ADAPT_CUSTOM") = -1;
m.attr("ARK_ADAPT_PID") = 0;
m.attr("ARK_ADAPT_PI") = 1;
m.attr("ARK_ADAPT_I") = 2;
m.attr("ARK_ADAPT_EXP_GUS") = 3;
m.attr("ARK_ADAPT_IMP_GUS") = 4;
m.attr("ARK_ADAPT_IMEX_GUS") = 5;
m.attr("ARK_FULLRHS_START") = 0;
m.attr("ARK_FULLRHS_END") = 1;
m.attr("ARK_FULLRHS_OTHER") = 2;
m.attr("ARK_INTERP_MAX_DEGREE") = 5;
m.attr("ARK_INTERP_NONE") = -1;
m.attr("ARK_INTERP_HERMITE") = 0;
m.attr("ARK_INTERP_LAGRANGE") = 1;
m.attr("ARK_SUCCESS") = 0;
m.attr("ARK_TSTOP_RETURN") = 1;
m.attr("ARK_ROOT_RETURN") = 2;
m.attr("ARK_WARNING") = 99;
m.attr("ARK_TOO_MUCH_WORK") = -1;
m.attr("ARK_TOO_MUCH_ACC") = -2;
m.attr("ARK_ERR_FAILURE") = -3;
m.attr("ARK_CONV_FAILURE") = -4;
m.attr("ARK_LINIT_FAIL") = -5;
m.attr("ARK_LSETUP_FAIL") = -6;
m.attr("ARK_LSOLVE_FAIL") = -7;
m.attr("ARK_RHSFUNC_FAIL") = -8;
m.attr("ARK_FIRST_RHSFUNC_ERR") = -9;
m.attr("ARK_REPTD_RHSFUNC_ERR") = -10;
m.attr("ARK_UNREC_RHSFUNC_ERR") = -11;
m.attr("ARK_RTFUNC_FAIL") = -12;
m.attr("ARK_LFREE_FAIL") = -13;
m.attr("ARK_MASSINIT_FAIL") = -14;
m.attr("ARK_MASSSETUP_FAIL") = -15;
m.attr("ARK_MASSSOLVE_FAIL") = -16;
m.attr("ARK_MASSFREE_FAIL") = -17;
m.attr("ARK_MASSMULT_FAIL") = -18;
m.attr("ARK_CONSTR_FAIL") = -19;
m.attr("ARK_MEM_FAIL") = -20;
m.attr("ARK_MEM_NULL") = -21;
m.attr("ARK_ILL_INPUT") = -22;
m.attr("ARK_NO_MALLOC") = -23;
m.attr("ARK_BAD_K") = -24;
m.attr("ARK_BAD_T") = -25;
m.attr("ARK_BAD_DKY") = -26;
m.attr("ARK_TOO_CLOSE") = -27;
m.attr("ARK_VECTOROP_ERR") = -28;
m.attr("ARK_NLS_INIT_FAIL") = -29;
m.attr("ARK_NLS_SETUP_FAIL") = -30;
m.attr("ARK_NLS_SETUP_RECVR") = -31;
m.attr("ARK_NLS_OP_ERR") = -32;
m.attr("ARK_INNERSTEP_ATTACH_ERR") = -33;
m.attr("ARK_INNERSTEP_FAIL") = -34;
m.attr("ARK_OUTERTOINNER_FAIL") = -35;
m.attr("ARK_INNERTOOUTER_FAIL") = -36;
m.attr("ARK_POSTPROCESS_FAIL") = -37;
m.attr("ARK_POSTPROCESS_STEP_FAIL") = -37;
m.attr("ARK_POSTPROCESS_STAGE_FAIL") = -38;
m.attr("ARK_USER_PREDICT_FAIL") = -39;
m.attr("ARK_INTERP_FAIL") = -40;
m.attr("ARK_INVALID_TABLE") = -41;
m.attr("ARK_CONTEXT_ERR") = -42;
m.attr("ARK_RELAX_FAIL") = -43;
m.attr("ARK_RELAX_MEM_NULL") = -44;
m.attr("ARK_RELAX_FUNC_FAIL") = -45;
m.attr("ARK_RELAX_JAC_FAIL") = -46;
m.attr("ARK_CONTROLLER_ERR") = -47;
m.attr("ARK_STEPPER_UNSUPPORTED") = -48;
m.attr("ARK_DOMEIG_FAIL") = -49;
m.attr("ARK_MAX_STAGE_LIMIT_FAIL") = -50;
m.attr("ARK_SUNSTEPPER_ERR") = -51;
m.attr("ARK_STEP_DIRECTION_ERR") = -52;
m.attr("ARK_ADJ_CHECKPOINT_FAIL") = -53;
m.attr("ARK_ADJ_RECOMPUTE_FAIL") = -54;
m.attr("ARK_SUNADJSTEPPER_ERR") = -55;
m.attr("ARK_UNRECOGNIZED_ERROR") = -99;

m.def("ARKodeResize",
    ARKodeResize, nb::arg("arkode_mem"), nb::arg("ynew"), nb::arg("hscale"), nb::arg("t0"), nb::arg("resize"), nb::arg("resize_data"));

m.def("ARKodeReset",
    ARKodeReset, nb::arg("arkode_mem"), nb::arg("tR"), nb::arg("yR"));

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
    [](void * arkode_mem, int rootdir) -> std::tuple<int, int>
    {
        auto ARKodeSetRootDirection_adapt_modifiable_immutable_to_return = [](void * arkode_mem, int rootdir) -> std::tuple<int, int>
        {
            int * rootdir_adapt_modifiable = & rootdir;

            int r = ARKodeSetRootDirection(arkode_mem, rootdir_adapt_modifiable);
            return std::make_tuple(r, rootdir);
        };

        return ARKodeSetRootDirection_adapt_modifiable_immutable_to_return(arkode_mem, rootdir);
    },     nb::arg("arkode_mem"), nb::arg("rootdir"));

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

m.def("ARKodeSetOwnUserData",
    ARKodeSetOwnUserData, nb::arg("ark_mem"), nb::arg("own_user_data"));

m.def("ARKodeSetNonlinearSolver",
    ARKodeSetNonlinearSolver, nb::arg("arkode_mem"), nb::arg("NLS"));

m.def("ARKodeSetLinear",
    ARKodeSetLinear, nb::arg("arkode_mem"), nb::arg("timedepend"));

m.def("ARKodeSetNonlinear",
    ARKodeSetNonlinear, nb::arg("arkode_mem"));

m.def("ARKodeSetAutonomous",
    ARKodeSetAutonomous, nb::arg("arkode_mem"), nb::arg("autonomous"));

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

m.def("ARKodeSetAdaptController",
    ARKodeSetAdaptController, nb::arg("arkode_mem"), nb::arg("C"));

m.def("ARKodeSetAdaptControllerByName",
    ARKodeSetAdaptControllerByName, nb::arg("arkode_mem"), nb::arg("cname"));

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

m.def("ARKodeSetAdjointCheckpointScheme",
    ARKodeSetAdjointCheckpointScheme, nb::arg("arkode_mem"), nb::arg("checkpoint_scheme"));

m.def("ARKodeSetAdjointCheckpointIndex",
    ARKodeSetAdjointCheckpointIndex, nb::arg("arkode_mem"), nb::arg("step_index"));

m.def("ARKodeSetUseCompensatedSums",
    ARKodeSetUseCompensatedSums, nb::arg("arkode_mem"), nb::arg("onoff"));

m.def("ARKodeSetAccumulatedErrorType",
    ARKodeSetAccumulatedErrorType, nb::arg("arkode_mem"), nb::arg("accum_type"));

m.def("ARKodeResetAccumulatedError",
    ARKodeResetAccumulatedError, nb::arg("arkode_mem"));

m.def("ARKodeEvolve",
    [](void * arkode_mem, double tout, N_Vector yout, double tret, int itask) -> std::tuple<int, double>
    {
        auto ARKodeEvolve_adapt_modifiable_immutable_to_return = [](void * arkode_mem, double tout, N_Vector yout, double tret, int itask) -> std::tuple<int, double>
        {
            double * tret_adapt_modifiable = & tret;

            int r = ARKodeEvolve(arkode_mem, tout, yout, tret_adapt_modifiable, itask);
            return std::make_tuple(r, tret);
        };

        return ARKodeEvolve_adapt_modifiable_immutable_to_return(arkode_mem, tout, yout, tret, itask);
    },     nb::arg("arkode_mem"), nb::arg("tout"), nb::arg("yout"), nb::arg("tret"), nb::arg("itask"));

m.def("ARKodeGetDky",
    ARKodeGetDky, nb::arg("arkode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dky"));

m.def("ARKodeComputeState",
    ARKodeComputeState, nb::arg("arkode_mem"), nb::arg("zcor"), nb::arg("z"));

m.def("ARKodeGetNumRhsEvals",
    [](void * arkode_mem, int partition_index, long num_rhs_evals) -> std::tuple<int, long>
    {
        auto ARKodeGetNumRhsEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, int partition_index, long num_rhs_evals) -> std::tuple<int, long>
        {
            long * num_rhs_evals_adapt_modifiable = & num_rhs_evals;

            int r = ARKodeGetNumRhsEvals(arkode_mem, partition_index, num_rhs_evals_adapt_modifiable);
            return std::make_tuple(r, num_rhs_evals);
        };

        return ARKodeGetNumRhsEvals_adapt_modifiable_immutable_to_return(arkode_mem, partition_index, num_rhs_evals);
    },     nb::arg("arkode_mem"), nb::arg("partition_index"), nb::arg("num_rhs_evals"));

m.def("ARKodeGetNumStepAttempts",
    [](void * arkode_mem, long step_attempts) -> std::tuple<int, long>
    {
        auto ARKodeGetNumStepAttempts_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long step_attempts) -> std::tuple<int, long>
        {
            long * step_attempts_adapt_modifiable = & step_attempts;

            int r = ARKodeGetNumStepAttempts(arkode_mem, step_attempts_adapt_modifiable);
            return std::make_tuple(r, step_attempts);
        };

        return ARKodeGetNumStepAttempts_adapt_modifiable_immutable_to_return(arkode_mem, step_attempts);
    },     nb::arg("arkode_mem"), nb::arg("step_attempts"));

m.def("ARKodeGetNumSteps",
    [](void * arkode_mem, long nsteps) -> std::tuple<int, long>
    {
        auto ARKodeGetNumSteps_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nsteps) -> std::tuple<int, long>
        {
            long * nsteps_adapt_modifiable = & nsteps;

            int r = ARKodeGetNumSteps(arkode_mem, nsteps_adapt_modifiable);
            return std::make_tuple(r, nsteps);
        };

        return ARKodeGetNumSteps_adapt_modifiable_immutable_to_return(arkode_mem, nsteps);
    },     nb::arg("arkode_mem"), nb::arg("nsteps"));

m.def("ARKodeGetLastStep",
    [](void * arkode_mem, double hlast) -> std::tuple<int, double>
    {
        auto ARKodeGetLastStep_adapt_modifiable_immutable_to_return = [](void * arkode_mem, double hlast) -> std::tuple<int, double>
        {
            double * hlast_adapt_modifiable = & hlast;

            int r = ARKodeGetLastStep(arkode_mem, hlast_adapt_modifiable);
            return std::make_tuple(r, hlast);
        };

        return ARKodeGetLastStep_adapt_modifiable_immutable_to_return(arkode_mem, hlast);
    },     nb::arg("arkode_mem"), nb::arg("hlast"));

m.def("ARKodeGetCurrentStep",
    [](void * arkode_mem, double hcur) -> std::tuple<int, double>
    {
        auto ARKodeGetCurrentStep_adapt_modifiable_immutable_to_return = [](void * arkode_mem, double hcur) -> std::tuple<int, double>
        {
            double * hcur_adapt_modifiable = & hcur;

            int r = ARKodeGetCurrentStep(arkode_mem, hcur_adapt_modifiable);
            return std::make_tuple(r, hcur);
        };

        return ARKodeGetCurrentStep_adapt_modifiable_immutable_to_return(arkode_mem, hcur);
    },     nb::arg("arkode_mem"), nb::arg("hcur"));

m.def("ARKodeGetStepDirection",
    [](void * arkode_mem, double stepdir) -> std::tuple<int, double>
    {
        auto ARKodeGetStepDirection_adapt_modifiable_immutable_to_return = [](void * arkode_mem, double stepdir) -> std::tuple<int, double>
        {
            double * stepdir_adapt_modifiable = & stepdir;

            int r = ARKodeGetStepDirection(arkode_mem, stepdir_adapt_modifiable);
            return std::make_tuple(r, stepdir);
        };

        return ARKodeGetStepDirection_adapt_modifiable_immutable_to_return(arkode_mem, stepdir);
    },     nb::arg("arkode_mem"), nb::arg("stepdir"));

m.def("ARKodeGetErrWeights",
    ARKodeGetErrWeights, nb::arg("arkode_mem"), nb::arg("eweight"));

m.def("ARKodeGetNumGEvals",
    [](void * arkode_mem, long ngevals) -> std::tuple<int, long>
    {
        auto ARKodeGetNumGEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long ngevals) -> std::tuple<int, long>
        {
            long * ngevals_adapt_modifiable = & ngevals;

            int r = ARKodeGetNumGEvals(arkode_mem, ngevals_adapt_modifiable);
            return std::make_tuple(r, ngevals);
        };

        return ARKodeGetNumGEvals_adapt_modifiable_immutable_to_return(arkode_mem, ngevals);
    },     nb::arg("arkode_mem"), nb::arg("ngevals"));

m.def("ARKodeGetRootInfo",
    [](void * arkode_mem, int rootsfound) -> std::tuple<int, int>
    {
        auto ARKodeGetRootInfo_adapt_modifiable_immutable_to_return = [](void * arkode_mem, int rootsfound) -> std::tuple<int, int>
        {
            int * rootsfound_adapt_modifiable = & rootsfound;

            int r = ARKodeGetRootInfo(arkode_mem, rootsfound_adapt_modifiable);
            return std::make_tuple(r, rootsfound);
        };

        return ARKodeGetRootInfo_adapt_modifiable_immutable_to_return(arkode_mem, rootsfound);
    },     nb::arg("arkode_mem"), nb::arg("rootsfound"));

m.def("ARKodePrintAllStats",
    ARKodePrintAllStats, nb::arg("arkode_mem"), nb::arg("outfile"), nb::arg("fmt"));

m.def("ARKodeGetReturnFlagName",
    ARKodeGetReturnFlagName, nb::arg("flag"));

m.def("ARKodeWriteParameters",
    ARKodeWriteParameters, nb::arg("arkode_mem"), nb::arg("fp"));

m.def("ARKodeGetNumExpSteps",
    [](void * arkode_mem, long expsteps) -> std::tuple<int, long>
    {
        auto ARKodeGetNumExpSteps_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long expsteps) -> std::tuple<int, long>
        {
            long * expsteps_adapt_modifiable = & expsteps;

            int r = ARKodeGetNumExpSteps(arkode_mem, expsteps_adapt_modifiable);
            return std::make_tuple(r, expsteps);
        };

        return ARKodeGetNumExpSteps_adapt_modifiable_immutable_to_return(arkode_mem, expsteps);
    },     nb::arg("arkode_mem"), nb::arg("expsteps"));

m.def("ARKodeGetNumAccSteps",
    [](void * arkode_mem, long accsteps) -> std::tuple<int, long>
    {
        auto ARKodeGetNumAccSteps_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long accsteps) -> std::tuple<int, long>
        {
            long * accsteps_adapt_modifiable = & accsteps;

            int r = ARKodeGetNumAccSteps(arkode_mem, accsteps_adapt_modifiable);
            return std::make_tuple(r, accsteps);
        };

        return ARKodeGetNumAccSteps_adapt_modifiable_immutable_to_return(arkode_mem, accsteps);
    },     nb::arg("arkode_mem"), nb::arg("accsteps"));

m.def("ARKodeGetNumErrTestFails",
    [](void * arkode_mem, long netfails) -> std::tuple<int, long>
    {
        auto ARKodeGetNumErrTestFails_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long netfails) -> std::tuple<int, long>
        {
            long * netfails_adapt_modifiable = & netfails;

            int r = ARKodeGetNumErrTestFails(arkode_mem, netfails_adapt_modifiable);
            return std::make_tuple(r, netfails);
        };

        return ARKodeGetNumErrTestFails_adapt_modifiable_immutable_to_return(arkode_mem, netfails);
    },     nb::arg("arkode_mem"), nb::arg("netfails"));

m.def("ARKodeGetEstLocalErrors",
    ARKodeGetEstLocalErrors, nb::arg("arkode_mem"), nb::arg("ele"));

m.def("ARKodeGetActualInitStep",
    [](void * arkode_mem, double hinused) -> std::tuple<int, double>
    {
        auto ARKodeGetActualInitStep_adapt_modifiable_immutable_to_return = [](void * arkode_mem, double hinused) -> std::tuple<int, double>
        {
            double * hinused_adapt_modifiable = & hinused;

            int r = ARKodeGetActualInitStep(arkode_mem, hinused_adapt_modifiable);
            return std::make_tuple(r, hinused);
        };

        return ARKodeGetActualInitStep_adapt_modifiable_immutable_to_return(arkode_mem, hinused);
    },     nb::arg("arkode_mem"), nb::arg("hinused"));

m.def("ARKodeGetTolScaleFactor",
    [](void * arkode_mem, double tolsfac) -> std::tuple<int, double>
    {
        auto ARKodeGetTolScaleFactor_adapt_modifiable_immutable_to_return = [](void * arkode_mem, double tolsfac) -> std::tuple<int, double>
        {
            double * tolsfac_adapt_modifiable = & tolsfac;

            int r = ARKodeGetTolScaleFactor(arkode_mem, tolsfac_adapt_modifiable);
            return std::make_tuple(r, tolsfac);
        };

        return ARKodeGetTolScaleFactor_adapt_modifiable_immutable_to_return(arkode_mem, tolsfac);
    },     nb::arg("arkode_mem"), nb::arg("tolsfac"));

m.def("ARKodeGetNumConstrFails",
    [](void * arkode_mem, long nconstrfails) -> std::tuple<int, long>
    {
        auto ARKodeGetNumConstrFails_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nconstrfails) -> std::tuple<int, long>
        {
            long * nconstrfails_adapt_modifiable = & nconstrfails;

            int r = ARKodeGetNumConstrFails(arkode_mem, nconstrfails_adapt_modifiable);
            return std::make_tuple(r, nconstrfails);
        };

        return ARKodeGetNumConstrFails_adapt_modifiable_immutable_to_return(arkode_mem, nconstrfails);
    },     nb::arg("arkode_mem"), nb::arg("nconstrfails"));

m.def("ARKodeGetStepStats",
    [](void * arkode_mem, long nsteps, double hinused, double hlast, double hcur, double tcur) -> std::tuple<int, long, double, double, double, double>
    {
        auto ARKodeGetStepStats_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nsteps, double hinused, double hlast, double hcur, double tcur) -> std::tuple<int, long, double, double, double, double>
        {
            long * nsteps_adapt_modifiable = & nsteps;
            double * hinused_adapt_modifiable = & hinused;
            double * hlast_adapt_modifiable = & hlast;
            double * hcur_adapt_modifiable = & hcur;
            double * tcur_adapt_modifiable = & tcur;

            int r = ARKodeGetStepStats(arkode_mem, nsteps_adapt_modifiable, hinused_adapt_modifiable, hlast_adapt_modifiable, hcur_adapt_modifiable, tcur_adapt_modifiable);
            return std::make_tuple(r, nsteps, hinused, hlast, hcur, tcur);
        };

        return ARKodeGetStepStats_adapt_modifiable_immutable_to_return(arkode_mem, nsteps, hinused, hlast, hcur, tcur);
    },     nb::arg("arkode_mem"), nb::arg("nsteps"), nb::arg("hinused"), nb::arg("hlast"), nb::arg("hcur"), nb::arg("tcur"));

m.def("ARKodeGetAccumulatedError",
    [](void * arkode_mem, double accum_error) -> std::tuple<int, double>
    {
        auto ARKodeGetAccumulatedError_adapt_modifiable_immutable_to_return = [](void * arkode_mem, double accum_error) -> std::tuple<int, double>
        {
            double * accum_error_adapt_modifiable = & accum_error;

            int r = ARKodeGetAccumulatedError(arkode_mem, accum_error_adapt_modifiable);
            return std::make_tuple(r, accum_error);
        };

        return ARKodeGetAccumulatedError_adapt_modifiable_immutable_to_return(arkode_mem, accum_error);
    },     nb::arg("arkode_mem"), nb::arg("accum_error"));

m.def("ARKodeGetNumLinSolvSetups",
    [](void * arkode_mem, long nlinsetups) -> std::tuple<int, long>
    {
        auto ARKodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nlinsetups) -> std::tuple<int, long>
        {
            long * nlinsetups_adapt_modifiable = & nlinsetups;

            int r = ARKodeGetNumLinSolvSetups(arkode_mem, nlinsetups_adapt_modifiable);
            return std::make_tuple(r, nlinsetups);
        };

        return ARKodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return(arkode_mem, nlinsetups);
    },     nb::arg("arkode_mem"), nb::arg("nlinsetups"));

m.def("ARKodeGetCurrentTime",
    [](void * arkode_mem, double tcur) -> std::tuple<int, double>
    {
        auto ARKodeGetCurrentTime_adapt_modifiable_immutable_to_return = [](void * arkode_mem, double tcur) -> std::tuple<int, double>
        {
            double * tcur_adapt_modifiable = & tcur;

            int r = ARKodeGetCurrentTime(arkode_mem, tcur_adapt_modifiable);
            return std::make_tuple(r, tcur);
        };

        return ARKodeGetCurrentTime_adapt_modifiable_immutable_to_return(arkode_mem, tcur);
    },     nb::arg("arkode_mem"), nb::arg("tcur"));

m.def("ARKodeGetCurrentGamma",
    [](void * arkode_mem, double gamma) -> std::tuple<int, double>
    {
        auto ARKodeGetCurrentGamma_adapt_modifiable_immutable_to_return = [](void * arkode_mem, double gamma) -> std::tuple<int, double>
        {
            double * gamma_adapt_modifiable = & gamma;

            int r = ARKodeGetCurrentGamma(arkode_mem, gamma_adapt_modifiable);
            return std::make_tuple(r, gamma);
        };

        return ARKodeGetCurrentGamma_adapt_modifiable_immutable_to_return(arkode_mem, gamma);
    },     nb::arg("arkode_mem"), nb::arg("gamma"));

m.def("ARKodeGetNumNonlinSolvIters",
    [](void * arkode_mem, long nniters) -> std::tuple<int, long>
    {
        auto ARKodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nniters) -> std::tuple<int, long>
        {
            long * nniters_adapt_modifiable = & nniters;

            int r = ARKodeGetNumNonlinSolvIters(arkode_mem, nniters_adapt_modifiable);
            return std::make_tuple(r, nniters);
        };

        return ARKodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(arkode_mem, nniters);
    },     nb::arg("arkode_mem"), nb::arg("nniters"));

m.def("ARKodeGetNumNonlinSolvConvFails",
    [](void * arkode_mem, long nnfails) -> std::tuple<int, long>
    {
        auto ARKodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nnfails) -> std::tuple<int, long>
        {
            long * nnfails_adapt_modifiable = & nnfails;

            int r = ARKodeGetNumNonlinSolvConvFails(arkode_mem, nnfails_adapt_modifiable);
            return std::make_tuple(r, nnfails);
        };

        return ARKodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(arkode_mem, nnfails);
    },     nb::arg("arkode_mem"), nb::arg("nnfails"));

m.def("ARKodeGetNonlinSolvStats",
    [](void * arkode_mem, long nniters, long nnfails) -> std::tuple<int, long, long>
    {
        auto ARKodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nniters, long nnfails) -> std::tuple<int, long, long>
        {
            long * nniters_adapt_modifiable = & nniters;
            long * nnfails_adapt_modifiable = & nnfails;

            int r = ARKodeGetNonlinSolvStats(arkode_mem, nniters_adapt_modifiable, nnfails_adapt_modifiable);
            return std::make_tuple(r, nniters, nnfails);
        };

        return ARKodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return(arkode_mem, nniters, nnfails);
    },     nb::arg("arkode_mem"), nb::arg("nniters"), nb::arg("nnfails"));

m.def("ARKodeGetNumStepSolveFails",
    [](void * arkode_mem, long nncfails) -> std::tuple<int, long>
    {
        auto ARKodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nncfails) -> std::tuple<int, long>
        {
            long * nncfails_adapt_modifiable = & nncfails;

            int r = ARKodeGetNumStepSolveFails(arkode_mem, nncfails_adapt_modifiable);
            return std::make_tuple(r, nncfails);
        };

        return ARKodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return(arkode_mem, nncfails);
    },     nb::arg("arkode_mem"), nb::arg("nncfails"));

m.def("ARKodeGetNumJacEvals",
    [](void * arkode_mem, long njevals) -> std::tuple<int, long>
    {
        auto ARKodeGetNumJacEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long njevals) -> std::tuple<int, long>
        {
            long * njevals_adapt_modifiable = & njevals;

            int r = ARKodeGetNumJacEvals(arkode_mem, njevals_adapt_modifiable);
            return std::make_tuple(r, njevals);
        };

        return ARKodeGetNumJacEvals_adapt_modifiable_immutable_to_return(arkode_mem, njevals);
    },     nb::arg("arkode_mem"), nb::arg("njevals"));

m.def("ARKodeGetNumPrecEvals",
    [](void * arkode_mem, long npevals) -> std::tuple<int, long>
    {
        auto ARKodeGetNumPrecEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long npevals) -> std::tuple<int, long>
        {
            long * npevals_adapt_modifiable = & npevals;

            int r = ARKodeGetNumPrecEvals(arkode_mem, npevals_adapt_modifiable);
            return std::make_tuple(r, npevals);
        };

        return ARKodeGetNumPrecEvals_adapt_modifiable_immutable_to_return(arkode_mem, npevals);
    },     nb::arg("arkode_mem"), nb::arg("npevals"));

m.def("ARKodeGetNumPrecSolves",
    [](void * arkode_mem, long npsolves) -> std::tuple<int, long>
    {
        auto ARKodeGetNumPrecSolves_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long npsolves) -> std::tuple<int, long>
        {
            long * npsolves_adapt_modifiable = & npsolves;

            int r = ARKodeGetNumPrecSolves(arkode_mem, npsolves_adapt_modifiable);
            return std::make_tuple(r, npsolves);
        };

        return ARKodeGetNumPrecSolves_adapt_modifiable_immutable_to_return(arkode_mem, npsolves);
    },     nb::arg("arkode_mem"), nb::arg("npsolves"));

m.def("ARKodeGetNumLinIters",
    [](void * arkode_mem, long nliters) -> std::tuple<int, long>
    {
        auto ARKodeGetNumLinIters_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nliters) -> std::tuple<int, long>
        {
            long * nliters_adapt_modifiable = & nliters;

            int r = ARKodeGetNumLinIters(arkode_mem, nliters_adapt_modifiable);
            return std::make_tuple(r, nliters);
        };

        return ARKodeGetNumLinIters_adapt_modifiable_immutable_to_return(arkode_mem, nliters);
    },     nb::arg("arkode_mem"), nb::arg("nliters"));

m.def("ARKodeGetNumLinConvFails",
    [](void * arkode_mem, long nlcfails) -> std::tuple<int, long>
    {
        auto ARKodeGetNumLinConvFails_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nlcfails) -> std::tuple<int, long>
        {
            long * nlcfails_adapt_modifiable = & nlcfails;

            int r = ARKodeGetNumLinConvFails(arkode_mem, nlcfails_adapt_modifiable);
            return std::make_tuple(r, nlcfails);
        };

        return ARKodeGetNumLinConvFails_adapt_modifiable_immutable_to_return(arkode_mem, nlcfails);
    },     nb::arg("arkode_mem"), nb::arg("nlcfails"));

m.def("ARKodeGetNumJTSetupEvals",
    [](void * arkode_mem, long njtsetups) -> std::tuple<int, long>
    {
        auto ARKodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long njtsetups) -> std::tuple<int, long>
        {
            long * njtsetups_adapt_modifiable = & njtsetups;

            int r = ARKodeGetNumJTSetupEvals(arkode_mem, njtsetups_adapt_modifiable);
            return std::make_tuple(r, njtsetups);
        };

        return ARKodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return(arkode_mem, njtsetups);
    },     nb::arg("arkode_mem"), nb::arg("njtsetups"));

m.def("ARKodeGetNumJtimesEvals",
    [](void * arkode_mem, long njvevals) -> std::tuple<int, long>
    {
        auto ARKodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long njvevals) -> std::tuple<int, long>
        {
            long * njvevals_adapt_modifiable = & njvevals;

            int r = ARKodeGetNumJtimesEvals(arkode_mem, njvevals_adapt_modifiable);
            return std::make_tuple(r, njvevals);
        };

        return ARKodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return(arkode_mem, njvevals);
    },     nb::arg("arkode_mem"), nb::arg("njvevals"));

m.def("ARKodeGetNumLinRhsEvals",
    [](void * arkode_mem, long nfevalsLS) -> std::tuple<int, long>
    {
        auto ARKodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nfevalsLS) -> std::tuple<int, long>
        {
            long * nfevalsLS_adapt_modifiable = & nfevalsLS;

            int r = ARKodeGetNumLinRhsEvals(arkode_mem, nfevalsLS_adapt_modifiable);
            return std::make_tuple(r, nfevalsLS);
        };

        return ARKodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return(arkode_mem, nfevalsLS);
    },     nb::arg("arkode_mem"), nb::arg("nfevalsLS"));

m.def("ARKodeGetLastLinFlag",
    [](void * arkode_mem, long flag) -> std::tuple<int, long>
    {
        auto ARKodeGetLastLinFlag_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long flag) -> std::tuple<int, long>
        {
            long * flag_adapt_modifiable = & flag;

            int r = ARKodeGetLastLinFlag(arkode_mem, flag_adapt_modifiable);
            return std::make_tuple(r, flag);
        };

        return ARKodeGetLastLinFlag_adapt_modifiable_immutable_to_return(arkode_mem, flag);
    },     nb::arg("arkode_mem"), nb::arg("flag"));

m.def("ARKodeGetLinReturnFlagName",
    ARKodeGetLinReturnFlagName, nb::arg("flag"));

m.def("ARKodeGetResWeights",
    ARKodeGetResWeights, nb::arg("arkode_mem"), nb::arg("rweight"));

m.def("ARKodeGetNumMassSetups",
    [](void * arkode_mem, long nmsetups) -> std::tuple<int, long>
    {
        auto ARKodeGetNumMassSetups_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nmsetups) -> std::tuple<int, long>
        {
            long * nmsetups_adapt_modifiable = & nmsetups;

            int r = ARKodeGetNumMassSetups(arkode_mem, nmsetups_adapt_modifiable);
            return std::make_tuple(r, nmsetups);
        };

        return ARKodeGetNumMassSetups_adapt_modifiable_immutable_to_return(arkode_mem, nmsetups);
    },     nb::arg("arkode_mem"), nb::arg("nmsetups"));

m.def("ARKodeGetNumMassMultSetups",
    [](void * arkode_mem, long nmvsetups) -> std::tuple<int, long>
    {
        auto ARKodeGetNumMassMultSetups_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nmvsetups) -> std::tuple<int, long>
        {
            long * nmvsetups_adapt_modifiable = & nmvsetups;

            int r = ARKodeGetNumMassMultSetups(arkode_mem, nmvsetups_adapt_modifiable);
            return std::make_tuple(r, nmvsetups);
        };

        return ARKodeGetNumMassMultSetups_adapt_modifiable_immutable_to_return(arkode_mem, nmvsetups);
    },     nb::arg("arkode_mem"), nb::arg("nmvsetups"));

m.def("ARKodeGetNumMassMult",
    [](void * arkode_mem, long nmvevals) -> std::tuple<int, long>
    {
        auto ARKodeGetNumMassMult_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nmvevals) -> std::tuple<int, long>
        {
            long * nmvevals_adapt_modifiable = & nmvevals;

            int r = ARKodeGetNumMassMult(arkode_mem, nmvevals_adapt_modifiable);
            return std::make_tuple(r, nmvevals);
        };

        return ARKodeGetNumMassMult_adapt_modifiable_immutable_to_return(arkode_mem, nmvevals);
    },     nb::arg("arkode_mem"), nb::arg("nmvevals"));

m.def("ARKodeGetNumMassSolves",
    [](void * arkode_mem, long nmsolves) -> std::tuple<int, long>
    {
        auto ARKodeGetNumMassSolves_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nmsolves) -> std::tuple<int, long>
        {
            long * nmsolves_adapt_modifiable = & nmsolves;

            int r = ARKodeGetNumMassSolves(arkode_mem, nmsolves_adapt_modifiable);
            return std::make_tuple(r, nmsolves);
        };

        return ARKodeGetNumMassSolves_adapt_modifiable_immutable_to_return(arkode_mem, nmsolves);
    },     nb::arg("arkode_mem"), nb::arg("nmsolves"));

m.def("ARKodeGetNumMassPrecEvals",
    [](void * arkode_mem, long nmpevals) -> std::tuple<int, long>
    {
        auto ARKodeGetNumMassPrecEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nmpevals) -> std::tuple<int, long>
        {
            long * nmpevals_adapt_modifiable = & nmpevals;

            int r = ARKodeGetNumMassPrecEvals(arkode_mem, nmpevals_adapt_modifiable);
            return std::make_tuple(r, nmpevals);
        };

        return ARKodeGetNumMassPrecEvals_adapt_modifiable_immutable_to_return(arkode_mem, nmpevals);
    },     nb::arg("arkode_mem"), nb::arg("nmpevals"));

m.def("ARKodeGetNumMassPrecSolves",
    [](void * arkode_mem, long nmpsolves) -> std::tuple<int, long>
    {
        auto ARKodeGetNumMassPrecSolves_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nmpsolves) -> std::tuple<int, long>
        {
            long * nmpsolves_adapt_modifiable = & nmpsolves;

            int r = ARKodeGetNumMassPrecSolves(arkode_mem, nmpsolves_adapt_modifiable);
            return std::make_tuple(r, nmpsolves);
        };

        return ARKodeGetNumMassPrecSolves_adapt_modifiable_immutable_to_return(arkode_mem, nmpsolves);
    },     nb::arg("arkode_mem"), nb::arg("nmpsolves"));

m.def("ARKodeGetNumMassIters",
    [](void * arkode_mem, long nmiters) -> std::tuple<int, long>
    {
        auto ARKodeGetNumMassIters_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nmiters) -> std::tuple<int, long>
        {
            long * nmiters_adapt_modifiable = & nmiters;

            int r = ARKodeGetNumMassIters(arkode_mem, nmiters_adapt_modifiable);
            return std::make_tuple(r, nmiters);
        };

        return ARKodeGetNumMassIters_adapt_modifiable_immutable_to_return(arkode_mem, nmiters);
    },     nb::arg("arkode_mem"), nb::arg("nmiters"));

m.def("ARKodeGetNumMassConvFails",
    [](void * arkode_mem, long nmcfails) -> std::tuple<int, long>
    {
        auto ARKodeGetNumMassConvFails_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nmcfails) -> std::tuple<int, long>
        {
            long * nmcfails_adapt_modifiable = & nmcfails;

            int r = ARKodeGetNumMassConvFails(arkode_mem, nmcfails_adapt_modifiable);
            return std::make_tuple(r, nmcfails);
        };

        return ARKodeGetNumMassConvFails_adapt_modifiable_immutable_to_return(arkode_mem, nmcfails);
    },     nb::arg("arkode_mem"), nb::arg("nmcfails"));

m.def("ARKodeGetNumMTSetups",
    [](void * arkode_mem, long nmtsetups) -> std::tuple<int, long>
    {
        auto ARKodeGetNumMTSetups_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long nmtsetups) -> std::tuple<int, long>
        {
            long * nmtsetups_adapt_modifiable = & nmtsetups;

            int r = ARKodeGetNumMTSetups(arkode_mem, nmtsetups_adapt_modifiable);
            return std::make_tuple(r, nmtsetups);
        };

        return ARKodeGetNumMTSetups_adapt_modifiable_immutable_to_return(arkode_mem, nmtsetups);
    },     nb::arg("arkode_mem"), nb::arg("nmtsetups"));

m.def("ARKodeGetLastMassFlag",
    [](void * arkode_mem, long flag) -> std::tuple<int, long>
    {
        auto ARKodeGetLastMassFlag_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long flag) -> std::tuple<int, long>
        {
            long * flag_adapt_modifiable = & flag;

            int r = ARKodeGetLastMassFlag(arkode_mem, flag_adapt_modifiable);
            return std::make_tuple(r, flag);
        };

        return ARKodeGetLastMassFlag_adapt_modifiable_immutable_to_return(arkode_mem, flag);
    },     nb::arg("arkode_mem"), nb::arg("flag"));

m.def("ARKodePrintMem",
    ARKodePrintMem, nb::arg("arkode_mem"), nb::arg("outfile"));

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
    [](void * arkode_mem, long r_evals) -> std::tuple<int, long>
    {
        auto ARKodeGetNumRelaxFnEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long r_evals) -> std::tuple<int, long>
        {
            long * r_evals_adapt_modifiable = & r_evals;

            int r = ARKodeGetNumRelaxFnEvals(arkode_mem, r_evals_adapt_modifiable);
            return std::make_tuple(r, r_evals);
        };

        return ARKodeGetNumRelaxFnEvals_adapt_modifiable_immutable_to_return(arkode_mem, r_evals);
    },     nb::arg("arkode_mem"), nb::arg("r_evals"));

m.def("ARKodeGetNumRelaxJacEvals",
    [](void * arkode_mem, long J_evals) -> std::tuple<int, long>
    {
        auto ARKodeGetNumRelaxJacEvals_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long J_evals) -> std::tuple<int, long>
        {
            long * J_evals_adapt_modifiable = & J_evals;

            int r = ARKodeGetNumRelaxJacEvals(arkode_mem, J_evals_adapt_modifiable);
            return std::make_tuple(r, J_evals);
        };

        return ARKodeGetNumRelaxJacEvals_adapt_modifiable_immutable_to_return(arkode_mem, J_evals);
    },     nb::arg("arkode_mem"), nb::arg("J_evals"));

m.def("ARKodeGetNumRelaxFails",
    [](void * arkode_mem, long relax_fails) -> std::tuple<int, long>
    {
        auto ARKodeGetNumRelaxFails_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long relax_fails) -> std::tuple<int, long>
        {
            long * relax_fails_adapt_modifiable = & relax_fails;

            int r = ARKodeGetNumRelaxFails(arkode_mem, relax_fails_adapt_modifiable);
            return std::make_tuple(r, relax_fails);
        };

        return ARKodeGetNumRelaxFails_adapt_modifiable_immutable_to_return(arkode_mem, relax_fails);
    },     nb::arg("arkode_mem"), nb::arg("relax_fails"));

m.def("ARKodeGetNumRelaxBoundFails",
    [](void * arkode_mem, long fails) -> std::tuple<int, long>
    {
        auto ARKodeGetNumRelaxBoundFails_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long fails) -> std::tuple<int, long>
        {
            long * fails_adapt_modifiable = & fails;

            int r = ARKodeGetNumRelaxBoundFails(arkode_mem, fails_adapt_modifiable);
            return std::make_tuple(r, fails);
        };

        return ARKodeGetNumRelaxBoundFails_adapt_modifiable_immutable_to_return(arkode_mem, fails);
    },     nb::arg("arkode_mem"), nb::arg("fails"));

m.def("ARKodeGetNumRelaxSolveFails",
    [](void * arkode_mem, long fails) -> std::tuple<int, long>
    {
        auto ARKodeGetNumRelaxSolveFails_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long fails) -> std::tuple<int, long>
        {
            long * fails_adapt_modifiable = & fails;

            int r = ARKodeGetNumRelaxSolveFails(arkode_mem, fails_adapt_modifiable);
            return std::make_tuple(r, fails);
        };

        return ARKodeGetNumRelaxSolveFails_adapt_modifiable_immutable_to_return(arkode_mem, fails);
    },     nb::arg("arkode_mem"), nb::arg("fails"));

m.def("ARKodeGetNumRelaxSolveIters",
    [](void * arkode_mem, long iters) -> std::tuple<int, long>
    {
        auto ARKodeGetNumRelaxSolveIters_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long iters) -> std::tuple<int, long>
        {
            long * iters_adapt_modifiable = & iters;

            int r = ARKodeGetNumRelaxSolveIters(arkode_mem, iters_adapt_modifiable);
            return std::make_tuple(r, iters);
        };

        return ARKodeGetNumRelaxSolveIters_adapt_modifiable_immutable_to_return(arkode_mem, iters);
    },     nb::arg("arkode_mem"), nb::arg("iters"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
// #ifndef _ARKLS_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("ARKodeSetLinearSolver",
    [](void * arkode_mem, SUNLinearSolver LS, std::optional<SUNMatrix> A = std::nullopt) -> int
    {
        auto ARKodeSetLinearSolver_adapt_optional_arg_with_default_null = [](void * arkode_mem, SUNLinearSolver LS, std::optional<SUNMatrix> A = std::nullopt) -> int
        {
            SUNMatrix A_adapt_default_null = nullptr;
            if (A.has_value())
                A_adapt_default_null = A.value();

            auto lambda_result = ARKodeSetLinearSolver(arkode_mem, LS, A_adapt_default_null);
            return lambda_result;
        };

        return ARKodeSetLinearSolver_adapt_optional_arg_with_default_null(arkode_mem, LS, A);
    },     nb::arg("arkode_mem"), nb::arg("LS"), nb::arg("A") = nb::none());

m.def("ARKodeSetJacEvalFrequency",
    ARKodeSetJacEvalFrequency, nb::arg("arkode_mem"), nb::arg("msbj"));

m.def("ARKodeSetLinearSolutionScaling",
    ARKodeSetLinearSolutionScaling, nb::arg("arkode_mem"), nb::arg("onoff"));

m.def("ARKodeSetEpsLin",
    ARKodeSetEpsLin, nb::arg("arkode_mem"), nb::arg("eplifac"));

m.def("ARKodeSetMassEpsLin",
    ARKodeSetMassEpsLin, nb::arg("arkode_mem"), nb::arg("eplifac"));

m.def("ARKodeSetLSNormFactor",
    ARKodeSetLSNormFactor, nb::arg("arkode_mem"), nb::arg("nrmfac"));

m.def("ARKodeSetMassLSNormFactor",
    ARKodeSetMassLSNormFactor, nb::arg("arkode_mem"), nb::arg("nrmfac"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
// #ifndef _ARKODE_BUTCHER_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyClassARKodeButcherTableMem =
    nb::class_<ARKodeButcherTableMem>
        (m, "ARKodeButcherTableMem", "")
    .def(nb::init<>()) // implicit default constructor
    ;


m.def("ARKodeButcherTable_Copy",
    ARKodeButcherTable_Copy, nb::arg("B"));

m.def("ARKodeButcherTable_Write",
    ARKodeButcherTable_Write, nb::arg("B"), nb::arg("outfile"));

m.def("ARKodeButcherTable_IsStifflyAccurate",
    ARKodeButcherTable_IsStifflyAccurate, nb::arg("B"));

m.def("ARKodeButcherTable_CheckOrder",
    [](ARKodeButcherTable B, int q, int p, FILE * outfile) -> std::tuple<int, int, int>
    {
        auto ARKodeButcherTable_CheckOrder_adapt_modifiable_immutable_to_return = [](ARKodeButcherTable B, int q, int p, FILE * outfile) -> std::tuple<int, int, int>
        {
            int * q_adapt_modifiable = & q;
            int * p_adapt_modifiable = & p;

            int r = ARKodeButcherTable_CheckOrder(B, q_adapt_modifiable, p_adapt_modifiable, outfile);
            return std::make_tuple(r, q, p);
        };

        return ARKodeButcherTable_CheckOrder_adapt_modifiable_immutable_to_return(B, q, p, outfile);
    },     nb::arg("B"), nb::arg("q"), nb::arg("p"), nb::arg("outfile"));

m.def("ARKodeButcherTable_CheckARKOrder",
    [](ARKodeButcherTable B1, ARKodeButcherTable B2, int q, int p, FILE * outfile) -> std::tuple<int, int, int>
    {
        auto ARKodeButcherTable_CheckARKOrder_adapt_modifiable_immutable_to_return = [](ARKodeButcherTable B1, ARKodeButcherTable B2, int q, int p, FILE * outfile) -> std::tuple<int, int, int>
        {
            int * q_adapt_modifiable = & q;
            int * p_adapt_modifiable = & p;

            int r = ARKodeButcherTable_CheckARKOrder(B1, B2, q_adapt_modifiable, p_adapt_modifiable, outfile);
            return std::make_tuple(r, q, p);
        };

        return ARKodeButcherTable_CheckARKOrder_adapt_modifiable_immutable_to_return(B1, B2, q, p, outfile);
    },     nb::arg("B1"), nb::arg("B2"), nb::arg("q"), nb::arg("p"), nb::arg("outfile"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
