// #ifndef _ARKODE_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("ARK_NORMAL")                 = 1;
m.attr("ARK_ONE_STEP")               = 2;
m.attr("ARK_ADAPT_CUSTOM")           = -1;
m.attr("ARK_ADAPT_PID")              = 0;
m.attr("ARK_ADAPT_PI")               = 1;
m.attr("ARK_ADAPT_I")                = 2;
m.attr("ARK_ADAPT_EXP_GUS")          = 3;
m.attr("ARK_ADAPT_IMP_GUS")          = 4;
m.attr("ARK_ADAPT_IMEX_GUS")         = 5;
m.attr("ARK_FULLRHS_START")          = 0;
m.attr("ARK_FULLRHS_END")            = 1;
m.attr("ARK_FULLRHS_OTHER")          = 2;
m.attr("ARK_INTERP_MAX_DEGREE")      = 5;
m.attr("ARK_INTERP_NONE")            = -1;
m.attr("ARK_INTERP_HERMITE")         = 0;
m.attr("ARK_INTERP_LAGRANGE")        = 1;
m.attr("ARK_SUCCESS")                = 0;
m.attr("ARK_TSTOP_RETURN")           = 1;
m.attr("ARK_ROOT_RETURN")            = 2;
m.attr("ARK_WARNING")                = 99;
m.attr("ARK_TOO_MUCH_WORK")          = -1;
m.attr("ARK_TOO_MUCH_ACC")           = -2;
m.attr("ARK_ERR_FAILURE")            = -3;
m.attr("ARK_CONV_FAILURE")           = -4;
m.attr("ARK_LINIT_FAIL")             = -5;
m.attr("ARK_LSETUP_FAIL")            = -6;
m.attr("ARK_LSOLVE_FAIL")            = -7;
m.attr("ARK_RHSFUNC_FAIL")           = -8;
m.attr("ARK_FIRST_RHSFUNC_ERR")      = -9;
m.attr("ARK_REPTD_RHSFUNC_ERR")      = -10;
m.attr("ARK_UNREC_RHSFUNC_ERR")      = -11;
m.attr("ARK_RTFUNC_FAIL")            = -12;
m.attr("ARK_LFREE_FAIL")             = -13;
m.attr("ARK_MASSINIT_FAIL")          = -14;
m.attr("ARK_MASSSETUP_FAIL")         = -15;
m.attr("ARK_MASSSOLVE_FAIL")         = -16;
m.attr("ARK_MASSFREE_FAIL")          = -17;
m.attr("ARK_MASSMULT_FAIL")          = -18;
m.attr("ARK_CONSTR_FAIL")            = -19;
m.attr("ARK_MEM_FAIL")               = -20;
m.attr("ARK_MEM_NULL")               = -21;
m.attr("ARK_ILL_INPUT")              = -22;
m.attr("ARK_NO_MALLOC")              = -23;
m.attr("ARK_BAD_K")                  = -24;
m.attr("ARK_BAD_T")                  = -25;
m.attr("ARK_BAD_DKY")                = -26;
m.attr("ARK_TOO_CLOSE")              = -27;
m.attr("ARK_VECTOROP_ERR")           = -28;
m.attr("ARK_NLS_INIT_FAIL")          = -29;
m.attr("ARK_NLS_SETUP_FAIL")         = -30;
m.attr("ARK_NLS_SETUP_RECVR")        = -31;
m.attr("ARK_NLS_OP_ERR")             = -32;
m.attr("ARK_INNERSTEP_ATTACH_ERR")   = -33;
m.attr("ARK_INNERSTEP_FAIL")         = -34;
m.attr("ARK_OUTERTOINNER_FAIL")      = -35;
m.attr("ARK_INNERTOOUTER_FAIL")      = -36;
m.attr("ARK_POSTPROCESS_FAIL")       = -37;
m.attr("ARK_POSTPROCESS_STEP_FAIL")  = -37;
m.attr("ARK_POSTPROCESS_STAGE_FAIL") = -38;
m.attr("ARK_PRESTEPFN_FAIL")         = -39;
m.attr("ARK_POSTSTEPFN_FAIL")        = -40;
m.attr("ARK_PRERHSFN_FAIL")          = -41;
m.attr("ARK_USER_PREDICT_FAIL")      = -42;
m.attr("ARK_INTERP_FAIL")            = -43;
m.attr("ARK_INVALID_TABLE")          = -44;
m.attr("ARK_CONTEXT_ERR")            = -45;
m.attr("ARK_RELAX_FAIL")             = -46;
m.attr("ARK_RELAX_MEM_NULL")         = -47;
m.attr("ARK_RELAX_FUNC_FAIL")        = -48;
m.attr("ARK_RELAX_JAC_FAIL")         = -49;
m.attr("ARK_CONTROLLER_ERR")         = -50;
m.attr("ARK_STEPPER_UNSUPPORTED")    = -51;
m.attr("ARK_DOMEIG_FAIL")            = -52;
m.attr("ARK_MAX_STAGE_LIMIT_FAIL")   = -53;
m.attr("ARK_SUNSTEPPER_ERR")         = -54;
m.attr("ARK_STEP_DIRECTION_ERR")     = -55;
m.attr("ARK_ADJ_CHECKPOINT_FAIL")    = -56;
m.attr("ARK_ADJ_RECOMPUTE_FAIL")     = -57;
m.attr("ARK_SUNADJSTEPPER_ERR")      = -58;
m.attr("ARK_DEE_FAIL")               = -59;
m.attr("ARK_STEP_H0_FAIL")           = -60;
m.attr("ARK_UNRECOGNIZED_ERROR")     = -99;

auto pyEnumARKRelaxSolver = nb::enum_<ARKRelaxSolver>(m, "ARKRelaxSolver",
                                                      nb::is_arithmetic(), "")
                              .value("ARK_RELAX_BRENT", ARK_RELAX_BRENT, "")
                              .value("ARK_RELAX_NEWTON", ARK_RELAX_NEWTON, "")
                              .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyEnumARKAccumError =
  nb::enum_<ARKAccumError>(m, "ARKAccumError", nb::is_arithmetic(), "")
    .value("ARK_ACCUMERROR_NONE", ARK_ACCUMERROR_NONE, "")
    .value("ARK_ACCUMERROR_MAX", ARK_ACCUMERROR_MAX, "")
    .value("ARK_ACCUMERROR_SUM", ARK_ACCUMERROR_SUM, "")
    .value("ARK_ACCUMERROR_AVG", ARK_ACCUMERROR_AVG, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

sundials4py::scoped_def(
  m, "ARKodeReset",
  [](const ARKodeView& arkode_mem, sunrealtype tR, N_Vector yR) -> int
  {
    auto ARKodeReset_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                              sunrealtype tR, N_Vector yR) -> int
    {
      auto lambda_result = ARKodeReset(arkode_mem.get(), tR, yR);
      return lambda_result;
    };

    return ARKodeReset_adapt_solver_memory(arkode_mem, tR, yR);
  },
  nb::arg("arkode_mem"), nb::arg("tR"), nb::arg("yR"));

sundials4py::scoped_def(
  m, "ARKodeInit",
  [](const ARKodeView& arkode_mem) -> int
  {
    auto ARKodeInit_adapt_solver_memory = [](const ARKodeView& arkode_mem) -> int
    {
      auto lambda_result = ARKodeInit(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeInit_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"), "Optional data allocation function");

sundials4py::scoped_def(
  m, "ARKodeCreateMRIStepInnerStepper",
  [](const ARKodeView& arkode_mem)
    -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<MRIStepInnerStepper>>>
  {
    auto ARKodeCreateMRIStepInnerStepper_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, MRIStepInnerStepper>
    {
      MRIStepInnerStepper stepper_adapt_modifiable;

      int r = ARKodeCreateMRIStepInnerStepper(arkode_mem,
                                              &stepper_adapt_modifiable);
      return std::make_tuple(r, stepper_adapt_modifiable);
    };
    auto ARKodeCreateMRIStepInnerStepper_adapt_return_type_to_shared_ptr =
      [&ARKodeCreateMRIStepInnerStepper_adapt_modifiable_immutable_to_return](
        void* arkode_mem)
      -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<MRIStepInnerStepper>>>
    {
      auto lambda_result =
        ARKodeCreateMRIStepInnerStepper_adapt_modifiable_immutable_to_return(
          arkode_mem);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<MRIStepInnerStepper>,
                                             MRIStepInnerStepperDeleter>(
                               std::get<1>(lambda_result)));
    };
    auto ARKodeCreateMRIStepInnerStepper_adapt_solver_memory =
      [&ARKodeCreateMRIStepInnerStepper_adapt_return_type_to_shared_ptr](
        const ARKodeView& arkode_mem)
      -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<MRIStepInnerStepper>>>
    {
      auto lambda_result =
        ARKodeCreateMRIStepInnerStepper_adapt_return_type_to_shared_ptr(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeCreateMRIStepInnerStepper_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"), "Utility to wrap ARKODE as an MRIStepInnerStepper",
  nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "ARKodeSStolerances",
  [](const ARKodeView& arkode_mem, sunrealtype reltol, sunrealtype abstol) -> int
  {
    auto ARKodeSStolerances_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype reltol,
         sunrealtype abstol) -> int
    {
      auto lambda_result = ARKodeSStolerances(arkode_mem.get(), reltol, abstol);
      return lambda_result;
    };

    return ARKodeSStolerances_adapt_solver_memory(arkode_mem, reltol, abstol);
  },
  nb::arg("arkode_mem"), nb::arg("reltol"), nb::arg("abstol"));

sundials4py::scoped_def(
  m, "ARKodeSVtolerances",
  [](const ARKodeView& arkode_mem, sunrealtype reltol, N_Vector abstol) -> int
  {
    auto ARKodeSVtolerances_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype reltol, N_Vector abstol) -> int
    {
      auto lambda_result = ARKodeSVtolerances(arkode_mem.get(), reltol, abstol);
      return lambda_result;
    };

    return ARKodeSVtolerances_adapt_solver_memory(arkode_mem, reltol, abstol);
  },
  nb::arg("arkode_mem"), nb::arg("reltol"), nb::arg("abstol"));

sundials4py::scoped_def(
  m, "ARKodeResStolerance",
  [](const ARKodeView& arkode_mem, sunrealtype rabstol) -> int
  {
    auto ARKodeResStolerance_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype rabstol) -> int
    {
      auto lambda_result = ARKodeResStolerance(arkode_mem.get(), rabstol);
      return lambda_result;
    };

    return ARKodeResStolerance_adapt_solver_memory(arkode_mem, rabstol);
  },
  nb::arg("arkode_mem"), nb::arg("rabstol"));

sundials4py::scoped_def(
  m, "ARKodeResVtolerance",
  [](const ARKodeView& arkode_mem, N_Vector rabstol) -> int
  {
    auto ARKodeResVtolerance_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, N_Vector rabstol) -> int
    {
      auto lambda_result = ARKodeResVtolerance(arkode_mem.get(), rabstol);
      return lambda_result;
    };

    return ARKodeResVtolerance_adapt_solver_memory(arkode_mem, rabstol);
  },
  nb::arg("arkode_mem"), nb::arg("rabstol"));

sundials4py::scoped_def(
  m, "ARKodeSetRootDirection",
  [](const ARKodeView& arkode_mem, std::vector<int> rootdir_1d) -> int
  {
    auto ARKodeSetRootDirection_adapt_arr_ptr_to_std_vector =
      [](void* arkode_mem, std::vector<int> rootdir_1d) -> int
    {
      int* rootdir_1d_ptr = rootdir_1d.empty() ? nullptr : rootdir_1d.data();

      auto lambda_result = ARKodeSetRootDirection(arkode_mem, rootdir_1d_ptr);
      return lambda_result;
    };
    auto ARKodeSetRootDirection_adapt_solver_memory =
      [&ARKodeSetRootDirection_adapt_arr_ptr_to_std_vector](const ARKodeView& arkode_mem,
                                                            std::vector<int>
                                                              rootdir_1d) -> int
    {
      auto lambda_result =
        ARKodeSetRootDirection_adapt_arr_ptr_to_std_vector(arkode_mem.get(),
                                                           rootdir_1d);
      return lambda_result;
    };

    return ARKodeSetRootDirection_adapt_solver_memory(arkode_mem, rootdir_1d);
  },
  nb::arg("arkode_mem"), nb::arg("rootdir_1d"));

sundials4py::scoped_def(
  m, "ARKodeSetNoInactiveRootWarn",
  [](const ARKodeView& arkode_mem) -> int
  {
    auto ARKodeSetNoInactiveRootWarn_adapt_solver_memory =
      [](const ARKodeView& arkode_mem) -> int
    {
      auto lambda_result = ARKodeSetNoInactiveRootWarn(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeSetNoInactiveRootWarn_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeSetDefaults",
  [](const ARKodeView& arkode_mem) -> int
  {
    auto ARKodeSetDefaults_adapt_solver_memory =
      [](const ARKodeView& arkode_mem) -> int
    {
      auto lambda_result = ARKodeSetDefaults(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeSetDefaults_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeSetOrder",
  [](const ARKodeView& arkode_mem, int maxord) -> int
  {
    auto ARKodeSetOrder_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                                 int maxord) -> int
    {
      auto lambda_result = ARKodeSetOrder(arkode_mem.get(), maxord);
      return lambda_result;
    };

    return ARKodeSetOrder_adapt_solver_memory(arkode_mem, maxord);
  },
  nb::arg("arkode_mem"), nb::arg("maxord"));

sundials4py::scoped_def(
  m, "ARKodeSetInterpolantType",
  [](const ARKodeView& arkode_mem, int itype) -> int
  {
    auto ARKodeSetInterpolantType_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int itype) -> int
    {
      auto lambda_result = ARKodeSetInterpolantType(arkode_mem.get(), itype);
      return lambda_result;
    };

    return ARKodeSetInterpolantType_adapt_solver_memory(arkode_mem, itype);
  },
  nb::arg("arkode_mem"), nb::arg("itype"));

sundials4py::scoped_def(
  m, "ARKodeSetInterpolantDegree",
  [](const ARKodeView& arkode_mem, int degree) -> int
  {
    auto ARKodeSetInterpolantDegree_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int degree) -> int
    {
      auto lambda_result = ARKodeSetInterpolantDegree(arkode_mem.get(), degree);
      return lambda_result;
    };

    return ARKodeSetInterpolantDegree_adapt_solver_memory(arkode_mem, degree);
  },
  nb::arg("arkode_mem"), nb::arg("degree"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxNumSteps",
  [](const ARKodeView& arkode_mem, long mxsteps) -> int
  {
    auto ARKodeSetMaxNumSteps_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, long mxsteps) -> int
    {
      auto lambda_result = ARKodeSetMaxNumSteps(arkode_mem.get(), mxsteps);
      return lambda_result;
    };

    return ARKodeSetMaxNumSteps_adapt_solver_memory(arkode_mem, mxsteps);
  },
  nb::arg("arkode_mem"), nb::arg("mxsteps"));

sundials4py::scoped_def(
  m, "ARKodeSetInterpolateStopTime",
  [](const ARKodeView& arkode_mem, sunbooleantype interp) -> int
  {
    auto ARKodeSetInterpolateStopTime_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunbooleantype interp) -> int
    {
      auto lambda_result = ARKodeSetInterpolateStopTime(arkode_mem.get(), interp);
      return lambda_result;
    };

    return ARKodeSetInterpolateStopTime_adapt_solver_memory(arkode_mem, interp);
  },
  nb::arg("arkode_mem"), nb::arg("interp"));

sundials4py::scoped_def(
  m, "ARKodeSetStopTime",
  [](const ARKodeView& arkode_mem, sunrealtype tstop) -> int
  {
    auto ARKodeSetStopTime_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                                    sunrealtype tstop) -> int
    {
      auto lambda_result = ARKodeSetStopTime(arkode_mem.get(), tstop);
      return lambda_result;
    };

    return ARKodeSetStopTime_adapt_solver_memory(arkode_mem, tstop);
  },
  nb::arg("arkode_mem"), nb::arg("tstop"));

sundials4py::scoped_def(
  m, "ARKodeClearStopTime",
  [](const ARKodeView& arkode_mem) -> int
  {
    auto ARKodeClearStopTime_adapt_solver_memory =
      [](const ARKodeView& arkode_mem) -> int
    {
      auto lambda_result = ARKodeClearStopTime(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeClearStopTime_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeSetFixedStep",
  [](const ARKodeView& arkode_mem, sunrealtype hfixed) -> int
  {
    auto ARKodeSetFixedStep_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype hfixed) -> int
    {
      auto lambda_result = ARKodeSetFixedStep(arkode_mem.get(), hfixed);
      return lambda_result;
    };

    return ARKodeSetFixedStep_adapt_solver_memory(arkode_mem, hfixed);
  },
  nb::arg("arkode_mem"), nb::arg("hfixed"));

sundials4py::scoped_def(
  m, "ARKodeSetStepDirection",
  [](const ARKodeView& arkode_mem, sunrealtype stepdir) -> int
  {
    auto ARKodeSetStepDirection_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype stepdir) -> int
    {
      auto lambda_result = ARKodeSetStepDirection(arkode_mem.get(), stepdir);
      return lambda_result;
    };

    return ARKodeSetStepDirection_adapt_solver_memory(arkode_mem, stepdir);
  },
  nb::arg("arkode_mem"), nb::arg("stepdir"));

sundials4py::scoped_def(
  m, "ARKodeSetNonlinearSolver",
  [](const ARKodeView& arkode_mem, SUNNonlinearSolver NLS) -> int
  {
    auto ARKodeSetNonlinearSolver_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, SUNNonlinearSolver NLS) -> int
    {
      auto lambda_result = ARKodeSetNonlinearSolver(arkode_mem.get(), NLS);
      return lambda_result;
    };

    return ARKodeSetNonlinearSolver_adapt_solver_memory(arkode_mem, NLS);
  },
  nb::arg("arkode_mem"), nb::arg("NLS"));

sundials4py::scoped_def(
  m, "ARKodeSetLinear",
  [](const ARKodeView& arkode_mem, int timedepend) -> int
  {
    auto ARKodeSetLinear_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                                  int timedepend) -> int
    {
      auto lambda_result = ARKodeSetLinear(arkode_mem.get(), timedepend);
      return lambda_result;
    };

    return ARKodeSetLinear_adapt_solver_memory(arkode_mem, timedepend);
  },
  nb::arg("arkode_mem"), nb::arg("timedepend"));

sundials4py::scoped_def(
  m, "ARKodeSetNonlinear",
  [](const ARKodeView& arkode_mem) -> int
  {
    auto ARKodeSetNonlinear_adapt_solver_memory =
      [](const ARKodeView& arkode_mem) -> int
    {
      auto lambda_result = ARKodeSetNonlinear(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeSetNonlinear_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeSetAutonomous",
  [](const ARKodeView& arkode_mem, sunbooleantype autonomous) -> int
  {
    auto ARKodeSetAutonomous_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunbooleantype autonomous) -> int
    {
      auto lambda_result = ARKodeSetAutonomous(arkode_mem.get(), autonomous);
      return lambda_result;
    };

    return ARKodeSetAutonomous_adapt_solver_memory(arkode_mem, autonomous);
  },
  nb::arg("arkode_mem"), nb::arg("autonomous"));

sundials4py::scoped_def(
  m, "ARKodeSetDeduceImplicitRhs",
  [](const ARKodeView& arkode_mem, sunbooleantype deduce) -> int
  {
    auto ARKodeSetDeduceImplicitRhs_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunbooleantype deduce) -> int
    {
      auto lambda_result = ARKodeSetDeduceImplicitRhs(arkode_mem.get(), deduce);
      return lambda_result;
    };

    return ARKodeSetDeduceImplicitRhs_adapt_solver_memory(arkode_mem, deduce);
  },
  nb::arg("arkode_mem"), nb::arg("deduce"));

sundials4py::scoped_def(
  m, "ARKodeSetNonlinCRDown",
  [](const ARKodeView& arkode_mem, sunrealtype crdown) -> int
  {
    auto ARKodeSetNonlinCRDown_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype crdown) -> int
    {
      auto lambda_result = ARKodeSetNonlinCRDown(arkode_mem.get(), crdown);
      return lambda_result;
    };

    return ARKodeSetNonlinCRDown_adapt_solver_memory(arkode_mem, crdown);
  },
  nb::arg("arkode_mem"), nb::arg("crdown"));

sundials4py::scoped_def(
  m, "ARKodeSetNonlinRDiv",
  [](const ARKodeView& arkode_mem, sunrealtype rdiv) -> int
  {
    auto ARKodeSetNonlinRDiv_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype rdiv) -> int
    {
      auto lambda_result = ARKodeSetNonlinRDiv(arkode_mem.get(), rdiv);
      return lambda_result;
    };

    return ARKodeSetNonlinRDiv_adapt_solver_memory(arkode_mem, rdiv);
  },
  nb::arg("arkode_mem"), nb::arg("rdiv"));

sundials4py::scoped_def(
  m, "ARKodeSetDeltaGammaMax",
  [](const ARKodeView& arkode_mem, sunrealtype dgmax) -> int
  {
    auto ARKodeSetDeltaGammaMax_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype dgmax) -> int
    {
      auto lambda_result = ARKodeSetDeltaGammaMax(arkode_mem.get(), dgmax);
      return lambda_result;
    };

    return ARKodeSetDeltaGammaMax_adapt_solver_memory(arkode_mem, dgmax);
  },
  nb::arg("arkode_mem"), nb::arg("dgmax"));

sundials4py::scoped_def(
  m, "ARKodeSetLSetupFrequency",
  [](const ARKodeView& arkode_mem, int msbp) -> int
  {
    auto ARKodeSetLSetupFrequency_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int msbp) -> int
    {
      auto lambda_result = ARKodeSetLSetupFrequency(arkode_mem.get(), msbp);
      return lambda_result;
    };

    return ARKodeSetLSetupFrequency_adapt_solver_memory(arkode_mem, msbp);
  },
  nb::arg("arkode_mem"), nb::arg("msbp"));

sundials4py::scoped_def(
  m, "ARKodeSetPredictorMethod",
  [](const ARKodeView& arkode_mem, int method) -> int
  {
    auto ARKodeSetPredictorMethod_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int method) -> int
    {
      auto lambda_result = ARKodeSetPredictorMethod(arkode_mem.get(), method);
      return lambda_result;
    };

    return ARKodeSetPredictorMethod_adapt_solver_memory(arkode_mem, method);
  },
  nb::arg("arkode_mem"), nb::arg("method"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxNonlinIters",
  [](const ARKodeView& arkode_mem, int maxcor) -> int
  {
    auto ARKodeSetMaxNonlinIters_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int maxcor) -> int
    {
      auto lambda_result = ARKodeSetMaxNonlinIters(arkode_mem.get(), maxcor);
      return lambda_result;
    };

    return ARKodeSetMaxNonlinIters_adapt_solver_memory(arkode_mem, maxcor);
  },
  nb::arg("arkode_mem"), nb::arg("maxcor"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxConvFails",
  [](const ARKodeView& arkode_mem, int maxncf) -> int
  {
    auto ARKodeSetMaxConvFails_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int maxncf) -> int
    {
      auto lambda_result = ARKodeSetMaxConvFails(arkode_mem.get(), maxncf);
      return lambda_result;
    };

    return ARKodeSetMaxConvFails_adapt_solver_memory(arkode_mem, maxncf);
  },
  nb::arg("arkode_mem"), nb::arg("maxncf"));

sundials4py::scoped_def(
  m, "ARKodeSetNonlinConvCoef",
  [](const ARKodeView& arkode_mem, sunrealtype nlscoef) -> int
  {
    auto ARKodeSetNonlinConvCoef_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype nlscoef) -> int
    {
      auto lambda_result = ARKodeSetNonlinConvCoef(arkode_mem.get(), nlscoef);
      return lambda_result;
    };

    return ARKodeSetNonlinConvCoef_adapt_solver_memory(arkode_mem, nlscoef);
  },
  nb::arg("arkode_mem"), nb::arg("nlscoef"));

sundials4py::scoped_def(
  m, "ARKodeSetAdaptController",
  [](const ARKodeView& arkode_mem, SUNAdaptController C) -> int
  {
    auto ARKodeSetAdaptController_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, SUNAdaptController C) -> int
    {
      auto lambda_result = ARKodeSetAdaptController(arkode_mem.get(), C);
      return lambda_result;
    };

    return ARKodeSetAdaptController_adapt_solver_memory(arkode_mem, C);
  },
  nb::arg("arkode_mem"), nb::arg("C"));

sundials4py::scoped_def(
  m, "ARKodeSetAdaptControllerByName",
  [](const ARKodeView& arkode_mem, const char* cname) -> int
  {
    auto ARKodeSetAdaptControllerByName_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, const char* cname) -> int
    {
      auto lambda_result = ARKodeSetAdaptControllerByName(arkode_mem.get(),
                                                          cname);
      return lambda_result;
    };

    return ARKodeSetAdaptControllerByName_adapt_solver_memory(arkode_mem, cname);
  },
  nb::arg("arkode_mem"), nb::arg("cname"));

sundials4py::scoped_def(
  m, "ARKodeSetAdaptivityAdjustment",
  [](const ARKodeView& arkode_mem, int adjust) -> int
  {
    auto ARKodeSetAdaptivityAdjustment_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int adjust) -> int
    {
      auto lambda_result = ARKodeSetAdaptivityAdjustment(arkode_mem.get(),
                                                         adjust);
      return lambda_result;
    };

    return ARKodeSetAdaptivityAdjustment_adapt_solver_memory(arkode_mem, adjust);
  },
  nb::arg("arkode_mem"), nb::arg("adjust"));

sundials4py::scoped_def(
  m, "ARKodeSetCFLFraction",
  [](const ARKodeView& arkode_mem, sunrealtype cfl_frac) -> int
  {
    auto ARKodeSetCFLFraction_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype cfl_frac) -> int
    {
      auto lambda_result = ARKodeSetCFLFraction(arkode_mem.get(), cfl_frac);
      return lambda_result;
    };

    return ARKodeSetCFLFraction_adapt_solver_memory(arkode_mem, cfl_frac);
  },
  nb::arg("arkode_mem"), nb::arg("cfl_frac"));

sundials4py::scoped_def(
  m, "ARKodeSetErrorBias",
  [](const ARKodeView& arkode_mem, sunrealtype bias) -> int
  {
    auto ARKodeSetErrorBias_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype bias) -> int
    {
      auto lambda_result = ARKodeSetErrorBias(arkode_mem.get(), bias);
      return lambda_result;
    };

    return ARKodeSetErrorBias_adapt_solver_memory(arkode_mem, bias);
  },
  nb::arg("arkode_mem"), nb::arg("bias"));

sundials4py::scoped_def(
  m, "ARKodeSetSafetyFactor",
  [](const ARKodeView& arkode_mem, sunrealtype safety) -> int
  {
    auto ARKodeSetSafetyFactor_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype safety) -> int
    {
      auto lambda_result = ARKodeSetSafetyFactor(arkode_mem.get(), safety);
      return lambda_result;
    };

    return ARKodeSetSafetyFactor_adapt_solver_memory(arkode_mem, safety);
  },
  nb::arg("arkode_mem"), nb::arg("safety"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxGrowth",
  [](const ARKodeView& arkode_mem, sunrealtype mx_growth) -> int
  {
    auto ARKodeSetMaxGrowth_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype mx_growth) -> int
    {
      auto lambda_result = ARKodeSetMaxGrowth(arkode_mem.get(), mx_growth);
      return lambda_result;
    };

    return ARKodeSetMaxGrowth_adapt_solver_memory(arkode_mem, mx_growth);
  },
  nb::arg("arkode_mem"), nb::arg("mx_growth"));

sundials4py::scoped_def(
  m, "ARKodeSetMinReduction",
  [](const ARKodeView& arkode_mem, sunrealtype eta_min) -> int
  {
    auto ARKodeSetMinReduction_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype eta_min) -> int
    {
      auto lambda_result = ARKodeSetMinReduction(arkode_mem.get(), eta_min);
      return lambda_result;
    };

    return ARKodeSetMinReduction_adapt_solver_memory(arkode_mem, eta_min);
  },
  nb::arg("arkode_mem"), nb::arg("eta_min"));

sundials4py::scoped_def(
  m, "ARKodeSetFixedStepBounds",
  [](const ARKodeView& arkode_mem, sunrealtype lb, sunrealtype ub) -> int
  {
    auto ARKodeSetFixedStepBounds_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype lb, sunrealtype ub) -> int
    {
      auto lambda_result = ARKodeSetFixedStepBounds(arkode_mem.get(), lb, ub);
      return lambda_result;
    };

    return ARKodeSetFixedStepBounds_adapt_solver_memory(arkode_mem, lb, ub);
  },
  nb::arg("arkode_mem"), nb::arg("lb"), nb::arg("ub"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxFirstGrowth",
  [](const ARKodeView& arkode_mem, sunrealtype etamx1) -> int
  {
    auto ARKodeSetMaxFirstGrowth_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype etamx1) -> int
    {
      auto lambda_result = ARKodeSetMaxFirstGrowth(arkode_mem.get(), etamx1);
      return lambda_result;
    };

    return ARKodeSetMaxFirstGrowth_adapt_solver_memory(arkode_mem, etamx1);
  },
  nb::arg("arkode_mem"), nb::arg("etamx1"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxEFailGrowth",
  [](const ARKodeView& arkode_mem, sunrealtype etamxf) -> int
  {
    auto ARKodeSetMaxEFailGrowth_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype etamxf) -> int
    {
      auto lambda_result = ARKodeSetMaxEFailGrowth(arkode_mem.get(), etamxf);
      return lambda_result;
    };

    return ARKodeSetMaxEFailGrowth_adapt_solver_memory(arkode_mem, etamxf);
  },
  nb::arg("arkode_mem"), nb::arg("etamxf"));

sundials4py::scoped_def(
  m, "ARKodeSetSmallNumEFails",
  [](const ARKodeView& arkode_mem, int small_nef) -> int
  {
    auto ARKodeSetSmallNumEFails_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int small_nef) -> int
    {
      auto lambda_result = ARKodeSetSmallNumEFails(arkode_mem.get(), small_nef);
      return lambda_result;
    };

    return ARKodeSetSmallNumEFails_adapt_solver_memory(arkode_mem, small_nef);
  },
  nb::arg("arkode_mem"), nb::arg("small_nef"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxCFailGrowth",
  [](const ARKodeView& arkode_mem, sunrealtype etacf) -> int
  {
    auto ARKodeSetMaxCFailGrowth_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype etacf) -> int
    {
      auto lambda_result = ARKodeSetMaxCFailGrowth(arkode_mem.get(), etacf);
      return lambda_result;
    };

    return ARKodeSetMaxCFailGrowth_adapt_solver_memory(arkode_mem, etacf);
  },
  nb::arg("arkode_mem"), nb::arg("etacf"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxErrTestFails",
  [](const ARKodeView& arkode_mem, int maxnef) -> int
  {
    auto ARKodeSetMaxErrTestFails_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int maxnef) -> int
    {
      auto lambda_result = ARKodeSetMaxErrTestFails(arkode_mem.get(), maxnef);
      return lambda_result;
    };

    return ARKodeSetMaxErrTestFails_adapt_solver_memory(arkode_mem, maxnef);
  },
  nb::arg("arkode_mem"), nb::arg("maxnef"));

sundials4py::scoped_def(
  m, "ARKodeSetConstraints",
  [](const ARKodeView& arkode_mem, N_Vector constraints) -> int
  {
    auto ARKodeSetConstraints_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, N_Vector constraints) -> int
    {
      auto lambda_result = ARKodeSetConstraints(arkode_mem.get(), constraints);
      return lambda_result;
    };

    return ARKodeSetConstraints_adapt_solver_memory(arkode_mem, constraints);
  },
  nb::arg("arkode_mem"), nb::arg("constraints"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxHnilWarns",
  [](const ARKodeView& arkode_mem, int mxhnil) -> int
  {
    auto ARKodeSetMaxHnilWarns_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int mxhnil) -> int
    {
      auto lambda_result = ARKodeSetMaxHnilWarns(arkode_mem.get(), mxhnil);
      return lambda_result;
    };

    return ARKodeSetMaxHnilWarns_adapt_solver_memory(arkode_mem, mxhnil);
  },
  nb::arg("arkode_mem"), nb::arg("mxhnil"));

sundials4py::scoped_def(
  m, "ARKodeSetInitStep",
  [](const ARKodeView& arkode_mem, sunrealtype hin) -> int
  {
    auto ARKodeSetInitStep_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                                    sunrealtype hin) -> int
    {
      auto lambda_result = ARKodeSetInitStep(arkode_mem.get(), hin);
      return lambda_result;
    };

    return ARKodeSetInitStep_adapt_solver_memory(arkode_mem, hin);
  },
  nb::arg("arkode_mem"), nb::arg("hin"));

sundials4py::scoped_def(
  m, "ARKodeSetMinStep",
  [](const ARKodeView& arkode_mem, sunrealtype hmin) -> int
  {
    auto ARKodeSetMinStep_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                                   sunrealtype hmin) -> int
    {
      auto lambda_result = ARKodeSetMinStep(arkode_mem.get(), hmin);
      return lambda_result;
    };

    return ARKodeSetMinStep_adapt_solver_memory(arkode_mem, hmin);
  },
  nb::arg("arkode_mem"), nb::arg("hmin"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxStep",
  [](const ARKodeView& arkode_mem, sunrealtype hmax) -> int
  {
    auto ARKodeSetMaxStep_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                                   sunrealtype hmax) -> int
    {
      auto lambda_result = ARKodeSetMaxStep(arkode_mem.get(), hmax);
      return lambda_result;
    };

    return ARKodeSetMaxStep_adapt_solver_memory(arkode_mem, hmax);
  },
  nb::arg("arkode_mem"), nb::arg("hmax"));

sundials4py::scoped_def(
  m, "ARKodeSetMaxNumConstrFails",
  [](const ARKodeView& arkode_mem, int maxfails) -> int
  {
    auto ARKodeSetMaxNumConstrFails_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int maxfails) -> int
    {
      auto lambda_result = ARKodeSetMaxNumConstrFails(arkode_mem.get(), maxfails);
      return lambda_result;
    };

    return ARKodeSetMaxNumConstrFails_adapt_solver_memory(arkode_mem, maxfails);
  },
  nb::arg("arkode_mem"), nb::arg("maxfails"));

sundials4py::scoped_def(
  m, "ARKodeSetAdjointCheckpointScheme",
  [](const ARKodeView& arkode_mem,
     SUNAdjointCheckpointScheme checkpoint_scheme) -> int
  {
    auto ARKodeSetAdjointCheckpointScheme_adapt_solver_memory =
      [](const ARKodeView& arkode_mem,
         SUNAdjointCheckpointScheme checkpoint_scheme) -> int
    {
      auto lambda_result = ARKodeSetAdjointCheckpointScheme(arkode_mem.get(),
                                                            checkpoint_scheme);
      return lambda_result;
    };

    return ARKodeSetAdjointCheckpointScheme_adapt_solver_memory(arkode_mem,
                                                                checkpoint_scheme);
  },
  nb::arg("arkode_mem"), nb::arg("checkpoint_scheme"));

sundials4py::scoped_def(
  m, "ARKodeSetAdjointCheckpointIndex",
  [](const ARKodeView& arkode_mem, suncountertype step_index) -> int
  {
    auto ARKodeSetAdjointCheckpointIndex_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, suncountertype step_index) -> int
    {
      auto lambda_result = ARKodeSetAdjointCheckpointIndex(arkode_mem.get(),
                                                           step_index);
      return lambda_result;
    };

    return ARKodeSetAdjointCheckpointIndex_adapt_solver_memory(arkode_mem,
                                                               step_index);
  },
  nb::arg("arkode_mem"), nb::arg("step_index"));

sundials4py::scoped_def(
  m, "ARKodeSetUseCompensatedSums",
  [](const ARKodeView& arkode_mem, sunbooleantype onoff) -> int
  {
    auto ARKodeSetUseCompensatedSums_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunbooleantype onoff) -> int
    {
      auto lambda_result = ARKodeSetUseCompensatedSums(arkode_mem.get(), onoff);
      return lambda_result;
    };

    return ARKodeSetUseCompensatedSums_adapt_solver_memory(arkode_mem, onoff);
  },
  nb::arg("arkode_mem"), nb::arg("onoff"));

sundials4py::scoped_def(
  m, "ARKodeSetAccumulatedErrorType",
  [](const ARKodeView& arkode_mem, ARKAccumError accum_type) -> int
  {
    auto ARKodeSetAccumulatedErrorType_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, ARKAccumError accum_type) -> int
    {
      auto lambda_result = ARKodeSetAccumulatedErrorType(arkode_mem.get(),
                                                         accum_type);
      return lambda_result;
    };

    return ARKodeSetAccumulatedErrorType_adapt_solver_memory(arkode_mem,
                                                             accum_type);
  },
  nb::arg("arkode_mem"), nb::arg("accum_type"));

sundials4py::scoped_def(
  m, "ARKodeResetAccumulatedError",
  [](const ARKodeView& arkode_mem) -> int
  {
    auto ARKodeResetAccumulatedError_adapt_solver_memory =
      [](const ARKodeView& arkode_mem) -> int
    {
      auto lambda_result = ARKodeResetAccumulatedError(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeResetAccumulatedError_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetDky",
  [](const ARKodeView& arkode_mem, sunrealtype t, int k, N_Vector dky) -> int
  {
    auto ARKodeGetDky_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype t, int k, N_Vector dky) -> int
    {
      auto lambda_result = ARKodeGetDky(arkode_mem.get(), t, k, dky);
      return lambda_result;
    };

    return ARKodeGetDky_adapt_solver_memory(arkode_mem, t, k, dky);
  },
  nb::arg("arkode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dky"),
  "Computes the kth derivative of the y function at time t");

sundials4py::scoped_def(
  m, "ARKodeComputeState",
  [](const ARKodeView& arkode_mem, N_Vector zcor, N_Vector z) -> int
  {
    auto ARKodeComputeState_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, N_Vector zcor, N_Vector z) -> int
    {
      auto lambda_result = ARKodeComputeState(arkode_mem.get(), zcor, z);
      return lambda_result;
    };

    return ARKodeComputeState_adapt_solver_memory(arkode_mem, zcor, z);
  },
  nb::arg("arkode_mem"), nb::arg("zcor"), nb::arg("z"),
  "Utility function to update/compute y based on zcor");

sundials4py::scoped_def(
  m, "ARKodeGetNumRhsEvals",
  [](const ARKodeView& arkode_mem, int partition_index) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, int partition_index) -> std::tuple<int, long>
    {
      long num_rhs_evals_adapt_modifiable;

      int r = ARKodeGetNumRhsEvals(arkode_mem, partition_index,
                                   &num_rhs_evals_adapt_modifiable);
      return std::make_tuple(r, num_rhs_evals_adapt_modifiable);
    };
    auto ARKodeGetNumRhsEvals_adapt_solver_memory =
      [&ARKodeGetNumRhsEvals_adapt_modifiable_immutable_to_return](const ARKodeView&
                                                                     arkode_mem,
                                                                   int partition_index)
      -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumRhsEvals_adapt_modifiable_immutable_to_return(arkode_mem.get(),
                                                                  partition_index);
      return lambda_result;
    };

    return ARKodeGetNumRhsEvals_adapt_solver_memory(arkode_mem, partition_index);
  },
  nb::arg("arkode_mem"), nb::arg("partition_index"));

sundials4py::scoped_def(
  m, "ARKodeGetNumStepAttempts",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumStepAttempts_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long step_attempts_adapt_modifiable;

      int r = ARKodeGetNumStepAttempts(arkode_mem,
                                       &step_attempts_adapt_modifiable);
      return std::make_tuple(r, step_attempts_adapt_modifiable);
    };
    auto ARKodeGetNumStepAttempts_adapt_solver_memory =
      [&ARKodeGetNumStepAttempts_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumStepAttempts_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumStepAttempts_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumSteps",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumSteps_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nsteps_adapt_modifiable;

      int r = ARKodeGetNumSteps(arkode_mem, &nsteps_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable);
    };
    auto ARKodeGetNumSteps_adapt_solver_memory =
      [&ARKodeGetNumSteps_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumSteps_adapt_modifiable_immutable_to_return(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumSteps_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetLastStep",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetLastStep_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hlast_adapt_modifiable;

      int r = ARKodeGetLastStep(arkode_mem, &hlast_adapt_modifiable);
      return std::make_tuple(r, hlast_adapt_modifiable);
    };
    auto ARKodeGetLastStep_adapt_solver_memory =
      [&ARKodeGetLastStep_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetLastStep_adapt_modifiable_immutable_to_return(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetLastStep_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetCurrentStep",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetCurrentStep_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hcur_adapt_modifiable;

      int r = ARKodeGetCurrentStep(arkode_mem, &hcur_adapt_modifiable);
      return std::make_tuple(r, hcur_adapt_modifiable);
    };
    auto ARKodeGetCurrentStep_adapt_solver_memory =
      [&ARKodeGetCurrentStep_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetCurrentStep_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetCurrentStep_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetStepDirection",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetStepDirection_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype stepdir_adapt_modifiable;

      int r = ARKodeGetStepDirection(arkode_mem, &stepdir_adapt_modifiable);
      return std::make_tuple(r, stepdir_adapt_modifiable);
    };
    auto ARKodeGetStepDirection_adapt_solver_memory =
      [&ARKodeGetStepDirection_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetStepDirection_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetStepDirection_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetErrWeights",
  [](const ARKodeView& arkode_mem, N_Vector eweight) -> int
  {
    auto ARKodeGetErrWeights_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, N_Vector eweight) -> int
    {
      auto lambda_result = ARKodeGetErrWeights(arkode_mem.get(), eweight);
      return lambda_result;
    };

    return ARKodeGetErrWeights_adapt_solver_memory(arkode_mem, eweight);
  },
  nb::arg("arkode_mem"), nb::arg("eweight"));

sundials4py::scoped_def(
  m, "ARKodeGetNumGEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumGEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long ngevals_adapt_modifiable;

      int r = ARKodeGetNumGEvals(arkode_mem, &ngevals_adapt_modifiable);
      return std::make_tuple(r, ngevals_adapt_modifiable);
    };
    auto ARKodeGetNumGEvals_adapt_solver_memory =
      [&ARKodeGetNumGEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumGEvals_adapt_modifiable_immutable_to_return(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumGEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetRootInfo",
  [](const ARKodeView& arkode_mem, sundials4py::IntArray1d rootsfound_1d) -> int
  {
    auto ARKodeGetRootInfo_adapt_arr_ptr_to_std_vector =
      [](void* arkode_mem, sundials4py::IntArray1d rootsfound_1d) -> int
    {
      int* rootsfound_1d_ptr = rootsfound_1d.size() == 0 ? nullptr
                                                         : rootsfound_1d.data();

      auto lambda_result = ARKodeGetRootInfo(arkode_mem, rootsfound_1d_ptr);
      return lambda_result;
    };
    auto ARKodeGetRootInfo_adapt_solver_memory =
      [&ARKodeGetRootInfo_adapt_arr_ptr_to_std_vector](const ARKodeView& arkode_mem,
                                                       sundials4py::IntArray1d
                                                         rootsfound_1d) -> int
    {
      auto lambda_result =
        ARKodeGetRootInfo_adapt_arr_ptr_to_std_vector(arkode_mem.get(),
                                                      rootsfound_1d);
      return lambda_result;
    };

    return ARKodeGetRootInfo_adapt_solver_memory(arkode_mem, rootsfound_1d);
  },
  nb::arg("arkode_mem"), nb::arg("rootsfound_1d"));

sundials4py::scoped_def(
  m, "ARKodePrintAllStats",
  [](const ARKodeView& arkode_mem, FILE* outfile, SUNOutputFormat fmt) -> int
  {
    auto ARKodePrintAllStats_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, FILE* outfile, SUNOutputFormat fmt) -> int
    {
      auto lambda_result = ARKodePrintAllStats(arkode_mem.get(), outfile, fmt);
      return lambda_result;
    };

    return ARKodePrintAllStats_adapt_solver_memory(arkode_mem, outfile, fmt);
  },
  nb::arg("arkode_mem"), nb::arg("outfile"), nb::arg("fmt"));

sundials4py::scoped_def(m, "ARKodeGetReturnFlagName", ARKodeGetReturnFlagName,
                        nb::arg("flag"));

sundials4py::scoped_def(
  m, "ARKodeWriteParameters",
  [](const ARKodeView& arkode_mem, FILE* fp) -> int
  {
    auto ARKodeWriteParameters_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, FILE* fp) -> int
    {
      auto lambda_result = ARKodeWriteParameters(arkode_mem.get(), fp);
      return lambda_result;
    };

    return ARKodeWriteParameters_adapt_solver_memory(arkode_mem, fp);
  },
  nb::arg("arkode_mem"), nb::arg("fp"));

sundials4py::scoped_def(
  m, "ARKodeGetStageIndex",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, int, int>
  {
    auto ARKodeGetStageIndex_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, int, int>
    {
      int stage_adapt_modifiable;
      int max_stages_adapt_modifiable;

      int r = ARKodeGetStageIndex(arkode_mem, &stage_adapt_modifiable,
                                  &max_stages_adapt_modifiable);
      return std::make_tuple(r, stage_adapt_modifiable,
                             max_stages_adapt_modifiable);
    };
    auto ARKodeGetStageIndex_adapt_solver_memory =
      [&ARKodeGetStageIndex_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, int, int>
    {
      auto lambda_result =
        ARKodeGetStageIndex_adapt_modifiable_immutable_to_return(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetStageIndex_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumExpSteps",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumExpSteps_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long expsteps_adapt_modifiable;

      int r = ARKodeGetNumExpSteps(arkode_mem, &expsteps_adapt_modifiable);
      return std::make_tuple(r, expsteps_adapt_modifiable);
    };
    auto ARKodeGetNumExpSteps_adapt_solver_memory =
      [&ARKodeGetNumExpSteps_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumExpSteps_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumExpSteps_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumAccSteps",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumAccSteps_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long accsteps_adapt_modifiable;

      int r = ARKodeGetNumAccSteps(arkode_mem, &accsteps_adapt_modifiable);
      return std::make_tuple(r, accsteps_adapt_modifiable);
    };
    auto ARKodeGetNumAccSteps_adapt_solver_memory =
      [&ARKodeGetNumAccSteps_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumAccSteps_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumAccSteps_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumErrTestFails",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long netfails_adapt_modifiable;

      int r = ARKodeGetNumErrTestFails(arkode_mem, &netfails_adapt_modifiable);
      return std::make_tuple(r, netfails_adapt_modifiable);
    };
    auto ARKodeGetNumErrTestFails_adapt_solver_memory =
      [&ARKodeGetNumErrTestFails_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumErrTestFails_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumErrTestFails_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetEstLocalErrors",
  [](const ARKodeView& arkode_mem, N_Vector ele) -> int
  {
    auto ARKodeGetEstLocalErrors_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, N_Vector ele) -> int
    {
      auto lambda_result = ARKodeGetEstLocalErrors(arkode_mem.get(), ele);
      return lambda_result;
    };

    return ARKodeGetEstLocalErrors_adapt_solver_memory(arkode_mem, ele);
  },
  nb::arg("arkode_mem"), nb::arg("ele"));

sundials4py::scoped_def(
  m, "ARKodeGetActualInitStep",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetActualInitStep_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hinused_adapt_modifiable;

      int r = ARKodeGetActualInitStep(arkode_mem, &hinused_adapt_modifiable);
      return std::make_tuple(r, hinused_adapt_modifiable);
    };
    auto ARKodeGetActualInitStep_adapt_solver_memory =
      [&ARKodeGetActualInitStep_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetActualInitStep_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetActualInitStep_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetTolScaleFactor",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetTolScaleFactor_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tolsfac_adapt_modifiable;

      int r = ARKodeGetTolScaleFactor(arkode_mem, &tolsfac_adapt_modifiable);
      return std::make_tuple(r, tolsfac_adapt_modifiable);
    };
    auto ARKodeGetTolScaleFactor_adapt_solver_memory =
      [&ARKodeGetTolScaleFactor_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetTolScaleFactor_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetTolScaleFactor_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumConstrFails",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumConstrFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nconstrfails_adapt_modifiable;

      int r = ARKodeGetNumConstrFails(arkode_mem, &nconstrfails_adapt_modifiable);
      return std::make_tuple(r, nconstrfails_adapt_modifiable);
    };
    auto ARKodeGetNumConstrFails_adapt_solver_memory =
      [&ARKodeGetNumConstrFails_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumConstrFails_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumConstrFails_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetStepStats",
  [](const ARKodeView& arkode_mem)
    -> std::tuple<int, long, sunrealtype, sunrealtype, sunrealtype, sunrealtype>
  {
    auto ARKodeGetStepStats_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem)
      -> std::tuple<int, long, sunrealtype, sunrealtype, sunrealtype, sunrealtype>
    {
      long nsteps_adapt_modifiable;
      sunrealtype hinused_adapt_modifiable;
      sunrealtype hlast_adapt_modifiable;
      sunrealtype hcur_adapt_modifiable;
      sunrealtype tcur_adapt_modifiable;

      int r = ARKodeGetStepStats(arkode_mem, &nsteps_adapt_modifiable,
                                 &hinused_adapt_modifiable,
                                 &hlast_adapt_modifiable,
                                 &hcur_adapt_modifiable, &tcur_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable,
                             hinused_adapt_modifiable, hlast_adapt_modifiable,
                             hcur_adapt_modifiable, tcur_adapt_modifiable);
    };
    auto ARKodeGetStepStats_adapt_solver_memory =
      [&ARKodeGetStepStats_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem)
      -> std::tuple<int, long, sunrealtype, sunrealtype, sunrealtype, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetStepStats_adapt_modifiable_immutable_to_return(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetStepStats_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetAccumulatedError",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetAccumulatedError_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype accum_error_adapt_modifiable;

      int r = ARKodeGetAccumulatedError(arkode_mem,
                                        &accum_error_adapt_modifiable);
      return std::make_tuple(r, accum_error_adapt_modifiable);
    };
    auto ARKodeGetAccumulatedError_adapt_solver_memory =
      [&ARKodeGetAccumulatedError_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetAccumulatedError_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetAccumulatedError_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumLinSolvSetups",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nlinsetups_adapt_modifiable;

      int r = ARKodeGetNumLinSolvSetups(arkode_mem, &nlinsetups_adapt_modifiable);
      return std::make_tuple(r, nlinsetups_adapt_modifiable);
    };
    auto ARKodeGetNumLinSolvSetups_adapt_solver_memory =
      [&ARKodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumLinSolvSetups_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetLastTime",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetLastTime_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tn_adapt_modifiable;

      int r = ARKodeGetLastTime(arkode_mem, &tn_adapt_modifiable);
      return std::make_tuple(r, tn_adapt_modifiable);
    };
    auto ARKodeGetLastTime_adapt_solver_memory =
      [&ARKodeGetLastTime_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetLastTime_adapt_modifiable_immutable_to_return(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetLastTime_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetLastState",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, N_Vector>
  {
    auto ARKodeGetLastState_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, N_Vector>
    {
      N_Vector state_adapt_modifiable;

      int r = ARKodeGetLastState(arkode_mem, &state_adapt_modifiable);
      return std::make_tuple(r, state_adapt_modifiable);
    };
    auto ARKodeGetLastState_adapt_solver_memory =
      [&ARKodeGetLastState_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, N_Vector>
    {
      auto lambda_result =
        ARKodeGetLastState_adapt_modifiable_immutable_to_return(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetLastState_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "ARKodeGetCurrentTime",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetCurrentTime_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tcur_adapt_modifiable;

      int r = ARKodeGetCurrentTime(arkode_mem, &tcur_adapt_modifiable);
      return std::make_tuple(r, tcur_adapt_modifiable);
    };
    auto ARKodeGetCurrentTime_adapt_solver_memory =
      [&ARKodeGetCurrentTime_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetCurrentTime_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetCurrentTime_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetCurrentState",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, N_Vector>
  {
    auto ARKodeGetCurrentState_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, N_Vector>
    {
      N_Vector state_adapt_modifiable;

      int r = ARKodeGetCurrentState(arkode_mem, &state_adapt_modifiable);
      return std::make_tuple(r, state_adapt_modifiable);
    };
    auto ARKodeGetCurrentState_adapt_solver_memory =
      [&ARKodeGetCurrentState_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, N_Vector>
    {
      auto lambda_result =
        ARKodeGetCurrentState_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetCurrentState_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "ARKodeGetCurrentGamma",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetCurrentGamma_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype gamma_adapt_modifiable;

      int r = ARKodeGetCurrentGamma(arkode_mem, &gamma_adapt_modifiable);
      return std::make_tuple(r, gamma_adapt_modifiable);
    };
    auto ARKodeGetCurrentGamma_adapt_solver_memory =
      [&ARKodeGetCurrentGamma_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetCurrentGamma_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetCurrentGamma_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumNonlinSolvIters",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nniters_adapt_modifiable;

      int r = ARKodeGetNumNonlinSolvIters(arkode_mem, &nniters_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable);
    };
    auto ARKodeGetNumNonlinSolvIters_adapt_solver_memory =
      [&ARKodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumNonlinSolvIters_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumNonlinSolvConvFails",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nnfails_adapt_modifiable;

      int r = ARKodeGetNumNonlinSolvConvFails(arkode_mem,
                                              &nnfails_adapt_modifiable);
      return std::make_tuple(r, nnfails_adapt_modifiable);
    };
    auto ARKodeGetNumNonlinSolvConvFails_adapt_solver_memory =
      [&ARKodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumNonlinSolvConvFails_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNonlinSolvStats",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long, long>
  {
    auto ARKodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long, long>
    {
      long nniters_adapt_modifiable;
      long nnfails_adapt_modifiable;

      int r = ARKodeGetNonlinSolvStats(arkode_mem, &nniters_adapt_modifiable,
                                       &nnfails_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable,
                             nnfails_adapt_modifiable);
    };
    auto ARKodeGetNonlinSolvStats_adapt_solver_memory =
      [&ARKodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long, long>
    {
      auto lambda_result =
        ARKodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNonlinSolvStats_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumStepSolveFails",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nncfails_adapt_modifiable;

      int r = ARKodeGetNumStepSolveFails(arkode_mem, &nncfails_adapt_modifiable);
      return std::make_tuple(r, nncfails_adapt_modifiable);
    };
    auto ARKodeGetNumStepSolveFails_adapt_solver_memory =
      [&ARKodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumStepSolveFails_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetJac",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, SUNMatrix>
  {
    auto ARKodeGetJac_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix J_adapt_modifiable;

      int r = ARKodeGetJac(arkode_mem, &J_adapt_modifiable);
      return std::make_tuple(r, J_adapt_modifiable);
    };
    auto ARKodeGetJac_adapt_solver_memory =
      [&ARKodeGetJac_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, SUNMatrix>
    {
      auto lambda_result =
        ARKodeGetJac_adapt_modifiable_immutable_to_return(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetJac_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "ARKodeGetJacTime",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetJacTime_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_J_adapt_modifiable;

      int r = ARKodeGetJacTime(arkode_mem, &t_J_adapt_modifiable);
      return std::make_tuple(r, t_J_adapt_modifiable);
    };
    auto ARKodeGetJacTime_adapt_solver_memory =
      [&ARKodeGetJacTime_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, sunrealtype>
    {
      auto lambda_result =
        ARKodeGetJacTime_adapt_modifiable_immutable_to_return(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetJacTime_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetJacNumSteps",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetJacNumSteps_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nst_J_adapt_modifiable;

      int r = ARKodeGetJacNumSteps(arkode_mem, &nst_J_adapt_modifiable);
      return std::make_tuple(r, nst_J_adapt_modifiable);
    };
    auto ARKodeGetJacNumSteps_adapt_solver_memory =
      [&ARKodeGetJacNumSteps_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetJacNumSteps_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetJacNumSteps_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumJacEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long njevals_adapt_modifiable;

      int r = ARKodeGetNumJacEvals(arkode_mem, &njevals_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable);
    };
    auto ARKodeGetNumJacEvals_adapt_solver_memory =
      [&ARKodeGetNumJacEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumJacEvals_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumJacEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumPrecEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long npevals_adapt_modifiable;

      int r = ARKodeGetNumPrecEvals(arkode_mem, &npevals_adapt_modifiable);
      return std::make_tuple(r, npevals_adapt_modifiable);
    };
    auto ARKodeGetNumPrecEvals_adapt_solver_memory =
      [&ARKodeGetNumPrecEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumPrecEvals_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumPrecEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumPrecSolves",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long npsolves_adapt_modifiable;

      int r = ARKodeGetNumPrecSolves(arkode_mem, &npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves_adapt_modifiable);
    };
    auto ARKodeGetNumPrecSolves_adapt_solver_memory =
      [&ARKodeGetNumPrecSolves_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumPrecSolves_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumPrecSolves_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumLinIters",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nliters_adapt_modifiable;

      int r = ARKodeGetNumLinIters(arkode_mem, &nliters_adapt_modifiable);
      return std::make_tuple(r, nliters_adapt_modifiable);
    };
    auto ARKodeGetNumLinIters_adapt_solver_memory =
      [&ARKodeGetNumLinIters_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumLinIters_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumLinIters_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumLinConvFails",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nlcfails_adapt_modifiable;

      int r = ARKodeGetNumLinConvFails(arkode_mem, &nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails_adapt_modifiable);
    };
    auto ARKodeGetNumLinConvFails_adapt_solver_memory =
      [&ARKodeGetNumLinConvFails_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumLinConvFails_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumLinConvFails_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumJTSetupEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long njtsetups_adapt_modifiable;

      int r = ARKodeGetNumJTSetupEvals(arkode_mem, &njtsetups_adapt_modifiable);
      return std::make_tuple(r, njtsetups_adapt_modifiable);
    };
    auto ARKodeGetNumJTSetupEvals_adapt_solver_memory =
      [&ARKodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumJTSetupEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumJtimesEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long njvevals_adapt_modifiable;

      int r = ARKodeGetNumJtimesEvals(arkode_mem, &njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals_adapt_modifiable);
    };
    auto ARKodeGetNumJtimesEvals_adapt_solver_memory =
      [&ARKodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumJtimesEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumLinRhsEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nfevalsLS_adapt_modifiable;

      int r = ARKodeGetNumLinRhsEvals(arkode_mem, &nfevalsLS_adapt_modifiable);
      return std::make_tuple(r, nfevalsLS_adapt_modifiable);
    };
    auto ARKodeGetNumLinRhsEvals_adapt_solver_memory =
      [&ARKodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumLinRhsEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetLastLinFlag",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = ARKodeGetLastLinFlag(arkode_mem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };
    auto ARKodeGetLastLinFlag_adapt_solver_memory =
      [&ARKodeGetLastLinFlag_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetLastLinFlag_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetLastLinFlag_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(m, "ARKodeGetLinReturnFlagName",
                        ARKodeGetLinReturnFlagName, nb::arg("flag"));

sundials4py::scoped_def(
  m, "ARKodeGetCurrentMassMatrix",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, SUNMatrix>
  {
    auto ARKodeGetCurrentMassMatrix_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix M_adapt_modifiable;

      int r = ARKodeGetCurrentMassMatrix(arkode_mem, &M_adapt_modifiable);
      return std::make_tuple(r, M_adapt_modifiable);
    };
    auto ARKodeGetCurrentMassMatrix_adapt_solver_memory =
      [&ARKodeGetCurrentMassMatrix_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, SUNMatrix>
    {
      auto lambda_result =
        ARKodeGetCurrentMassMatrix_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetCurrentMassMatrix_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"), " Optional output functions (non-identity mass matrices)\n\n nb::rv_policy::reference",
  nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "ARKodeGetResWeights",
  [](const ARKodeView& arkode_mem, N_Vector rweight) -> int
  {
    auto ARKodeGetResWeights_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, N_Vector rweight) -> int
    {
      auto lambda_result = ARKodeGetResWeights(arkode_mem.get(), rweight);
      return lambda_result;
    };

    return ARKodeGetResWeights_adapt_solver_memory(arkode_mem, rweight);
  },
  nb::arg("arkode_mem"), nb::arg("rweight"));

sundials4py::scoped_def(
  m, "ARKodeGetNumMassSetups",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassSetups_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmsetups_adapt_modifiable;

      int r = ARKodeGetNumMassSetups(arkode_mem, &nmsetups_adapt_modifiable);
      return std::make_tuple(r, nmsetups_adapt_modifiable);
    };
    auto ARKodeGetNumMassSetups_adapt_solver_memory =
      [&ARKodeGetNumMassSetups_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumMassSetups_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumMassSetups_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumMassMultSetups",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassMultSetups_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmvsetups_adapt_modifiable;

      int r = ARKodeGetNumMassMultSetups(arkode_mem, &nmvsetups_adapt_modifiable);
      return std::make_tuple(r, nmvsetups_adapt_modifiable);
    };
    auto ARKodeGetNumMassMultSetups_adapt_solver_memory =
      [&ARKodeGetNumMassMultSetups_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumMassMultSetups_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumMassMultSetups_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumMassMult",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassMult_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmvevals_adapt_modifiable;

      int r = ARKodeGetNumMassMult(arkode_mem, &nmvevals_adapt_modifiable);
      return std::make_tuple(r, nmvevals_adapt_modifiable);
    };
    auto ARKodeGetNumMassMult_adapt_solver_memory =
      [&ARKodeGetNumMassMult_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumMassMult_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumMassMult_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumMassSolves",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassSolves_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmsolves_adapt_modifiable;

      int r = ARKodeGetNumMassSolves(arkode_mem, &nmsolves_adapt_modifiable);
      return std::make_tuple(r, nmsolves_adapt_modifiable);
    };
    auto ARKodeGetNumMassSolves_adapt_solver_memory =
      [&ARKodeGetNumMassSolves_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumMassSolves_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumMassSolves_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumMassPrecEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmpevals_adapt_modifiable;

      int r = ARKodeGetNumMassPrecEvals(arkode_mem, &nmpevals_adapt_modifiable);
      return std::make_tuple(r, nmpevals_adapt_modifiable);
    };
    auto ARKodeGetNumMassPrecEvals_adapt_solver_memory =
      [&ARKodeGetNumMassPrecEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumMassPrecEvals_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumMassPrecEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumMassPrecSolves",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmpsolves_adapt_modifiable;

      int r = ARKodeGetNumMassPrecSolves(arkode_mem, &nmpsolves_adapt_modifiable);
      return std::make_tuple(r, nmpsolves_adapt_modifiable);
    };
    auto ARKodeGetNumMassPrecSolves_adapt_solver_memory =
      [&ARKodeGetNumMassPrecSolves_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumMassPrecSolves_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumMassPrecSolves_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumMassIters",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmiters_adapt_modifiable;

      int r = ARKodeGetNumMassIters(arkode_mem, &nmiters_adapt_modifiable);
      return std::make_tuple(r, nmiters_adapt_modifiable);
    };
    auto ARKodeGetNumMassIters_adapt_solver_memory =
      [&ARKodeGetNumMassIters_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumMassIters_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumMassIters_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumMassConvFails",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassConvFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmcfails_adapt_modifiable;

      int r = ARKodeGetNumMassConvFails(arkode_mem, &nmcfails_adapt_modifiable);
      return std::make_tuple(r, nmcfails_adapt_modifiable);
    };
    auto ARKodeGetNumMassConvFails_adapt_solver_memory =
      [&ARKodeGetNumMassConvFails_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumMassConvFails_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumMassConvFails_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumMTSetups",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMTSetups_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmtsetups_adapt_modifiable;

      int r = ARKodeGetNumMTSetups(arkode_mem, &nmtsetups_adapt_modifiable);
      return std::make_tuple(r, nmtsetups_adapt_modifiable);
    };
    auto ARKodeGetNumMTSetups_adapt_solver_memory =
      [&ARKodeGetNumMTSetups_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumMTSetups_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumMTSetups_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetLastMassFlag",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetLastMassFlag_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = ARKodeGetLastMassFlag(arkode_mem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };
    auto ARKodeGetLastMassFlag_adapt_solver_memory =
      [&ARKodeGetLastMassFlag_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetLastMassFlag_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetLastMassFlag_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodePrintMem",
  [](const ARKodeView& arkode_mem, FILE* outfile)
  {
    auto ARKodePrintMem_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, FILE* outfile)
    { ARKodePrintMem(arkode_mem.get(), outfile); };

    ARKodePrintMem_adapt_solver_memory(arkode_mem, outfile);
  },
  nb::arg("arkode_mem"), nb::arg("outfile"),
  "Output the ARKODE memory structure (useful when debugging)");

sundials4py::scoped_def(
  m, "ARKodeSetRelaxEtaFail",
  [](const ARKodeView& arkode_mem, sunrealtype eta_rf) -> int
  {
    auto ARKodeSetRelaxEtaFail_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype eta_rf) -> int
    {
      auto lambda_result = ARKodeSetRelaxEtaFail(arkode_mem.get(), eta_rf);
      return lambda_result;
    };

    return ARKodeSetRelaxEtaFail_adapt_solver_memory(arkode_mem, eta_rf);
  },
  nb::arg("arkode_mem"), nb::arg("eta_rf"));

sundials4py::scoped_def(
  m, "ARKodeSetRelaxLowerBound",
  [](const ARKodeView& arkode_mem, sunrealtype lower) -> int
  {
    auto ARKodeSetRelaxLowerBound_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype lower) -> int
    {
      auto lambda_result = ARKodeSetRelaxLowerBound(arkode_mem.get(), lower);
      return lambda_result;
    };

    return ARKodeSetRelaxLowerBound_adapt_solver_memory(arkode_mem, lower);
  },
  nb::arg("arkode_mem"), nb::arg("lower"));

sundials4py::scoped_def(
  m, "ARKodeSetRelaxMaxFails",
  [](const ARKodeView& arkode_mem, int max_fails) -> int
  {
    auto ARKodeSetRelaxMaxFails_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int max_fails) -> int
    {
      auto lambda_result = ARKodeSetRelaxMaxFails(arkode_mem.get(), max_fails);
      return lambda_result;
    };

    return ARKodeSetRelaxMaxFails_adapt_solver_memory(arkode_mem, max_fails);
  },
  nb::arg("arkode_mem"), nb::arg("max_fails"));

sundials4py::scoped_def(
  m, "ARKodeSetRelaxMaxIters",
  [](const ARKodeView& arkode_mem, int max_iters) -> int
  {
    auto ARKodeSetRelaxMaxIters_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, int max_iters) -> int
    {
      auto lambda_result = ARKodeSetRelaxMaxIters(arkode_mem.get(), max_iters);
      return lambda_result;
    };

    return ARKodeSetRelaxMaxIters_adapt_solver_memory(arkode_mem, max_iters);
  },
  nb::arg("arkode_mem"), nb::arg("max_iters"));

sundials4py::scoped_def(
  m, "ARKodeSetRelaxSolver",
  [](const ARKodeView& arkode_mem, ARKRelaxSolver solver) -> int
  {
    auto ARKodeSetRelaxSolver_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, ARKRelaxSolver solver) -> int
    {
      auto lambda_result = ARKodeSetRelaxSolver(arkode_mem.get(), solver);
      return lambda_result;
    };

    return ARKodeSetRelaxSolver_adapt_solver_memory(arkode_mem, solver);
  },
  nb::arg("arkode_mem"), nb::arg("solver"));

sundials4py::scoped_def(
  m, "ARKodeSetRelaxResTol",
  [](const ARKodeView& arkode_mem, sunrealtype res_tol) -> int
  {
    auto ARKodeSetRelaxResTol_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype res_tol) -> int
    {
      auto lambda_result = ARKodeSetRelaxResTol(arkode_mem.get(), res_tol);
      return lambda_result;
    };

    return ARKodeSetRelaxResTol_adapt_solver_memory(arkode_mem, res_tol);
  },
  nb::arg("arkode_mem"), nb::arg("res_tol"));

sundials4py::scoped_def(
  m, "ARKodeSetRelaxTol",
  [](const ARKodeView& arkode_mem, sunrealtype rel_tol, sunrealtype abs_tol) -> int
  {
    auto ARKodeSetRelaxTol_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                                    sunrealtype rel_tol,
                                                    sunrealtype abs_tol) -> int
    {
      auto lambda_result = ARKodeSetRelaxTol(arkode_mem.get(), rel_tol, abs_tol);
      return lambda_result;
    };

    return ARKodeSetRelaxTol_adapt_solver_memory(arkode_mem, rel_tol, abs_tol);
  },
  nb::arg("arkode_mem"), nb::arg("rel_tol"), nb::arg("abs_tol"));

sundials4py::scoped_def(
  m, "ARKodeSetRelaxUpperBound",
  [](const ARKodeView& arkode_mem, sunrealtype upper) -> int
  {
    auto ARKodeSetRelaxUpperBound_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype upper) -> int
    {
      auto lambda_result = ARKodeSetRelaxUpperBound(arkode_mem.get(), upper);
      return lambda_result;
    };

    return ARKodeSetRelaxUpperBound_adapt_solver_memory(arkode_mem, upper);
  },
  nb::arg("arkode_mem"), nb::arg("upper"));

sundials4py::scoped_def(
  m, "ARKodeGetNumRelaxFnEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxFnEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long r_evals_adapt_modifiable;

      int r = ARKodeGetNumRelaxFnEvals(arkode_mem, &r_evals_adapt_modifiable);
      return std::make_tuple(r, r_evals_adapt_modifiable);
    };
    auto ARKodeGetNumRelaxFnEvals_adapt_solver_memory =
      [&ARKodeGetNumRelaxFnEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumRelaxFnEvals_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumRelaxFnEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumRelaxJacEvals",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxJacEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long J_evals_adapt_modifiable;

      int r = ARKodeGetNumRelaxJacEvals(arkode_mem, &J_evals_adapt_modifiable);
      return std::make_tuple(r, J_evals_adapt_modifiable);
    };
    auto ARKodeGetNumRelaxJacEvals_adapt_solver_memory =
      [&ARKodeGetNumRelaxJacEvals_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumRelaxJacEvals_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumRelaxJacEvals_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumRelaxFails",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long relax_fails_adapt_modifiable;

      int r = ARKodeGetNumRelaxFails(arkode_mem, &relax_fails_adapt_modifiable);
      return std::make_tuple(r, relax_fails_adapt_modifiable);
    };
    auto ARKodeGetNumRelaxFails_adapt_solver_memory =
      [&ARKodeGetNumRelaxFails_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumRelaxFails_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumRelaxFails_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumRelaxBoundFails",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxBoundFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long fails_adapt_modifiable;

      int r = ARKodeGetNumRelaxBoundFails(arkode_mem, &fails_adapt_modifiable);
      return std::make_tuple(r, fails_adapt_modifiable);
    };
    auto ARKodeGetNumRelaxBoundFails_adapt_solver_memory =
      [&ARKodeGetNumRelaxBoundFails_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumRelaxBoundFails_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumRelaxBoundFails_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumRelaxSolveFails",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxSolveFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long fails_adapt_modifiable;

      int r = ARKodeGetNumRelaxSolveFails(arkode_mem, &fails_adapt_modifiable);
      return std::make_tuple(r, fails_adapt_modifiable);
    };
    auto ARKodeGetNumRelaxSolveFails_adapt_solver_memory =
      [&ARKodeGetNumRelaxSolveFails_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumRelaxSolveFails_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumRelaxSolveFails_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeGetNumRelaxSolveIters",
  [](const ARKodeView& arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxSolveIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long iters_adapt_modifiable;

      int r = ARKodeGetNumRelaxSolveIters(arkode_mem, &iters_adapt_modifiable);
      return std::make_tuple(r, iters_adapt_modifiable);
    };
    auto ARKodeGetNumRelaxSolveIters_adapt_solver_memory =
      [&ARKodeGetNumRelaxSolveIters_adapt_modifiable_immutable_to_return](
        const ARKodeView& arkode_mem) -> std::tuple<int, long>
    {
      auto lambda_result =
        ARKodeGetNumRelaxSolveIters_adapt_modifiable_immutable_to_return(
          arkode_mem.get());
      return lambda_result;
    };

    return ARKodeGetNumRelaxSolveIters_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKodeCreateSUNStepper",
  [](const ARKodeView& arkode_mem)
    -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<SUNStepper>>>
  {
    auto ARKodeCreateSUNStepper_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, SUNStepper>
    {
      SUNStepper stepper_adapt_modifiable;

      int r = ARKodeCreateSUNStepper(arkode_mem, &stepper_adapt_modifiable);
      return std::make_tuple(r, stepper_adapt_modifiable);
    };
    auto ARKodeCreateSUNStepper_adapt_return_type_to_shared_ptr =
      [&ARKodeCreateSUNStepper_adapt_modifiable_immutable_to_return](
        void* arkode_mem)
      -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<SUNStepper>>>
    {
      auto lambda_result =
        ARKodeCreateSUNStepper_adapt_modifiable_immutable_to_return(arkode_mem);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<SUNStepper>,
                                             SUNStepperDeleter>(
                               std::get<1>(lambda_result)));
    };
    auto ARKodeCreateSUNStepper_adapt_solver_memory =
      [&ARKodeCreateSUNStepper_adapt_return_type_to_shared_ptr](
        const ARKodeView& arkode_mem)
      -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<SUNStepper>>>
    {
      auto lambda_result =
        ARKodeCreateSUNStepper_adapt_return_type_to_shared_ptr(arkode_mem.get());
      return lambda_result;
    };

    return ARKodeCreateSUNStepper_adapt_solver_memory(arkode_mem);
  },
  nb::arg("arkode_mem"), "SUNStepper functions", nb::rv_policy::reference);
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
m.attr("ARKLS_SUCCESS")          = 0;
m.attr("ARKLS_MEM_NULL")         = -1;
m.attr("ARKLS_LMEM_NULL")        = -2;
m.attr("ARKLS_ILL_INPUT")        = -3;
m.attr("ARKLS_MEM_FAIL")         = -4;
m.attr("ARKLS_PMEM_NULL")        = -5;
m.attr("ARKLS_MASSMEM_NULL")     = -6;
m.attr("ARKLS_JACFUNC_UNRECVR")  = -7;
m.attr("ARKLS_JACFUNC_RECVR")    = -8;
m.attr("ARKLS_MASSFUNC_UNRECVR") = -9;
m.attr("ARKLS_MASSFUNC_RECVR")   = -10;
m.attr("ARKLS_SUNMAT_FAIL")      = -11;
m.attr("ARKLS_SUNLS_FAIL")       = -12;

sundials4py::scoped_def(
  m, "ARKodeSetLinearSolver",
  [](const ARKodeView& arkode_mem, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto ARKodeSetLinearSolver_adapt_optional_arg_with_default_null =
      [](void* arkode_mem, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = ARKodeSetLinearSolver(arkode_mem, LS,
                                                 A_adapt_default_null);
      return lambda_result;
    };
    auto ARKodeSetLinearSolver_adapt_solver_memory =
      [&ARKodeSetLinearSolver_adapt_optional_arg_with_default_null](const ARKodeView&
                                                                      arkode_mem,
                                                                    SUNLinearSolver LS,
                                                                    std::optional<SUNMatrix>
                                                                      A = std::nullopt)
      -> int
    {
      auto lambda_result =
        ARKodeSetLinearSolver_adapt_optional_arg_with_default_null(arkode_mem.get(),
                                                                   LS, A);
      return lambda_result;
    };

    return ARKodeSetLinearSolver_adapt_solver_memory(arkode_mem, LS, A);
  },
  nb::arg("arkode_mem"), nb::arg("LS"), nb::arg("A").none() = nb::none());

sundials4py::scoped_def(
  m, "ARKodeSetJacEvalFrequency",
  [](const ARKodeView& arkode_mem, long msbj) -> int
  {
    auto ARKodeSetJacEvalFrequency_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, long msbj) -> int
    {
      auto lambda_result = ARKodeSetJacEvalFrequency(arkode_mem.get(), msbj);
      return lambda_result;
    };

    return ARKodeSetJacEvalFrequency_adapt_solver_memory(arkode_mem, msbj);
  },
  nb::arg("arkode_mem"), nb::arg("msbj"));

sundials4py::scoped_def(
  m, "ARKodeSetLinearSolutionScaling",
  [](const ARKodeView& arkode_mem, sunbooleantype onoff) -> int
  {
    auto ARKodeSetLinearSolutionScaling_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunbooleantype onoff) -> int
    {
      auto lambda_result = ARKodeSetLinearSolutionScaling(arkode_mem.get(),
                                                          onoff);
      return lambda_result;
    };

    return ARKodeSetLinearSolutionScaling_adapt_solver_memory(arkode_mem, onoff);
  },
  nb::arg("arkode_mem"), nb::arg("onoff"));

sundials4py::scoped_def(
  m, "ARKodeSetEpsLin",
  [](const ARKodeView& arkode_mem, sunrealtype eplifac) -> int
  {
    auto ARKodeSetEpsLin_adapt_solver_memory = [](const ARKodeView& arkode_mem,
                                                  sunrealtype eplifac) -> int
    {
      auto lambda_result = ARKodeSetEpsLin(arkode_mem.get(), eplifac);
      return lambda_result;
    };

    return ARKodeSetEpsLin_adapt_solver_memory(arkode_mem, eplifac);
  },
  nb::arg("arkode_mem"), nb::arg("eplifac"));

sundials4py::scoped_def(
  m, "ARKodeSetMassEpsLin",
  [](const ARKodeView& arkode_mem, sunrealtype eplifac) -> int
  {
    auto ARKodeSetMassEpsLin_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype eplifac) -> int
    {
      auto lambda_result = ARKodeSetMassEpsLin(arkode_mem.get(), eplifac);
      return lambda_result;
    };

    return ARKodeSetMassEpsLin_adapt_solver_memory(arkode_mem, eplifac);
  },
  nb::arg("arkode_mem"), nb::arg("eplifac"));

sundials4py::scoped_def(
  m, "ARKodeSetLSNormFactor",
  [](const ARKodeView& arkode_mem, sunrealtype nrmfac) -> int
  {
    auto ARKodeSetLSNormFactor_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype nrmfac) -> int
    {
      auto lambda_result = ARKodeSetLSNormFactor(arkode_mem.get(), nrmfac);
      return lambda_result;
    };

    return ARKodeSetLSNormFactor_adapt_solver_memory(arkode_mem, nrmfac);
  },
  nb::arg("arkode_mem"), nb::arg("nrmfac"));

sundials4py::scoped_def(
  m, "ARKodeSetMassLSNormFactor",
  [](const ARKodeView& arkode_mem, sunrealtype nrmfac) -> int
  {
    auto ARKodeSetMassLSNormFactor_adapt_solver_memory =
      [](const ARKodeView& arkode_mem, sunrealtype nrmfac) -> int
    {
      auto lambda_result = ARKodeSetMassLSNormFactor(arkode_mem.get(), nrmfac);
      return lambda_result;
    };

    return ARKodeSetMassLSNormFactor_adapt_solver_memory(arkode_mem, nrmfac);
  },
  nb::arg("arkode_mem"), nb::arg("nrmfac"));
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
  nb::class_<ARKodeButcherTableMem>(m, "ARKodeButcherTableMem", "")
    .def(nb::init<>()) // implicit default constructor
  ;

sundials4py::scoped_def(
  m, "ARKodeButcherTable_Copy",
  [](ARKodeButcherTable B)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_Copy_adapt_return_type_to_shared_ptr =
      [](ARKodeButcherTable B)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_Copy(B);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_Copy_adapt_return_type_to_shared_ptr(B);
  },
  nb::arg("B"));

sundials4py::scoped_def(m, "ARKodeButcherTable_Write", ARKodeButcherTable_Write,
                        nb::arg("B"), nb::arg("outfile"));

sundials4py::scoped_def(m, "ARKodeButcherTable_IsStifflyAccurate",
                        ARKodeButcherTable_IsStifflyAccurate, nb::arg("B"));

sundials4py::scoped_def(
  m, "ARKodeButcherTable_CheckOrder",
  [](ARKodeButcherTable B, FILE* outfile) -> std::tuple<int, int, int>
  {
    auto ARKodeButcherTable_CheckOrder_adapt_modifiable_immutable_to_return =
      [](ARKodeButcherTable B, FILE* outfile) -> std::tuple<int, int, int>
    {
      int q_adapt_modifiable;
      int p_adapt_modifiable;

      int r = ARKodeButcherTable_CheckOrder(B, &q_adapt_modifiable,
                                            &p_adapt_modifiable, outfile);
      return std::make_tuple(r, q_adapt_modifiable, p_adapt_modifiable);
    };

    return ARKodeButcherTable_CheckOrder_adapt_modifiable_immutable_to_return(B,
                                                                              outfile);
  },
  nb::arg("B"), nb::arg("outfile"));

sundials4py::scoped_def(
  m, "ARKodeButcherTable_CheckARKOrder",
  [](ARKodeButcherTable B1, ARKodeButcherTable B2,
     FILE* outfile) -> std::tuple<int, int, int>
  {
    auto ARKodeButcherTable_CheckARKOrder_adapt_modifiable_immutable_to_return =
      [](ARKodeButcherTable B1, ARKodeButcherTable B2,
         FILE* outfile) -> std::tuple<int, int, int>
    {
      int q_adapt_modifiable;
      int p_adapt_modifiable;

      int r = ARKodeButcherTable_CheckARKOrder(B1, B2, &q_adapt_modifiable,
                                               &p_adapt_modifiable, outfile);
      return std::make_tuple(r, q_adapt_modifiable, p_adapt_modifiable);
    };

    return ARKodeButcherTable_CheckARKOrder_adapt_modifiable_immutable_to_return(B1,
                                                                                 B2,
                                                                                 outfile);
  },
  nb::arg("B1"), nb::arg("B2"), nb::arg("outfile"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _ARKODE_ERK_TABLES_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumARKODE_ERKTableID =
  nb::enum_<ARKODE_ERKTableID>(m, "ARKODE_ERKTableID", nb::is_arithmetic(), "")
    .value("ARKODE_ERK_NONE", ARKODE_ERK_NONE, "ensure enum is signed int")
    .value("ARKODE_HEUN_EULER_2_1_2", ARKODE_HEUN_EULER_2_1_2, "")
    .value("ARKODE_MIN_ERK_NUM", ARKODE_MIN_ERK_NUM, "")
    .value("ARKODE_BOGACKI_SHAMPINE_4_2_3", ARKODE_BOGACKI_SHAMPINE_4_2_3, "")
    .value("ARKODE_ARK324L2SA_ERK_4_2_3", ARKODE_ARK324L2SA_ERK_4_2_3, "")
    .value("ARKODE_ZONNEVELD_5_3_4", ARKODE_ZONNEVELD_5_3_4, "")
    .value("ARKODE_ARK436L2SA_ERK_6_3_4", ARKODE_ARK436L2SA_ERK_6_3_4, "")
    .value("ARKODE_SAYFY_ABURUB_6_3_4", ARKODE_SAYFY_ABURUB_6_3_4, "")
    .value("ARKODE_CASH_KARP_6_4_5", ARKODE_CASH_KARP_6_4_5, "")
    .value("ARKODE_FEHLBERG_6_4_5", ARKODE_FEHLBERG_6_4_5, "")
    .value("ARKODE_DORMAND_PRINCE_7_4_5", ARKODE_DORMAND_PRINCE_7_4_5, "")
    .value("ARKODE_ARK548L2SA_ERK_8_4_5", ARKODE_ARK548L2SA_ERK_8_4_5, "")
    .value("ARKODE_VERNER_8_5_6", ARKODE_VERNER_8_5_6, "")
    .value("ARKODE_FEHLBERG_13_7_8", ARKODE_FEHLBERG_13_7_8, "")
    .value("ARKODE_KNOTH_WOLKE_3_3", ARKODE_KNOTH_WOLKE_3_3, "")
    .value("ARKODE_ARK437L2SA_ERK_7_3_4", ARKODE_ARK437L2SA_ERK_7_3_4, "")
    .value("ARKODE_ARK548L2SAb_ERK_8_4_5", ARKODE_ARK548L2SAb_ERK_8_4_5, "")
    .value("ARKODE_GKC21_ERK_3_1_2", ARKODE_GKC21_ERK_3_1_2, "")
    .value("ARKODE_SOFRONIOU_SPALETTA_5_3_4", ARKODE_SOFRONIOU_SPALETTA_5_3_4, "")
    .value("ARKODE_SHU_OSHER_3_2_3", ARKODE_SHU_OSHER_3_2_3, "")
    .value("ARKODE_VERNER_9_5_6", ARKODE_VERNER_9_5_6, "")
    .value("ARKODE_VERNER_10_6_7", ARKODE_VERNER_10_6_7, "")
    .value("ARKODE_VERNER_13_7_8", ARKODE_VERNER_13_7_8, "")
    .value("ARKODE_VERNER_16_8_9", ARKODE_VERNER_16_8_9, "")
    .value("ARKODE_FORWARD_EULER_1_1", ARKODE_FORWARD_EULER_1_1, "")
    .value("ARKODE_RALSTON_EULER_2_1_2", ARKODE_RALSTON_EULER_2_1_2, "")
    .value("ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2",
           ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2, "")
    .value("ARKODE_RALSTON_3_1_2", ARKODE_RALSTON_3_1_2, "")
    .value("ARKODE_TSITOURAS_7_4_5", ARKODE_TSITOURAS_7_4_5, "")
    .value("ARKODE_SSP_ERK_3_1_2", ARKODE_SSP_ERK_3_1_2, "")
    .value("ARKODE_SSP_ERK_4_1_2", ARKODE_SSP_ERK_4_1_2, "")
    .value("ARKODE_SSP_ERK_4_2_3", ARKODE_SSP_ERK_4_2_3, "")
    .value("ARKODE_SSP_ERK_10_3_4", ARKODE_SSP_ERK_10_3_4, "")
    .value("ARKODE_SSP_LSPUM_ERK_3_1_2", ARKODE_SSP_LSPUM_ERK_3_1_2, "")
    .value("ARKODE_ARS222_ERK_3_1_2", ARKODE_ARS222_ERK_3_1_2, "")
    .value("ARKODE_KUTTA_RK4a_4_4", ARKODE_KUTTA_RK4a_4_4, "")
    .value("ARKODE_KUTTA_RK4b_4_4", ARKODE_KUTTA_RK4b_4_4, "")
    .value("ARKODE_MAX_ERK_NUM", ARKODE_MAX_ERK_NUM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

sundials4py::scoped_def(
  m, "ARKodeButcherTable_LoadERK",
  [](ARKODE_ERKTableID emethod)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_LoadERK_adapt_return_type_to_shared_ptr =
      [](ARKODE_ERKTableID emethod)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_LoadERK(emethod);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_LoadERK_adapt_return_type_to_shared_ptr(emethod);
  },
  nb::arg("emethod"), "Accessor routine to load built-in ERK table");

sundials4py::scoped_def(
  m, "ARKodeButcherTable_LoadERKByName",
  [](const char* emethod)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_LoadERKByName_adapt_return_type_to_shared_ptr =
      [](const char* emethod)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_LoadERKByName(emethod);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_LoadERKByName_adapt_return_type_to_shared_ptr(
      emethod);
  },
  nb::arg("emethod"));

sundials4py::scoped_def(m, "ARKodeButcherTable_ERKIDToName",
                        ARKodeButcherTable_ERKIDToName, nb::arg("emethod"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _ARKODE_DIRK_TABLES_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumARKODE_DIRKTableID =
  nb::enum_<ARKODE_DIRKTableID>(m, "ARKODE_DIRKTableID", nb::is_arithmetic(), "")
    .value("ARKODE_DIRK_NONE", ARKODE_DIRK_NONE, "ensure enum is signed int")
    .value("ARKODE_SDIRK_2_1_2", ARKODE_SDIRK_2_1_2, "")
    .value("ARKODE_MIN_DIRK_NUM", ARKODE_MIN_DIRK_NUM, "")
    .value("ARKODE_BILLINGTON_3_3_2", ARKODE_BILLINGTON_3_3_2, "")
    .value("ARKODE_TRBDF2_3_3_2", ARKODE_TRBDF2_3_3_2, "")
    .value("ARKODE_KVAERNO_4_2_3", ARKODE_KVAERNO_4_2_3, "")
    .value("ARKODE_ARK324L2SA_DIRK_4_2_3", ARKODE_ARK324L2SA_DIRK_4_2_3, "")
    .value("ARKODE_CASH_5_2_4", ARKODE_CASH_5_2_4, "")
    .value("ARKODE_CASH_5_3_4", ARKODE_CASH_5_3_4, "")
    .value("ARKODE_SDIRK_5_3_4", ARKODE_SDIRK_5_3_4, "")
    .value("ARKODE_KVAERNO_5_3_4", ARKODE_KVAERNO_5_3_4, "")
    .value("ARKODE_ARK436L2SA_DIRK_6_3_4", ARKODE_ARK436L2SA_DIRK_6_3_4, "")
    .value("ARKODE_KVAERNO_7_4_5", ARKODE_KVAERNO_7_4_5, "")
    .value("ARKODE_ARK548L2SA_DIRK_8_4_5", ARKODE_ARK548L2SA_DIRK_8_4_5, "")
    .value("ARKODE_ARK437L2SA_DIRK_7_3_4", ARKODE_ARK437L2SA_DIRK_7_3_4, "")
    .value("ARKODE_ARK548L2SAb_DIRK_8_4_5", ARKODE_ARK548L2SAb_DIRK_8_4_5, "")
    .value("ARKODE_ESDIRK324L2SA_4_2_3", ARKODE_ESDIRK324L2SA_4_2_3, "")
    .value("ARKODE_ESDIRK325L2SA_5_2_3", ARKODE_ESDIRK325L2SA_5_2_3, "")
    .value("ARKODE_ESDIRK32I5L2SA_5_2_3", ARKODE_ESDIRK32I5L2SA_5_2_3, "")
    .value("ARKODE_ESDIRK436L2SA_6_3_4", ARKODE_ESDIRK436L2SA_6_3_4, "")
    .value("ARKODE_ESDIRK43I6L2SA_6_3_4", ARKODE_ESDIRK43I6L2SA_6_3_4, "")
    .value("ARKODE_QESDIRK436L2SA_6_3_4", ARKODE_QESDIRK436L2SA_6_3_4, "")
    .value("ARKODE_ESDIRK437L2SA_7_3_4", ARKODE_ESDIRK437L2SA_7_3_4, "")
    .value("ARKODE_ESDIRK547L2SA_7_4_5", ARKODE_ESDIRK547L2SA_7_4_5, "")
    .value("ARKODE_ESDIRK547L2SA2_7_4_5", ARKODE_ESDIRK547L2SA2_7_4_5, "")
    .value("ARKODE_GKC21_ESDIRK_3_1_2", ARKODE_GKC21_ESDIRK_3_1_2, "")
    .value("ARKODE_BACKWARD_EULER_1_1", ARKODE_BACKWARD_EULER_1_1, "")
    .value("ARKODE_IMPLICIT_MIDPOINT_1_2", ARKODE_IMPLICIT_MIDPOINT_1_2, "")
    .value("ARKODE_IMPLICIT_TRAPEZOIDAL_2_2", ARKODE_IMPLICIT_TRAPEZOIDAL_2_2, "")
    .value("ARKODE_SSP_DIRK_3_1_2", ARKODE_SSP_DIRK_3_1_2, "")
    .value("ARKODE_SSP_LSPUM_SDIRK_3_1_2", ARKODE_SSP_LSPUM_SDIRK_3_1_2, "")
    .value("ARKODE_ESDIRK_4_2_3", ARKODE_ESDIRK_4_2_3, "")
    .value("ARKODE_ARS222_ESDIRK_3_1_2", ARKODE_ARS222_ESDIRK_3_1_2, "")
    .value("ARKODE_MAX_DIRK_NUM", ARKODE_MAX_DIRK_NUM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

sundials4py::scoped_def(
  m, "ARKodeButcherTable_LoadDIRK",
  [](ARKODE_DIRKTableID imethod)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_LoadDIRK_adapt_return_type_to_shared_ptr =
      [](ARKODE_DIRKTableID imethod)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_LoadDIRK(imethod);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_LoadDIRK_adapt_return_type_to_shared_ptr(imethod);
  },
  nb::arg("imethod"), "Accessor routine to load built-in DIRK table");

sundials4py::scoped_def(
  m, "ARKodeButcherTable_LoadDIRKByName",
  [](const char* imethod)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_LoadDIRKByName_adapt_return_type_to_shared_ptr =
      [](const char* imethod)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_LoadDIRKByName(imethod);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_LoadDIRKByName_adapt_return_type_to_shared_ptr(
      imethod);
  },
  nb::arg("imethod"), "Accessor routine to load built-in DIRK table");

sundials4py::scoped_def(m, "ARKodeButcherTable_DIRKIDToName",
                        ARKodeButcherTable_DIRKIDToName, nb::arg("imethod"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _ARKODE_SPRKTABLE_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumARKODE_SPRKMethodID =
  nb::enum_<ARKODE_SPRKMethodID>(m, "ARKODE_SPRKMethodID", nb::is_arithmetic(), "")
    .value("ARKODE_SPRK_NONE", ARKODE_SPRK_NONE, "ensure enum is signed int")
    .value("ARKODE_SPRK_EULER_1_1", ARKODE_SPRK_EULER_1_1, "")
    .value("ARKODE_MIN_SPRK_NUM", ARKODE_MIN_SPRK_NUM, "")
    .value("ARKODE_SPRK_LEAPFROG_2_2", ARKODE_SPRK_LEAPFROG_2_2, "")
    .value("ARKODE_SPRK_PSEUDO_LEAPFROG_2_2", ARKODE_SPRK_PSEUDO_LEAPFROG_2_2, "")
    .value("ARKODE_SPRK_RUTH_3_3", ARKODE_SPRK_RUTH_3_3, "")
    .value("ARKODE_SPRK_MCLACHLAN_2_2", ARKODE_SPRK_MCLACHLAN_2_2, "")
    .value("ARKODE_SPRK_MCLACHLAN_3_3", ARKODE_SPRK_MCLACHLAN_3_3, "")
    .value("ARKODE_SPRK_CANDY_ROZMUS_4_4", ARKODE_SPRK_CANDY_ROZMUS_4_4, "")
    .value("ARKODE_SPRK_MCLACHLAN_4_4", ARKODE_SPRK_MCLACHLAN_4_4, "")
    .value("ARKODE_SPRK_MCLACHLAN_5_6", ARKODE_SPRK_MCLACHLAN_5_6, "")
    .value("ARKODE_SPRK_YOSHIDA_6_8", ARKODE_SPRK_YOSHIDA_6_8, "")
    .value("ARKODE_SPRK_SUZUKI_UMENO_8_16", ARKODE_SPRK_SUZUKI_UMENO_8_16, "")
    .value("ARKODE_SPRK_SOFRONIOU_10_36", ARKODE_SPRK_SOFRONIOU_10_36, "")
    .value("ARKODE_MAX_SPRK_NUM", ARKODE_MAX_SPRK_NUM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClassARKodeSPRKTableMem =
  nb::class_<ARKodeSPRKTableMem>(m, "ARKodeSPRKTableMem", "")
    .def(nb::init<>()) // implicit default constructor
  ;

sundials4py::scoped_def(
  m, "ARKodeSPRKTable_Create",
  [](int s, int q, sundials4py::Array1d a_1d, sundials4py::Array1d ahat_1d)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
  {
    auto ARKodeSPRKTable_Create_adapt_arr_ptr_to_std_vector =
      [](int s, int q, sundials4py::Array1d a_1d,
         sundials4py::Array1d ahat_1d) -> ARKodeSPRKTable
    {
      sunrealtype* a_1d_ptr    = a_1d.size() == 0 ? nullptr : a_1d.data();
      sunrealtype* ahat_1d_ptr = ahat_1d.size() == 0 ? nullptr : ahat_1d.data();

      auto lambda_result = ARKodeSPRKTable_Create(s, q, a_1d_ptr, ahat_1d_ptr);
      return lambda_result;
    };
    auto ARKodeSPRKTable_Create_adapt_return_type_to_shared_ptr =
      [&ARKodeSPRKTable_Create_adapt_arr_ptr_to_std_vector](int s, int q,
                                                            sundials4py::Array1d a_1d,
                                                            sundials4py::Array1d ahat_1d)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
    {
      auto lambda_result =
        ARKodeSPRKTable_Create_adapt_arr_ptr_to_std_vector(s, q, a_1d, ahat_1d);

      return our_make_shared<std::remove_pointer_t<ARKodeSPRKTable>,
                             ARKodeSPRKTableDeleter>(lambda_result);
    };

    return ARKodeSPRKTable_Create_adapt_return_type_to_shared_ptr(s, q, a_1d,
                                                                  ahat_1d);
  },
  nb::arg("s"), nb::arg("q"), nb::arg("a_1d"), nb::arg("ahat_1d"));

sundials4py::scoped_def(
  m, "ARKodeSPRKTable_Load",
  [](ARKODE_SPRKMethodID id)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
  {
    auto ARKodeSPRKTable_Load_adapt_return_type_to_shared_ptr =
      [](ARKODE_SPRKMethodID id)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
    {
      auto lambda_result = ARKodeSPRKTable_Load(id);

      return our_make_shared<std::remove_pointer_t<ARKodeSPRKTable>,
                             ARKodeSPRKTableDeleter>(lambda_result);
    };

    return ARKodeSPRKTable_Load_adapt_return_type_to_shared_ptr(id);
  },
  nb::arg("id"));

sundials4py::scoped_def(
  m, "ARKodeSPRKTable_LoadByName",
  [](const char* method) -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
  {
    auto ARKodeSPRKTable_LoadByName_adapt_return_type_to_shared_ptr =
      [](const char* method)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
    {
      auto lambda_result = ARKodeSPRKTable_LoadByName(method);

      return our_make_shared<std::remove_pointer_t<ARKodeSPRKTable>,
                             ARKodeSPRKTableDeleter>(lambda_result);
    };

    return ARKodeSPRKTable_LoadByName_adapt_return_type_to_shared_ptr(method);
  },
  nb::arg("method"));

sundials4py::scoped_def(
  m, "ARKodeSPRKTable_Copy",
  [](ARKodeSPRKTable that_sprk_storage)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
  {
    auto ARKodeSPRKTable_Copy_adapt_return_type_to_shared_ptr =
      [](ARKodeSPRKTable that_sprk_storage)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
    {
      auto lambda_result = ARKodeSPRKTable_Copy(that_sprk_storage);

      return our_make_shared<std::remove_pointer_t<ARKodeSPRKTable>,
                             ARKodeSPRKTableDeleter>(lambda_result);
    };

    return ARKodeSPRKTable_Copy_adapt_return_type_to_shared_ptr(that_sprk_storage);
  },
  nb::arg("that_sprk_storage"));

sundials4py::scoped_def(m, "ARKodeSPRKTable_Write", ARKodeSPRKTable_Write,
                        nb::arg("sprk_table"), nb::arg("outfile"));

sundials4py::scoped_def(
  m, "ARKodeSPRKTable_ToButcher",
  [](ARKodeSPRKTable sprk_storage)
    -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>,
                  std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>>
  {
    auto ARKodeSPRKTable_ToButcher_adapt_modifiable_immutable_to_return =
      [](ARKodeSPRKTable sprk_storage)
      -> std::tuple<int, ARKodeButcherTable, ARKodeButcherTable>
    {
      ARKodeButcherTable a_ptr_adapt_modifiable;
      ARKodeButcherTable b_ptr_adapt_modifiable;

      int r = ARKodeSPRKTable_ToButcher(sprk_storage, &a_ptr_adapt_modifiable,
                                        &b_ptr_adapt_modifiable);
      return std::make_tuple(r, a_ptr_adapt_modifiable, b_ptr_adapt_modifiable);
    };
    auto ARKodeSPRKTable_ToButcher_adapt_return_type_to_shared_ptr =
      [&ARKodeSPRKTable_ToButcher_adapt_modifiable_immutable_to_return](
        ARKodeSPRKTable sprk_storage)
      -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>,
                    std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>>
    {
      auto lambda_result =
        ARKodeSPRKTable_ToButcher_adapt_modifiable_immutable_to_return(
          sprk_storage);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                                             ARKodeButcherTableDeleter>(
                               std::get<1>(lambda_result)),
                             our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                                             ARKodeButcherTableDeleter>(
                               std::get<2>(lambda_result)));
    };

    return ARKodeSPRKTable_ToButcher_adapt_return_type_to_shared_ptr(sprk_storage);
  },
  nb::arg("sprk_storage"), nb::rv_policy::reference);
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
