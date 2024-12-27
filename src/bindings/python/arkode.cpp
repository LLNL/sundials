#include <memory>
#include <iostream>
#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

#include <sundials/sundials_core.h>
#include <arkode/arkode.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode.hpp>

namespace nb = nanobind;

// TODO(CJB: we will need these wrappers for every callback function
using erk_rhsfn_type = int(sunrealtype, N_Vector, N_Vector, void*);
int erk_rhsfn_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data) {
  auto* object = static_cast<nb::object*>(user_data);
  auto erk_rhs = nb::cast<std::function<erk_rhsfn_type>>(*object);
  return erk_rhs(t, y, ydot, user_data);
}

void bind_arkode(nb::module_ &m) {
  //
  // ARKODE Constants
  //

  m.attr("ARK_NORMAL") = ARK_NORMAL;
  m.attr("ARK_ONE_STEP") = ARK_ONE_STEP;
  m.attr("ARK_ADAPT_CUSTOM") = ARK_ADAPT_CUSTOM;
  m.attr("ARK_ADAPT_PID") = ARK_ADAPT_PID;
  m.attr("ARK_ADAPT_PI") = ARK_ADAPT_PI;
  m.attr("ARK_ADAPT_I") = ARK_ADAPT_I;
  m.attr("ARK_ADAPT_EXP_GUS") = ARK_ADAPT_EXP_GUS;
  m.attr("ARK_ADAPT_IMP_GUS") = ARK_ADAPT_IMP_GUS;
  m.attr("ARK_ADAPT_IMEX_GUS") = ARK_ADAPT_IMEX_GUS;
  m.attr("ARK_FULLRHS_START") = ARK_FULLRHS_START;
  m.attr("ARK_FULLRHS_END") = ARK_FULLRHS_END;
  m.attr("ARK_FULLRHS_OTHER") = ARK_FULLRHS_OTHER;
  m.attr("ARK_INTERP_MAX_DEGREE") = ARK_INTERP_MAX_DEGREE;
  m.attr("ARK_INTERP_NONE") = ARK_INTERP_NONE;
  m.attr("ARK_INTERP_HERMITE") = ARK_INTERP_HERMITE;
  m.attr("ARK_INTERP_LAGRANGE") = ARK_INTERP_LAGRANGE;
  m.attr("ARK_SUCCESS") = ARK_SUCCESS;
  m.attr("ARK_TSTOP_RETURN") = ARK_TSTOP_RETURN;
  m.attr("ARK_ROOT_RETURN") = ARK_ROOT_RETURN;
  m.attr("ARK_WARNING") = ARK_WARNING;
  m.attr("ARK_TOO_MUCH_WORK") = ARK_TOO_MUCH_WORK;
  m.attr("ARK_TOO_MUCH_ACC") = ARK_TOO_MUCH_ACC;
  m.attr("ARK_ERR_FAILURE") = ARK_ERR_FAILURE;
  m.attr("ARK_CONV_FAILURE") = ARK_CONV_FAILURE;
  m.attr("ARK_LINIT_FAIL") = ARK_LINIT_FAIL;
  m.attr("ARK_LSETUP_FAIL") = ARK_LSETUP_FAIL;
  m.attr("ARK_LSOLVE_FAIL") = ARK_LSOLVE_FAIL;
  m.attr("ARK_RHSFUNC_FAIL") = ARK_RHSFUNC_FAIL;
  m.attr("ARK_FIRST_RHSFUNC_ERR") = ARK_FIRST_RHSFUNC_ERR;
  m.attr("ARK_REPTD_RHSFUNC_ERR") = ARK_REPTD_RHSFUNC_ERR;
  m.attr("ARK_UNREC_RHSFUNC_ERR") = ARK_UNREC_RHSFUNC_ERR;
  m.attr("ARK_RTFUNC_FAIL") = ARK_RTFUNC_FAIL;
  m.attr("ARK_LFREE_FAIL") = ARK_LFREE_FAIL;
  m.attr("ARK_MASSINIT_FAIL") = ARK_MASSINIT_FAIL;
  m.attr("ARK_MASSSETUP_FAIL") = ARK_MASSSETUP_FAIL;
  m.attr("ARK_MASSSOLVE_FAIL") = ARK_MASSSOLVE_FAIL;
  m.attr("ARK_MASSFREE_FAIL") = ARK_MASSFREE_FAIL;
  m.attr("ARK_MASSMULT_FAIL") = ARK_MASSMULT_FAIL;
  m.attr("ARK_CONSTR_FAIL") = ARK_CONSTR_FAIL;
  m.attr("ARK_MEM_FAIL") = ARK_MEM_FAIL;
  m.attr("ARK_MEM_NULL") = ARK_MEM_NULL;
  m.attr("ARK_ILL_INPUT") = ARK_ILL_INPUT;
  m.attr("ARK_NO_MALLOC") = ARK_NO_MALLOC;
  m.attr("ARK_BAD_K") = ARK_BAD_K;
  m.attr("ARK_BAD_T") = ARK_BAD_T;
  m.attr("ARK_BAD_DKY") = ARK_BAD_DKY;
  m.attr("ARK_TOO_CLOSE") = ARK_TOO_CLOSE;
  m.attr("ARK_VECTOROP_ERR") = ARK_VECTOROP_ERR;
  m.attr("ARK_NLS_INIT_FAIL") = ARK_NLS_INIT_FAIL;
  m.attr("ARK_NLS_SETUP_FAIL") = ARK_NLS_SETUP_FAIL;
  m.attr("ARK_NLS_SETUP_RECVR") = ARK_NLS_SETUP_RECVR;
  m.attr("ARK_NLS_OP_ERR") = ARK_NLS_OP_ERR;
  m.attr("ARK_INNERSTEP_ATTACH_ERR") = ARK_INNERSTEP_ATTACH_ERR;
  m.attr("ARK_INNERSTEP_FAIL") = ARK_INNERSTEP_FAIL;
  m.attr("ARK_OUTERTOINNER_FAIL") = ARK_OUTERTOINNER_FAIL;
  m.attr("ARK_INNERTOOUTER_FAIL") = ARK_INNERTOOUTER_FAIL;
  m.attr("ARK_POSTPROCESS_FAIL") = ARK_POSTPROCESS_FAIL;
  m.attr("ARK_POSTPROCESS_STEP_FAIL") = ARK_POSTPROCESS_STEP_FAIL;
  m.attr("ARK_POSTPROCESS_STAGE_FAIL") = ARK_POSTPROCESS_STAGE_FAIL;
  m.attr("ARK_USER_PREDICT_FAIL") = ARK_USER_PREDICT_FAIL;
  m.attr("ARK_INTERP_FAIL") = ARK_INTERP_FAIL;
  m.attr("ARK_INVALID_TABLE") = ARK_INVALID_TABLE;
  m.attr("ARK_CONTEXT_ERR") = ARK_CONTEXT_ERR;
  m.attr("ARK_RELAX_FAIL") = ARK_RELAX_FAIL;
  m.attr("ARK_RELAX_MEM_NULL") = ARK_RELAX_MEM_NULL;
  m.attr("ARK_RELAX_FUNC_FAIL") = ARK_RELAX_FUNC_FAIL;
  m.attr("ARK_RELAX_JAC_FAIL") = ARK_RELAX_JAC_FAIL;
  m.attr("ARK_CONTROLLER_ERR") = ARK_CONTROLLER_ERR;
  m.attr("ARK_STEPPER_UNSUPPORTED") = ARK_STEPPER_UNSUPPORTED;
  m.attr("ARK_DOMEIG_FAIL") = ARK_DOMEIG_FAIL;
  m.attr("ARK_MAX_STAGE_LIMIT_FAIL") = ARK_MAX_STAGE_LIMIT_FAIL;
  m.attr("ARK_SUNSTEPPER_ERR") = ARK_SUNSTEPPER_ERR;
  m.attr("ARK_STEP_DIRECTION_ERR") = ARK_STEP_DIRECTION_ERR;
  m.attr("ARK_UNRECOGNIZED_ERROR") = ARK_UNRECOGNIZED_ERROR;

  //
  // ARKODE functions
  //

  nb::class_<sundials::experimental::ARKodeView>(m, "ARKodeView")
    .def(nb::init<>())
    .def(nb::init<void*>())
    .def(nb::init_implicit<void*>())
    .def("Convert", nb::overload_cast<>(&sundials::experimental::ARKodeView::Convert, nb::const_), nb::rv_policy::reference);

  m.def("ARKodeResize", &ARKodeResize);
  m.def("ARKodeReset", &ARKodeReset);
  // m.def("ARKodeCreateMRIStepInnerStepper", &ARKodeCreateMRIStepInnerStepper);
  m.def("ARKodeSStolerances", &ARKodeSStolerances);
  m.def("ARKodeSVtolerances", &ARKodeSVtolerances);
  m.def("ARKodeWFtolerances", &ARKodeWFtolerances);
  m.def("ARKodeResStolerance", &ARKodeResStolerance);
  m.def("ARKodeResVtolerance", &ARKodeResVtolerance);
  m.def("ARKodeResFtolerance", &ARKodeResFtolerance);
  m.def("ARKodeRootInit", &ARKodeRootInit);
  m.def("ARKodeSetRootDirection", &ARKodeSetRootDirection);
  m.def("ARKodeSetNoInactiveRootWarn", &ARKodeSetNoInactiveRootWarn);
  m.def("ARKodeSetDefaults", &ARKodeSetDefaults);
  m.def("ARKodeSetOrder", &ARKodeSetOrder);
  m.def("ARKodeSetInterpolantType", &ARKodeSetInterpolantType);
  m.def("ARKodeSetInterpolantDegree", &ARKodeSetInterpolantDegree);
  m.def("ARKodeSetMaxNumSteps", &ARKodeSetMaxNumSteps);
  m.def("ARKodeSetInterpolateStopTime", &ARKodeSetInterpolateStopTime);
  m.def("ARKodeSetStopTime", &ARKodeSetStopTime);
  m.def("ARKodeClearStopTime", &ARKodeClearStopTime);
  m.def("ARKodeSetFixedStep", &ARKodeSetFixedStep);
  m.def("ARKodeSetStepDirection", &ARKodeSetStepDirection);
  m.def("ARKodeSetUserData", &ARKodeSetUserData);
  // m.def("ARKodeSetPostprocessStepFn", &ARKodeSetPostprocessStepFn);
  // m.def("ARKodeSetPostprocessStageFn", &ARKodeSetPostprocessStageFn);
  m.def("ARKodeSetNonlinearSolver", &ARKodeSetNonlinearSolver);
  m.def("ARKodeSetLinear", &ARKodeSetLinear);
  m.def("ARKodeSetNonlinear", &ARKodeSetNonlinear);
  m.def("ARKodeSetAutonomous", &ARKodeSetAutonomous);
  // m.def("ARKodeSetNlsRhsFn", &ARKodeSetNlsRhsFn);
  m.def("ARKodeSetDeduceImplicitRhs", &ARKodeSetDeduceImplicitRhs);
  m.def("ARKodeSetNonlinCRDown", &ARKodeSetNonlinCRDown);
  m.def("ARKodeSetNonlinRDiv", &ARKodeSetNonlinRDiv);
  m.def("ARKodeSetDeltaGammaMax", &ARKodeSetDeltaGammaMax);
  m.def("ARKodeSetLSetupFrequency", &ARKodeSetLSetupFrequency);
  m.def("ARKodeSetPredictorMethod", &ARKodeSetPredictorMethod);
  m.def("ARKodeSetMaxNonlinIters", &ARKodeSetMaxNonlinIters);
  m.def("ARKodeSetMaxConvFails", &ARKodeSetMaxConvFails);
  m.def("ARKodeSetNonlinConvCoef", &ARKodeSetNonlinConvCoef);
  // m.def("ARKodeSetStagePredictFn", &ARKodeSetStagePredictFn);
  m.def("ARKodeSetAdaptController", &ARKodeSetAdaptController);
  m.def("ARKodeSetAdaptivityAdjustment", &ARKodeSetAdaptivityAdjustment);
  m.def("ARKodeSetCFLFraction", &ARKodeSetCFLFraction);
  m.def("ARKodeSetErrorBias", &ARKodeSetErrorBias);
  m.def("ARKodeSetSafetyFactor", &ARKodeSetSafetyFactor);
  m.def("ARKodeSetMaxGrowth", &ARKodeSetMaxGrowth);
  m.def("ARKodeSetMinReduction", &ARKodeSetMinReduction);
  m.def("ARKodeSetFixedStepBounds", &ARKodeSetFixedStepBounds);
  m.def("ARKodeSetMaxFirstGrowth", &ARKodeSetMaxFirstGrowth);
  m.def("ARKodeSetMaxEFailGrowth", &ARKodeSetMaxEFailGrowth);
  m.def("ARKodeSetSmallNumEFails", &ARKodeSetSmallNumEFails);
  m.def("ARKodeSetMaxCFailGrowth", &ARKodeSetMaxCFailGrowth);
  // m.def("ARKodeSetStabilityFn", &ARKodeSetStabilityFn);
  m.def("ARKodeSetMaxErrTestFails", &ARKodeSetMaxErrTestFails);
  m.def("ARKodeSetConstraints", &ARKodeSetConstraints);
  m.def("ARKodeSetMaxHnilWarns", &ARKodeSetMaxHnilWarns);
  m.def("ARKodeSetInitStep", &ARKodeSetInitStep);
  m.def("ARKodeSetMinStep", &ARKodeSetMinStep);
  m.def("ARKodeSetMaxStep", &ARKodeSetMaxStep);
  m.def("ARKodeSetMaxNumConstrFails", &ARKodeSetMaxNumConstrFails);
  m.def("ARKodeSetAccumulatedErrorType", &ARKodeSetAccumulatedErrorType);
  m.def("ARKodeResetAccumulatedError", &ARKodeResetAccumulatedError);
  m.def("ARKodeEvolve", &ARKodeEvolve);
  m.def("ARKodeGetDky", &ARKodeGetDky);
  m.def("ARKodeComputeState", &ARKodeComputeState);
  m.def("ARKodeGetNumRhsEvals", &ARKodeGetNumRhsEvals);
  m.def("ARKodeGetNumStepAttempts", &ARKodeGetNumStepAttempts);
  m.def("ARKodeGetWorkSpace", &ARKodeGetWorkSpace);
  m.def("ARKodeGetNumSteps", &ARKodeGetNumSteps);
  m.def("ARKodeGetLastStep", &ARKodeGetLastStep);
  m.def("ARKodeGetCurrentStep", &ARKodeGetCurrentStep);
  m.def("ARKodeGetStepDirection", &ARKodeGetStepDirection);
  m.def("ARKodeGetErrWeights", &ARKodeGetErrWeights);
  m.def("ARKodeGetNumGEvals", &ARKodeGetNumGEvals);
  m.def("ARKodeGetRootInfo", &ARKodeGetRootInfo);
  // m.def("ARKodeGetUserData", &ARKodeGetUserData);
  m.def("ARKodePrintAllStats", &ARKodePrintAllStats);
  m.def("ARKodeGetReturnFlagName", &ARKodeGetReturnFlagName);
  m.def("ARKodeWriteParameters", &ARKodeWriteParameters);
  m.def("ARKodeGetNumExpSteps", &ARKodeGetNumExpSteps);
  m.def("ARKodeGetNumAccSteps", &ARKodeGetNumAccSteps);
  m.def("ARKodeGetNumErrTestFails", &ARKodeGetNumErrTestFails);
  m.def("ARKodeGetEstLocalErrors", &ARKodeGetEstLocalErrors);
  m.def("ARKodeGetActualInitStep", &ARKodeGetActualInitStep);
  m.def("ARKodeGetTolScaleFactor", &ARKodeGetTolScaleFactor);
  m.def("ARKodeGetNumConstrFails", &ARKodeGetNumConstrFails);
  m.def("ARKodeGetStepStats", &ARKodeGetStepStats);
  m.def("ARKodeGetAccumulatedError", &ARKodeGetAccumulatedError);
  m.def("ARKodeGetNumLinSolvSetups", &ARKodeGetNumLinSolvSetups);
  m.def("ARKodeGetCurrentTime", &ARKodeGetCurrentTime);
  // m.def("ARKodeGetCurrentState", &ARKodeGetCurrentState);
  m.def("ARKodeGetCurrentGamma", &ARKodeGetCurrentGamma);
  // m.def("ARKodeGetNonlinearSystemData", &ARKodeGetNonlinearSystemData);
  m.def("ARKodeGetNumNonlinSolvIters", &ARKodeGetNumNonlinSolvIters);
  m.def("ARKodeGetNumNonlinSolvConvFails", &ARKodeGetNumNonlinSolvConvFails);
  m.def("ARKodeGetNonlinSolvStats", &ARKodeGetNonlinSolvStats);
  m.def("ARKodeGetNumStepSolveFails", &ARKodeGetNumStepSolveFails);
  // m.def("ARKodeGetJac", &ARKodeGetJac);
  m.def("ARKodeGetJacTime", &ARKodeGetJacTime);
  m.def("ARKodeGetJacNumSteps", &ARKodeGetJacNumSteps);
  m.def("ARKodeGetLinWorkSpace", &ARKodeGetLinWorkSpace);
  m.def("ARKodeGetNumJacEvals", &ARKodeGetNumJacEvals);
  m.def("ARKodeGetNumPrecEvals", &ARKodeGetNumPrecEvals);
  m.def("ARKodeGetNumPrecSolves", &ARKodeGetNumPrecSolves);
  m.def("ARKodeGetNumLinIters", &ARKodeGetNumLinIters);
  m.def("ARKodeGetNumLinConvFails", &ARKodeGetNumLinConvFails);
  m.def("ARKodeGetNumJTSetupEvals", &ARKodeGetNumJTSetupEvals);
  m.def("ARKodeGetNumJtimesEvals", &ARKodeGetNumJtimesEvals);
  m.def("ARKodeGetNumLinRhsEvals", &ARKodeGetNumLinRhsEvals);
  m.def("ARKodeGetLastLinFlag", &ARKodeGetLastLinFlag);
  m.def("ARKodeGetLinReturnFlagName", &ARKodeGetLinReturnFlagName);
  // m.def("ARKodeGetCurrentMassMatrix", &ARKodeGetCurrentMassMatrix);
  m.def("ARKodeGetResWeights", &ARKodeGetResWeights);
  m.def("ARKodeGetMassWorkSpace", &ARKodeGetMassWorkSpace);
  m.def("ARKodeGetNumMassSetups", &ARKodeGetNumMassSetups);
  m.def("ARKodeGetNumMassMultSetups", &ARKodeGetNumMassMultSetups);
  m.def("ARKodeGetNumMassMult", &ARKodeGetNumMassMult);
  m.def("ARKodeGetNumMassSolves", &ARKodeGetNumMassSolves);
  m.def("ARKodeGetNumMassPrecEvals", &ARKodeGetNumMassPrecEvals);
  m.def("ARKodeGetNumMassPrecSolves", &ARKodeGetNumMassPrecSolves);
  m.def("ARKodeGetNumMassIters", &ARKodeGetNumMassIters);
  m.def("ARKodeGetNumMassConvFails", &ARKodeGetNumMassConvFails);
  m.def("ARKodeGetNumMTSetups", &ARKodeGetNumMTSetups);
  m.def("ARKodeGetLastMassFlag", &ARKodeGetLastMassFlag);
  // m.def("ARKodeFree", &ARKodeFree);
  m.def("ARKodePrintMem", &ARKodePrintMem);
  // m.def("ARKodeSetRelaxFn", &ARKodeSetRelaxFn);
  m.def("ARKodeSetRelaxEtaFail", &ARKodeSetRelaxEtaFail);
  m.def("ARKodeSetRelaxLowerBound", &ARKodeSetRelaxLowerBound);
  m.def("ARKodeSetRelaxMaxFails", &ARKodeSetRelaxMaxFails);
  m.def("ARKodeSetRelaxMaxIters", &ARKodeSetRelaxMaxIters);
  m.def("ARKodeSetRelaxSolver", &ARKodeSetRelaxSolver);
  m.def("ARKodeSetRelaxResTol", &ARKodeSetRelaxResTol);
  m.def("ARKodeSetRelaxTol", &ARKodeSetRelaxTol);
  m.def("ARKodeSetRelaxUpperBound", &ARKodeSetRelaxUpperBound);
  m.def("ARKodeGetNumRelaxFnEvals", &ARKodeGetNumRelaxFnEvals);
  m.def("ARKodeGetNumRelaxJacEvals", &ARKodeGetNumRelaxJacEvals);
  m.def("ARKodeGetNumRelaxFails", &ARKodeGetNumRelaxFails);
  m.def("ARKodeGetNumRelaxBoundFails", &ARKodeGetNumRelaxBoundFails);
  m.def("ARKodeGetNumRelaxSolveFails", &ARKodeGetNumRelaxSolveFails);
  m.def("ARKodeGetNumRelaxSolveIters", &ARKodeGetNumRelaxSolveIters);
  // m.def("ARKodeCreateSUNStepper", &ARKodeCreateSUNStepper);

  //
  // ERKStep functions
  //
  m.def("ERKStepCreate", [](std::function<erk_rhsfn_type> rhs, sunrealtype t0, N_Vector y0, SUNContext sunctx) {
    static nb::object the_static_rhs = nb::steal(nb::cast(rhs));
    auto ark_mem = ERKStepCreate(erk_rhsfn_wrapper, t0, y0, sunctx);
    ARKodeSetUserData(ark_mem, &the_static_rhs);
    return ark_mem;
  });
  m.def("ERKStepReInit", &ERKStepReInit);
  m.def("ERKStepSetTable", &ERKStepSetTable);
  m.def("ERKStepSetTableNum", &ERKStepSetTableNum);
  m.def("ERKStepSetTableName", &ERKStepSetTableName);
  // m.def("ERKStepGetCurrentButcherTable", &ERKStepGetCurrentButcherTable);
  m.def("ERKStepGetTimestepperStats", &ERKStepGetTimestepperStats);
}
