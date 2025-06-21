#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

#include <sundials/sundials_core.hpp>

#include <arkode/arkode.h>
#include <arkode/arkode.hpp>
#include <arkode/arkode_ls.h>

#include "sundials_adjointcheckpointscheme_impl.h"

namespace nb = nanobind;

// Forward declarations of functions defined in other translation units
void bind_arkode_erkstep(nb::module_& m);
void bind_arkode_arkstep(nb::module_& m);

void bind_arkode(nb::module_& m)
{
#include "pysundials_arkode_generated.hpp"

  nb::class_<sundials::experimental::ARKodeView>(m, "ARKodeView")
    .def(nb::init<>())
    .def(nb::init<void*>())
    .def("get",
         nb::overload_cast<>(&sundials::experimental::ARKodeView::get, nb::const_),
         nb::rv_policy::reference);

  //
  // arkode.h definitions
  //

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
  // arkode_ls.h definitions
  //

  m.def("ARKodeSetLinearSolver", &ARKodeSetLinearSolver, nb::arg("arkode_mem"),
        nb::arg("LS"), nb::arg("A").none());
  m.def("ARKodeSetJacFn", &ARKodeSetJacFn);
  m.def("ARKodeSetJacEvalFrequency", &ARKodeSetJacEvalFrequency);
  m.def("ARKodeSetLinSysFn", &ARKodeSetLinSysFn);
  m.def("ARKodeSetPreconditioner", &ARKodeSetPreconditioner);
  m.def("ARKodeSetLSNormFactor", &ARKodeSetLSNormFactor);
  m.def("ARKodeSetEpsLin", &ARKodeSetEpsLin);
  m.def("ARKodeGetNumLinIters", &ARKodeGetNumLinIters);
  m.def("ARKodeGetNumLinConvFails", &ARKodeGetNumLinConvFails);
  m.def("ARKodeGetNumPrecEvals", &ARKodeGetNumPrecEvals);
  m.def("ARKodeGetNumPrecSolves", &ARKodeGetNumPrecSolves);
  m.def("ARKodeGetNumLinRhsEvals", &ARKodeGetNumLinRhsEvals);
  m.def("ARKodeGetLastLinFlag", &ARKodeGetLastLinFlag);
  m.def("ARKodeGetLinReturnFlagName", &ARKodeGetLinReturnFlagName);

  bind_arkode_erkstep(m);
  bind_arkode_arkstep(m);
}
