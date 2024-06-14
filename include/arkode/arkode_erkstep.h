/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the ARKODE ERKStep module.
 * -----------------------------------------------------------------*/

#ifndef _ERKSTEP_H
#define _ERKSTEP_H

#include <arkode/arkode.h>
#include <arkode/arkode_butcher_erk.h>
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * ERKStep Constants
 * ----------------- */

/* Default Butcher tables for each order */

static const int ERKSTEP_DEFAULT_1 = ARKODE_FORWARD_EULER_1_1;
static const int ERKSTEP_DEFAULT_2 = ARKODE_HEUN_EULER_2_1_2;
static const int ERKSTEP_DEFAULT_3 = ARKODE_BOGACKI_SHAMPINE_4_2_3;
static const int ERKSTEP_DEFAULT_4 = ARKODE_ZONNEVELD_5_3_4;
static const int ERKSTEP_DEFAULT_5 = ARKODE_CASH_KARP_6_4_5;
static const int ERKSTEP_DEFAULT_6 = ARKODE_VERNER_8_5_6;
static const int ERKSTEP_DEFAULT_7 = ARKODE_VERNER_10_6_7;
static const int ERKSTEP_DEFAULT_8 = ARKODE_FEHLBERG_13_7_8;
static const int ERKSTEP_DEFAULT_9 = ARKODE_VERNER_16_8_9;

/* -------------------
 * Exported Functions
 * ------------------- */

/* Creation and Reinitialization functions */
SUNDIALS_EXPORT void* ERKStepCreate(ARKRhsFn f, sunrealtype t0, N_Vector y0,
                                    SUNContext sunctx);
SUNDIALS_EXPORT int ERKStepReInit(void* arkode_mem, ARKRhsFn f, sunrealtype t0,
                                  N_Vector y0);

/* Optional input functions -- must be called AFTER ERKStepCreate */
SUNDIALS_EXPORT int ERKStepSetTable(void* arkode_mem, ARKodeButcherTable B);
SUNDIALS_EXPORT int ERKStepSetTableNum(void* arkode_mem,
                                       ARKODE_ERKTableID etable);
SUNDIALS_EXPORT int ERKStepSetTableName(void* arkode_mem, const char* etable);

/* Optional output functions */
SUNDIALS_EXPORT int ERKStepGetNumRhsEvals(void* arkode_mem, long int* nfevals);
SUNDIALS_EXPORT int ERKStepGetCurrentButcherTable(void* arkode_mem,
                                                  ARKodeButcherTable* B);

/* Grouped optional output functions */
SUNDIALS_EXPORT int ERKStepGetTimestepperStats(
  void* arkode_mem, long int* expsteps, long int* accsteps,
  long int* step_attempts, long int* nfevals, long int* netfails);

/* --------------------------------------------------------------------------
 * Deprecated Functions -- all are superseded by shared ARKODE-level routines
 * -------------------------------------------------------------------------- */

SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeResize instead")
int ERKStepResize(void* arkode_mem, N_Vector ynew, sunrealtype hscale,
                  sunrealtype t0, ARKVecResizeFn resize, void* resize_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeReset instead")
int ERKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSStolerances instead")
int ERKStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSVtolerances instead")
int ERKStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeWFtolerances instead")
int ERKStepWFtolerances(void* arkode_mem, ARKEwtFn efun);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeRootInit instead")
int ERKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDefaults instead")
int ERKStepSetDefaults(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetOrder instead")
int ERKStepSetOrder(void* arkode_mem, int maxord);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantType instead")
int ERKStepSetInterpolantType(void* arkode_mem, int itype);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantDegree instead")
int ERKStepSetInterpolantDegree(void* arkode_mem, int degree);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantDegree instead")
int ERKStepSetDenseOrder(void* arkode_mem, int dord);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetAdaptController instead")
int ERKStepSetAdaptController(void* arkode_mem, SUNAdaptController C);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetAdaptivityAdjustment instead")
int ERKStepSetAdaptivityAdjustment(void* arkode_mem, int adjust);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetCFLFraction instead")
int ERKStepSetCFLFraction(void* arkode_mem, sunrealtype cfl_frac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetSafetyFactor instead")
int ERKStepSetSafetyFactor(void* arkode_mem, sunrealtype safety);
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController instead")
int ERKStepSetErrorBias(void* arkode_mem, sunrealtype bias);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxGrowth instead")
int ERKStepSetMaxGrowth(void* arkode_mem, sunrealtype mx_growth);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMinReduction instead")
int ERKStepSetMinReduction(void* arkode_mem, sunrealtype eta_min);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepSetFiARKodeBounds instead")
int ERKStepSetFixedStepBounds(void* arkode_mem, sunrealtype lb, sunrealtype ub);
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController instead")
int ERKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault,
                               int pq, sunrealtype adapt_params[3]);
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController instead")
int ERKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxFirstGrowth instead")
int ERKStepSetMaxFirstGrowth(void* arkode_mem, sunrealtype etamx1);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxEFailGrowth instead")
int ERKStepSetMaxEFailGrowth(void* arkode_mem, sunrealtype etamxf);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetSmallNumEFails instead")
int ERKStepSetSmallNumEFails(void* arkode_mem, int small_nef);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStabilityFn instead")
int ERKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxErrTestFails instead")
int ERKStepSetMaxErrTestFails(void* arkode_mem, int maxnef);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetConstraints instead")
int ERKStepSetConstraints(void* arkode_mem, N_Vector constraints);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepSetMaxARKodes instead")
int ERKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxHnilWarns instead")
int ERKStepSetMaxHnilWarns(void* arkode_mem, int mxhnil);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepSetIARKode instead")
int ERKStepSetInitStep(void* arkode_mem, sunrealtype hin);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepSetARKode instead")
int ERKStepSetMinStep(void* arkode_mem, sunrealtype hmin);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepSetARKode instead")
int ERKStepSetMaxStep(void* arkode_mem, sunrealtype hmax);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolateStopTime instead")
int ERKStepSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStopTime instead")
int ERKStepSetStopTime(void* arkode_mem, sunrealtype tstop);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeClearStopTime instead")
int ERKStepClearStopTime(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepSetFiARKode instead")
int ERKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxNumConstrFails instead")
int ERKStepSetMaxNumConstrFails(void* arkode_mem, int maxfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRootDirection instead")
int ERKStepSetRootDirection(void* arkode_mem, int* rootdir);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNoInactiveRootWarn instead")
int ERKStepSetNoInactiveRootWarn(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetUserData instead")
int ERKStepSetUserData(void* arkode_mem, void* user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepSetPostprocARKodeFn instead")
int ERKStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPostprocessStageFn instead")
int ERKStepSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeEvolve instead")
int ERKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret, int itask);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetDky instead")
int ERKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepGetNumARKodes instead")
int ERKStepGetNumExpSteps(void* arkode_mem, long int* expsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepGetNumARKodes instead")
int ERKStepGetNumAccSteps(void* arkode_mem, long int* accsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepGetARKodeAttempts instead")
int ERKStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumErrTestFails instead")
int ERKStepGetNumErrTestFails(void* arkode_mem, long int* netfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetEstLocalErrors instead")
int ERKStepGetEstLocalErrors(void* arkode_mem, N_Vector ele);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetWorkSpace instead")
int ERKStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepGetARKodes instead")
int ERKStepGetNumSteps(void* arkode_mem, long int* nsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepGetActualIARKode instead")
int ERKStepGetActualInitStep(void* arkode_mem, sunrealtype* hinused);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepGetLARKode instead")
int ERKStepGetLastStep(void* arkode_mem, sunrealtype* hlast);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepGetCurrARKode instead")
int ERKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentTime instead")
int ERKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetTolScaleFactor instead")
int ERKStepGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetErrWeights instead")
int ERKStepGetErrWeights(void* arkode_mem, N_Vector eweight);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumGEvals instead")
int ERKStepGetNumGEvals(void* arkode_mem, long int* ngevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetRootInfo instead")
int ERKStepGetRootInfo(void* arkode_mem, int* rootsfound);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumConstrFails instead")
int ERKStepGetNumConstrFails(void* arkode_mem, long int* nconstrfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetUserData instead")
int ERKStepGetUserData(void* arkode_mem, void** user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodePrintAllStats instead")
int ERKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetReturnFlagName instead")
char* ERKStepGetReturnFlagName(long int flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeWriteParameters instead")
int ERKStepWriteParameters(void* arkode_mem, FILE* fp);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use ERKStepGetCurrentButcherTable and ARKodeButcherTable_Write instead")
int ERKStepWriteButcher(void* arkode_mem, FILE* fp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ERKStepARKodeStats instead")
int ERKStepGetStepStats(void* arkode_mem, long int* nsteps, sunrealtype* hinused,
                        sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeFree instead")
void ERKStepFree(void** arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodePrintMem instead")
void ERKStepPrintMem(void* arkode_mem, FILE* outfile);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxFn instead")
int ERKStepSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxEtaFail instead")
int ERKStepSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_rf);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxLowerBound instead")
int ERKStepSetRelaxLowerBound(void* arkode_mem, sunrealtype lower);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxMaxFails instead")
int ERKStepSetRelaxMaxFails(void* arkode_mem, int max_fails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxMaxIters instead")
int ERKStepSetRelaxMaxIters(void* arkode_mem, int max_iters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxSolver instead")
int ERKStepSetRelaxSolver(void* arkode_mem, ARKRelaxSolver solver);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxResTol instead")
int ERKStepSetRelaxResTol(void* arkode_mem, sunrealtype res_tol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxTol instead")
int ERKStepSetRelaxTol(void* arkode_mem, sunrealtype rel_tol,
                       sunrealtype abs_tol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxUpperBound instead")
int ERKStepSetRelaxUpperBound(void* arkode_mem, sunrealtype upper);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxFnEvals instead")
int ERKStepGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxJacEvals instead")
int ERKStepGetNumRelaxJacEvals(void* arkode_mem, long int* J_evals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxFails instead")
int ERKStepGetNumRelaxFails(void* arkode_mem, long int* relax_fails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxBoundFails instead")
int ERKStepGetNumRelaxBoundFails(void* arkode_mem, long int* fails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxSolveFails instead")
int ERKStepGetNumRelaxSolveFails(void* arkode_mem, long int* fails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxSolveIters instead")
int ERKStepGetNumRelaxSolveIters(void* arkode_mem, long int* iters);

#ifdef __cplusplus
}
#endif

#endif
