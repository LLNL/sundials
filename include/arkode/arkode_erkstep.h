/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the ARKode ERKStep module.
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

/* Create, Resize, and Reinitialization functions */
SUNDIALS_EXPORT void* ERKStepCreate(ARKRhsFn f, sunrealtype t0, N_Vector y0,
                                    SUNContext sunctx);

SUNDIALS_EXPORT int ERKStepResize(void* arkode_mem, N_Vector ynew,
                                  sunrealtype hscale, sunrealtype t0,
                                  ARKVecResizeFn resize, void* resize_data);

SUNDIALS_EXPORT int ERKStepReInit(void* arkode_mem, ARKRhsFn f, sunrealtype t0,
                                  N_Vector y0);

SUNDIALS_EXPORT int ERKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR);

/* Tolerance input functions */
SUNDIALS_EXPORT int ERKStepSStolerances(void* arkode_mem, sunrealtype reltol,
                                        sunrealtype abstol);
SUNDIALS_EXPORT int ERKStepSVtolerances(void* arkode_mem, sunrealtype reltol,
                                        N_Vector abstol);
SUNDIALS_EXPORT int ERKStepWFtolerances(void* arkode_mem, ARKEwtFn efun);

/* Rootfinding initialization */
SUNDIALS_EXPORT int ERKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g);

/* Optional input functions -- must be called AFTER ERKStepCreate */
SUNDIALS_EXPORT int ERKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int ERKStepSetOrder(void* arkode_mem, int maxord);
SUNDIALS_EXPORT int ERKStepSetInterpolantType(void* arkode_mem, int itype);
SUNDIALS_EXPORT int ERKStepSetInterpolantDegree(void* arkode_mem, int degree);
SUNDIALS_EXPORT int ERKStepSetDenseOrder(void* arkode_mem, int dord);
SUNDIALS_EXPORT int ERKStepSetTable(void* arkode_mem, ARKodeButcherTable B);
SUNDIALS_EXPORT int ERKStepSetTableNum(void* arkode_mem,
                                       ARKODE_ERKTableID etable);
SUNDIALS_EXPORT int ERKStepSetTableName(void* arkode_mem, const char* etable);
SUNDIALS_EXPORT int ERKStepSetAdaptController(void* arkode_mem,
                                              SUNAdaptController C);
SUNDIALS_EXPORT int ERKStepSetAdaptivityAdjustment(void* arkode_mem, int adjust);
SUNDIALS_EXPORT int ERKStepSetCFLFraction(void* arkode_mem, sunrealtype cfl_frac);
SUNDIALS_EXPORT int ERKStepSetSafetyFactor(void* arkode_mem, sunrealtype safety);
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController instead")
int ERKStepSetErrorBias(void* arkode_mem, sunrealtype bias);
SUNDIALS_EXPORT int ERKStepSetMaxGrowth(void* arkode_mem, sunrealtype mx_growth);
SUNDIALS_EXPORT int ERKStepSetMinReduction(void* arkode_mem, sunrealtype eta_min);
SUNDIALS_EXPORT int ERKStepSetFixedStepBounds(void* arkode_mem, sunrealtype lb,
                                              sunrealtype ub);
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController instead")
int ERKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault,
                               int pq, sunrealtype adapt_params[3]);
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController instead")
int ERKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data);
SUNDIALS_EXPORT int ERKStepSetMaxFirstGrowth(void* arkode_mem,
                                             sunrealtype etamx1);
SUNDIALS_EXPORT int ERKStepSetMaxEFailGrowth(void* arkode_mem,
                                             sunrealtype etamxf);
SUNDIALS_EXPORT int ERKStepSetSmallNumEFails(void* arkode_mem, int small_nef);
SUNDIALS_EXPORT int ERKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab,
                                          void* estab_data);
SUNDIALS_EXPORT int ERKStepSetMaxErrTestFails(void* arkode_mem, int maxnef);
SUNDIALS_EXPORT int ERKStepSetConstraints(void* arkode_mem, N_Vector constraints);
SUNDIALS_EXPORT int ERKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps);
SUNDIALS_EXPORT int ERKStepSetMaxHnilWarns(void* arkode_mem, int mxhnil);
SUNDIALS_EXPORT int ERKStepSetInitStep(void* arkode_mem, sunrealtype hin);
SUNDIALS_EXPORT int ERKStepSetMinStep(void* arkode_mem, sunrealtype hmin);
SUNDIALS_EXPORT int ERKStepSetMaxStep(void* arkode_mem, sunrealtype hmax);
SUNDIALS_EXPORT int ERKStepSetInterpolateStopTime(void* arkode_mem,
                                                  sunbooleantype interp);
SUNDIALS_EXPORT int ERKStepSetStopTime(void* arkode_mem, sunrealtype tstop);
SUNDIALS_EXPORT int ERKStepClearStopTime(void* arkode_mem);
SUNDIALS_EXPORT int ERKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed);
SUNDIALS_EXPORT int ERKStepSetMaxNumConstrFails(void* arkode_mem, int maxfails);

SUNDIALS_EXPORT int ERKStepSetRootDirection(void* arkode_mem, int* rootdir);
SUNDIALS_EXPORT int ERKStepSetNoInactiveRootWarn(void* arkode_mem);

SUNDIALS_EXPORT int ERKStepSetErrHandlerFn(void* arkode_mem,
                                           ARKErrHandlerFn ehfun, void* eh_data);
SUNDIALS_EXPORT int ERKStepSetErrFile(void* arkode_mem, FILE* errfp);
SUNDIALS_EXPORT int ERKStepSetUserData(void* arkode_mem, void* user_data);

SUNDIALS_EXPORT int ERKStepSetPostprocessStepFn(void* arkode_mem,
                                                ARKPostProcessFn ProcessStep);
SUNDIALS_EXPORT int ERKStepSetPostprocessStageFn(void* arkode_mem,
                                                 ARKPostProcessFn ProcessStage);

/* Integrate the ODE over an interval in t */
SUNDIALS_EXPORT int ERKStepEvolve(void* arkode_mem, sunrealtype tout,
                                  N_Vector yout, sunrealtype* tret, int itask);

/* Computes the kth derivative of the y function at time t */
SUNDIALS_EXPORT int ERKStepGetDky(void* arkode_mem, sunrealtype t, int k,
                                  N_Vector dky);

/* Optional output functions */
SUNDIALS_EXPORT int ERKStepGetNumExpSteps(void* arkode_mem, long int* expsteps);
SUNDIALS_EXPORT int ERKStepGetNumAccSteps(void* arkode_mem, long int* accsteps);
SUNDIALS_EXPORT int ERKStepGetNumStepAttempts(void* arkode_mem,
                                              long int* step_attempts);
SUNDIALS_EXPORT int ERKStepGetNumRhsEvals(void* arkode_mem, long int* nfevals);
SUNDIALS_EXPORT int ERKStepGetNumErrTestFails(void* arkode_mem,
                                              long int* netfails);
SUNDIALS_EXPORT int ERKStepGetCurrentButcherTable(void* arkode_mem,
                                                  ARKodeButcherTable* B);
SUNDIALS_EXPORT int ERKStepGetEstLocalErrors(void* arkode_mem, N_Vector ele);
SUNDIALS_EXPORT int ERKStepGetWorkSpace(void* arkode_mem, long int* lenrw,
                                        long int* leniw);
SUNDIALS_EXPORT int ERKStepGetNumSteps(void* arkode_mem, long int* nsteps);
SUNDIALS_EXPORT int ERKStepGetActualInitStep(void* arkode_mem,
                                             sunrealtype* hinused);
SUNDIALS_EXPORT int ERKStepGetLastStep(void* arkode_mem, sunrealtype* hlast);
SUNDIALS_EXPORT int ERKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur);
SUNDIALS_EXPORT int ERKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur);
SUNDIALS_EXPORT int ERKStepGetTolScaleFactor(void* arkode_mem,
                                             sunrealtype* tolsfac);
SUNDIALS_EXPORT int ERKStepGetErrWeights(void* arkode_mem, N_Vector eweight);
SUNDIALS_EXPORT int ERKStepGetNumGEvals(void* arkode_mem, long int* ngevals);
SUNDIALS_EXPORT int ERKStepGetRootInfo(void* arkode_mem, int* rootsfound);
SUNDIALS_EXPORT int ERKStepGetNumConstrFails(void* arkode_mem,
                                             long int* nconstrfails);
SUNDIALS_EXPORT int ERKStepGetUserData(void* arkode_mem, void** user_data);
SUNDIALS_EXPORT int ERKStepPrintAllStats(void* arkode_mem, FILE* outfile,
                                         SUNOutputFormat fmt);
SUNDIALS_EXPORT char* ERKStepGetReturnFlagName(long int flag);

SUNDIALS_EXPORT int ERKStepWriteParameters(void* arkode_mem, FILE* fp);

SUNDIALS_EXPORT int ERKStepWriteButcher(void* arkode_mem, FILE* fp);

/* Grouped optional output functions */
SUNDIALS_EXPORT int ERKStepGetTimestepperStats(
  void* arkode_mem, long int* expsteps, long int* accsteps,
  long int* step_attempts, long int* nfevals, long int* netfails);
SUNDIALS_EXPORT int ERKStepGetStepStats(void* arkode_mem, long int* nsteps,
                                        sunrealtype* hinused, sunrealtype* hlast,
                                        sunrealtype* hcur, sunrealtype* tcur);

/* Free function */
SUNDIALS_EXPORT void ERKStepFree(void** arkode_mem);

/* Output the ERKStep memory structure (useful when debugging) */
SUNDIALS_EXPORT void ERKStepPrintMem(void* arkode_mem, FILE* outfile);

/* Relaxation functions */
SUNDIALS_EXPORT int ERKStepSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn,
                                      ARKRelaxJacFn rjac);
SUNDIALS_EXPORT int ERKStepSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_rf);
SUNDIALS_EXPORT int ERKStepSetRelaxLowerBound(void* arkode_mem,
                                              sunrealtype lower);
SUNDIALS_EXPORT int ERKStepSetRelaxMaxFails(void* arkode_mem, int max_fails);
SUNDIALS_EXPORT int ERKStepSetRelaxMaxIters(void* arkode_mem, int max_iters);
SUNDIALS_EXPORT int ERKStepSetRelaxSolver(void* arkode_mem,
                                          ARKRelaxSolver solver);
SUNDIALS_EXPORT int ERKStepSetRelaxResTol(void* arkode_mem, sunrealtype res_tol);
SUNDIALS_EXPORT int ERKStepSetRelaxTol(void* arkode_mem, sunrealtype rel_tol,
                                       sunrealtype abs_tol);
SUNDIALS_EXPORT int ERKStepSetRelaxUpperBound(void* arkode_mem,
                                              sunrealtype upper);
SUNDIALS_EXPORT int ERKStepGetNumRelaxFnEvals(void* arkode_mem,
                                              long int* r_evals);
SUNDIALS_EXPORT int ERKStepGetNumRelaxJacEvals(void* arkode_mem,
                                               long int* J_evals);
SUNDIALS_EXPORT int ERKStepGetNumRelaxFails(void* arkode_mem,
                                            long int* relax_fails);
SUNDIALS_EXPORT int ERKStepGetNumRelaxBoundFails(void* arkode_mem,
                                                 long int* fails);
SUNDIALS_EXPORT int ERKStepGetNumRelaxSolveFails(void* arkode_mem,
                                                 long int* fails);
SUNDIALS_EXPORT int ERKStepGetNumRelaxSolveIters(void* arkode_mem,
                                                 long int* iters);

#ifdef __cplusplus
}
#endif

#endif
