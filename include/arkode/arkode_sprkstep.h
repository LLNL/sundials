/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the ARKode SPRKStep module.
 * -----------------------------------------------------------------*/

#ifndef _SPRKSTEP_H
#define _SPRKSTEP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <arkode/arkode.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_sprk.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * SPRKStep Constants
 * ----------------- */

static const int SPRKSTEP_DEFAULT_1 = ARKODE_SYMPLECTIC_EULER_1;
static const int SPRKSTEP_DEFAULT_2 = ARKODE_SYMPLECTIC_MCLACHLAN_2;
static const int SPRKSTEP_DEFAULT_3 = ARKODE_SYMPLECTIC_MCLACHLAN_3;
static const int SPRKSTEP_DEFAULT_4 = ARKODE_SYMPLECTIC_MCLACHLAN_4;
static const int SPRKSTEP_DEFAULT_5 = ARKODE_SYMPLECTIC_MCLACHLAN_5;
static const int SPRKSTEP_DEFAULT_6 = ARKODE_SYMPLECTIC_YOSHIDA_6;
static const int SPRKSTEP_DEFAULT_8 = ARKODE_SYMPLECTIC_MCLACHLAN_8;
static const int SPRKSTEP_DEFAULT_10 = ARKODE_SYMPLECTIC_SOFRONIOU_10;

/* -------------------
 * Exported Functions
 * ------------------- */

/* Create, Resize, and Reinitialization functions */
SUNDIALS_EXPORT void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2,
                                     realtype t0, N_Vector y0,
                                     SUNContext sunctx);

SUNDIALS_EXPORT int SPRKStepResize(void *arkode_mem, N_Vector ynew,
                                   realtype hscale, realtype t0,
                                   ARKVecResizeFn resize,
                                   void *resize_data);

SUNDIALS_EXPORT int SPRKStepReInit(void* arkode_mem, ARKRhsFn f1, ARKRhsFn f2, realtype t0, N_Vector y0);

SUNDIALS_EXPORT int SPRKStepReset(void* arkode_mem, realtype tR, N_Vector yR);

/* Tolerance input functions */
SUNDIALS_EXPORT int SPRKStepSStolerances(void *arkode_mem,
                                        realtype reltol,
                                        realtype abstol);
SUNDIALS_EXPORT int SPRKStepSVtolerances(void *arkode_mem,
                                        realtype reltol,
                                        N_Vector abstol);
SUNDIALS_EXPORT int SPRKStepWFtolerances(void *arkode_mem,
                                        ARKEwtFn efun);

/* Residual tolerance input functions */
SUNDIALS_EXPORT int SPRKStepResStolerance(void *arkode_mem,
                                         realtype rabstol);
SUNDIALS_EXPORT int SPRKStepResVtolerance(void *arkode_mem,
                                         N_Vector rabstol);
SUNDIALS_EXPORT int SPRKStepResFtolerance(void *arkode_mem,
                                         ARKRwtFn rfun);

// SUNDIALS_EXPORT int SPRKStepSetMassLinearSolver(void *arkode_mem,
//                                                SUNLinearSolver LS,
//                                                SUNMatrix M,
//                                                booleantype time_dep);

/* Optional input functions -- must be called AFTER SPRKStepCreate */
SUNDIALS_EXPORT int SPRKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int SPRKStepSetOptimalParams(void *arkode_mem);
SUNDIALS_EXPORT int SPRKStepSetUseCompSums(void *arkode_mem, sunbooleantype onoff);
SUNDIALS_EXPORT int SPRKStepSetMethod(void *arkode_mem, ARKODE_SPRKMethodID id);
SUNDIALS_EXPORT int SPRKStepSetOrder(void *arkode_mem, int maxord);
SUNDIALS_EXPORT int SPRKStepSetInterpolantType(void *arkode_mem, int itype);
SUNDIALS_EXPORT int SPRKStepSetInterpolantDegree(void *arkode_mem, int degree);
SUNDIALS_EXPORT int SPRKStepSetDenseOrder(void *arkode_mem, int dord);

SUNDIALS_EXPORT int SPRKStepSetSafetyFactor(void *arkode_mem,
                                            realtype safety);
SUNDIALS_EXPORT int SPRKStepSetErrorBias(void *arkode_mem,
                                         realtype bias);
SUNDIALS_EXPORT int SPRKStepSetMaxGrowth(void *arkode_mem,
                                         realtype mx_growth);
SUNDIALS_EXPORT int SPRKStepSetMinReduction(void *arkode_mem,
                                            realtype eta_min);
SUNDIALS_EXPORT int SPRKStepSetFixedStepBounds(void *arkode_mem,
                                               realtype lb, realtype ub);
SUNDIALS_EXPORT int SPRKStepSetAdaptivityMethod(void *arkode_mem,
                                                int imethod,
                                                int idefault, int pq,
                                                realtype adapt_params[3]);
SUNDIALS_EXPORT int SPRKStepSetAdaptivityFn(void *arkode_mem,
                                            ARKAdaptFn hfun,
                                            void *h_data);
SUNDIALS_EXPORT int SPRKStepSetMaxFirstGrowth(void *arkode_mem,
                                             realtype etamx1);
SUNDIALS_EXPORT int SPRKStepSetMaxEFailGrowth(void *arkode_mem,
                                             realtype etamxf);
SUNDIALS_EXPORT int SPRKStepSetSmallNumEFails(void *arkode_mem,
                                             int small_nef);
SUNDIALS_EXPORT int SPRKStepSetMaxCFailGrowth(void *arkode_mem,
                                             realtype etacf);
SUNDIALS_EXPORT int SPRKStepSetNonlinCRDown(void *arkode_mem,
                                           realtype crdown);
SUNDIALS_EXPORT int SPRKStepSetNonlinRDiv(void *arkode_mem,
                                         realtype rdiv);
SUNDIALS_EXPORT int SPRKStepSetStabilityFn(void *arkode_mem,
                                          ARKExpStabFn EStab,
                                          void *estab_data);
SUNDIALS_EXPORT int SPRKStepSetMaxErrTestFails(void *arkode_mem,
                                              int maxnef);
SUNDIALS_EXPORT int SPRKStepSetMaxNumSteps(void *arkode_mem,
                                          long int mxsteps);
SUNDIALS_EXPORT int SPRKStepSetMaxHnilWarns(void *arkode_mem,
                                           int mxhnil);
SUNDIALS_EXPORT int SPRKStepSetInitStep(void *arkode_mem,
                                       realtype hin);
SUNDIALS_EXPORT int SPRKStepSetMinStep(void *arkode_mem,
                                      realtype hmin);
SUNDIALS_EXPORT int SPRKStepSetMaxStep(void *arkode_mem,
                                      realtype hmax);
SUNDIALS_EXPORT int SPRKStepSetStopTime(void *arkode_mem,
                                       realtype tstop);
SUNDIALS_EXPORT int SPRKStepSetFixedStep(void *arkode_mem,
                                        realtype hfixed);
SUNDIALS_EXPORT int SPRKStepSetErrHandlerFn(void *arkode_mem,
                                           ARKErrHandlerFn ehfun,
                                           void *eh_data);
SUNDIALS_EXPORT int SPRKStepSetErrFile(void *arkode_mem,
                                      FILE *errfp);
SUNDIALS_EXPORT int SPRKStepSetUserData(void *arkode_mem,
                                       void *user_data);

SUNDIALS_EXPORT int SPRKStepSetPostprocessStepFn(void *arkode_mem,
                                                ARKPostProcessFn ProcessStep);
SUNDIALS_EXPORT int SPRKStepSetPostprocessStageFn(void *arkode_mem,
                                                 ARKPostProcessFn ProcessStage);

/* Integrate the ODE over an interval in t */
SUNDIALS_EXPORT int SPRKStepEvolve(void *arkode_mem, realtype tout,
                                  N_Vector yout, realtype *tret,
                                  int itask);

/* Computes the kth derivative of the y function at time t */
SUNDIALS_EXPORT int SPRKStepGetDky(void *arkode_mem, realtype t,
                                  int k, N_Vector dky);

/* Optional output functions */
SUNDIALS_EXPORT int SPRKStepGetNumExpSteps(void *arkode_mem,
                                          long int *expsteps);
SUNDIALS_EXPORT int SPRKStepGetNumAccSteps(void *arkode_mem,
                                          long int *accsteps);
SUNDIALS_EXPORT int SPRKStepGetNumStepAttempts(void *arkode_mem,
                                              long int *step_attempts);
SUNDIALS_EXPORT int SPRKStepGetNumRhsEvals(void *arkode_mem, long int* nf1, long int* nf2);
SUNDIALS_EXPORT int SPRKStepGetNumErrTestFails(void *arkode_mem,
                                               long int *netfails);
SUNDIALS_EXPORT int SPRKStepGetEstLocalErrors(void *arkode_mem,
                                              N_Vector ele);
SUNDIALS_EXPORT int SPRKStepGetWorkSpace(void *arkode_mem,
                                        long int *lenrw,
                                        long int *leniw);
SUNDIALS_EXPORT int SPRKStepGetNumSteps(void *arkode_mem,
                                       long int *nsteps);
SUNDIALS_EXPORT int SPRKStepGetActualInitStep(void *arkode_mem,
                                              realtype *hinused);
SUNDIALS_EXPORT int SPRKStepGetLastStep(void *arkode_mem,
                                        realtype *hlast);
SUNDIALS_EXPORT int SPRKStepGetCurrentStep(void *arkode_mem,
                                           realtype *hcur);
SUNDIALS_EXPORT int SPRKStepGetCurrentTime(void *arkode_mem,
                                           realtype *tcur);
SUNDIALS_EXPORT int SPRKStepGetCurrentState(void *arkode_mem,
                                            N_Vector *state);
SUNDIALS_EXPORT int SPRKStepGetCurrentGamma(void *arkode_mem,
                                            realtype *gamma);
SUNDIALS_EXPORT int SPRKStepGetCurrentMassMatrix(void *arkode_mem,
                                                 SUNMatrix *M);
SUNDIALS_EXPORT int SPRKStepGetTolScaleFactor(void *arkode_mem,
                                              realtype *tolsfac);
SUNDIALS_EXPORT int SPRKStepGetErrWeights(void *arkode_mem,
                                          N_Vector eweight);
SUNDIALS_EXPORT int SPRKStepGetResWeights(void *arkode_mem,
                                          N_Vector rweight);
SUNDIALS_EXPORT int SPRKStepGetUserData(void *arkode_mem,
                                        void **user_data);
SUNDIALS_EXPORT int SPRKStepPrintAllStats(void *arkode_mem, FILE *outfile,
                                          SUNOutputFormat fmt);
SUNDIALS_EXPORT char *SPRKStepGetReturnFlagName(long int flag);
SUNDIALS_EXPORT int SPRKStepWriteParameters(void *arkode_mem, FILE *fp);


/* Grouped optional output functions */
SUNDIALS_EXPORT int SPRKStepGetTimestepperStats(void *arkode_mem,
                                               long int *expsteps,
                                               long int *accsteps,
                                               long int *step_attempts,
                                               long int *nf1, long int *nf2,
                                               long int *nlinsetups,
                                               long int *netfails);
SUNDIALS_EXPORT int SPRKStepGetStepStats(void *arkode_mem,
                                        long int *nsteps,
                                        realtype *hinused,
                                        realtype *hlast,
                                        realtype *hcur,
                                        realtype *tcur);

SUNDIALS_EXPORT int SPRKStepGetNumStepSolveFails(void *arkode_mem,
                                                long int *nncfails);

/* Free function */
SUNDIALS_EXPORT void SPRKStepFree(void **arkode_mem);

/* Output the SPRKStep memory structure (useful when debugging) */
SUNDIALS_EXPORT void SPRKStepPrintMem(void* arkode_mem, FILE* outfile);

/* MRIStep interface functions */
SUNDIALS_EXPORT int SPRKStepCreateMRIStepInnerStepper(void *arkode_mem,
                                                      MRIStepInnerStepper *stepper);

#ifdef __cplusplus
}
#endif

#endif
