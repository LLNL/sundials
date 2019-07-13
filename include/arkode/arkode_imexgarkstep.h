/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the ARKode IMEXGARKStep module.
 * -----------------------------------------------------------------*/

#ifndef _IMEXGARKSTEP_H
#define _IMEXGARKSTEP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <arkode/arkode.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_butcher_dirk.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------
 * IMEXGARKStep Constants
 * ---------------------- */

/* -------------------
 * Exported Functions
 * ------------------- */

/* Create, Resize, and Reinitializeation functions */
SUNDIALS_EXPORT void* IMEXGARKStepCreate(ARKRhsFn fe, ARKRhsFn fi,
                                         realtype t0, N_Vector y0);

SUNDIALS_EXPORT int IMEXGARKStepResize(void *arkode_mem, N_Vector ynew,
                                       realtype hscale, realtype t0,
                                       ARKVecResizeFn resize,
                                       void *resize_data);

SUNDIALS_EXPORT int IMEXGARKStepReInit(void* arkode_mem, ARKRhsFn fe,
                                       ARKRhsFn fi, realtype t0, N_Vector y0);

/* Tolerance input functions */
SUNDIALS_EXPORT int IMEXGARKStepSStolerances(void *arkode_mem,
                                             realtype reltol,
                                             realtype abstol);
SUNDIALS_EXPORT int IMEXGARKStepSVtolerances(void *arkode_mem,
                                             realtype reltol,
                                             N_Vector abstol);
SUNDIALS_EXPORT int IMEXGARKStepWFtolerances(void *arkode_mem,
                                             ARKEwtFn efun);

/* Resudal tolerance input functions */
SUNDIALS_EXPORT int IMEXGARKStepResStolerance(void *arkode_mem,
                                              realtype rabstol);
SUNDIALS_EXPORT int IMEXGARKStepResVtolerance(void *arkode_mem,
                                              N_Vector rabstol);
SUNDIALS_EXPORT int IMEXGARKStepResFtolerance(void *arkode_mem,
                                              ARKRwtFn rfun);


/* Linear solver set functions */
SUNDIALS_EXPORT int IMEXGARKStepSetLinearSolver(void *arkode_mem,
                                                SUNLinearSolver LS,
                                                SUNMatrix A);
SUNDIALS_EXPORT int IMEXGARKStepSetMassLinearSolver(void *arkode_mem,
                                                    SUNLinearSolver LS,
                                                    SUNMatrix M,
                                                    booleantype time_dep);

/* Rootfinding initialization */
SUNDIALS_EXPORT int IMEXGARKStepRootInit(void *arkode_mem, int nrtfn,
                                         ARKRootFn g);

/* Optional input functions -- must be called AFTER IMEXGARKStepCreate */
SUNDIALS_EXPORT int IMEXGARKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int IMEXGARKStepSetDenseOrder(void *arkode_mem, int dord);
SUNDIALS_EXPORT int IMEXGARKStepSetNonlinearSolver(void *arkode_mem,
                                                   SUNNonlinearSolver NLS);
SUNDIALS_EXPORT int IMEXGARKStepSetLinear(void *arkode_mem, int timedepend);
SUNDIALS_EXPORT int IMEXGARKStepSetNonlinear(void *arkode_mem);
SUNDIALS_EXPORT int IMEXGARKStepSetTables(void *arkode_mem, int s,
                                          int q, int p,
                                          realtype *ce, realtype *ci,
                                          realtype *Aee, realtype *Aei,
                                          realtype *Aie, realtype *Aii,
                                          realtype *be, realtype *bi,
                                          realtype *de, realtype *di);
SUNDIALS_EXPORT int IMEXGARKStepSetCFLFraction(void *arkode_mem,
                                               realtype cfl_frac);
SUNDIALS_EXPORT int IMEXGARKStepSetSafetyFactor(void *arkode_mem,
                                                realtype safety);
SUNDIALS_EXPORT int IMEXGARKStepSetErrorBias(void *arkode_mem,
                                             realtype bias);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxGrowth(void *arkode_mem,
                                             realtype mx_growth);
SUNDIALS_EXPORT int IMEXGARKStepSetFixedStepBounds(void *arkode_mem,
                                                   realtype lb, realtype ub);
SUNDIALS_EXPORT int IMEXGARKStepSetAdaptivityMethod(void *arkode_mem,
                                                    int imethod,
                                                    int idefault, int pq,
                                                    realtype *adapt_params);
SUNDIALS_EXPORT int IMEXGARKStepSetAdaptivityFn(void *arkode_mem,
                                                ARKAdaptFn hfun,
                                                void *h_data);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxFirstGrowth(void *arkode_mem,
                                                  realtype etamx1);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxEFailGrowth(void *arkode_mem,
                                                  realtype etamxf);
SUNDIALS_EXPORT int IMEXGARKStepSetSmallNumEFails(void *arkode_mem,
                                                  int small_nef);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxCFailGrowth(void *arkode_mem,
                                                  realtype etacf);
SUNDIALS_EXPORT int IMEXGARKStepSetNonlinCRDown(void *arkode_mem,
                                                realtype crdown);
SUNDIALS_EXPORT int IMEXGARKStepSetNonlinRDiv(void *arkode_mem,
                                              realtype rdiv);
SUNDIALS_EXPORT int IMEXGARKStepSetDeltaGammaMax(void *arkode_mem,
                                                 realtype dgmax);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxStepsBetweenLSet(void *arkode_mem,
                                                       int msbp);
SUNDIALS_EXPORT int IMEXGARKStepSetPredictorMethod(void *arkode_mem,
                                                   int method);
SUNDIALS_EXPORT int IMEXGARKStepSetStabilityFn(void *arkode_mem,
                                               ARKExpStabFn EStab,
                                               void *estab_data);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxErrTestFails(void *arkode_mem,
                                                   int maxnef);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxNonlinIters(void *arkode_mem,
                                                  int maxcor);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxConvFails(void *arkode_mem,
                                                int maxncf);
SUNDIALS_EXPORT int IMEXGARKStepSetNonlinConvCoef(void *arkode_mem,
                                                  realtype nlscoef);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxNumSteps(void *arkode_mem,
                                               long int mxsteps);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxHnilWarns(void *arkode_mem,
                                                int mxhnil);
SUNDIALS_EXPORT int IMEXGARKStepSetInitStep(void *arkode_mem,
                                            realtype hin);
SUNDIALS_EXPORT int IMEXGARKStepSetMinStep(void *arkode_mem,
                                           realtype hmin);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxStep(void *arkode_mem,
                                           realtype hmax);
SUNDIALS_EXPORT int IMEXGARKStepSetStopTime(void *arkode_mem,
                                            realtype tstop);
SUNDIALS_EXPORT int IMEXGARKStepSetFixedStep(void *arkode_mem,
                                             realtype hfixed);

SUNDIALS_EXPORT int IMEXGARKStepSetRootDirection(void *arkode_mem,
                                                 int *rootdir);
SUNDIALS_EXPORT int IMEXGARKStepSetNoInactiveRootWarn(void *arkode_mem);

SUNDIALS_EXPORT int IMEXGARKStepSetErrHandlerFn(void *arkode_mem,
                                                ARKErrHandlerFn ehfun,
                                                void *eh_data);
SUNDIALS_EXPORT int IMEXGARKStepSetErrFile(void *arkode_mem,
                                           FILE *errfp);
SUNDIALS_EXPORT int IMEXGARKStepSetUserData(void *arkode_mem,
                                            void *user_data);
SUNDIALS_EXPORT int IMEXGARKStepSetDiagnostics(void *arkode_mem,
                                               FILE *diagfp);

SUNDIALS_EXPORT int IMEXGARKStepSetPostprocessStepFn(void *arkode_mem,
                                                     ARKPostProcessStepFn ProcessStep);

/* Linear solver interface optional input functions -- must be called
   AFTER IMEXGARKStepSetLinearSolver and/or IMEXGARKStepSetMassLinearSolver */
SUNDIALS_EXPORT int IMEXGARKStepSetJacFn(void *arkode_mem, ARKLsJacFn jac);
SUNDIALS_EXPORT int IMEXGARKStepSetMassFn(void *arkode_mem, ARKLsMassFn mass);
SUNDIALS_EXPORT int IMEXGARKStepSetMaxStepsBetweenJac(void *arkode_mem,
                                                      long int msbj);
SUNDIALS_EXPORT int IMEXGARKStepSetEpsLin(void *arkode_mem, realtype eplifac);
SUNDIALS_EXPORT int IMEXGARKStepSetMassEpsLin(void *arkode_mem, realtype eplifac);
SUNDIALS_EXPORT int IMEXGARKStepSetPreconditioner(void *arkode_mem,
                                                  ARKLsPrecSetupFn psetup,
                                                  ARKLsPrecSolveFn psolve);
SUNDIALS_EXPORT int IMEXGARKStepSetMassPreconditioner(void *arkode_mem,
                                                      ARKLsMassPrecSetupFn psetup,
                                                      ARKLsMassPrecSolveFn psolve);
SUNDIALS_EXPORT int IMEXGARKStepSetJacTimes(void *arkode_mem,
                                            ARKLsJacTimesSetupFn jtsetup,
                                            ARKLsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int IMEXGARKStepSetMassTimes(void *arkode_mem,
                                             ARKLsMassTimesSetupFn msetup,
                                             ARKLsMassTimesVecFn mtimes,
                                             void *mtimes_data);
SUNDIALS_EXPORT int IMEXGARKStepSetLinSysFn(void *arkode_mem, ARKLsLinSysFn linsys);

/* Integrate the ODE over an interval in t */
SUNDIALS_EXPORT int IMEXGARKStepEvolve(void *arkode_mem, realtype tout,
                                       N_Vector yout, realtype *tret,
                                       int itask);

/* Computes the kth derivative of the y function at time t */
SUNDIALS_EXPORT int IMEXGARKStepGetDky(void *arkode_mem, realtype t,
                                       int k, N_Vector dky);

/* Optional output functions */
SUNDIALS_EXPORT int IMEXGARKStepGetNumExpSteps(void *arkode_mem,
                                               long int *expsteps);
SUNDIALS_EXPORT int IMEXGARKStepGetNumAccSteps(void *arkode_mem,
                                               long int *accsteps);
SUNDIALS_EXPORT int IMEXGARKStepGetNumStepAttempts(void *arkode_mem,
                                                   long int *step_attempts);
SUNDIALS_EXPORT int IMEXGARKStepGetNumRhsEvals(void *arkode_mem,
                                               long int *nfe_evals,
                                               long int *nfi_evals);
SUNDIALS_EXPORT int IMEXGARKStepGetNumLinSolvSetups(void *arkode_mem,
                                                    long int *nlinsetups);
SUNDIALS_EXPORT int IMEXGARKStepGetNumErrTestFails(void *arkode_mem,
                                                   long int *netfails);
SUNDIALS_EXPORT int IMEXGARKStepGetCurrentButcherTables(void *arkode_mem,
                                                        ARKodeButcherTable *Bee,
                                                        ARKodeButcherTable *Bei,
                                                        ARKodeButcherTable *Bie,
                                                        ARKodeButcherTable *Bii);
SUNDIALS_EXPORT int IMEXGARKStepGetEstLocalErrors(void *arkode_mem,
                                                  N_Vector ele);
SUNDIALS_EXPORT int IMEXGARKStepGetWorkSpace(void *arkode_mem,
                                             long int *lenrw,
                                             long int *leniw);
SUNDIALS_EXPORT int IMEXGARKStepGetNumSteps(void *arkode_mem,
                                            long int *nsteps);
SUNDIALS_EXPORT int IMEXGARKStepGetActualInitStep(void *arkode_mem,
                                                  realtype *hinused);
SUNDIALS_EXPORT int IMEXGARKStepGetLastStep(void *arkode_mem,
                                            realtype *hlast);
SUNDIALS_EXPORT int IMEXGARKStepGetCurrentStep(void *arkode_mem,
                                               realtype *hcur);
SUNDIALS_EXPORT int IMEXGARKStepGetCurrentTime(void *arkode_mem,
                                               realtype *tcur);
SUNDIALS_EXPORT int IMEXGARKStepGetTolScaleFactor(void *arkode_mem,
                                                  realtype *tolsfac);
SUNDIALS_EXPORT int IMEXGARKStepGetErrWeights(void *arkode_mem,
                                              N_Vector eweight);
SUNDIALS_EXPORT int IMEXGARKStepGetResWeights(void *arkode_mem,
                                              N_Vector rweight);
SUNDIALS_EXPORT int IMEXGARKStepGetNumGEvals(void *arkode_mem,
                                             long int *ngevals);
SUNDIALS_EXPORT int IMEXGARKStepGetRootInfo(void *arkode_mem,
                                            int *rootsfound);
SUNDIALS_EXPORT char *IMEXGARKStepGetReturnFlagName(long int flag);

SUNDIALS_EXPORT int IMEXGARKStepWriteParameters(void *arkode_mem, FILE *fp);

SUNDIALS_EXPORT int IMEXGARKStepWriteButcher(void *arkode_mem, FILE *fp);


/* Grouped optional output functions */
SUNDIALS_EXPORT int IMEXGARKStepGetTimestepperStats(void *arkode_mem,
                                                    long int *expsteps,
                                                    long int *accsteps,
                                                    long int *step_attempts,
                                                    long int *nfe_evals,
                                                    long int *nfi_evals,
                                                    long int *nlinsetups,
                                                    long int *netfails);
SUNDIALS_EXPORT int IMEXGARKStepGetStepStats(void *arkode_mem,
                                            long int *nsteps,
                                            realtype *hinused,
                                            realtype *hlast,
                                            realtype *hcur,
                                            realtype *tcur);

/* Nonlinear solver optional output functions */
SUNDIALS_EXPORT int IMEXGARKStepGetNumNonlinSolvIters(void *arkode_mem,
                                                      long int *nniters);
SUNDIALS_EXPORT int IMEXGARKStepGetNumNonlinSolvConvFails(void *arkode_mem,
                                                          long int *nncfails);
SUNDIALS_EXPORT int IMEXGARKStepGetNonlinSolvStats(void *arkode_mem,
                                                   long int *nniters,
                                                   long int *nncfails);

/* Linear solver optional output functions */
SUNDIALS_EXPORT int IMEXGARKStepGetLinWorkSpace(void *arkode_mem,
                                                long int *lenrwLS,
                                                long int *leniwLS);
SUNDIALS_EXPORT int IMEXGARKStepGetNumJacEvals(void *arkode_mem,
                                               long int *njevals);
SUNDIALS_EXPORT int IMEXGARKStepGetNumPrecEvals(void *arkode_mem,
                                                long int *npevals);
SUNDIALS_EXPORT int IMEXGARKStepGetNumPrecSolves(void *arkode_mem,
                                                 long int *npsolves);
SUNDIALS_EXPORT int IMEXGARKStepGetNumLinIters(void *arkode_mem,
                                               long int *nliters);
SUNDIALS_EXPORT int IMEXGARKStepGetNumLinConvFails(void *arkode_mem,
                                                   long int *nlcfails);
SUNDIALS_EXPORT int IMEXGARKStepGetNumJTSetupEvals(void *arkode_mem,
                                                   long int *njtsetups);
SUNDIALS_EXPORT int IMEXGARKStepGetNumJtimesEvals(void *arkode_mem,
                                                  long int *njvevals);
SUNDIALS_EXPORT int IMEXGARKStepGetNumLinRhsEvals(void *arkode_mem,
                                                  long int *nfevalsLS);
SUNDIALS_EXPORT int IMEXGARKStepGetLastLinFlag(void *arkode_mem,
                                               long int *flag);

SUNDIALS_EXPORT int IMEXGARKStepGetMassWorkSpace(void *arkode_mem,
                                                 long int *lenrwMLS,
                                                 long int *leniwMLS);
SUNDIALS_EXPORT int IMEXGARKStepGetNumMassSetups(void *arkode_mem,
                                                 long int *nmsetups);
SUNDIALS_EXPORT int IMEXGARKStepGetNumMassMultSetups(void *arkode_mem,
                                                     long int *nmvsetups);
SUNDIALS_EXPORT int IMEXGARKStepGetNumMassMult(void *arkode_mem,
                                               long int *nmvevals);
SUNDIALS_EXPORT int IMEXGARKStepGetNumMassSolves(void *arkode_mem,
                                                 long int *nmsolves);
SUNDIALS_EXPORT int IMEXGARKStepGetNumMassPrecEvals(void *arkode_mem,
                                                    long int *nmpevals);
SUNDIALS_EXPORT int IMEXGARKStepGetNumMassPrecSolves(void *arkode_mem,
                                                     long int *nmpsolves);
SUNDIALS_EXPORT int IMEXGARKStepGetNumMassIters(void *arkode_mem,
                                                long int *nmiters);
SUNDIALS_EXPORT int IMEXGARKStepGetNumMassConvFails(void *arkode_mem,
                                                    long int *nmcfails);
SUNDIALS_EXPORT int IMEXGARKStepGetNumMTSetups(void *arkode_mem,
                                               long int *nmtsetups);
SUNDIALS_EXPORT int IMEXGARKStepGetLastMassFlag(void *arkode_mem,
                                                long int *flag);

SUNDIALS_EXPORT char *IMEXGARKStepGetLinReturnFlagName(long int flag);


/* Free function */
SUNDIALS_EXPORT void IMEXGARKStepFree(void **arkode_mem);

/* Output the IMEXGARKStep memory structure (useful when debugging) */
SUNDIALS_EXPORT void IMEXGARKStepPrintMem(void* arkode_mem, FILE* outfile);


#ifdef __cplusplus
}
#endif

#endif
