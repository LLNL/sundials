/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * Based on arkode_arkstep.h written by Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------------------
 * This is the header file for ARKode's Implicit-Explicit (IMEX) Generalized
 * Additive Runge Kutta (GARK) time step module.
 * ---------------------------------------------------------------------------*/

#ifndef _IMEXGARKSTEP_H
#define _IMEXGARKSTEP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* =============================================================================
 * IMEXGARKSTEP Constants
 * ===========================================================================*/



/* =============================================================================
 * IMEXGARKSTEP Exported functions
 * ===========================================================================*/

/* Create the IMEX IMEXGARK step module and attaches it to ARKode */
SUNDIALS_EXPORT int IMEXGARKStepCreate(void* arkode_mem, ARKRhsFn fe,
                                   ARKRhsFn fi, realtype t0, N_Vector y0);

/* Re-initialize the ARK time step module */
SUNDIALS_EXPORT int IMEXGARKStepReInit(void* arkode_mem, ARKRhsFn fe,
                                   ARKRhsFn fi, realtype t0, N_Vector y0);

/* Optional input functions (must be called AFTER IMEXGARKStepCreate) */
SUNDIALS_EXPORT int IMEXGARKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int IMEXGARKStepSetLinear(void *arkode_mem, int timedepend);
SUNDIALS_EXPORT int IMEXGARKStepSetNonlinear(void *arkode_mem);
SUNDIALS_EXPORT int IMEXGARKStepSetFixedPoint(void *arkode_mem, long int fp_m);
SUNDIALS_EXPORT int IMEXGARKStepSetNewton(void *arkode_mem);
SUNDIALS_EXPORT int IMEXGARKStepSetButcherTables(void *arkode_mem, int s, 
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
SUNDIALS_EXPORT int IMEXGARKStepGetTimestepperStats(void *arkode_mem, 
                                                long int *expsteps, 
                                                long int *accsteps, 
                                                long int *step_attempts, 
                                                long int *nfe_evals, 
                                                long int *nfi_evals, 
                                                long int *nlinsetups, 
                                                long int *netfails);

/* Nonlinear solver optional output functions */
SUNDIALS_EXPORT int IMEXGARKStepGetNumNonlinSolvIters(void *arkode_mem, 
                                                  long int *nniters);
SUNDIALS_EXPORT int IMEXGARKStepGetNumNonlinSolvConvFails(void *arkode_mem, 
                                                      long int *nncfails);
SUNDIALS_EXPORT int IMEXGARKStepGetNonlinSolvStats(void *arkode_mem, 
                                               long int *nniters,
                                               long int *nncfails);

/* Output all timestepper module parameters to the provided file pointer */
SUNDIALS_EXPORT int IMEXGARKStepWriteParameters(void *arkode_mem, FILE *fp);

/* Output the Butcher tables to the provided file pointer */
SUNDIALS_EXPORT int IMEXGARKStepWriteButcher(void *arkode_mem, FILE *fp);


#ifdef __cplusplus
}
#endif

#endif
