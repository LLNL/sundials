/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This is the header file for the ARKode SPRKStep module.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_SPRKSTEP_H
#define _ARKODE_SPRKSTEP_H

#include <arkode/arkode.h>
#include <arkode/arkode_sprk.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * SPRKStep Constants
 * ----------------- */

static const int SPRKSTEP_DEFAULT_1  = ARKODE_SPRK_EULER_1_1;
static const int SPRKSTEP_DEFAULT_2  = ARKODE_SPRK_LEAPFROG_2_2;
static const int SPRKSTEP_DEFAULT_3  = ARKODE_SPRK_MCLACHLAN_3_3;
static const int SPRKSTEP_DEFAULT_4  = ARKODE_SPRK_MCLACHLAN_4_4;
static const int SPRKSTEP_DEFAULT_5  = ARKODE_SPRK_MCLACHLAN_5_6;
static const int SPRKSTEP_DEFAULT_6  = ARKODE_SPRK_YOSHIDA_6_8;
static const int SPRKSTEP_DEFAULT_8  = ARKODE_SPRK_SUZUKI_UMENO_8_16;
static const int SPRKSTEP_DEFAULT_10 = ARKODE_SPRK_SOFRONIOU_10_36;

/* -------------------
 * Exported Functions
 * ------------------- */

/* Create, Resize, and Reinitialization functions */
SUNDIALS_EXPORT void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0,
                                     N_Vector y0, SUNContext sunctx);

SUNDIALS_EXPORT int SPRKStepReInit(void* arkode_mem, ARKRhsFn f1, ARKRhsFn f2,
                                   sunrealtype t0, N_Vector y0);

SUNDIALS_EXPORT int SPRKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR);

/* Rootfinding functions */

/* Rootfinding initialization */
SUNDIALS_EXPORT int SPRKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g);

/* Optional input functions -- must be called AFTER SPRKStepCreate */
SUNDIALS_EXPORT int SPRKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int SPRKStepSetUseCompensatedSums(void* arkode_mem,
                                                  sunbooleantype onoff);
SUNDIALS_EXPORT int SPRKStepSetMethod(void* arkode_mem,
                                      ARKodeSPRKTable sprk_storage);
SUNDIALS_EXPORT int SPRKStepSetMethodName(void* arkode_mem, const char* method);
SUNDIALS_EXPORT int SPRKStepSetOrder(void* arkode_mem, int maxord);
SUNDIALS_EXPORT int SPRKStepSetInterpolantType(void* arkode_mem, int itype);
SUNDIALS_EXPORT int SPRKStepSetInterpolantDegree(void* arkode_mem, int degree);
SUNDIALS_EXPORT int SPRKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps);
SUNDIALS_EXPORT int SPRKStepSetStopTime(void* arkode_mem, sunrealtype tstop);
SUNDIALS_EXPORT int SPRKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed);
SUNDIALS_EXPORT int SPRKStepSetErrHandlerFn(void* arkode_mem,
                                            ARKErrHandlerFn ehfun, void* eh_data);
SUNDIALS_EXPORT int SPRKStepSetErrFile(void* arkode_mem, FILE* errfp);
SUNDIALS_EXPORT int SPRKStepSetUserData(void* arkode_mem, void* user_data);

SUNDIALS_EXPORT int SPRKStepSetPostprocessStepFn(void* arkode_mem,
                                                 ARKPostProcessFn ProcessStep);
SUNDIALS_EXPORT int SPRKStepSetPostprocessStageFn(void* arkode_mem,
                                                  ARKPostProcessFn ProcessStage);

/* Integrate the ODE over an interval in t */
SUNDIALS_EXPORT int SPRKStepEvolve(void* arkode_mem, sunrealtype tout,
                                   N_Vector yout, sunrealtype* tret, int itask);

/* Computes the kth derivative of the y function at time t */
SUNDIALS_EXPORT int SPRKStepGetDky(void* arkode_mem, sunrealtype t, int k,
                                   N_Vector dky);

/* Optional output functions */
SUNDIALS_EXPORT char* SPRKStepGetReturnFlagName(long int flag);
SUNDIALS_EXPORT int SPRKStepGetCurrentMethod(void* arkode_mem,
                                             ARKodeSPRKTable* sprk_storage);
SUNDIALS_EXPORT int SPRKStepGetCurrentState(void* arkode_mem, N_Vector* state);
SUNDIALS_EXPORT int SPRKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur);
SUNDIALS_EXPORT int SPRKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur);
SUNDIALS_EXPORT int SPRKStepGetLastStep(void* arkode_mem, sunrealtype* hlast);
SUNDIALS_EXPORT int SPRKStepGetNumRhsEvals(void* arkode_mem, long int* nf1,
                                           long int* nf2);
SUNDIALS_EXPORT int SPRKStepGetNumStepAttempts(void* arkode_mem,
                                               long int* step_attempts);
SUNDIALS_EXPORT int SPRKStepGetNumSteps(void* arkode_mem, long int* nsteps);
SUNDIALS_EXPORT int SPRKStepGetRootInfo(void* arkode_mem, int* rootsfound);
SUNDIALS_EXPORT int SPRKStepGetUserData(void* arkode_mem, void** user_data);
SUNDIALS_EXPORT int SPRKStepPrintAllStats(void* arkode_mem, FILE* outfile,
                                          SUNOutputFormat fmt);
SUNDIALS_EXPORT int SPRKStepWriteParameters(void* arkode_mem, FILE* fp);

/* Grouped optional output functions */
SUNDIALS_EXPORT int SPRKStepGetStepStats(void* arkode_mem, long int* nsteps,
                                         sunrealtype* hinused, sunrealtype* hlast,
                                         sunrealtype* hcur, sunrealtype* tcur);

/* Free function */
SUNDIALS_EXPORT void SPRKStepFree(void** arkode_mem);

#ifdef __cplusplus
}
#endif

#endif
