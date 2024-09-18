/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds and Mustafa Aggul @ SMU
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
 * This is the header file for the ARKODE LSRKStep module.
 * -----------------------------------------------------------------*/

#ifndef _LSRKSTEP_H
#define _LSRKSTEP_H

#include <arkode/arkode.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef int (*ARKDomEigFn)(sunrealtype t, N_Vector y, sunrealtype* lambdaR,
                           sunrealtype* lambdaI, void* user_data);

/* ------------------
 * LSRKStep Constants
 * ------------------ */

typedef enum
{
  ARKODE_LSRK_RKC_2   = 1, /* ensure enum is int */
  ARKODE_LSRK_RKL_2   = 2,
  ARKODE_LSRK_RKG_2   = 3,
  ARKODE_LSRK_SSPs_2  = 4,
  ARKODE_LSRK_SSPs_3  = 5,
  ARKODE_LSRK_SSP10_4 = 6
} ARKODE_LSRKMethodType;

/* -------------------
 * Exported Functions
 * ------------------- */

/* Creation and Reinitialization functions */

SUNDIALS_EXPORT void* LSRKStepCreate(ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0,
                                     N_Vector y0, SUNContext sunctx);

SUNDIALS_EXPORT int LSRKStepReInit(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi,
                                   sunrealtype t0, N_Vector y0);

/* Optional input functions -- must be called AFTER a creation routine above */

SUNDIALS_EXPORT int LSRKStepSetMethod(void* arkode_mem,
                                      ARKODE_LSRKMethodType method);

SUNDIALS_EXPORT int LSRKStepSetDomEigFn(void* arkode_mem, ARKDomEigFn dom_eig);

SUNDIALS_EXPORT int LSRKStepSetDomEigFrequency(void* arkode_mem, int nsteps);

SUNDIALS_EXPORT int LSRKStepSetMaxNumStages(void* arkode_mem, int stagemaxlimit);

SUNDIALS_EXPORT int LSRKStepSetDomEigSafetyFactor(void* arkode_mem,
                                                  sunrealtype domeigsfty);

SUNDIALS_EXPORT int LSRKStepSetSSPStageNum(void* arkode_mem, int numofstages);

/* Optional output functions */

SUNDIALS_EXPORT int LSRKStepGetNumRhsEvals(void* arkode_mem, long int* fe_evals,
                                           long int* fi_evals);

SUNDIALS_EXPORT int LSRKStepGetNumDomEigUpdates(void* arkode_mem,
                                                long int* ndomeigupdates);

SUNDIALS_EXPORT int LSRKStepGetMaxNumStages(void* arkode_mem, int* stagemax);

SUNDIALS_EXPORT int LSRKStepGetAverageStageNum(void* arkode_mem,
                                               sunrealtype* averstage);

#ifdef __cplusplus
}
#endif

#endif
