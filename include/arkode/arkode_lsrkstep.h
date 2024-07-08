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
 * This is the header file for the ARKODE LSRKStep module.
 * -----------------------------------------------------------------*/

#ifndef _LSRKSTEP_H
#define _LSRKSTEP_H

#include <arkode/arkode.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef int (*ARKSprFn)(sunrealtype t, sunrealtype* extsprad, 
                                   void* user_data);
                                   
/* ------------------
 * LSRKStep Constants
 * ------------------ */

/* -------------------
 * Exported Functions
 * ------------------- */

/* Creation and Reinitialization functions */

SUNDIALS_EXPORT void* LSRKStepCreate(ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0, N_Vector y0,
                                     SUNContext sunctx);

SUNDIALS_EXPORT int LSRKodeSetSprRadFn(void* arkode_mem, ARKSprFn spr);

SUNDIALS_EXPORT int LSRKodeSetConstJac(void* arkode_mem);

SUNDIALS_EXPORT int LSRKodeSetSprRadFrequency(void* arkode_mem, int nsteps);

SUNDIALS_EXPORT int LSRKodeSetMaxStageNum(void* arkode_mem, int stagemaxlimit);

SUNDIALS_EXPORT int LSRKodeSetMaxStepNum(void* arkode_mem, int stepmaxlimit);

SUNDIALS_EXPORT int LSRKodeSetSprRadSafetyFactor(void* arkode_mem, sunrealtype sprsfty);

SUNDIALS_EXPORT int LSRKStepReInit(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0,
                                   N_Vector y0);

/* Optional input functions -- must be called AFTER a creation routine above */

/* SUNDIALS_EXPORT int LSRKStepSetTableName(void* arkode_mem, const char* etable); */

/* Optional output functions */

SUNDIALS_EXPORT int LSRKStepGetNumRhsEvals(void* arkode_mem, long int* nfevals);

/* Grouped optional output functions */
SUNDIALS_EXPORT int LSRKStepGetTimestepperStats(
  void* arkode_mem, long int* expsteps, long int* accsteps,
  long int* step_attempts, long int* nfevals, long int* netfails);

#ifdef __cplusplus
}
#endif

#endif
