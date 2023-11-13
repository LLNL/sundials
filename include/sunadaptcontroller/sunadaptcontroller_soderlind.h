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
 * This is the header file for the SUNAdaptController_Soderlind
 * module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_SODERLIND_H
#define _SUNADAPTCONTROLLER_SODERLIND_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------------------
 * Soderlind implementation of SUNAdaptController
 * ---------------------------------------------------- */

struct _SUNAdaptControllerContent_Soderlind {
  sunrealtype k1;        /* internal controller parameters */
  sunrealtype k2;
  sunrealtype k3;
  sunrealtype k4;
  sunrealtype k5;
  sunrealtype bias;      /* error bias factor */
  sunrealtype ep;        /* error from previous step */
  sunrealtype epp;       /* error from 2 steps ago */
  sunrealtype hp;        /* previous step size */
  sunrealtype hpp;       /* step size from 2 steps ago */
  int firststeps;        /* flag to handle first few steps */
};

typedef struct _SUNAdaptControllerContent_Soderlind *SUNAdaptControllerContent_Soderlind;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_Soderlind(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_Soderlind(SUNAdaptController C,
                                           sunrealtype k1, sunrealtype k2,
                                           sunrealtype k3, sunrealtype k4,
                                           sunrealtype k5);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_Soderlind(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStep_Soderlind(SUNAdaptController C, sunrealtype h,
                                              int p, sunrealtype dsm, sunrealtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_Reset_Soderlind(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_Soderlind(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_Soderlind(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_Soderlind(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_UpdateH_Soderlind(SUNAdaptController C, sunrealtype h, sunrealtype dsm);
SUNDIALS_EXPORT
int SUNAdaptController_Space_Soderlind(SUNAdaptController C, long int *lenrw,
                                       long int *leniw);

/* Convenience routines to construct subsidiary controllers */
SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_PID(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_PID(SUNAdaptController C,
                                     sunrealtype k1, sunrealtype k2,
                                     sunrealtype k3);
SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_PI(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_PI(SUNAdaptController C,
                                    sunrealtype k1, sunrealtype k2);
SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_I(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_I(SUNAdaptController C, sunrealtype k1);

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_ExpGus(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_ExpGus(SUNAdaptController C,
                                        sunrealtype k1, sunrealtype k2);

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_ImpGus(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_ImpGus(SUNAdaptController C,
                                        sunrealtype k1, sunrealtype k2);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_SODERLIND_H */
