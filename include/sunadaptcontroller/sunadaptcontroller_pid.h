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
 * This is the header file for the SUNAdaptController_PID module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_PID_H
#define _SUNADAPTCONTROLLER_PID_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------
 * PID implementation of SUNAdaptController
 * ---------------------------------------- */

struct _SUNAdaptControllerContent_PID {
  sunrealtype k1;    /* internal controller parameters */
  sunrealtype k2;
  sunrealtype k3;
  sunrealtype bias;  /* error bias factor */
  sunrealtype ep;    /* error from previous step */
  sunrealtype epp;   /* error from 2 steps ago */
  int p;             /* method/embedding order of accuracy */
  int adj;           /* order of accuracy adjustment to use for controller */
};

typedef struct _SUNAdaptControllerContent_PID *SUNAdaptControllerContent_PID;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_PID(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_PID(SUNAdaptController C, sunrealtype k1,
                                     sunrealtype k2, sunrealtype k3);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_PID(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStep_PID(SUNAdaptController C, sunrealtype h,
                                        sunrealtype dsm, sunrealtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_Reset_PID(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_PID(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_PID(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetMethodOrder_PID(SUNAdaptController C, int p);
SUNDIALS_EXPORT
int SUNAdaptController_AdjustControllerOrder_PID(SUNAdaptController C, int adj);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_PID(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_Update_PID(SUNAdaptController C, sunrealtype h, sunrealtype dsm);
SUNDIALS_EXPORT
int SUNAdaptController_Space_PID(SUNAdaptController C, long int *lenrw,
                                 long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_PID_H */
