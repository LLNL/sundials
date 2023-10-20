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
 * This is the header file for the SUNAdaptController_PI module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_PI_H
#define _SUNADAPTCONTROLLER_PI_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------
 * PI implementation of SUNAdaptController
 * --------------------------------------- */

struct _SUNAdaptControllerContent_PI {
  sunrealtype k1;     /* internal controller parameters */
  sunrealtype k2;
  sunrealtype bias;   /* error bias factor */
  sunrealtype ep;     /* error from previous step */
};

typedef struct _SUNAdaptControllerContent_PI *SUNAdaptControllerContent_PI;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_PI(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_PI(SUNAdaptController C,
                                    sunrealtype k1, sunrealtype k2);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_PI(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStep_PI(SUNAdaptController C, sunrealtype h,
                                       int p, sunrealtype dsm, sunrealtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_Reset_PI(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_PI(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_PI(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_PI(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_Update_PI(SUNAdaptController C, sunrealtype h, sunrealtype dsm);
SUNDIALS_EXPORT
int SUNAdaptController_Space_PI(SUNAdaptController C, long int *lenrw,
                                long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_PI_H */
