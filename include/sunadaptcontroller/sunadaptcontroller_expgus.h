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
 * This is the header file for the SUNAdaptController_ExpGus module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_EXPGUS_H
#define _SUNADAPTCONTROLLER_EXPGUS_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------------------
 * Explicit Gustafsson implementation of SUNAdaptController
 * -------------------------------------------------------- */

struct _SUNAdaptControllerContent_ExpGus {
  sunrealtype k1;           /* internal controller parameters */
  sunrealtype k2;
  sunrealtype bias;         /* error bias factor */
  sunrealtype ep;           /* error from previous step */
  sunbooleantype firststep; /* flag indicating first step */
};

typedef struct _SUNAdaptControllerContent_ExpGus *SUNAdaptControllerContent_ExpGus;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_ExpGus(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_ExpGus(SUNAdaptController C,
                                        sunrealtype k1, sunrealtype k2);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_ExpGus(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStep_ExpGus(SUNAdaptController C, sunrealtype h,
                                           int p, sunrealtype dsm, sunrealtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_Reset_ExpGus(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_ExpGus(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_ExpGus(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_ExpGus(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_Update_ExpGus(SUNAdaptController C, sunrealtype h, sunrealtype dsm);
SUNDIALS_EXPORT
int SUNAdaptController_Space_ExpGus(SUNAdaptController C, long int *lenrw,
                                    long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_EXPGUS_H */
