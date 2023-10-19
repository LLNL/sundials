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
 * This is the header file for the SUNAdaptController_I module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_I_H
#define _SUNADAPTCONTROLLER_I_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------
 * I implementation of SUNAdaptController
 * -------------------------------------- */

struct _SUNAdaptControllerContent_I {
  sunrealtype k1;    /* internal controller parameters */
  sunrealtype bias;  /* error bias factor */
  int p;             /* method/embedding order of accuracy */
};

typedef struct _SUNAdaptControllerContent_I *SUNAdaptControllerContent_I;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_I(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_I(SUNAdaptController C, sunrealtype k1);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_I(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStep_I(SUNAdaptController C, sunrealtype h,
                                      sunrealtype dsm, sunrealtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_I(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_I(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetMethodOrder_I(SUNAdaptController C, int p);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_I(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_Space_I(SUNAdaptController C, long int *lenrw,
                               long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_I_H */
