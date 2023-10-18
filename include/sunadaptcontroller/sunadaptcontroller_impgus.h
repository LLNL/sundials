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
 * This is the header file for the SUNAdaptController_ImpGus module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_IMPGUS_H
#define _SUNADAPTCONTROLLER_IMPGUS_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------------------
 * Implicit Gustafsson implementation of SUNAdaptController
 * -------------------------------------------------------- */

struct _SUNAdaptControllerContent_ImpGus {
  sunrealtype k1;           /* internal controller parameters */
  sunrealtype k2;
  sunrealtype bias;         /* error bias factor */
  sunrealtype ep;           /* error from previous step */
  sunrealtype hp;           /* previous step size */
  int p;                    /* method/embedding order of accuracy */
  int adj;                  /* order of accuracy adjustment to use for controller */
  int pq;                   /* p is order of embedding (0), method (1), or minimum (-1) */
  sunbooleantype firststep; /* flag indicating first step */
};

typedef struct _SUNAdaptControllerContent_ImpGus *SUNAdaptControllerContent_ImpGus;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_ImpGus(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_ImpGus(SUNAdaptController C, int pq,
                                        sunrealtype k1, sunrealtype k2);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_ImpGus(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStep_ImpGus(SUNAdaptController C, sunrealtype h,
                                           sunrealtype dsm, sunrealtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_Reset_ImpGus(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_ImpGus(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_ImpGus(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetMethodOrder_ImpGus(SUNAdaptController C, int p, int q);
SUNDIALS_EXPORT
int SUNAdaptController_AdjustControllerOrder_ImpGus(SUNAdaptController C, int adj);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_ImpGus(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_Update_ImpGus(SUNAdaptController C, sunrealtype h, sunrealtype dsm);
SUNDIALS_EXPORT
int SUNAdaptController_Space_ImpGus(SUNAdaptController C, long int *lenrw,
                                    long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_IMPGUS_H */
