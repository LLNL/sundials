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

struct SUNAdaptControllerContent_PID_ {
  realtype k1;        /* internal controller parameters */
  realtype k2;
  realtype k3;
  realtype bias;      /* error bias factor */
  realtype ep;        /* error from previous step */
  realtype epp;       /* error from 2 steps ago */
  int p;              /* order of accuracy to use for controller */
  sunbooleantype pq;  /* p is embedding order (FALSE) or method order (TRUE) */
};

typedef struct SUNAdaptControllerContent_PID_ *SUNAdaptControllerContent_PID;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_PID(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_PID(SUNAdaptController C, sunbooleantype pq,
                                     realtype k1, realtype k2, realtype k3);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_PID(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStep_PID(SUNAdaptController C, realtype h,
                                        realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_Reset_PID(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_PID(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_PID(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetMethodOrder_PID(SUNAdaptController C, int q);
SUNDIALS_EXPORT
int SUNAdaptController_SetEmbeddingOrder_PID(SUNAdaptController C, int p);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_PID(SUNAdaptController C, realtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_Update_PID(SUNAdaptController C, realtype h, realtype dsm);
SUNDIALS_EXPORT
int SUNAdaptController_Space_PID(SUNAdaptController C, long int *lenrw,
                                 long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_PID_H */
