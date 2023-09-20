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

struct SUNAdaptControllerContent_PI_ {
  realtype k1;        /* internal controller parameters */
  realtype k2;
  realtype bias;      /* error bias factor */
  realtype ep;        /* error from previous step */
  int p;              /* order of accuracy to use for controller */
  sunbooleantype pq;  /* p is embedding order (FALSE) or method order (TRUE) */
};

typedef struct SUNAdaptControllerContent_PI_ *SUNAdaptControllerContent_PI;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptControllerPI(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptControllerPI_SetParams(SUNAdaptController C, sunbooleantype pq,
                                   realtype k1, realtype k2);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptControllerGetType_PI(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptControllerEstimateStep_PI(SUNAdaptController C, realtype h,
                                      realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptControllerReset_PI(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptControllerSetDefaults_PI(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptControllerWrite_PI(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptControllerSetMethodOrder_PI(SUNAdaptController C, int q);
SUNDIALS_EXPORT
int SUNAdaptControllerSetEmbeddingOrder_PI(SUNAdaptController C, int p);
SUNDIALS_EXPORT
int SUNAdaptControllerSetErrorBias_PI(SUNAdaptController C, realtype bias);
SUNDIALS_EXPORT
int SUNAdaptControllerUpdate_PI(SUNAdaptController C, realtype h, realtype dsm);
SUNDIALS_EXPORT
int SUNAdaptControllerSpace_PI(SUNAdaptController C, long int *lenrw,
                               long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_PI_H */
