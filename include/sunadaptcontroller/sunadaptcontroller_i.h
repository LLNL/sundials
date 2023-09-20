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

struct SUNAdaptControllerContent_I_ {
  realtype k1;        /* internal controller parameters */
  realtype bias;      /* error bias factor */
  int p;              /* order of accuracy to use for controller */
  sunbooleantype pq;  /* p is embedding order (FALSE) or method order (TRUE) */
};

typedef struct SUNAdaptControllerContent_I_ *SUNAdaptControllerContent_I;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_I(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_I(SUNAdaptController C, sunbooleantype pq,
                                   realtype k1);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_I(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStep_I(SUNAdaptController C, realtype h,
                                      realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_I(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_I(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetMethodOrder_I(SUNAdaptController C, int q);
SUNDIALS_EXPORT
int SUNAdaptController_SetEmbeddingOrder_I(SUNAdaptController C, int p);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_I(SUNAdaptController C, realtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_Space_I(SUNAdaptController C, long int *lenrw,
                               long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_I_H */
