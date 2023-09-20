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
 * This is the header file for the SUNAdaptController_ImExGus module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_IMEXGUS_H
#define _SUNADAPTCONTROLLER_IMEXGUS_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------------------
 * ImEx Gustafsson implementation of SUNAdaptController
 * ---------------------------------------------------- */

struct SUNAdaptControllerContent_ImExGus_ {
  realtype k1i;             /* internal controller parameters */
  realtype k2i;
  realtype k1e;
  realtype k2e;
  realtype bias;            /* error bias factor */
  realtype ep;              /* error from previous step */
  realtype hp;              /* previous step size */
  int p;                    /* order of accuracy to use for controller */
  sunbooleantype pq;        /* p is embedding order (FALSE) or method order (TRUE) */
  sunbooleantype firststep; /* flag indicating first step */
};

typedef struct SUNAdaptControllerContent_ImExGus_ *SUNAdaptControllerContent_ImExGus;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptControllerImExGus(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNAdaptControllerImExGus_SetParams(SUNAdaptController C, sunbooleantype pq,
                                        realtype k1e, realtype k2e,
                                        realtype k1i, realtype k2i);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptControllerGetType_ImExGus(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptControllerEstimateStep_ImExGus(SUNAdaptController C, realtype h,
                                           realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptControllerReset_ImExGus(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptControllerSetDefaults_ImExGus(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptControllerWrite_ImExGus(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptControllerSetMethodOrder_ImExGus(SUNAdaptController C, int q);
SUNDIALS_EXPORT
int SUNAdaptControllerSetEmbeddingOrder_ImExGus(SUNAdaptController C, int p);
SUNDIALS_EXPORT
int SUNAdaptControllerSetErrorBias_ImExGus(SUNAdaptController C, realtype bias);
SUNDIALS_EXPORT
int SUNAdaptControllerUpdate_ImExGus(SUNAdaptController C, realtype h, realtype dsm);
SUNDIALS_EXPORT
int SUNAdaptControllerSpace_ImExGus(SUNAdaptController C, long int *lenrw,
                                    long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_IMPGUS_H */
