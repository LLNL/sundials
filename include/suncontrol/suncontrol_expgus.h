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
 * This is the header file for the SUNControl_ExpGus module.
 * -----------------------------------------------------------------*/

#ifndef _SUNCONTROL_EXPGUS_H
#define _SUNCONTROL_EXPGUS_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ------------------------------------------------
 * Explicit Gustafsson implementation of SUNControl
 * ------------------------------------------------ */

struct _SUNControlContent_ExpGus {
  realtype k1;              /* internal controller parameters */
  realtype k2;
  realtype bias;            /* error bias factor */
  realtype ep;              /* error from previous step */
  int p;                    /* order of accuracy to use for controller */
  sunbooleantype pq;        /* p is embedding order (FALSE) or method order (TRUE) */
  sunbooleantype firststep; /* flag indicating first step */
};

typedef struct _SUNControlContent_ExpGus *SUNControlContent_ExpGus;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNControl SUNControlExpGus(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNControlExpGus_SetParams(SUNControl C, sunbooleantype pq,
                               realtype k1, realtype k2);
SUNDIALS_EXPORT
SUNControl_Type SUNControlGetType_ExpGus(SUNControl C);
SUNDIALS_EXPORT
int SUNControlEstimateStep_ExpGus(SUNControl C, realtype h,
                                  realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNControlReset_ExpGus(SUNControl C);
SUNDIALS_EXPORT
int SUNControlSetDefaults_ExpGus(SUNControl C);
SUNDIALS_EXPORT
int SUNControlWrite_ExpGus(SUNControl C, FILE* fptr);
SUNDIALS_EXPORT
int SUNControlSetMethodOrder_ExpGus(SUNControl C, int q);
SUNDIALS_EXPORT
int SUNControlSetEmbeddingOrder_ExpGus(SUNControl C, int p);
SUNDIALS_EXPORT
int SUNControlSetErrorBias_ExpGus(SUNControl C, realtype bias);
SUNDIALS_EXPORT
int SUNControlUpdate_ExpGus(SUNControl C, realtype h, realtype dsm);
SUNDIALS_EXPORT
int SUNControlSpace_ExpGus(SUNControl C, long int *lenrw,
                           long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNCONTROL_EXPGUS_H */
