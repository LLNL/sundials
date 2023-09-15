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
 * This is the header file for the SUNControl_I module.
 * -----------------------------------------------------------------*/

#ifndef _SUNCONTROL_I_H
#define _SUNCONTROL_I_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ------------------------------
 * I implementation of SUNControl
 * ------------------------------ */

struct _SUNControlContent_I {
  realtype k1;        /* internal controller parameters */
  realtype bias;      /* error bias factor */
  int p;              /* order of accuracy to use for controller */
  sunbooleantype pq;  /* p is embedding order (FALSE) or method order (TRUE) */
};

typedef struct _SUNControlContent_I *SUNControlContent_I;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNControl SUNControlI(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNControlI_SetParams(SUNControl C, sunbooleantype pq,
                          realtype k1);
SUNDIALS_EXPORT
SUNControl_Type SUNControlGetType_I(SUNControl C);
SUNDIALS_EXPORT
int SUNControlEstimateStep_I(SUNControl C, realtype h,
                             realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNControlSetDefaults_I(SUNControl C);
SUNDIALS_EXPORT
int SUNControlWrite_I(SUNControl C, FILE* fptr);
SUNDIALS_EXPORT
int SUNControlSetMethodOrder_I(SUNControl C, int q);
SUNDIALS_EXPORT
int SUNControlSetEmbeddingOrder_I(SUNControl C, int p);
SUNDIALS_EXPORT
int SUNControlSetErrorBias_I(SUNControl C, realtype bias);
SUNDIALS_EXPORT
int SUNControlSpace_I(SUNControl C, long int *lenrw,
                      long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNCONTROL_I_H */
