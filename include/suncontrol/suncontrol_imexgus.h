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
 * This is the header file for the SUNControl_ImExGus module.
 * -----------------------------------------------------------------*/

#ifndef _SUNCONTROL_IMEXGUS_H
#define _SUNCONTROL_IMEXGUS_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------
 * ImEx Gustafsson implementation of SUNControl
 * -------------------------------------------- */

struct _SUNControlContent_ImExGus {
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

typedef struct _SUNControlContent_ImExGus *SUNControlContent_ImExGus;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNControl SUNControlImExGus(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNControlImExGus_SetParams(SUNControl C, sunbooleantype pq,
                                realtype k1e, realtype k2e,
                                realtype k1i, realtype k2i);
SUNDIALS_EXPORT
SUNControl_Type SUNControlGetType_ImExGus(SUNControl C);
SUNDIALS_EXPORT
int SUNControlEstimateStep_ImExGus(SUNControl C, realtype h,
                                   realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNControlReset_ImExGus(SUNControl C);
SUNDIALS_EXPORT
int SUNControlSetDefaults_ImExGus(SUNControl C);
SUNDIALS_EXPORT
int SUNControlWrite_ImExGus(SUNControl C, FILE* fptr);
SUNDIALS_EXPORT
int SUNControlSetMethodOrder_ImExGus(SUNControl C, int q);
SUNDIALS_EXPORT
int SUNControlSetEmbeddingOrder_ImExGus(SUNControl C, int p);
SUNDIALS_EXPORT
int SUNControlSetErrorBias_ImExGus(SUNControl C, realtype bias);
SUNDIALS_EXPORT
int SUNControlUpdate_ImExGus(SUNControl C, realtype h, realtype dsm);
SUNDIALS_EXPORT
int SUNControlSpace_ImExGus(SUNControl C, long int *lenrw,
                            long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNCONTROL_IMPGUS_H */
