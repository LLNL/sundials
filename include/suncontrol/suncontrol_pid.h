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
 * This is the header file for the SUNControl_PID module.
 * -----------------------------------------------------------------*/

#ifndef _SUNCONTROL_PID_H
#define _SUNCONTROL_PID_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------
 * PID implementation of SUNControl
 * -------------------------------- */

struct _SUNControlContent_PID {
  realtype k1;        /* internal controller parameters */
  realtype k2;
  realtype k3;
  realtype bias;      /* error bias factor */
  realtype ep;        /* error from previous step */
  realtype epp;       /* error from 2 steps ago */
  int p;              /* order of accuracy to use for controller */
  sunbooleantype pq;  /* p is embedding order (FALSE) or method order (TRUE) */
};

typedef struct _SUNControlContent_PID *SUNControlContent_PID;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNControl SUNControlPID(SUNContext sunctx);
SUNDIALS_EXPORT
int SUNControlPID_SetParams(SUNControl C, sunbooleantype pq,
                            realtype k1, realtype k2, realtype k3);
SUNDIALS_EXPORT
SUNControl_Type SUNControlGetType_PID(SUNControl C);
SUNDIALS_EXPORT
int SUNControlEstimateStep_PID(SUNControl C, realtype h,
                               realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNControlReset_PID(SUNControl C);
SUNDIALS_EXPORT
int SUNControlSetDefaults_PID(SUNControl C);
SUNDIALS_EXPORT
int SUNControlWrite_PID(SUNControl C, FILE* fptr);
SUNDIALS_EXPORT
int SUNControlSetMethodOrder_PID(SUNControl C, int q);
SUNDIALS_EXPORT
int SUNControlSetEmbeddingOrder_PID(SUNControl C, int p);
SUNDIALS_EXPORT
int SUNControlSetErrorBias_PID(SUNControl C, realtype bias);
SUNDIALS_EXPORT
int SUNControlUpdate_PID(SUNControl C, realtype h, realtype dsm);
SUNDIALS_EXPORT
int SUNControlSpace_PID(SUNControl C, long int *lenrw,
                        long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNCONTROL_PID_H */
