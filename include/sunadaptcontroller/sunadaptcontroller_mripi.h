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
 * This is the header file for the SUNControl_MRIPI module.
 * -----------------------------------------------------------------*/

#ifndef _SUNCONTROL_MRIPI_H
#define _SUNCONTROL_MRIPI_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------
 * MRI PI implementation of SUNControl
 * ----------------------------------- */

struct _SUNControlContent_MRIPI {
  realtype k11;    /* internal controller parameters */
  realtype k12;
  realtype k21;
  realtype k22;
  realtype bias;   /* error bias factor */
  realtype esp;    /* slow error from previous step */
  realtype efp;    /* fast error from previous step */
  int P;           /* slow order of accuracy to use for controller */
  int p;           /* fast order of accuracy to use for controller */
};

typedef struct _SUNControlContent_MRIPI *SUNControlContent_MRIPI;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNControl SUNControlMRIPI(SUNContext sunctx, int P, int p);
SUNDIALS_EXPORT
int SUNControlMRIPI_SetParams(SUNControl C, realtype k11, realtype k12,
                              realtype k21, realtype k22);
SUNDIALS_EXPORT
SUNControl_ID SUNControlGetID_MRIPI(SUNControl C);
SUNDIALS_EXPORT
int SUNControlEstimateMRISteps_MRIPI(SUNControl C, realtype H, realtype h,
                                     realtype DSM, realtype dsm,
                                     realtype* Hnew, realtype *hnew);
SUNDIALS_EXPORT
int SUNControlReset_MRIPI(SUNControl C);
SUNDIALS_EXPORT
int SUNControlSetDefaults_MRIPI(SUNControl C);
SUNDIALS_EXPORT
int SUNControlWrite_MRIPI(SUNControl C, FILE* fptr);
SUNDIALS_EXPORT
int SUNControlSetErrorBias_MRIPI(SUNControl C, realtype bias);
SUNDIALS_EXPORT
int SUNControlUpdateMRIH_MRIPI(SUNControl C, realtype H, realtype h,
                               realtype DSM, realtype dsm);
SUNDIALS_EXPORT
int SUNControlSpace_MRIPI(SUNControl C, long int *lenrw,
                          long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNCONTROL_MRIPI_H */
