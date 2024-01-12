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
 * This is the header file for the SUNControl_MRICC module.
 * -----------------------------------------------------------------*/

#ifndef _SUNCONTROL_MRICC_H
#define _SUNCONTROL_MRICC_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------------
 * MRI constant-constant implementation of SUNControl
 * -------------------------------------------------- */

struct _SUNControlContent_MRICC {
  realtype k1;     /* internal controller parameters */
  realtype k2;
  realtype bias;   /* error bias factor */
  int P;           /* slow order of accuracy to use for controller */
  int p;           /* fast order of accuracy to use for controller */
};

typedef struct _SUNControlContent_MRICC *SUNControlContent_MRICC;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNControl SUNControlMRICC(SUNContext sunctx, int P, int p);
SUNDIALS_EXPORT
int SUNControlMRICC_SetParams(SUNControl C, realtype k1, realtype k2);
SUNDIALS_EXPORT
SUNControl_ID SUNControlGetID_MRICC(SUNControl C);
SUNDIALS_EXPORT
int SUNControlEstimateMRISteps_MRICC(SUNControl C, realtype H, realtype h,
                                     realtype DSM, realtype dsm,
                                     realtype* Hnew, realtype *hnew);
SUNDIALS_EXPORT
int SUNControlSetDefaults_MRICC(SUNControl C);
SUNDIALS_EXPORT
int SUNControlWrite_MRICC(SUNControl C, FILE* fptr);
SUNDIALS_EXPORT
int SUNControlSetErrorBias_MRICC(SUNControl C, realtype bias);
SUNDIALS_EXPORT
int SUNControlSpace_MRICC(SUNControl C, long int *lenrw,
                          long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNCONTROL_MRICC_H */
