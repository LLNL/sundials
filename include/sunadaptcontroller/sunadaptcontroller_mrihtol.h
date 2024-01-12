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
 * This is the header file for the SUNControl_MRIHTol module.
 * -----------------------------------------------------------------*/

#ifndef _SUNCONTROL_MRIHTOL_H
#define _SUNCONTROL_MRIHTOL_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------
 * MRI H+tolerance implementation of SUNControl
 * -------------------------------------------- */

struct _SUNControlContent_MRIHTol {
  SUNControl HControl;
  SUNControl TolControl;
};

typedef struct _SUNControlContent_MRIHTol *SUNControlContent_MRIHTol;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNControl SUNControlMRIHTol(SUNContext sunctx, SUNControl HControl,
                             SUNControl TolControl);
SUNDIALS_EXPORT
SUNControl SUNControlMRIHTol_GetSlowController(SUNControl C);
SUNDIALS_EXPORT
SUNControl SUNControlMRIHTol_GetFastController(SUNControl C);
SUNDIALS_EXPORT
SUNControl_ID SUNControlGetID_MRIHTol(SUNControl C);
SUNDIALS_EXPORT
int SUNControlEstimateStepTol_MRIHTol(SUNControl C, realtype H,
                                      realtype tolfac, realtype DSM,
                                      realtype dsm, realtype *Hnew,
                                      realtype* tolfacnew);
SUNDIALS_EXPORT
int SUNControlReset_MRIHTol(SUNControl C);
SUNDIALS_EXPORT
int SUNControlSetDefaults_MRIHTol(SUNControl C);
SUNDIALS_EXPORT
int SUNControlWrite_MRIHTol(SUNControl C, FILE* fptr);
SUNDIALS_EXPORT
int SUNControlSetMethodOrder_MRIHTol(SUNControl C, int q);
SUNDIALS_EXPORT
int SUNControlSetEmbeddingOrder_MRIHTol(SUNControl C, int p);
SUNDIALS_EXPORT
int SUNControlSetErrorBias_MRIHTol(SUNControl C, realtype bias);
SUNDIALS_EXPORT
int SUNControlUpdateMRITol_MRIHTol(SUNControl C, realtype H, realtype tolfac,
                                   realtype DSM, realtype dsm);
SUNDIALS_EXPORT
int SUNControlSpace_MRIHTol(SUNControl C, long int *lenrw,
                            long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNCONTROL_MRIHTOL_H */
