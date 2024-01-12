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
 * This is the header file for the SUNControl_MRILL module.
 * -----------------------------------------------------------------*/

#ifndef _SUNCONTROL_MRILL_H
#define _SUNCONTROL_MRILL_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------------
 * MRI linear-linear implementation of SUNControl
 * ---------------------------------------------- */

struct _SUNControlContent_MRILL {
  realtype k11;    /* internal controller parameters */
  realtype k12;
  realtype k21;
  realtype k22;
  realtype bias;   /* error bias factor */
  realtype esp;    /* slow error from previous step */
  realtype efp;    /* fast error from previous step */
  realtype hsp;    /* slow previous step size */
  realtype hfp;    /* fast previous step size */
  int P;           /* slow order of accuracy to use for controller */
  int p;           /* fast order of accuracy to use for controller */
  sunbooleantype firststep; /* flag indicating first step */
};

typedef struct _SUNControlContent_MRILL *SUNControlContent_MRILL;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNControl SUNControlMRILL(SUNContext sunctx, int P, int p);
SUNDIALS_EXPORT
int SUNControlMRILL_SetParams(SUNControl C, realtype k11, realtype k12,
                              realtype k21, realtype k22);
SUNDIALS_EXPORT
SUNControl_ID SUNControlGetID_MRILL(SUNControl C);
SUNDIALS_EXPORT
int SUNControlEstimateMRISteps_MRILL(SUNControl C, realtype H, realtype h,
                                     realtype DSM, realtype dsm,
                                     realtype* Hnew, realtype *hnew);
SUNDIALS_EXPORT
int SUNControlReset_MRILL(SUNControl C);
SUNDIALS_EXPORT
int SUNControlSetDefaults_MRILL(SUNControl C);
SUNDIALS_EXPORT
int SUNControlWrite_MRILL(SUNControl C, FILE* fptr);
SUNDIALS_EXPORT
int SUNControlSetErrorBias_MRILL(SUNControl C, realtype bias);
SUNDIALS_EXPORT
int SUNControlUpdateMRIH_MRILL(SUNControl C, realtype H, realtype h,
                               realtype DSM, realtype dsm);
SUNDIALS_EXPORT
int SUNControlSpace_MRILL(SUNControl C, long int *lenrw,
                          long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNCONTROL_MRILL_H */
