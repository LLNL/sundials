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
 * This is the header file for the fixed step implementation of the
 * SUNControl module.
 * -----------------------------------------------------------------*/

#ifndef _SUNCONTROL_FIXED_H
#define _SUNCONTROL_FIXED_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------
 * Fixed-step implementation of SUNControl
 * --------------------------------------- */

struct _SUNControlContent_Fixed {
  realtype hfixed;
};

typedef struct _SUNControlContent_Fixed *SUNControlContent_Fixed;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNControl SUNControlFixed(realtype hfixed, SUNContext sunctx);
SUNDIALS_EXPORT
SUNControl_ID SUNControlGetID_Fixed(SUNControl C);
SUNDIALS_EXPORT
void SUNControlDestroy_Fixed(SUNControl C);
SUNDIALS_EXPORT
int SUNControlEstimateStep_Fixed(SUNControl C, realtype h,
                                 realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNControlWrite_Fixed(SUNControl C, FILE* fptr);
SUNDIALS_EXPORT
int SUNControlSpace_Fixed(SUNControl C, long int *lenrw,
                          long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNCONTROL_FIXED_H */
