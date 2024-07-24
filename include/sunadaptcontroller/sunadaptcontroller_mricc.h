/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the SUNAdaptController_MRICC module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_MRICC_H
#define _SUNADAPTCONTROLLER_MRICC_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------------
 * MRI constant-constant implementation of SUNAdaptController
 * -------------------------------------------------- */

struct _SUNAdaptControllerContent_MRICC
{
  sunrealtype k1; /* internal controller parameters */
  sunrealtype k2;
  sunrealtype bias; /* error bias factor */
  int p;            /* fast order of accuracy to use */
};

typedef struct _SUNAdaptControllerContent_MRICC* SUNAdaptControllerContent_MRICC;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_MRICC(SUNContext sunctx, int p);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_MRICC(SUNAdaptController C, sunrealtype k1,
                                       sunrealtype k2);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_MRICC(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateMRISteps_MRICC(SUNAdaptController C,
                                              sunrealtype H, sunrealtype h,
                                              int P, sunrealtype DSM,
                                              sunrealtype dsm, sunrealtype* Hnew,
                                              sunrealtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_MRICC(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_MRICC(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_MRICC(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_Space_MRICC(SUNAdaptController C, long int* lenrw,
                                   long int* leniw);

#ifdef __cplusplus
}
#endif

#endif /* _SUNADAPTCONTROLLER_MRICC_H */
