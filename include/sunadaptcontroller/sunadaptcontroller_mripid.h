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
 * This is the header file for the SUNAdaptController_MRIPID module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_MRIPID_H
#define _SUNADAPTCONTROLLER_MRIPID_H

#include <stdio.h>
#include <sundials/sundials_control.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------
 * MRI PI implementation of SUNAdaptController
 * ----------------------------------- */

struct _SUNAdaptControllerContent_MRIPID {
  sunrealtype k11;    /* internal controller parameters */
  sunrealtype k12;
  sunrealtype k13;
  sunrealtype k21;
  sunrealtype k22;
  sunrealtype k23;
  sunrealtype bias;   /* error bias factor */
  sunrealtype esp;    /* slow error from previous step */
  sunrealtype efp;    /* fast error from previous step */
  sunrealtype espp;   /* slow error from two previous steps ago */
  sunrealtype efpp;   /* fast error from two previous steps ago */
  int p;              /* fast order of accuracy to use */
};

typedef struct _SUNAdaptControllerContent_MRIPID* SUNAdaptControllerContent_MRIPID;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_MRIPID(SUNContext sunctx, int p);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_MRIPID(SUNAdaptController C, sunrealtype k11,
                                        sunrealtype k12, sunrealtype k13,
                                        sunrealtype k21, sunrealtype k22,
                                        sunrealtype k23);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_MRIPID(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateMRISteps_MRIPID(SUNAdaptController C, sunrealtype H,
                                               sunrealtype h, int P, sunrealtype DSM,
                                               sunrealtype dsm, sunrealtype* Hnew,
                                               sunrealtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_Reset_MRIPID(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_MRIPID(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_MRIPID(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_MRIPID(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_UpdateMRIH_MRIPID(SUNAdaptController C, sunrealtype H,
                                         sunrealtype h, sunrealtype DSM,
                                         sunrealtype dsm);
SUNDIALS_EXPORT
int SUNAdaptController_Space_MRIPID(SUNAdaptController C, long int* lenrw,
                                    long int* leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_MRIPID_H */
