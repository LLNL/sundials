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
 * This is the header file for the SUNAdaptController_MRIHTol module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_MRIHTOL_H
#define _SUNADAPTCONTROLLER_MRIHTOL_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------
 * MRI H+tolerance implementation of SUNAdaptController
 * -------------------------------------------- */

struct _SUNAdaptControllerContent_MRIHTol {
  SUNAdaptController HControl;
  SUNAdaptController TolControl;
};

typedef struct _SUNAdaptControllerContent_MRIHTol* SUNAdaptControllerContent_MRIHTol;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_MRIHTol(SUNContext sunctx,
                                              SUNAdaptController HControl,
                                              SUNAdaptController TolControl);
SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_GetSlowController_MRIHTol(SUNAdaptController C);
SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_GetFastController_MRIHTol(SUNAdaptController C);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_MRIHTol(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStepTol_MRIHTol(SUNAdaptController C, sunrealtype H,
                                               sunrealtype tolfac, int P,
                                               sunrealtype DSM, sunrealtype dsm,
                                               sunrealtype* Hnew,
                                               sunrealtype* tolfacnew);
SUNDIALS_EXPORT
int SUNAdaptController_Reset_MRIHTol(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_MRIHTol(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_MRIHTol(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_MRIHTol(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_UpdateMRITol_MRIHTol(SUNAdaptController C, sunrealtype H,
                                            sunrealtype tolfac, sunrealtype DSM,
                                            sunrealtype dsm);
SUNDIALS_EXPORT
int SUNAdaptController_Space_MRIHTol(SUNAdaptController C, long int* lenrw,
                                     long int* leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_MRIHTOL_H */