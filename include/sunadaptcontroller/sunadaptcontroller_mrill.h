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
 * This is the header file for the SUNAdaptController_MRILL module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_MRILL_H
#define _SUNADAPTCONTROLLER_MRILL_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------------
 * MRI linear-linear implementation of SUNAdaptController
 * ---------------------------------------------- */

struct _SUNAdaptControllerContent_MRILL
{
  sunrealtype k11; /* internal controller parameters */
  sunrealtype k12;
  sunrealtype k21;
  sunrealtype k22;
  sunrealtype bias;         /* error bias factor */
  sunrealtype esp;          /* slow error from previous step */
  sunrealtype efp;          /* fast error from previous step */
  sunrealtype hsp;          /* slow previous step size */
  sunrealtype hfp;          /* fast previous step size */
  int p;                    /* fast order of accuracy to use */
  sunbooleantype firststep; /* flag indicating first step */
};

typedef struct _SUNAdaptControllerContent_MRILL* SUNAdaptControllerContent_MRILL;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_MRILL(SUNContext sunctx, int p);
SUNDIALS_EXPORT
int SUNAdaptController_SetParams_MRILL(SUNAdaptController C, sunrealtype k11,
                                       sunrealtype k12, sunrealtype k21,
                                       sunrealtype k22);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_MRILL(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_EstimateMRISteps_MRILL(SUNAdaptController C,
                                              sunrealtype H, sunrealtype h,
                                              int P, sunrealtype DSM,
                                              sunrealtype dsm, sunrealtype* Hnew,
                                              sunrealtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptController_Reset_MRILL(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_MRILL(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptController_Write_MRILL(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_MRILL(SUNAdaptController C, sunrealtype bias);
SUNDIALS_EXPORT
int SUNAdaptController_UpdateMRIH_MRILL(SUNAdaptController C, sunrealtype H,
                                        sunrealtype h, sunrealtype DSM,
                                        sunrealtype dsm);
SUNDIALS_EXPORT
int SUNAdaptController_Space_MRILL(SUNAdaptController C, long int* lenrw,
                                   long int* leniw);

#ifdef __cplusplus
}
#endif

#endif /* _SUNADAPTCONTROLLER_MRILL_H */
