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
 * This is the header file for a user-provided controller function
 * SUNAdaptController module implementation.  This provides backwards-
 * compatibility for ARKODE's previous "ARKAdaptFn"
 * -----------------------------------------------------------------*/

#ifndef _ARK_USERCONTROL_H
#define _ARK_USERCONTROL_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>
#include "arkode_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------------------
 * ARKUserControl implementation of SUNAdaptController
 * --------------------------------------------------- */

struct ARKUserControlContent_ {
  realtype hp;        /* h from previous step */
  realtype hpp;       /* h from 2 steps ago */
  realtype ep;        /* error from previous step */
  realtype epp;       /* error from 2 steps ago */
  int p;              /* embedding order of accuracy */
  int q;              /* method order of accuracy */
  ARKodeMem ark_mem;  /* main ARKODE memory structure */
  ARKAdaptFn hadapt;  /* user-provided adaptivity fn */
  void* hadapt_data;  /* user-provided data pointer */
};

typedef struct ARKUserControlContent_ *ARKUserControlContent;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController ARKUserControl(SUNContext sunctx, void* arkode_mem,
                                  ARKAdaptFn hadapt, void* hadapt_data);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptControllerGetType_ARKUserControl(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptControllerEstimateStep_ARKUserControl(SUNAdaptController C, realtype h,
                                                  realtype dsm, realtype* hnew);
SUNDIALS_EXPORT
int SUNAdaptControllerReset_ARKUserControl(SUNAdaptController C);
SUNDIALS_EXPORT
int SUNAdaptControllerWrite_ARKUserControl(SUNAdaptController C, FILE* fptr);
SUNDIALS_EXPORT
int SUNAdaptControllerSetMethodOrder_ARKUserControl(SUNAdaptController C, int q);
SUNDIALS_EXPORT
int SUNAdaptControllerSetEmbeddingOrder_ARKUserControl(SUNAdaptController C, int p);
SUNDIALS_EXPORT
int SUNAdaptControllerUpdate_ARKUserControl(SUNAdaptController C, realtype h, realtype dsm);
SUNDIALS_EXPORT
int SUNAdaptControllerSpace_ARKUserControl(SUNAdaptController C, long int *lenrw,
                                           long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _ARK_USERCONTROL_H */
