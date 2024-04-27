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
 * This is the implementation file for the
 * SUNAdaptController_MRIHTol module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_mrihtol.h>
#include <sundials/sundials_core.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"

/* ---------------
 * Macro accessors
 * --------------- */

#define MRIHTOL_CONTENT(C) ((SUNAdaptControllerContent_MRIHTol)(C->content))
#define MRIHTOL_CSLOW(C)   (MRIHTOL_CONTENT(C)->HControl)
#define MRIHTOL_CFAST(C)   (MRIHTOL_CONTENT(C)->TolControl)

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRIHTol controller
 */

SUNAdaptController SUNAdaptController_MRIHTol(SUNContext sunctx,
                                              SUNAdaptController HControl,
                                              SUNAdaptController TolControl)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  SUNAdaptControllerContent_MRIHTol content;

  /* Verify that input controllers have the appropriate type */
  SUNAssertNull(SUNAdaptController_GetType(HControl) == SUN_ADAPTCONTROLLER_H,
                SUN_ERR_ARG_INCOMPATIBLE);
  SUNAssertNull(SUNAdaptController_GetType(TolControl) == SUN_ADAPTCONTROLLER_H,
                SUN_ERR_ARG_INCOMPATIBLE);

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  C->ops->gettype         = SUNAdaptController_GetType_MRIHTol;
  C->ops->estimatesteptol = SUNAdaptController_EstimateStepTol_MRIHTol;
  C->ops->reset           = SUNAdaptController_Reset_MRIHTol;
  C->ops->setdefaults     = SUNAdaptController_SetDefaults_MRIHTol;
  C->ops->write           = SUNAdaptController_Write_MRIHTol;
  C->ops->seterrorbias    = SUNAdaptController_SetErrorBias_MRIHTol;
  C->ops->updatemritol    = SUNAdaptController_UpdateMRITol_MRIHTol;
  C->ops->space           = SUNAdaptController_Space_MRIHTol;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_MRIHTol)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach input controllers */
  content->HControl   = HControl;
  content->TolControl = TolControl;

  /* Attach content */
  C->content = content;

  return (C);
}

/* -----------------------------------------------------------------
 * Function to get slow and fast sub-controllers
 */

SUNAdaptController SUNAdaptController_GetSlowController_MRIHTol(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  return MRIHTOL_CSLOW(C);
}

SUNAdaptController SUNAdaptController_GetFastController_MRIHTol(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  return MRIHTOL_CFAST(C);
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_MRIHTol(SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_MRI_TOL;
}

SUNErrCode SUNAdaptController_EstimateStepTol_MRIHTol(
  SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
  sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Hnew, SUN_ERR_ARG_CORRUPT);
  SUNAssert(tolfacnew, SUN_ERR_ARG_CORRUPT);

  /* Call slow time scale sub-controller to fill Hnew */
  SUNCheckCall(SUNAdaptController_EstimateStep(MRIHTOL_CSLOW(C), H, P, DSM, Hnew));

  /* Call fast time sacle sub-controller with order=1: no matter the integrator
     order, we expect its error to be proportional to the tolerance factor */
  SUNCheckCall(SUNAdaptController_EstimateStep(MRIHTOL_CFAST(C), tolfac, 1, dsm,
                                               tolfacnew));

  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Reset_MRIHTol(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SUNCheckCall(SUNAdaptController_Reset(MRIHTOL_CSLOW(C)));
  SUNCheckCall(SUNAdaptController_Reset(MRIHTOL_CFAST(C)));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetDefaults_MRIHTol(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SUNCheckCall(SUNAdaptController_SetDefaults(MRIHTOL_CSLOW(C)));
  SUNCheckCall(SUNAdaptController_SetDefaults(MRIHTOL_CFAST(C)));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Write_MRIHTol(SUNAdaptController C, FILE* fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "Multirate H-Tol SUNAdaptController module:\n");
  fprintf(fptr, "\nSlow step controller:\n");
  SUNCheckCall(SUNAdaptController_Write(MRIHTOL_CSLOW(C), fptr));
  fprintf(fptr, "\nFast tolerance controller:\n");
  SUNCheckCall(SUNAdaptController_Write(MRIHTOL_CFAST(C), fptr));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetErrorBias_MRIHTol(SUNAdaptController C,
                                                   sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);
  SUNCheckCall(SUNAdaptController_SetErrorBias(MRIHTOL_CSLOW(C), bias));
  SUNCheckCall(SUNAdaptController_SetErrorBias(MRIHTOL_CFAST(C), bias));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_UpdateMRITol_MRIHTol(SUNAdaptController C,
                                                   sunrealtype H,
                                                   sunrealtype tolfac,
                                                   sunrealtype DSM,
                                                   sunrealtype dsm)
{
  SUNFunctionBegin(C->sunctx);
  SUNCheckCall(SUNAdaptController_UpdateH(MRIHTOL_CSLOW(C), H, DSM));
  SUNCheckCall(SUNAdaptController_UpdateH(MRIHTOL_CFAST(C), tolfac, dsm));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Space_MRIHTol(SUNAdaptController C,
                                            long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  long int lrw, liw;
  SUNCheckCall(SUNAdaptController_Space(MRIHTOL_CSLOW(C), lenrw, leniw));
  SUNCheckCall(SUNAdaptController_Space(MRIHTOL_CFAST(C), &lrw, &liw));
  *lenrw += lrw;
  *leniw += liw;
  return SUN_SUCCESS;
}
