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

#include "sundials_macros.h"

/* ------------------
 * Default parameters
 * ------------------ */

/*   maximum relative change for inner tolerance factor */
#define INNER_MAX_RELCH SUN_RCONST(20.0)
/*   minimum tolerance factor for inner solver */
#define INNER_MIN_TOLFAC SUN_RCONST(1.e-5)
/*   maximum tolerance factor for inner solver */
#define INNER_MAX_TOLFAC SUN_RCONST(1.0)

/* ---------------
 * Macro accessors
 * --------------- */

#define MRIHTOL_CONTENT(C)          ((SUNAdaptControllerContent_MRIHTol)(C->content))
#define MRIHTOL_CSLOW(C)            (MRIHTOL_CONTENT(C)->HControl)
#define MRIHTOL_CFAST(C)            (MRIHTOL_CONTENT(C)->TolControl)
#define MRIHTOL_INNER_MAX_RELCH(C)  (MRIHTOL_CONTENT(C)->inner_max_relch)
#define MRIHTOL_INNER_MIN_TOLFAC(C) (MRIHTOL_CONTENT(C)->inner_min_tolfac)
#define MRIHTOL_INNER_MAX_TOLFAC(C) (MRIHTOL_CONTENT(C)->inner_max_tolfac)

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRIHTol controller
 */

SUNAdaptController SUNAdaptController_MRIHTol(SUNAdaptController HControl,
                                              SUNAdaptController TolControl,
                                              SUNContext sunctx)
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
  C->ops->updatemrihtol   = SUNAdaptController_UpdateMRIHTol_MRIHTol;
  C->ops->space           = SUNAdaptController_Space_MRIHTol;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_MRIHTol)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach input controllers */
  content->HControl   = HControl;
  content->TolControl = TolControl;

  /* Set parameters to default values */
  content->inner_max_relch  = INNER_MAX_RELCH;
  content->inner_min_tolfac = INNER_MIN_TOLFAC;
  content->inner_max_tolfac = INNER_MAX_TOLFAC;

  /* Attach content */
  C->content = content;

  return C;
}

/* -----------------------------------------------------------------
 * Function to set MRIHTol parameters
 */

SUNErrCode SUNAdaptController_SetParams_MRIHTol(SUNAdaptController C,
                                                sunrealtype inner_max_relch,
                                                sunrealtype inner_min_tolfac,
                                                sunrealtype inner_max_tolfac)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(inner_max_tolfac > inner_min_tolfac, SUN_ERR_ARG_OUTOFRANGE);
  if (inner_max_relch < SUN_RCONST(1.0))
  {
    MRIHTOL_INNER_MAX_RELCH(C) = INNER_MAX_RELCH;
  }
  else { MRIHTOL_INNER_MAX_RELCH(C) = inner_max_relch; }
  if (inner_min_tolfac <= SUN_RCONST(0.0))
  {
    MRIHTOL_INNER_MIN_TOLFAC(C) = INNER_MIN_TOLFAC;
  }
  else { MRIHTOL_INNER_MIN_TOLFAC(C) = inner_min_tolfac; }
  if ((inner_max_tolfac <= SUN_RCONST(0.0)) ||
      (inner_max_tolfac > SUN_RCONST(1.0)))
  {
    MRIHTOL_INNER_MAX_TOLFAC(C) = INNER_MAX_TOLFAC;
  }
  else { MRIHTOL_INNER_MAX_TOLFAC(C) = inner_max_tolfac; }
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to get slow and fast sub-controllers
 */

SUNErrCode SUNAdaptController_GetSlowController_MRIHTol(SUNAdaptController C,
                                                        SUNAdaptController* Cslow)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Cslow, SUN_ERR_ARG_CORRUPT);
  *Cslow = MRIHTOL_CSLOW(C);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_GetFastController_MRIHTol(SUNAdaptController C,
                                                        SUNAdaptController* Cfast)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Cfast, SUN_ERR_ARG_CORRUPT);
  *Cfast = MRIHTOL_CFAST(C);
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_MRIHTol(
  SUNDIALS_MAYBE_UNUSED SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_MRI_H_TOL;
}

SUNErrCode SUNAdaptController_EstimateStepTol_MRIHTol(
  SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
  sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Hnew, SUN_ERR_ARG_CORRUPT);
  SUNAssert(tolfacnew, SUN_ERR_ARG_CORRUPT);
  sunrealtype tolfacest;

  /* Call slow time scale sub-controller to fill Hnew -- note that all heuristics
     bounds on Hnew will be enforced by the time integrator itself */
  SUNCheckCall(SUNAdaptController_EstimateStep(MRIHTOL_CSLOW(C), H, P, DSM, Hnew));

  /* Call fast time scale sub-controller with order=1: no matter the integrator
     order, we expect its error to be proportional to the tolerance factor */
  SUNCheckCall(SUNAdaptController_EstimateStep(MRIHTOL_CFAST(C), tolfac, 0, dsm,
                                               &tolfacest));

  /* Enforce bounds on estimated tolerance factor */
  /*     keep relative change within bounds */
  tolfacest = SUNMAX(tolfacest, tolfac / MRIHTOL_INNER_MAX_RELCH(C));
  tolfacest = SUNMIN(tolfacest, tolfac * MRIHTOL_INNER_MAX_RELCH(C));
  /*     enforce absolute min/max bounds */
  tolfacest = SUNMAX(tolfacest, MRIHTOL_INNER_MIN_TOLFAC(C));
  tolfacest = SUNMIN(tolfacest, MRIHTOL_INNER_MAX_TOLFAC(C));

  /* Set result and return */
  *tolfacnew = tolfacest;
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
  MRIHTOL_INNER_MAX_RELCH(C)  = INNER_MAX_RELCH;
  MRIHTOL_INNER_MIN_TOLFAC(C) = INNER_MIN_TOLFAC;
  MRIHTOL_INNER_MAX_TOLFAC(C) = INNER_MAX_TOLFAC;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Write_MRIHTol(SUNAdaptController C, FILE* fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "Multirate H-Tol SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  inner_max_relch  = %32Lg\n", MRIHTOL_INNER_MAX_RELCH(C));
  fprintf(fptr, "  inner_min_tolfac = %32Lg\n", MRIHTOL_INNER_MIN_TOLFAC(C));
  fprintf(fptr, "  inner_max_tolfac = %32Lg\n", MRIHTOL_INNER_MAX_TOLFAC(C));
#else
  fprintf(fptr, "  inner_max_relch  = %16g\n", MRIHTOL_INNER_MAX_RELCH(C));
  fprintf(fptr, "  inner_min_tolfac = %16g\n", MRIHTOL_INNER_MIN_TOLFAC(C));
  fprintf(fptr, "  inner_max_tolfac = %16g\n", MRIHTOL_INNER_MAX_TOLFAC(C));
#endif
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

SUNErrCode SUNAdaptController_UpdateMRIHTol_MRIHTol(SUNAdaptController C,
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
