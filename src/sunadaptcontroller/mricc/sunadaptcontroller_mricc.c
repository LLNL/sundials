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
 * SUNAdaptController_MRICC module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_mricc.h>
#include <sundials/sundials_core.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"

/* ---------------
 * Macro accessors
 * --------------- */

#define MRICC_CONTENT(C) ((SUNAdaptControllerContent_MRICC)(C->content))
#define MRICC_K1(C)      (MRICC_CONTENT(C)->k1)
#define MRICC_K2(C)      (MRICC_CONTENT(C)->k2)
#define MRICC_BIAS(C)    (MRICC_CONTENT(C)->bias)
#define MRICC_PFAST(C)   (MRICC_CONTENT(C)->p)

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1   SUN_RCONST(0.42)
#define DEFAULT_K2   SUN_RCONST(0.44)
#define DEFAULT_BIAS SUN_RCONST(1.5)
#define ONE          SUN_RCONST(1.0)
#define TINY         (SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF)

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRICC controller
 */

SUNAdaptController SUNAdaptController_MRICC(SUNContext sunctx, int p)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  SUNAdaptControllerContent_MRICC content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  C->ops->gettype          = SUNAdaptController_GetType_MRICC;
  C->ops->estimatemristeps = SUNAdaptController_EstimateMRISteps_MRICC;
  C->ops->setdefaults      = SUNAdaptController_SetDefaults_MRICC;
  C->ops->write            = SUNAdaptController_Write_MRICC;
  C->ops->seterrorbias     = SUNAdaptController_SetErrorBias_MRICC;
  C->ops->space            = SUNAdaptController_Space_MRICC;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_MRICC)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  C->content = content;

  /* Set fast method order */
  content->p = p;

  /* Fill content with default/reset values */
  SUNCheckCallNull(SUNAdaptController_SetDefaults_MRICC(C));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set MRICC parameters
 */

SUNErrCode SUNAdaptController_SetParams_MRICC(SUNAdaptController C,
                                              sunrealtype k1, sunrealtype k2)
{
  SUNFunctionBegin(C->sunctx);
  MRICC_K1(C) = k1;
  MRICC_K2(C) = k2;
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_MRICC(SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_MRI_H;
}

SUNErrCode SUNAdaptController_EstimateMRISteps_MRICC(
  SUNAdaptController C, sunrealtype H, sunrealtype h, int P, sunrealtype DSM,
  sunrealtype dsm, sunrealtype* Hnew, sunrealtype* hnew)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Hnew, SUN_ERR_ARG_CORRUPT);
  SUNAssert(hnew, SUN_ERR_ARG_CORRUPT);

  /* set usable time-step adaptivity parameters */
  const int p          = MRICC_PFAST(C);
  const sunrealtype k1 = MRICC_K1(C);
  const sunrealtype k2 = MRICC_K2(C);
  const sunrealtype al = k1 / P;
  const sunrealtype b1 = (p + 1) * k1 / (P * p);
  const sunrealtype b2 = -k2 / p;
  const sunrealtype es = ONE / SUNMAX(MRICC_BIAS(C) * DSM, TINY);
  const sunrealtype ef = ONE / SUNMAX(MRICC_BIAS(C) * dsm, TINY);
  const sunrealtype M  = SUNRceil(H / h);

  /* compute estimated optimal time step size */
  *Hnew                  = H * SUNRpowerR(es, al);
  const sunrealtype Mnew = M * SUNRpowerR(es, b1) * SUNRpowerR(ef, b2);
  *hnew                  = (*Hnew) / Mnew;

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetDefaults_MRICC(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  MRICC_K1(C)   = DEFAULT_K1;
  MRICC_K2(C)   = DEFAULT_K2;
  MRICC_BIAS(C) = DEFAULT_BIAS;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Write_MRICC(SUNAdaptController C, FILE* fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "Multirate constant-constant SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", MRICC_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", MRICC_K2(C));
  fprintf(fptr, "  bias factor = %32Lg\n", MRICC_BIAS(C));
#else
  fprintf(fptr, "  k1 = %16g\n", MRICC_K1(C));
  fprintf(fptr, "  k2 = %16g\n", MRICC_K2(C));
  fprintf(fptr, "  bias factor = %16g\n", MRICC_BIAS(C));
#endif
  fprintf(fptr, "  p = %i (fast method order)\n", MRICC_PFAST(C));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetErrorBias_MRICC(SUNAdaptController C,
                                                 sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);

  /* set allowed value, otherwise set default */
  if (bias <= SUN_RCONST(0.0)) { MRICC_BIAS(C) = DEFAULT_BIAS; }
  else { MRICC_BIAS(C) = bias; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Space_MRICC(SUNAdaptController C, long int* lenrw,
                                          long int* leniw)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = 3;
  *leniw = 1;
  return SUN_SUCCESS;
}
