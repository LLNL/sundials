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
 * This is the implementation file for the SUNAdaptController_ImExGus
 * module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>

#include "sundials_macros.h"

/* ---------------
 * Macro accessors
 * --------------- */

#define SACIMEXGUS_CONTENT(C)   ((SUNAdaptControllerContent_ImExGus)(C->content))
#define SACIMEXGUS_K1E(C)       (SACIMEXGUS_CONTENT(C)->k1e)
#define SACIMEXGUS_K2E(C)       (SACIMEXGUS_CONTENT(C)->k2e)
#define SACIMEXGUS_K1I(C)       (SACIMEXGUS_CONTENT(C)->k1i)
#define SACIMEXGUS_K2I(C)       (SACIMEXGUS_CONTENT(C)->k2i)
#define SACIMEXGUS_BIAS(C)      (SACIMEXGUS_CONTENT(C)->bias)
#define SACIMEXGUS_EP(C)        (SACIMEXGUS_CONTENT(C)->ep)
#define SACIMEXGUS_HP(C)        (SACIMEXGUS_CONTENT(C)->hp)
#define SACIMEXGUS_FIRSTSTEP(C) (SACIMEXGUS_CONTENT(C)->firststep)

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1E  SUN_RCONST(0.367)
#define DEFAULT_K2E  SUN_RCONST(0.268)
#define DEFAULT_K1I  SUN_RCONST(0.95)
#define DEFAULT_K2I  SUN_RCONST(0.95)
#define DEFAULT_BIAS SUN_RCONST(1.5)
#define TINY         SUN_RCONST(1.0e-10)

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new ImExGus controller
 */

SUNAdaptController SUNAdaptController_ImExGus(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  SUNAdaptControllerContent_ImExGus content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  C->ops->gettype      = SUNAdaptController_GetType_ImExGus;
  C->ops->estimatestep = SUNAdaptController_EstimateStep_ImExGus;
  C->ops->reset        = SUNAdaptController_Reset_ImExGus;
  C->ops->setdefaults  = SUNAdaptController_SetDefaults_ImExGus;
  C->ops->write        = SUNAdaptController_Write_ImExGus;
  C->ops->seterrorbias = SUNAdaptController_SetErrorBias_ImExGus;
  C->ops->updateh      = SUNAdaptController_UpdateH_ImExGus;
  C->ops->space        = SUNAdaptController_Space_ImExGus;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_ImExGus)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  C->content = content;

  /* Fill content with default/reset values */
  SUNCheckCallNull(SUNAdaptController_SetDefaults_ImExGus(C));
  SUNCheckCallNull(SUNAdaptController_Reset_ImExGus(C));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set ImExGus parameters
 */

SUNErrCode SUNAdaptController_SetParams_ImExGus(SUNAdaptController C,
                                                sunrealtype k1e, sunrealtype k2e,
                                                sunrealtype k1i, sunrealtype k2i)
{
  SUNFunctionBegin(C->sunctx);
  SACIMEXGUS_K1E(C) = k1e;
  SACIMEXGUS_K2E(C) = k2e;
  SACIMEXGUS_K1I(C) = k1i;
  SACIMEXGUS_K2I(C) = k2i;
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_ImExGus(
  SUNDIALS_MAYBE_UNUSED SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_H;
}

SUNErrCode SUNAdaptController_EstimateStep_ImExGus(SUNAdaptController C,
                                                   sunrealtype h, int p,
                                                   sunrealtype dsm,
                                                   sunrealtype* hnew)
{
  SUNFunctionBegin(C->sunctx);

  SUNAssert(hnew, SUN_ERR_ARG_CORRUPT);

  /* order parameter to use */
  const int ord = p + 1;

  /* set usable time-step adaptivity parameters -- first step */
  const sunrealtype k = -SUN_RCONST(1.0) / ord;
  const sunrealtype e = SUNMAX(SACIMEXGUS_BIAS(C) * dsm, TINY);

  /* set usable time-step adaptivity parameters -- subsequent steps */
  const sunrealtype k1e  = -SACIMEXGUS_K1E(C) / ord;
  const sunrealtype k2e  = -SACIMEXGUS_K2E(C) / ord;
  const sunrealtype k1i  = -SACIMEXGUS_K1I(C) / ord;
  const sunrealtype k2i  = -SACIMEXGUS_K2I(C) / ord;
  const sunrealtype e1   = SUNMAX(SACIMEXGUS_BIAS(C) * dsm, TINY);
  const sunrealtype e2   = e1 / SUNMAX(SACIMEXGUS_EP(C), TINY);
  const sunrealtype hrat = h / SACIMEXGUS_HP(C);

  /* compute estimated time step size, modifying the first step formula */
  if (SACIMEXGUS_FIRSTSTEP(C)) { *hnew = h * SUNRpowerR(e, k); }
  else
  {
    *hnew = h * SUNMIN(hrat * SUNRpowerR(e1, k1i) * SUNRpowerR(e2, k2i),
                       SUNRpowerR(e1, k1e) * SUNRpowerR(e2, k2e));
  }

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Reset_ImExGus(SUNAdaptController C)
{
  SACIMEXGUS_EP(C)        = SUN_RCONST(1.0);
  SACIMEXGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetDefaults_ImExGus(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SACIMEXGUS_K1E(C)  = DEFAULT_K1E;
  SACIMEXGUS_K2E(C)  = DEFAULT_K2E;
  SACIMEXGUS_K1I(C)  = DEFAULT_K1I;
  SACIMEXGUS_K2I(C)  = DEFAULT_K2I;
  SACIMEXGUS_BIAS(C) = DEFAULT_BIAS;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Write_ImExGus(SUNAdaptController C, FILE* fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "ImEx Gustafsson SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1e = %32Lg\n", SACIMEXGUS_K1E(C));
  fprintf(fptr, "  k2e = %32Lg\n", SACIMEXGUS_K2E(C));
  fprintf(fptr, "  k1i = %32Lg\n", SACIMEXGUS_K1I(C));
  fprintf(fptr, "  k2i = %32Lg\n", SACIMEXGUS_K2I(C));
  fprintf(fptr, "  bias factor = %22Lg\n", SACIMEXGUS_BIAS(C));
  fprintf(fptr, "  previous error = %22Lg\n", SACIMEXGUS_EP(C));
  fprintf(fptr, "  previous step = %22Lg\n", SACIMEXGUS_HP(C));
#else
  fprintf(fptr, "  k1e = %16g\n", SACIMEXGUS_K1E(C));
  fprintf(fptr, "  k2e = %16g\n", SACIMEXGUS_K2E(C));
  fprintf(fptr, "  k1i = %16g\n", SACIMEXGUS_K1I(C));
  fprintf(fptr, "  k2i = %16g\n", SACIMEXGUS_K2I(C));
  fprintf(fptr, "  bias factor = %16g\n", SACIMEXGUS_BIAS(C));
  fprintf(fptr, "  previous error = %16g\n", SACIMEXGUS_EP(C));
  fprintf(fptr, "  previous step = %16g\n", SACIMEXGUS_HP(C));
#endif
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetErrorBias_ImExGus(SUNAdaptController C,
                                                   sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);
  /* set allowed value, otherwise set default */
  if (bias <= SUN_RCONST(0.0)) { SACIMEXGUS_BIAS(C) = DEFAULT_BIAS; }
  else { SACIMEXGUS_BIAS(C) = bias; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_UpdateH_ImExGus(SUNAdaptController C,
                                              sunrealtype h, sunrealtype dsm)
{
  SUNFunctionBegin(C->sunctx);
  SACIMEXGUS_EP(C)        = SACIMEXGUS_BIAS(C) * dsm;
  SACIMEXGUS_HP(C)        = h;
  SACIMEXGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Space_ImExGus(SUNAdaptController C,
                                            long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = 7;
  *leniw = 1;
  return SUN_SUCCESS;
}
