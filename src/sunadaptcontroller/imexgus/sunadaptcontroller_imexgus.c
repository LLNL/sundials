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
 * This is the implementation file for the SUNAdaptController_ImExGus
 * module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sundials/sundials_math.h>

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
  SUNAdaptController C;
  SUNAdaptControllerContent_ImExGus content;

  if (sunctx == NULL) { return NULL; }

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

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
  if (content == NULL)
  {
    (void)SUNAdaptController_Destroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Fill content with default/reset values */
  SUNAdaptController_SetDefaults_ImExGus(C);
  SUNAdaptController_Reset_ImExGus(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set ImExGus parameters
 */

int SUNAdaptController_SetParams_ImExGus(SUNAdaptController C, sunrealtype k1e,
                                         sunrealtype k2e, sunrealtype k1i,
                                         sunrealtype k2i)
{
  if (C == NULL) { return SUNADAPTCONTROLLER_ILL_INPUT; }
  SACIMEXGUS_K1E(C) = k1e;
  SACIMEXGUS_K2E(C) = k2e;
  SACIMEXGUS_K1I(C) = k1i;
  SACIMEXGUS_K2I(C) = k2i;
  return SUNADAPTCONTROLLER_SUCCESS;
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_ImExGus(SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_H;
}

int SUNAdaptController_EstimateStep_ImExGus(SUNAdaptController C, sunrealtype h,
                                            int p, sunrealtype dsm,
                                            sunrealtype* hnew)
{
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

  if (C == NULL || hnew == NULL) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* compute estimated time step size, modifying the first step formula */
  if (SACIMEXGUS_FIRSTSTEP(C)) { *hnew = h * SUNRpowerR(e, k); }
  else
  {
    *hnew = h * SUNMIN(hrat * SUNRpowerR(e1, k1i) * SUNRpowerR(e2, k2i),
                       SUNRpowerR(e1, k1e) * SUNRpowerR(e2, k2e));
  }

  /* return with success */
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Reset_ImExGus(SUNAdaptController C)
{
  SACIMEXGUS_EP(C)        = SUN_RCONST(1.0);
  SACIMEXGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetDefaults_ImExGus(SUNAdaptController C)
{
  if (C == NULL) { return SUNADAPTCONTROLLER_ILL_INPUT; }
  SACIMEXGUS_K1E(C)  = DEFAULT_K1E;
  SACIMEXGUS_K2E(C)  = DEFAULT_K2E;
  SACIMEXGUS_K1I(C)  = DEFAULT_K1I;
  SACIMEXGUS_K2I(C)  = DEFAULT_K2I;
  SACIMEXGUS_BIAS(C) = DEFAULT_BIAS;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Write_ImExGus(SUNAdaptController C, FILE* fptr)
{
  if (C == NULL || fptr == NULL) { return SUNADAPTCONTROLLER_ILL_INPUT; }
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
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetErrorBias_ImExGus(SUNAdaptController C, sunrealtype bias)
{
  if (C == NULL) { return SUNADAPTCONTROLLER_ILL_INPUT; }
  /* set allowed value, otherwise set default */
  if (bias <= SUN_RCONST(0.0)) { SACIMEXGUS_BIAS(C) = DEFAULT_BIAS; }
  else { SACIMEXGUS_BIAS(C) = bias; }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_UpdateH_ImExGus(SUNAdaptController C, sunrealtype h,
                                       sunrealtype dsm)
{
  if (C == NULL) { return SUNADAPTCONTROLLER_ILL_INPUT; }
  SACIMEXGUS_EP(C)        = SACIMEXGUS_BIAS(C) * dsm;
  SACIMEXGUS_HP(C)        = h;
  SACIMEXGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Space_ImExGus(SUNAdaptController C, long int* lenrw,
                                     long int* leniw)
{
  if (C == NULL || lenrw == NULL || leniw == NULL)
  {
    return SUNADAPTCONTROLLER_ILL_INPUT;
  }
  *lenrw = 7;
  *leniw = 1;
  return SUNADAPTCONTROLLER_SUCCESS;
}
