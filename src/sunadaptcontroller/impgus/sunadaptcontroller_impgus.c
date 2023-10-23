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
 * This is the implementation file for the SUNAdaptController_ImpGus
 * module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_impgus.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SACIMPGUS_CONTENT(C)     ( (SUNAdaptControllerContent_ImpGus)(C->content) )
#define SACIMPGUS_K1(C)          ( SACIMPGUS_CONTENT(C)->k1 )
#define SACIMPGUS_K2(C)          ( SACIMPGUS_CONTENT(C)->k2 )
#define SACIMPGUS_BIAS(C)        ( SACIMPGUS_CONTENT(C)->bias )
#define SACIMPGUS_EP(C)          ( SACIMPGUS_CONTENT(C)->ep )
#define SACIMPGUS_HP(C)          ( SACIMPGUS_CONTENT(C)->hp )
#define SACIMPGUS_FIRSTSTEP(C)   ( SACIMPGUS_CONTENT(C)->firststep )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1     RCONST(0.98)
#define DEFAULT_K2     RCONST(0.95)
#define DEFAULT_BIAS   RCONST(1.5)
#define TINY           RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new ImpGus controller
 */

SUNAdaptController SUNAdaptController_ImpGus(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptControllerContent_ImpGus content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype        = SUNAdaptController_GetType_ImpGus;
  C->ops->estimatestep   = SUNAdaptController_EstimateStep_ImpGus;
  C->ops->reset          = SUNAdaptController_Reset_ImpGus;
  C->ops->setdefaults    = SUNAdaptController_SetDefaults_ImpGus;
  C->ops->write          = SUNAdaptController_Write_ImpGus;
  C->ops->seterrorbias   = SUNAdaptController_SetErrorBias_ImpGus;
  C->ops->update         = SUNAdaptController_Update_ImpGus;
  C->ops->space          = SUNAdaptController_Space_ImpGus;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_ImpGus)malloc(sizeof *content);
  if (content == NULL)
  {
    (void) SUNAdaptController_Destroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Fill content with default/reset values */
  SUNAdaptController_SetDefaults_ImpGus(C);
  SUNAdaptController_Reset_ImpGus(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set ImpGus parameters
 */

int SUNAdaptController_SetParams_ImpGus(SUNAdaptController C,
                                        sunrealtype k1, sunrealtype k2)
{
  SACIMPGUS_K1(C) = k1;
  SACIMPGUS_K2(C) = k2;
  return SUNADAPTCONTROLLER_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_ImpGus(SUNAdaptController C)
{ return SUN_ADAPTCONTROLLER_H; }

int SUNAdaptController_EstimateStep_ImpGus(SUNAdaptController C, sunrealtype h,
                                           int p, sunrealtype dsm, sunrealtype* hnew)
{
  /* order parameter to use */
  const int ord = p + 1;

  /* modified method for first step */
  if (SACIMPGUS_FIRSTSTEP(C))
  {
    /* set usable time-step adaptivity parameters */
    const sunrealtype k = -RCONST(1.0) / ord;
    const sunrealtype e = SUNMAX(SACIMPGUS_BIAS(C) * dsm, TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e,k);
  }
  else
  {
    /* set usable time-step adaptivity parameters */
    const sunrealtype k1 = -SACIMPGUS_K1(C) / ord;
    const sunrealtype k2 = -SACIMPGUS_K2(C) / ord;
    const sunrealtype ecur = SACIMPGUS_BIAS(C) * dsm;
    const sunrealtype e1 = SUNMAX(ecur, TINY);
    const sunrealtype e2 = e1 / SUNMAX(SACIMPGUS_EP(C), TINY);
    const sunrealtype hrat = h / SACIMPGUS_HP(C);

    /* compute estimated optimal time step size */
    *hnew = h * hrat * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  }

  /* return with success */
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Reset_ImpGus(SUNAdaptController C)
{
  SACIMPGUS_EP(C) = RCONST(1.0);
  SACIMPGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetDefaults_ImpGus(SUNAdaptController C)
{
  SACIMPGUS_K1(C)   = DEFAULT_K1;
  SACIMPGUS_K2(C)   = DEFAULT_K2;
  SACIMPGUS_BIAS(C) = DEFAULT_BIAS;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Write_ImpGus(SUNAdaptController C, FILE *fptr)
{
  fprintf(fptr, "Implicit Gustafsson SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SACIMPGUS_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", SACIMPGUS_K2(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SACIMPGUS_BIAS(C));
  fprintf(fptr, "  previous error = %32Lg\n", SACIMPGUS_EP(C));
  fprintf(fptr, "  previous step = %32Lg\n", SACIMPGUS_HP(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SACIMPGUS_K1(C));
  fprintf(fptr, "  k2 = %16g\n", SACIMPGUS_K2(C));
  fprintf(fptr, "  bias factor = %16g\n", SACIMPGUS_BIAS(C));
  fprintf(fptr, "  previous error = %16g\n", SACIMPGUS_EP(C));
  fprintf(fptr, "  previous step = %16g\n", SACIMPGUS_HP(C));
#endif
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetErrorBias_ImpGus(SUNAdaptController C, sunrealtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SACIMPGUS_BIAS(C) = DEFAULT_BIAS;
  } else {
    SACIMPGUS_BIAS(C) = bias;
  }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Update_ImpGus(SUNAdaptController C, sunrealtype h, sunrealtype dsm)
{
  SACIMPGUS_EP(C) = SACIMPGUS_BIAS(C) * dsm;
  SACIMPGUS_HP(C) = h;
  SACIMPGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Space_ImpGus(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 5;
  *leniw = 2;
  return SUNADAPTCONTROLLER_SUCCESS;
}
