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
 * This is the implementation file for the SUNAdaptController_PI
 * module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_pi.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SACPI_CONTENT(C)     ( (SUNAdaptControllerContent_PI)(C->content) )
#define SACPI_K1(C)          ( SACPI_CONTENT(C)->k1 )
#define SACPI_K2(C)          ( SACPI_CONTENT(C)->k2 )
#define SACPI_BIAS(C)        ( SACPI_CONTENT(C)->bias )
#define SACPI_EP(C)          ( SACPI_CONTENT(C)->ep )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1     RCONST(0.8)
#define DEFAULT_K2     RCONST(0.31)
#define DEFAULT_BIAS   RCONST(1.5)
#define TINY           RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new PI controller
 */

SUNAdaptController SUNAdaptController_PI(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptControllerContent_PI content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype        = SUNAdaptController_GetType_PI;
  C->ops->estimatestep   = SUNAdaptController_EstimateStep_PI;
  C->ops->reset          = SUNAdaptController_Reset_PI;
  C->ops->setdefaults    = SUNAdaptController_SetDefaults_PI;
  C->ops->write          = SUNAdaptController_Write_PI;
  C->ops->seterrorbias   = SUNAdaptController_SetErrorBias_PI;
  C->ops->update         = SUNAdaptController_Update_PI;
  C->ops->space          = SUNAdaptController_Space_PI;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_PI)malloc(sizeof *content);
  if (content == NULL)
  {
    (void) SUNAdaptController_Destroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Fill content with default/reset values */
  SUNAdaptController_SetDefaults_PI(C);
  SUNAdaptController_Reset_PI(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PI parameters
 */

int SUNAdaptController_SetParams_PI(SUNAdaptController C,
                                    sunrealtype k1, sunrealtype k2)
{
  SACPI_K1(C) = k1;
  SACPI_K2(C) = k2;
  return SUNADAPTCONTROLLER_SUCCESS;
}


/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_PI(SUNAdaptController C)
{ return SUN_ADAPTCONTROLLER_H; }

int SUNAdaptController_EstimateStep_PI(SUNAdaptController C, sunrealtype h,
                                       int p, sunrealtype dsm, sunrealtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const int ord = p + 1;
  const sunrealtype k1 = -SACPI_K1(C) / ord;
  const sunrealtype k2 =  SACPI_K2(C) / ord;
  const sunrealtype ecur = SACPI_BIAS(C) * dsm;
  const sunrealtype e1 = SUNMAX(ecur, TINY);
  const sunrealtype e2 = SUNMAX(SACPI_EP(C), TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = h * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Reset_PI(SUNAdaptController C)
{
  SACPI_EP(C) = RCONST(1.0);
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetDefaults_PI(SUNAdaptController C)
{
  SACPI_K1(C)   = DEFAULT_K1;
  SACPI_K2(C)   = DEFAULT_K2;
  SACPI_BIAS(C) = DEFAULT_BIAS;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Write_PI(SUNAdaptController C, FILE *fptr)
{
  fprintf(fptr, "PI SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SACPI_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", SACPI_K2(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SACPI_BIAS(C));
  fprintf(fptr, "  previous error = %32Lg\n", SACPI_EP(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SACPI_K1(C));
  fprintf(fptr, "  k2 = %16g\n", SACPI_K2(C));
  fprintf(fptr, "  bias factor = %16g\n", SACPI_BIAS(C));
  fprintf(fptr, "  previous error = %16g\n", SACPI_EP(C));
#endif
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetErrorBias_PI(SUNAdaptController C, sunrealtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SACPI_BIAS(C) = DEFAULT_BIAS;
  } else {
    SACPI_BIAS(C) = bias;
  }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Update_PI(SUNAdaptController C, sunrealtype h, sunrealtype dsm)
{
  SACPI_EP(C) = SACPI_BIAS(C) * dsm;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Space_PI(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 4;
  *leniw = 1;
  return SUNADAPTCONTROLLER_SUCCESS;
}
