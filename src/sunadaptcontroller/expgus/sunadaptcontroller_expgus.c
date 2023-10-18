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
 * This is the implementation file for the SUNAdaptController_ExpGus
 * module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_expgus.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SACEXPGUS_CONTENT(C)     ( (SUNAdaptControllerContent_ExpGus)(C->content) )
#define SACEXPGUS_K1(C)          ( SACEXPGUS_CONTENT(C)->k1 )
#define SACEXPGUS_K2(C)          ( SACEXPGUS_CONTENT(C)->k2 )
#define SACEXPGUS_BIAS(C)        ( SACEXPGUS_CONTENT(C)->bias )
#define SACEXPGUS_EP(C)          ( SACEXPGUS_CONTENT(C)->ep )
#define SACEXPGUS_P(C)           ( SACEXPGUS_CONTENT(C)->p )
#define SACEXPGUS_PQ(C)          ( SACEXPGUS_CONTENT(C)->pq )
#define SACEXPGUS_ADJ(C)         ( SACEXPGUS_CONTENT(C)->adj )
#define SACEXPGUS_FIRSTSTEP(C)   ( SACEXPGUS_CONTENT(C)->firststep )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1     RCONST(0.367)
#define DEFAULT_K2     RCONST(0.268)
#define DEFAULT_BIAS   RCONST(1.5)
#define DEFAULT_PQ     0
#define DEFAULT_ADJ    -1
/* #define DEFAULT_PQ     -1 */
/* #define DEFAULT_ADJ    0 */
#define TINY           RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new ExpGus controller
 */

SUNAdaptController SUNAdaptController_ExpGus(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptControllerContent_ExpGus content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype               = SUNAdaptController_GetType_ExpGus;
  C->ops->estimatestep          = SUNAdaptController_EstimateStep_ExpGus;
  C->ops->reset                 = SUNAdaptController_Reset_ExpGus;
  C->ops->setdefaults           = SUNAdaptController_SetDefaults_ExpGus;
  C->ops->write                 = SUNAdaptController_Write_ExpGus;
  C->ops->setmethodorder        = SUNAdaptController_SetMethodOrder_ExpGus;
  C->ops->adjustcontrollerorder = SUNAdaptController_AdjustControllerOrder_ExpGus;
  C->ops->seterrorbias          = SUNAdaptController_SetErrorBias_ExpGus;
  C->ops->update                = SUNAdaptController_Update_ExpGus;
  C->ops->space                 = SUNAdaptController_Space_ExpGus;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_ExpGus)malloc(sizeof *content);
  if (content == NULL)
  {
    (void) SUNAdaptController_Destroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Initialize method order */
  content->p = 1;

  /* Fill content with default/reset values */
  SUNAdaptController_SetDefaults_ExpGus(C);
  SUNAdaptController_Reset_ExpGus(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set ExpGus parameters
 */

int SUNAdaptController_SetParams_ExpGus(SUNAdaptController C, int pq,
                                        sunrealtype k1, sunrealtype k2)
{
  /* store legal inputs, and return with success */
  SACEXPGUS_PQ(C) = pq;
  if (k1 >= RCONST(0.0)) { SACEXPGUS_K1(C) = k1; }
  if (k2 >= RCONST(0.0)) { SACEXPGUS_K2(C) = k2; }
  return SUNADAPTCONTROLLER_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_ExpGus(SUNAdaptController C)
{ return SUN_ADAPTCONTROLLER_H; }

int SUNAdaptController_EstimateStep_ExpGus(SUNAdaptController C, sunrealtype h,
                                           sunrealtype dsm, sunrealtype* hnew)
{
  /* order parameter to use in controller */
  const int ord = SACEXPGUS_P(C) + SACEXPGUS_ADJ(C) + 1;

  /* modified method for first step */
  if (SACEXPGUS_FIRSTSTEP(C))
  {
    /* set usable time-step adaptivity parameters */
    const sunrealtype k = -RCONST(1.0) / ord;
    const sunrealtype e = SUNMAX(SACEXPGUS_BIAS(C) * dsm, TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e,k);
  }
  else
  {
    /* set usable time-step adaptivity parameters */
    const sunrealtype k1 = -SACEXPGUS_K1(C) / ord;
    const sunrealtype k2 = -SACEXPGUS_K2(C) / ord;
    const sunrealtype ecur = SACEXPGUS_BIAS(C) * dsm;
    const sunrealtype e1 = SUNMAX(ecur, TINY);
    const sunrealtype e2 = e1 / SUNMAX(SACEXPGUS_EP(C), TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  }

  /* return with success */
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Reset_ExpGus(SUNAdaptController C)
{
  SACEXPGUS_EP(C) = RCONST(1.0);
  SACEXPGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetDefaults_ExpGus(SUNAdaptController C)
{
  SACEXPGUS_K1(C)     = DEFAULT_K1;
  SACEXPGUS_K2(C)     = DEFAULT_K2;
  SACEXPGUS_BIAS(C)   = DEFAULT_BIAS;
  SACEXPGUS_PQ(C)     = DEFAULT_PQ;
  SACEXPGUS_ADJ(C)    = DEFAULT_ADJ;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Write_ExpGus(SUNAdaptController C, FILE *fptr)
{
  fprintf(fptr, "Explicit Gustafsson SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SACEXPGUS_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", SACEXPGUS_K2(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SACEXPGUS_BIAS(C));
  fprintf(fptr, "  previous error = %32g\n", SACEXPGUS_EP(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SACEXPGUS_K1(C));
  fprintf(fptr, "  k2 = %16g\n", SACEXPGUS_K2(C));
  fprintf(fptr, "  bias factor = %16g\n", SACEXPGUS_BIAS(C));
  fprintf(fptr, "  previous error = %16g\n", SACEXPGUS_EP(C));
#endif
  if (SACEXPGUS_PQ(C) == 1)
  {
    fprintf(fptr, "  p = %i (method order)\n", SACEXPGUS_P(C));
  }
  else if (SACEXPGUS_PQ(C) == 0)
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SACEXPGUS_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (minimum of method & embedding order)\n", SACEXPGUS_P(C));
  }
  fprintf(fptr, "  adj = %i\n", SACEXPGUS_ADJ(C));
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetMethodOrder_ExpGus(SUNAdaptController C, int p, int q)
{
  /* check for legal inputs */
  if ((p <= 0) || (q <= 0)) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store appropriate value based on "pq" */
  if (SACEXPGUS_PQ(C) == 1) {
    SACEXPGUS_P(C) = p;
  } else if (SACEXPGUS_PQ(C) == 0) {
    SACEXPGUS_P(C) = q;
  } else {
    SACEXPGUS_P(C) = SUNMIN(p,q);
  }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_AdjustControllerOrder_ExpGus(SUNAdaptController C, int adj)
{
  SACEXPGUS_ADJ(C) = adj;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetErrorBias_ExpGus(SUNAdaptController C, sunrealtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SACEXPGUS_BIAS(C) = DEFAULT_BIAS;
  } else {
    SACEXPGUS_BIAS(C) = bias;
  }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Update_ExpGus(SUNAdaptController C, sunrealtype h, sunrealtype dsm)
{
  SACEXPGUS_EP(C) = SACEXPGUS_BIAS(C) * dsm;
  SACEXPGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Space_ExpGus(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 5;
  *leniw = 4;
  return SUNADAPTCONTROLLER_SUCCESS;
}
