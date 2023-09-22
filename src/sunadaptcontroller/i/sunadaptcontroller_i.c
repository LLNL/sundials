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
 * This is the implementation file for the SUNAdaptController_I
 * module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_i.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SACI_CONTENT(C)     ( (SUNAdaptControllerContent_I)(C->content) )
#define SACI_K1(C)          ( SACI_CONTENT(C)->k1 )
#define SACI_BIAS(C)        ( SACI_CONTENT(C)->bias )
#define SACI_P(C)           ( SACI_CONTENT(C)->p )
#define SACI_PQ(C)          ( SACI_CONTENT(C)->pq )
#define SACI_ADJ(C)         ( SACI_CONTENT(C)->adj )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1     RCONST(1.0)
#define DEFAULT_BIAS   RCONST(1.5)
#define DEFAULT_PQ     0
#define DEFAULT_ADJ    -1
#define TINY           RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new I controller
 */

SUNAdaptController SUNAdaptController_I(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptControllerContent_I content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype               = SUNAdaptController_GetType_I;
  C->ops->estimatestep          = SUNAdaptController_EstimateStep_I;
  C->ops->setdefaults           = SUNAdaptController_SetDefaults_I;
  C->ops->write                 = SUNAdaptController_Write_I;
  C->ops->setmethodorder        = SUNAdaptController_SetMethodOrder_I;
  C->ops->adjustcontrollerorder = SUNAdaptController_AdjustControllerOrder_I;
  C->ops->seterrorbias          = SUNAdaptController_SetErrorBias_I;
  C->ops->space                 = SUNAdaptController_Space_I;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_I)malloc(sizeof *content);
  if (content == NULL)
  {
    (void) SUNAdaptController_Destroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Initialize method order */
  content->p = 1;

  /* Fill content with default values */
  SUNAdaptController_SetDefaults_I(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set I parameters
 */

int SUNAdaptController_SetParams_I(SUNAdaptController C, int pq,
                                   sunrealtype k1)
{
  /* store legal inputs, and return with success */
  SACI_PQ(C) = pq;
  if (k1 >= RCONST(0.0)) { SACI_K1(C) = k1; }
  return SUNADAPTCONTROLLER_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_I(SUNAdaptController C)
{ return SUN_ADAPTCONTROLLER_H; }

int SUNAdaptController_EstimateStep_I(SUNAdaptController C, sunrealtype h,
                                      sunrealtype dsm, sunrealtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const int ord = SACI_P(C) + SACI_ADJ(C);
  const sunrealtype k1 = -SACI_K1(C) / ord;
  const sunrealtype ecur = SACI_BIAS(C) * dsm;
  const sunrealtype e1 = SUNMAX(ecur, TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = h * SUNRpowerR(e1,k1);
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetDefaults_I(SUNAdaptController C)
{
  SACI_K1(C)     = DEFAULT_K1;
  SACI_BIAS(C)   = DEFAULT_BIAS;
  SACI_PQ(C)     = DEFAULT_PQ;
  SACI_ADJ(C)    = DEFAULT_ADJ;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Write_I(SUNAdaptController C, FILE *fptr)
{
  fprintf(fptr, "I SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SACI_K1(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SACI_BIAS(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SACI_K1(C));
  fprintf(fptr, "  bias factor = %16g\n", SACI_BIAS(C));
#endif
  if (SACI_PQ(C) == 1)
  {
    fprintf(fptr, "  p = %i (method order)\n", SACI_P(C));
  }
  else if (SACI_PQ(C) == 0)
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SACI_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (minimum of method & embedding order)\n", SACI_P(C));
  }
  fprintf(fptr, "  adj = %i\n", SACI_ADJ(C));
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetMethodOrder_I(SUNAdaptController C, int p, int q)
{
  /* check for legal inputs */
  if ((p <= 0) || (q <= 0)) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store appropriate value based on "pq" */
  if (SACI_PQ(C) == 1) {
    SACI_P(C) = p;
  } else if (SACI_PQ(C) == 0) {
    SACI_P(C) = q;
  } else {
    SACI_P(C) = SUNMIN(p,q);
  }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_AdjustControllerOrder_I(SUNAdaptController C, int adj)
{
  SACI_ADJ(C) = adj;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetErrorBias_I(SUNAdaptController C, sunrealtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SACI_BIAS(C) = DEFAULT_BIAS;
  } else {
    SACI_BIAS(C) = bias;
  }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Space_I(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 3;
  *leniw = 3;
  return SUNADAPTCONTROLLER_SUCCESS;
}
