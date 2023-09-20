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

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1     RCONST(1.0)
#define DEFAULT_BIAS   RCONST(1.5)
#define DEFAULT_PQ     SUNFALSE
#define TINY           RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new I controller
 */

SUNAdaptController SUNAdaptControllerI(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptControllerContent_I content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype           = SUNAdaptControllerGetType_I;
  C->ops->estimatestep      = SUNAdaptControllerEstimateStep_I;
  C->ops->setdefaults       = SUNAdaptControllerSetDefaults_I;
  C->ops->write             = SUNAdaptControllerWrite_I;
  C->ops->setmethodorder    = SUNAdaptControllerSetMethodOrder_I;
  C->ops->setembeddingorder = SUNAdaptControllerSetEmbeddingOrder_I;
  C->ops->seterrorbias      = SUNAdaptControllerSetErrorBias_I;
  C->ops->space             = SUNAdaptControllerSpace_I;

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
  SUNAdaptControllerSetDefaults_I(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set I parameters
 */

int SUNAdaptControllerI_SetParams(SUNAdaptController C, sunbooleantype pq,
                                  realtype k1)
{
  /* store legal inputs, and return with success */
  SACI_PQ(C) = pq;
  if (k1 >= RCONST(0.0)) { SACI_K1(C) = k1; }
  return SUNADAPTCONTROLLER_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptControllerGetType_I(SUNAdaptController C)
{ return SUN_ADAPTCONTROLLER_H; }

int SUNAdaptControllerEstimateStep_I(SUNAdaptController C, realtype h,
                                     realtype dsm, realtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype k1 = -SACI_K1(C) / SACI_P(C);
  const realtype ecur = SACI_BIAS(C) * dsm;
  const realtype e1 = SUNMAX(ecur, TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = h * SUNRpowerR(e1,k1);
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetDefaults_I(SUNAdaptController C)
{
  SACI_K1(C)     = DEFAULT_K1;
  SACI_BIAS(C)   = DEFAULT_BIAS;
  SACI_PQ(C)     = DEFAULT_PQ;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerWrite_I(SUNAdaptController C, FILE *fptr)
{
  fprintf(fptr, "I SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SACI_K1(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SACI_BIAS(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SACI_K1(C));
  fprintf(fptr, "  bias factor = %16g\n", SACI_BIAS(C));
#endif
  if (SACI_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SACI_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SACI_P(C));
  }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetMethodOrder_I(SUNAdaptController C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SACI_PQ(C)) { SACI_P(C) = q; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetEmbeddingOrder_I(SUNAdaptController C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SACI_PQ(C)) { SACI_P(C) = p; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetErrorBias_I(SUNAdaptController C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SACI_BIAS(C) = DEFAULT_BIAS;
  } else {
    SACI_BIAS(C) = bias;
  }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSpace_I(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 3;
  *leniw = 2;
  return SUNADAPTCONTROLLER_SUCCESS;
}
