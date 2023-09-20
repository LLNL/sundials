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
#define SACEXPGUS_FIRSTSTEP(C)   ( SACEXPGUS_CONTENT(C)->firststep )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1     RCONST(0.367)
#define DEFAULT_K2     RCONST(0.268)
#define DEFAULT_BIAS   RCONST(1.5)
#define DEFAULT_PQ     SUNFALSE
#define TINY           RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new ExpGus controller
 */

SUNAdaptController SUNAdaptControllerExpGus(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptControllerContent_ExpGus content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype           = SUNAdaptControllerGetType_ExpGus;
  C->ops->estimatestep      = SUNAdaptControllerEstimateStep_ExpGus;
  C->ops->reset             = SUNAdaptControllerReset_ExpGus;
  C->ops->setdefaults       = SUNAdaptControllerSetDefaults_ExpGus;
  C->ops->write             = SUNAdaptControllerWrite_ExpGus;
  C->ops->setmethodorder    = SUNAdaptControllerSetMethodOrder_ExpGus;
  C->ops->setembeddingorder = SUNAdaptControllerSetEmbeddingOrder_ExpGus;
  C->ops->seterrorbias      = SUNAdaptControllerSetErrorBias_ExpGus;
  C->ops->update            = SUNAdaptControllerUpdate_ExpGus;
  C->ops->space             = SUNAdaptControllerSpace_ExpGus;

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
  SUNAdaptControllerSetDefaults_ExpGus(C);
  SUNAdaptControllerReset_ExpGus(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set ExpGus parameters
 */

int SUNAdaptControllerExpGus_SetParams(SUNAdaptController C, sunbooleantype pq,
                                       realtype k1, realtype k2)
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

SUNAdaptController_Type SUNAdaptControllerGetType_ExpGus(SUNAdaptController C)
{ return SUN_ADAPTCONTROLLER_H; }

int SUNAdaptControllerEstimateStep_ExpGus(SUNAdaptController C, realtype h,
                                          realtype dsm, realtype* hnew)
{
  /* modified method for first step */
  if (SACEXPGUS_FIRSTSTEP(C))
  {
    /* set usable time-step adaptivity parameters */
    const realtype k = -RCONST(1.0) / SACEXPGUS_P(C);
    const realtype e = SUNMAX(SACEXPGUS_BIAS(C) * dsm, TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e,k);
  }
  else
  {
    /* set usable time-step adaptivity parameters */
    const realtype k1 = -SACEXPGUS_K1(C) / SACEXPGUS_P(C);
    const realtype k2 = -SACEXPGUS_K2(C) / SACEXPGUS_P(C);
    const realtype ecur = SACEXPGUS_BIAS(C) * dsm;
    const realtype e1 = SUNMAX(ecur, TINY);
    const realtype e2 = e1 / SUNMAX(SACEXPGUS_EP(C), TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  }

  /* return with success */
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerReset_ExpGus(SUNAdaptController C)
{
  SACEXPGUS_EP(C) = RCONST(1.0);
  SACEXPGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetDefaults_ExpGus(SUNAdaptController C)
{
  SACEXPGUS_K1(C)     = DEFAULT_K1;
  SACEXPGUS_K2(C)     = DEFAULT_K2;
  SACEXPGUS_BIAS(C)   = DEFAULT_BIAS;
  SACEXPGUS_PQ(C)     = DEFAULT_PQ;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerWrite_ExpGus(SUNAdaptController C, FILE *fptr)
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
  if (SACEXPGUS_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SACEXPGUS_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SACEXPGUS_P(C));
  }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetMethodOrder_ExpGus(SUNAdaptController C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SACEXPGUS_PQ(C)) { SACEXPGUS_P(C) = q; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetEmbeddingOrder_ExpGus(SUNAdaptController C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SACEXPGUS_PQ(C)) { SACEXPGUS_P(C) = p; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetErrorBias_ExpGus(SUNAdaptController C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SACEXPGUS_BIAS(C) = DEFAULT_BIAS;
  } else {
    SACEXPGUS_BIAS(C) = bias;
  }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerUpdate_ExpGus(SUNAdaptController C, realtype h, realtype dsm)
{
  SACEXPGUS_EP(C) = SACEXPGUS_BIAS(C) * dsm;
  SACEXPGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSpace_ExpGus(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 5;
  *leniw = 3;
  return SUNADAPTCONTROLLER_SUCCESS;
}
