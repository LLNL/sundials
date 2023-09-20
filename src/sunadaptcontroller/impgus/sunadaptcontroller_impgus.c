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
#define SACIMPGUS_P(C)           ( SACIMPGUS_CONTENT(C)->p )
#define SACIMPGUS_PQ(C)          ( SACIMPGUS_CONTENT(C)->pq )
#define SACIMPGUS_FIRSTSTEP(C)   ( SACIMPGUS_CONTENT(C)->firststep )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1     RCONST(0.98)
#define DEFAULT_K2     RCONST(0.95)
#define DEFAULT_BIAS   RCONST(1.5)
#define DEFAULT_PQ     SUNFALSE
#define TINY           RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new ImpGus controller
 */

SUNAdaptController SUNAdaptControllerImpGus(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptControllerContent_ImpGus content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype           = SUNAdaptControllerGetType_ImpGus;
  C->ops->estimatestep      = SUNAdaptControllerEstimateStep_ImpGus;
  C->ops->reset             = SUNAdaptControllerReset_ImpGus;
  C->ops->setdefaults       = SUNAdaptControllerSetDefaults_ImpGus;
  C->ops->write             = SUNAdaptControllerWrite_ImpGus;
  C->ops->setmethodorder    = SUNAdaptControllerSetMethodOrder_ImpGus;
  C->ops->setembeddingorder = SUNAdaptControllerSetEmbeddingOrder_ImpGus;
  C->ops->seterrorbias      = SUNAdaptControllerSetErrorBias_ImpGus;
  C->ops->update            = SUNAdaptControllerUpdate_ImpGus;
  C->ops->space             = SUNAdaptControllerSpace_ImpGus;

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

  /* Initialize method order */
  content->p = 1;

  /* Fill content with default/reset values */
  SUNAdaptControllerSetDefaults_ImpGus(C);
  SUNAdaptControllerReset_ImpGus(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set ImpGus parameters
 */

int SUNAdaptControllerImpGus_SetParams(SUNAdaptController C, sunbooleantype pq,
                                       realtype k1, realtype k2)
{
  /* store legal inputs, and return with success */
  SACIMPGUS_PQ(C) = pq;
  if (k1 >= RCONST(0.0)) { SACIMPGUS_K1(C) = k1; }
  if (k2 >= RCONST(0.0)) { SACIMPGUS_K2(C) = k2; }
  return SUNADAPTCONTROLLER_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptControllerGetType_ImpGus(SUNAdaptController C)
{ return SUN_ADAPTCONTROLLER_H; }

int SUNAdaptControllerEstimateStep_ImpGus(SUNAdaptController C, realtype h,
                                          realtype dsm, realtype* hnew)
{
  /* modified method for first step */
  if (SACIMPGUS_FIRSTSTEP(C))
  {
    /* set usable time-step adaptivity parameters */
    const realtype k = -RCONST(1.0) / SACIMPGUS_P(C);
    const realtype e = SUNMAX(SACIMPGUS_BIAS(C) * dsm, TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e,k);
  }
  else
  {
    /* set usable time-step adaptivity parameters */
    const realtype k1 = -SACIMPGUS_K1(C) / SACIMPGUS_P(C);
    const realtype k2 = -SACIMPGUS_K2(C) / SACIMPGUS_P(C);
    const realtype ecur = SACIMPGUS_BIAS(C) * dsm;
    const realtype e1 = SUNMAX(ecur, TINY);
    const realtype e2 = e1 / SUNMAX(SACIMPGUS_EP(C), TINY);
    const realtype hrat = h / SACIMPGUS_HP(C);

    /* compute estimated optimal time step size */
    *hnew = h * hrat * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  }

  /* return with success */
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerReset_ImpGus(SUNAdaptController C)
{
  SACIMPGUS_EP(C) = RCONST(1.0);
  SACIMPGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetDefaults_ImpGus(SUNAdaptController C)
{
  SACIMPGUS_K1(C)     = DEFAULT_K1;
  SACIMPGUS_K2(C)     = DEFAULT_K2;
  SACIMPGUS_BIAS(C)   = DEFAULT_BIAS;
  SACIMPGUS_PQ(C)     = DEFAULT_PQ;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerWrite_ImpGus(SUNAdaptController C, FILE *fptr)
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
  if (SACIMPGUS_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SACIMPGUS_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SACIMPGUS_P(C));
  }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetMethodOrder_ImpGus(SUNAdaptController C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SACIMPGUS_PQ(C)) { SACIMPGUS_P(C) = q; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetEmbeddingOrder_ImpGus(SUNAdaptController C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SACIMPGUS_PQ(C)) { SACIMPGUS_P(C) = p; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetErrorBias_ImpGus(SUNAdaptController C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SACIMPGUS_BIAS(C) = DEFAULT_BIAS;
  } else {
    SACIMPGUS_BIAS(C) = bias;
  }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerUpdate_ImpGus(SUNAdaptController C, realtype h, realtype dsm)
{
  SACIMPGUS_EP(C) = SACIMPGUS_BIAS(C) * dsm;
  SACIMPGUS_HP(C) = h;
  SACIMPGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSpace_ImpGus(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 6;
  *leniw = 3;
  return SUNADAPTCONTROLLER_SUCCESS;
}
