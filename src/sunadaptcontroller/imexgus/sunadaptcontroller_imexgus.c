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

#define SACIMEXGUS_CONTENT(C)     ( (SUNAdaptControllerContent_ImExGus)(C->content) )
#define SACIMEXGUS_K1E(C)         ( SACIMEXGUS_CONTENT(C)->k1e )
#define SACIMEXGUS_K2E(C)         ( SACIMEXGUS_CONTENT(C)->k2e )
#define SACIMEXGUS_K1I(C)         ( SACIMEXGUS_CONTENT(C)->k1i )
#define SACIMEXGUS_K2I(C)         ( SACIMEXGUS_CONTENT(C)->k2i )
#define SACIMEXGUS_BIAS(C)        ( SACIMEXGUS_CONTENT(C)->bias )
#define SACIMEXGUS_EP(C)          ( SACIMEXGUS_CONTENT(C)->ep )
#define SACIMEXGUS_HP(C)          ( SACIMEXGUS_CONTENT(C)->hp )
#define SACIMEXGUS_P(C)           ( SACIMEXGUS_CONTENT(C)->p )
#define SACIMEXGUS_PQ(C)          ( SACIMEXGUS_CONTENT(C)->pq )
#define SACIMEXGUS_FIRSTSTEP(C)   ( SACIMEXGUS_CONTENT(C)->firststep )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1E    RCONST(0.367)
#define DEFAULT_K2E    RCONST(0.268)
#define DEFAULT_K1I    RCONST(0.95)
#define DEFAULT_K2I    RCONST(0.95)
#define DEFAULT_BIAS   RCONST(1.5)
#define DEFAULT_PQ     SUNFALSE
#define TINY           RCONST(1.0e-10)


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

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype           = SUNAdaptController_GetType_ImExGus;
  C->ops->estimatestep      = SUNAdaptController_EstimateStep_ImExGus;
  C->ops->reset             = SUNAdaptController_Reset_ImExGus;
  C->ops->setdefaults       = SUNAdaptController_SetDefaults_ImExGus;
  C->ops->write             = SUNAdaptController_Write_ImExGus;
  C->ops->setmethodorder    = SUNAdaptController_SetMethodOrder_ImExGus;
  C->ops->setembeddingorder = SUNAdaptController_SetEmbeddingOrder_ImExGus;
  C->ops->seterrorbias      = SUNAdaptController_SetErrorBias_ImExGus;
  C->ops->update            = SUNAdaptController_Update_ImExGus;
  C->ops->space             = SUNAdaptController_Space_ImExGus;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_ImExGus)malloc(sizeof *content);
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
  SUNAdaptController_SetDefaults_ImExGus(C);
  SUNAdaptController_Reset_ImExGus(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set ImExGus parameters
 */

int SUNAdaptController_SetParams_ImExGus(SUNAdaptController C, sunbooleantype pq,
                                         realtype k1e, realtype k2e,
                                         realtype k1i, realtype k2i)
{
  /* store legal inputs, and return with success */
  SACIMEXGUS_PQ(C) = pq;
  if (k1e >= RCONST(0.0)) { SACIMEXGUS_K1E(C) = k1e; }
  if (k2e >= RCONST(0.0)) { SACIMEXGUS_K2E(C) = k2e; }
  if (k1i >= RCONST(0.0)) { SACIMEXGUS_K1I(C) = k1i; }
  if (k2i >= RCONST(0.0)) { SACIMEXGUS_K2I(C) = k2i; }
  return SUNADAPTCONTROLLER_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_ImExGus(SUNAdaptController C)
{ return SUN_ADAPTCONTROLLER_H; }

int SUNAdaptController_EstimateStep_ImExGus(SUNAdaptController C, realtype h,
                                            realtype dsm, realtype* hnew)
{
  /* modified method for first step */
  if (SACIMEXGUS_FIRSTSTEP(C))
  {
    /* set usable time-step adaptivity parameters */
    const realtype k = -RCONST(1.0) / SACIMEXGUS_P(C);
    const realtype e = SUNMAX(SACIMEXGUS_BIAS(C) * dsm, TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e,k);
  }
  else
  {
    /* set usable time-step adaptivity parameters */
    const realtype k1e = -SACIMEXGUS_K1E(C) / SACIMEXGUS_P(C);
    const realtype k2e = -SACIMEXGUS_K2E(C) / SACIMEXGUS_P(C);
    const realtype k1i = -SACIMEXGUS_K1I(C) / SACIMEXGUS_P(C);
    const realtype k2i = -SACIMEXGUS_K2I(C) / SACIMEXGUS_P(C);
    const realtype e1 = SUNMAX(SACIMEXGUS_BIAS(C) * dsm, TINY);
    const realtype e2 = e1 / SUNMAX(SACIMEXGUS_EP(C), TINY);
    const realtype hrat = h / SACIMEXGUS_HP(C);

    /* compute estimated optimal time step size */
    *hnew = h * SUNMIN(hrat * SUNRpowerR(e1,k1i) * SUNRpowerR(e2,k2i),
                       SUNRpowerR(e1,k1e) * SUNRpowerR(e2,k2e));
  }

  /* return with success */
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Reset_ImExGus(SUNAdaptController C)
{
  SACIMEXGUS_EP(C) = RCONST(1.0);
  SACIMEXGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetDefaults_ImExGus(SUNAdaptController C)
{
  SACIMEXGUS_K1E(C)    = DEFAULT_K1E;
  SACIMEXGUS_K2E(C)    = DEFAULT_K2E;
  SACIMEXGUS_K1I(C)    = DEFAULT_K1I;
  SACIMEXGUS_K2I(C)    = DEFAULT_K2I;
  SACIMEXGUS_BIAS(C)   = DEFAULT_BIAS;
  SACIMEXGUS_PQ(C)     = DEFAULT_PQ;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Write_ImExGus(SUNAdaptController C, FILE *fptr)
{
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
  if (SACIMEXGUS_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SACIMEXGUS_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SACIMEXGUS_P(C));
  }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetMethodOrder_ImExGus(SUNAdaptController C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SACIMEXGUS_PQ(C)) { SACIMEXGUS_P(C) = q; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetEmbeddingOrder_ImExGus(SUNAdaptController C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SACIMEXGUS_PQ(C)) { SACIMEXGUS_P(C) = p; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_SetErrorBias_ImExGus(SUNAdaptController C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SACIMEXGUS_BIAS(C) = DEFAULT_BIAS;
  } else {
    SACIMEXGUS_BIAS(C) = bias;
  }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Update_ImExGus(SUNAdaptController C, realtype h, realtype dsm)
{
  SACIMEXGUS_EP(C) = SACIMEXGUS_BIAS(C) * dsm;
  SACIMEXGUS_HP(C) = h;
  SACIMEXGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptController_Space_ImExGus(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 8;
  *leniw = 3;
  return SUNADAPTCONTROLLER_SUCCESS;
}
