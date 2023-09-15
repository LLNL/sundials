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
 * This is the implementation file for the SUNControl_ImExGus module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_imexgus.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_IMEXGUS_CONTENT(C)     ( (SUNControlContent_ImExGus)(C->content) )
#define SC_IMEXGUS_K1E(C)         ( SC_IMEXGUS_CONTENT(C)->k1e )
#define SC_IMEXGUS_K2E(C)         ( SC_IMEXGUS_CONTENT(C)->k2e )
#define SC_IMEXGUS_K1I(C)         ( SC_IMEXGUS_CONTENT(C)->k1i )
#define SC_IMEXGUS_K2I(C)         ( SC_IMEXGUS_CONTENT(C)->k2i )
#define SC_IMEXGUS_BIAS(C)        ( SC_IMEXGUS_CONTENT(C)->bias )
#define SC_IMEXGUS_EP(C)          ( SC_IMEXGUS_CONTENT(C)->ep )
#define SC_IMEXGUS_HP(C)          ( SC_IMEXGUS_CONTENT(C)->hp )
#define SC_IMEXGUS_P(C)           ( SC_IMEXGUS_CONTENT(C)->p )
#define SC_IMEXGUS_PQ(C)          ( SC_IMEXGUS_CONTENT(C)->pq )
#define SC_IMEXGUS_FIRSTSTEP(C)   ( SC_IMEXGUS_CONTENT(C)->firststep )

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

SUNControl SUNControlImExGus(SUNContext sunctx)
{
  SUNControl C;
  SUNControlContent_ImExGus content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControl_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype           = SUNControlGetType_ImExGus;
  C->ops->estimatestep      = SUNControlEstimateStep_ImExGus;
  C->ops->reset             = SUNControlReset_ImExGus;
  C->ops->setdefaults       = SUNControlSetDefaults_ImExGus;
  C->ops->write             = SUNControlWrite_ImExGus;
  C->ops->setmethodorder    = SUNControlSetMethodOrder_ImExGus;
  C->ops->setembeddingorder = SUNControlSetEmbeddingOrder_ImExGus;
  C->ops->seterrorbias      = SUNControlSetErrorBias_ImExGus;
  C->ops->update            = SUNControlUpdate_ImExGus;
  C->ops->space             = SUNControlSpace_ImExGus;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_ImExGus)malloc(sizeof *content);
  if (content == NULL)
  {
    (void) SUNControl_Destroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Initialize method order */
  content->p = 1;

  /* Fill content with default/reset values */
  SUNControlSetDefaults_ImExGus(C);
  SUNControlReset_ImExGus(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set ImExGus parameters
 */

int SUNControlImExGus_SetParams(SUNControl C, sunbooleantype pq,
                                realtype k1e, realtype k2e,
                                realtype k1i, realtype k2i)
{
  /* store legal inputs, and return with success */
  SC_IMEXGUS_PQ(C) = pq;
  if (k1e >= RCONST(0.0)) { SC_IMEXGUS_K1E(C) = k1e; }
  if (k2e >= RCONST(0.0)) { SC_IMEXGUS_K2E(C) = k2e; }
  if (k1i >= RCONST(0.0)) { SC_IMEXGUS_K1I(C) = k1i; }
  if (k2i >= RCONST(0.0)) { SC_IMEXGUS_K2I(C) = k2i; }
  return SUNCONTROL_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_Type SUNControlGetType_ImExGus(SUNControl C) { return SUNDIALS_CONTROL_H; }

int SUNControlEstimateStep_ImExGus(SUNControl C, realtype h,
                                   realtype dsm, realtype* hnew)
{
  /* modified method for first step */
  if (SC_IMEXGUS_FIRSTSTEP(C))
  {
    /* set usable time-step adaptivity parameters */
    const realtype k = -RCONST(1.0) / SC_IMEXGUS_P(C);
    const realtype e = SUNMAX(SC_IMEXGUS_BIAS(C) * dsm, TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e,k);
  }
  else
  {
    /* set usable time-step adaptivity parameters */
    const realtype k1e = -SC_IMEXGUS_K1E(C) / SC_IMEXGUS_P(C);
    const realtype k2e = -SC_IMEXGUS_K2E(C) / SC_IMEXGUS_P(C);
    const realtype k1i = -SC_IMEXGUS_K1I(C) / SC_IMEXGUS_P(C);
    const realtype k2i = -SC_IMEXGUS_K2I(C) / SC_IMEXGUS_P(C);
    const realtype e1 = SUNMAX(SC_IMEXGUS_BIAS(C) * dsm, TINY);
    const realtype e2 = e1 / SUNMAX(SC_IMEXGUS_EP(C), TINY);
    const realtype hrat = h / SC_IMEXGUS_HP(C);

    /* compute estimated optimal time step size */
    *hnew = h * SUNMIN(hrat * SUNRpowerR(e1,k1i) * SUNRpowerR(e2,k2i),
                       SUNRpowerR(e1,k1e) * SUNRpowerR(e2,k2e));
  }

  /* return with success */
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_ImExGus(SUNControl C)
{
  SC_IMEXGUS_EP(C) = RCONST(1.0);
  SC_IMEXGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_ImExGus(SUNControl C)
{
  SC_IMEXGUS_K1E(C)    = DEFAULT_K1E;
  SC_IMEXGUS_K2E(C)    = DEFAULT_K2E;
  SC_IMEXGUS_K1I(C)    = DEFAULT_K1I;
  SC_IMEXGUS_K2I(C)    = DEFAULT_K2I;
  SC_IMEXGUS_BIAS(C)   = DEFAULT_BIAS;
  SC_IMEXGUS_PQ(C)     = DEFAULT_PQ;
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_ImExGus(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "ImEx Gustafsson SUNControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1e = %32Lg\n", SC_IMEXGUS_K1E(C));
  fprintf(fptr, "  k2e = %32Lg\n", SC_IMEXGUS_K2E(C));
  fprintf(fptr, "  k1i = %32Lg\n", SC_IMEXGUS_K1I(C));
  fprintf(fptr, "  k2i = %32Lg\n", SC_IMEXGUS_K2I(C));
  fprintf(fptr, "  bias factor = %22Lg\n", SC_IMEXGUS_BIAS(C));
  fprintf(fptr, "  previous error = %22Lg\n", SC_IMEXGUS_EP(C));
  fprintf(fptr, "  previous step = %22Lg\n", SC_IMEXGUS_HP(C));
#else
  fprintf(fptr, "  k1e = %16g\n", SC_IMEXGUS_K1E(C));
  fprintf(fptr, "  k2e = %16g\n", SC_IMEXGUS_K2E(C));
  fprintf(fptr, "  k1i = %16g\n", SC_IMEXGUS_K1I(C));
  fprintf(fptr, "  k2i = %16g\n", SC_IMEXGUS_K2I(C));
  fprintf(fptr, "  bias factor = %16g\n", SC_IMEXGUS_BIAS(C));
  fprintf(fptr, "  previous error = %16g\n", SC_IMEXGUS_EP(C));
  fprintf(fptr, "  previous step = %16g\n", SC_IMEXGUS_HP(C));
#endif
  if (SC_IMEXGUS_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SC_IMEXGUS_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SC_IMEXGUS_P(C));
  }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_ImExGus(SUNControl C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SC_IMEXGUS_PQ(C)) { SC_IMEXGUS_P(C) = q; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetEmbeddingOrder_ImExGus(SUNControl C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SC_IMEXGUS_PQ(C)) { SC_IMEXGUS_P(C) = p; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_ImExGus(SUNControl C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SC_IMEXGUS_BIAS(C) = DEFAULT_BIAS;
  } else {
    SC_IMEXGUS_BIAS(C) = bias;
  }

  return SUNCONTROL_SUCCESS;
}

int SUNControlUpdate_ImExGus(SUNControl C, realtype h, realtype dsm)
{
  SC_IMEXGUS_EP(C) = SC_IMEXGUS_BIAS(C) * dsm;
  SC_IMEXGUS_HP(C) = h;
  SC_IMEXGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_ImExGus(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 8;
  *leniw = 3;
  return SUNCONTROL_SUCCESS;
}
