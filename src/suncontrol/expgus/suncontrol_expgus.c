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
 * This is the implementation file for the SUNControl_ExpGus module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_expgus.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_EXPGUS_CONTENT(C)     ( (SUNControlContent_ExpGus)(C->content) )
#define SC_EXPGUS_K1(C)          ( SC_EXPGUS_CONTENT(C)->k1 )
#define SC_EXPGUS_K2(C)          ( SC_EXPGUS_CONTENT(C)->k2 )
#define SC_EXPGUS_BIAS(C)        ( SC_EXPGUS_CONTENT(C)->bias )
#define SC_EXPGUS_EP(C)          ( SC_EXPGUS_CONTENT(C)->ep )
#define SC_EXPGUS_P(C)           ( SC_EXPGUS_CONTENT(C)->p )
#define SC_EXPGUS_PQ(C)          ( SC_EXPGUS_CONTENT(C)->pq )
#define SC_EXPGUS_FIRSTSTEP(C)   ( SC_EXPGUS_CONTENT(C)->firststep )

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

SUNControl SUNControlExpGus(SUNContext sunctx)
{
  SUNControl C;
  SUNControlContent_ExpGus content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControl_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype           = SUNControlGetType_ExpGus;
  C->ops->estimatestep      = SUNControlEstimateStep_ExpGus;
  C->ops->reset             = SUNControlReset_ExpGus;
  C->ops->setdefaults       = SUNControlSetDefaults_ExpGus;
  C->ops->write             = SUNControlWrite_ExpGus;
  C->ops->setmethodorder    = SUNControlSetMethodOrder_ExpGus;
  C->ops->setembeddingorder = SUNControlSetEmbeddingOrder_ExpGus;
  C->ops->seterrorbias      = SUNControlSetErrorBias_ExpGus;
  C->ops->update            = SUNControlUpdate_ExpGus;
  C->ops->space             = SUNControlSpace_ExpGus;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_ExpGus)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNControl_Destroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Initialize method order */
  content->p = 1;

  /* Fill content with default/reset values */
  SUNControlSetDefaults_ExpGus(C);
  SUNControlReset_ExpGus(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set ExpGus parameters
 */

int SUNControlExpGus_SetParams(SUNControl C, sunbooleantype pq,
                               realtype k1, realtype k2)
{
  /* store legal inputs, and return with success */
  SC_EXPGUS_PQ(C) = pq;
  if (k1 >= RCONST(0.0)) { SC_EXPGUS_K1(C) = k1; }
  if (k2 >= RCONST(0.0)) { SC_EXPGUS_K2(C) = k2; }
  return SUNCONTROL_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_Type SUNControlGetType_ExpGus(SUNControl C) { return SUNDIALS_CONTROL_H; }

int SUNControlEstimateStep_ExpGus(SUNControl C, realtype h,
                               realtype dsm, realtype* hnew)
{
  /* modified method for first step */
  if (SC_EXPGUS_FIRSTSTEP(C))
  {
    /* set usable time-step adaptivity parameters */
    const realtype k = -RCONST(1.0) / SC_EXPGUS_P(C);
    const realtype e = SUNMAX(SC_EXPGUS_BIAS(C) * dsm, TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e,k);
  }
  else
  {
    /* set usable time-step adaptivity parameters */
    const realtype k1 = -SC_EXPGUS_K1(C) / SC_EXPGUS_P(C);
    const realtype k2 = -SC_EXPGUS_K2(C) / SC_EXPGUS_P(C);
    const realtype ecur = SC_EXPGUS_BIAS(C) * dsm;
    const realtype e1 = SUNMAX(ecur, TINY);
    const realtype e2 = e1 / SUNMAX(SC_EXPGUS_EP(C), TINY);

    /* compute estimated optimal time step size */
    *hnew = h * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  }

  /* return with success */
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_ExpGus(SUNControl C)
{
  SC_EXPGUS_EP(C) = RCONST(1.0);
  SC_EXPGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_ExpGus(SUNControl C)
{
  SC_EXPGUS_K1(C)     = DEFAULT_K1;
  SC_EXPGUS_K2(C)     = DEFAULT_K2;
  SC_EXPGUS_BIAS(C)   = DEFAULT_BIAS;
  SC_EXPGUS_PQ(C)     = DEFAULT_PQ;
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_ExpGus(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "Explicit Gustafsson SUNControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SC_EXPGUS_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", SC_EXPGUS_K2(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SC_EXPGUS_BIAS(C));
  fprintf(fptr, "  previous error = %32g\n", SC_EXPGUS_EP(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SC_EXPGUS_K1(C));
  fprintf(fptr, "  k2 = %16g\n", SC_EXPGUS_K2(C));
  fprintf(fptr, "  bias factor = %16g\n", SC_EXPGUS_BIAS(C));
  fprintf(fptr, "  previous error = %16g\n", SC_EXPGUS_EP(C));
#endif
  if (SC_EXPGUS_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SC_EXPGUS_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SC_EXPGUS_P(C));
  }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_ExpGus(SUNControl C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SC_EXPGUS_PQ(C)) { SC_EXPGUS_P(C) = q; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetEmbeddingOrder_ExpGus(SUNControl C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SC_EXPGUS_PQ(C)) { SC_EXPGUS_P(C) = p; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_ExpGus(SUNControl C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SC_EXPGUS_BIAS(C) = DEFAULT_BIAS;
  } else {
    SC_EXPGUS_BIAS(C) = bias;
  }

  return SUNCONTROL_SUCCESS;
}

int SUNControlUpdate_ExpGus(SUNControl C, realtype h, realtype dsm)
{
  SC_EXPGUS_EP(C) = SC_EXPGUS_BIAS(C) * dsm;
  SC_EXPGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_ExpGus(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 5;
  *leniw = 3;
  return SUNCONTROL_SUCCESS;
}
