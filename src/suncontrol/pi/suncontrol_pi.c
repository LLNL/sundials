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
 * This is the implementation file for the SUNControl_PI module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_pi.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_PI_CONTENT(C)     ( (SUNControlContent_PI)(C->content) )
#define SC_PI_K1(C)          ( SC_PI_CONTENT(C)->k1 )
#define SC_PI_K2(C)          ( SC_PI_CONTENT(C)->k2 )
#define SC_PI_BIAS(C)        ( SC_PI_CONTENT(C)->bias )
#define SC_PI_EP(C)          ( SC_PI_CONTENT(C)->ep )
#define SC_PI_P(C)           ( SC_PI_CONTENT(C)->p )
#define SC_PI_PQ(C)          ( SC_PI_CONTENT(C)->pq )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1     RCONST(0.8)
#define DEFAULT_K2     RCONST(0.31)
#define DEFAULT_BIAS   RCONST(1.5)
#define DEFAULT_PQ     SUNFALSE
#define TINY           RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new PI controller
 */

SUNControl SUNControlPI(SUNContext sunctx)
{
  SUNControl C;
  SUNControlContent_PI content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid             = SUNControlGetType_PI;
  C->ops->estimatestep      = SUNControlEstimateStep_PI;
  C->ops->reset             = SUNControlReset_PI;
  C->ops->setdefaults       = SUNControlSetDefaults_PI;
  C->ops->write             = SUNControlWrite_PI;
  C->ops->setmethodorder    = SUNControlSetMethodOrder_PI;
  C->ops->setembeddingorder = SUNControlSetEmbeddingOrder_PI;
  C->ops->seterrorbias      = SUNControlSetErrorBias_PI;
  C->ops->update            = SUNControlUpdate_PI;
  C->ops->space             = SUNControlSpace_PI;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_PI)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNControlDestroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Initialize method order */
  content->p = 1;

  /* Fill content with default/reset values */
  SUNControlSetDefaults_PI(C);
  SUNControlReset_PI(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PI parameters
 */

int SUNControlPI_SetParams(SUNControl C, sunbooleantype pq,
                           realtype k1, realtype k2)
{
  /* store legal inputs, and return with success */
  SC_PI_PQ(C) = pq;
  if (k1 >= RCONST(0.0)) { SC_PI_K1(C) = k1; }
  if (k2 >= RCONST(0.0)) { SC_PI_K2(C) = k2; }
  return SUNCONTROL_SUCCESS;
}


/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_Type SUNControlGetType_PI(SUNControl C) { return SUNDIALS_CONTROL_H; }

int SUNControlEstimateStep_PI(SUNControl C, realtype h,
                               realtype dsm, realtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype k1 = -SC_PI_K1(C) / SC_PI_P(C);
  const realtype k2 =  SC_PI_K2(C) / SC_PI_P(C);
  const realtype ecur = SC_PI_BIAS(C) * dsm;
  const realtype e1 = SUNMAX(ecur, TINY);
  const realtype e2 = SUNMAX(SC_PI_EP(C), TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = h * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_PI(SUNControl C)
{
  SC_PI_EP(C) = RCONST(1.0);
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_PI(SUNControl C)
{
  SC_PI_K1(C)     = DEFAULT_K1;
  SC_PI_K2(C)     = DEFAULT_K2;
  SC_PI_BIAS(C)   = DEFAULT_BIAS;
  SC_PI_PQ(C)     = DEFAULT_PQ;
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_PI(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "PI-controller SUNControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SC_PI_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", SC_PI_K2(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SC_PI_BIAS(C));
  fprintf(fptr, "  previous error = %32Lg\n", SC_PI_EP(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SC_PI_K1(C));
  fprintf(fptr, "  k2 = %16g\n", SC_PI_K2(C));
  fprintf(fptr, "  bias factor = %16g\n", SC_PI_BIAS(C));
  fprintf(fptr, "  previous error = %16g\n", SC_PI_EP(C));
#endif
  if (SC_PI_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SC_PI_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SC_PI_P(C));
  }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_PI(SUNControl C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SC_PI_PQ(C)) { SC_PI_P(C) = q; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetEmbeddingOrder_PI(SUNControl C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SC_PI_PQ(C)) { SC_PI_P(C) = p; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_PI(SUNControl C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SC_PI_BIAS(C) = DEFAULT_BIAS;
  } else {
    SC_PI_BIAS(C) = bias;
  }

  return SUNCONTROL_SUCCESS;
}

int SUNControlUpdate_PI(SUNControl C, realtype h, realtype dsm)
{
  SC_PI_EP(C) = SC_PI_BIAS(C) * dsm;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_PI(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 5;
  *leniw = 2;
  return SUNCONTROL_SUCCESS;
}
