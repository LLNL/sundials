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
 * This is the implementation file for the SUNControl_I module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_i.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_I_CONTENT(C)     ( (SUNControlContent_I)(C->content) )
#define SC_I_K1(C)          ( SC_I_CONTENT(C)->k1 )
#define SC_I_BIAS(C)        ( SC_I_CONTENT(C)->bias )
#define SC_I_P(C)           ( SC_I_CONTENT(C)->p )
#define SC_I_PQ(C)          ( SC_I_CONTENT(C)->pq )

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

SUNControl SUNControlI(SUNContext sunctx)
{
  SUNControl C;
  SUNControlContent_I content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid             = SUNControlGetType_I;
  C->ops->estimatestep      = SUNControlEstimateStep_I;
  C->ops->setdefaults       = SUNControlSetDefaults_I;
  C->ops->write             = SUNControlWrite_I;
  C->ops->setmethodorder    = SUNControlSetMethodOrder_I;
  C->ops->setembeddingorder = SUNControlSetEmbeddingOrder_I;
  C->ops->seterrorbias      = SUNControlSetErrorBias_I;
  C->ops->space             = SUNControlSpace_I;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_I)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNControlDestroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Initialize method order */
  content->p = 1;

  /* Fill content with default values */
  SUNControlSetDefaults_I(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set I parameters
 */

int SUNControlI_SetParams(SUNControl C, sunbooleantype pq,
                          realtype k1)
{
  /* store legal inputs, and return with success */
  SC_I_PQ(C) = pq;
  if (k1 >= RCONST(0.0)) { SC_I_K1(C) = k1; }
  return SUNCONTROL_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_Type SUNControlGetType_I(SUNControl C) { return SUNDIALS_CONTROL_H; }

int SUNControlEstimateStep_I(SUNControl C, realtype h,
                               realtype dsm, realtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype k1 = -SC_I_K1(C) / SC_I_P(C);
  const realtype ecur = SC_I_BIAS(C) * dsm;
  const realtype e1 = SUNMAX(ecur, TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = h * SUNRpowerR(e1,k1);
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_I(SUNControl C)
{
  SC_I_K1(C)     = DEFAULT_K1;
  SC_I_BIAS(C)   = DEFAULT_BIAS;
  SC_I_PQ(C)     = DEFAULT_PQ;
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_I(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "I-controller SUNControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SC_I_K1(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SC_I_BIAS(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SC_I_K1(C));
  fprintf(fptr, "  bias factor = %16g\n", SC_I_BIAS(C));
#endif
  if (SC_I_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SC_I_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SC_I_P(C));
  }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_I(SUNControl C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SC_I_PQ(C)) { SC_I_P(C) = q; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetEmbeddingOrder_I(SUNControl C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SC_I_PQ(C)) { SC_I_P(C) = p; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_I(SUNControl C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SC_I_BIAS(C) = DEFAULT_BIAS;
  } else {
    SC_I_BIAS(C) = bias;
  }

  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_I(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 3;
  *leniw = 2;
  return SUNCONTROL_SUCCESS;
}
