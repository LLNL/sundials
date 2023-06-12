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
#define SC_I_K2(C)          ( SC_I_CONTENT(C)->k2 )
#define SC_I_BIAS(C)        ( SC_I_CONTENT(C)->bias )
#define SC_I_SAFETY(C)      ( SC_I_CONTENT(C)->safety )
#define SC_I_EP(C)          ( SC_I_CONTENT(C)->ep )
#define SC_I_P(C)           ( SC_I_CONTENT(C)->p )

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
  C->ops->getid           = SUNControlGetID_I;
  C->ops->destroy         = SUNControlDestroy_I;
  C->ops->estimatestep    = SUNControlEstimateStep_I;
  C->ops->setdefaults     = SUNControlSetDefaults_I;
  C->ops->write           = SUNControlWrite_I;
  C->ops->setmethodorder  = SUNControlSetMethodOrder_I;
  C->ops->setsafetyfactor = SUNControlSetSafetyFactor_I;
  C->ops->seterrorbias    = SUNControlSetErrorBias_I;
  C->ops->space           = SUNControlSpace_I;

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

  /* Attach content and return */
  return (C);
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_I(SUNControl C) { return SUNDIALS_CONTROL_H; }

void SUNControlDestroy_I(SUNControl C)
{
  if (C == NULL) { return; }

  /* free content */
  if (C->content != NULL) {
    free(C->content);
    C->content = NULL;
  }

  /* free ops and controller */
  if (C->ops) {
    free(C->ops);
    C->ops = NULL;
  }
  free(C);
  C = NULL;

  return;
}

int SUNControlEstimateStep_I(SUNControl C, realtype h,
                               realtype dsm, realtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype k1 = -SC_I_K1(C) / SC_I_P(C);
  const realtype ecur = SC_I_BIAS(C) * dsm;
  const realtype TINY = RCONST(1.0e-10);
  const realtype e1 = SUNMAX(ecur, TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = SC_I_SAFETY(C) * h * SUNRpowerR(e1,k1);
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_I(SUNControl C)
{
  SC_I_K1(C)     = RCONST(1.0);
  SC_I_BIAS(C)   = RCONST(1.5);
  SC_I_SAFETY(C) = RCONST(0.96);
}

int SUNControlWrite_I(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "SUNControl_I module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %12Lg", SC_I_K1(C));
  fprintf(fptr, "  bias = %12Lg", SC_I_BIAS(C));
  fprintf(fptr, "  safety = %12Lg", SC_I_SAFETY(C));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  fprintf(fptr, "  k1 = %12g", SC_I_K1(C));
  fprintf(fptr, "  bias = %12g", SC_I_BIAS(C));
  fprintf(fptr, "  safety = %12g", SC_I_SAFETY(C));
#else
  fprintf(fptr, "  k1 = %12g", SC_I_K1(C));
  fprintf(fptr, "  bias = %12g", SC_I_BIAS(C));
  fprintf(fptr, "  safety = %12g", SC_I_SAFETY(C));
#endif
  fprintf(fptr, "  p = %i", SC_I_P(C));
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_I(SUNControl C, int p)
{
  SC_I_P(C) = p;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetSafetyFactor_I(SUNControl C, realtype safety)
{
  SC_I_SAFETY(C) = safety;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_I(SUNControl C, realtype bias)
{
  SC_I_BIAS(C) = bias;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_I(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 3;
  *leniw = 1;
  return SUNCONTROL_SUCCESS;
}
