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
#define SC_PI_SAFETY(C)      ( SC_PI_CONTENT(C)->safety )
#define SC_PI_EP(C)          ( SC_PI_CONTENT(C)->ep )
#define SC_PI_P(C)           ( SC_PI_CONTENT(C)->p )

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
  C->ops->getid           = SUNControlGetID_PI;
  C->ops->destroy         = SUNControlDestroy_PI;
  C->ops->estimatestep    = SUNControlEstimateStep_PI;
  C->ops->reset           = SUNControlReset_PI;
  C->ops->setdefaults     = SUNControlSetDefaults_PI;
  C->ops->write           = SUNControlWrite_PI;
  C->ops->setmethodorder  = SUNControlSetMethodOrder_PI;
  C->ops->setsafetyfactor = SUNControlSetSafetyFactor_PI;
  C->ops->seterrorbias    = SUNControlSetErrorBias_PI;
  C->ops->update          = SUNControlUpdate_PI;
  C->ops->space           = SUNControlSpace_PI;

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

  /* Attach content and return */
  return (C);
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_PI(SUNControl C) { return SUNDIALS_CONTROL_H; }

void SUNControlDestroy_PI(SUNControl C)
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

int SUNControlEstimateStep_PI(SUNControl C, realtype h,
                               realtype dsm, realtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype k1 = -SC_PI_K1(C) / SC_PI_P(C);
  const realtype k2 =  SC_PI_K2(C) / SC_PI_P(C);
  const realtype ecur = SC_PI_BIAS(C) * dsm;
  const realtype TINY = RCONST(1.0e-10);
  const realtype e1 = SUNMAX(ecur, TINY);
  const realtype e2 = SUNMAX(SC_PI_EP(C), TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = SC_PI_SAFETY(C) * h * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_PI(SUNControl C)
{
  SC_PI_K1(C)     = RCONST(0.8);
  SC_PI_K2(C)     = RCONST(0.31);
  SC_PI_BIAS(C)   = RCONST(1.5);
  SC_PI_SAFETY(C) = RCONST(0.96);
}

int SUNControlWrite_PI(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "SUNControl_PI module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %12Lg", SC_PI_K1(C));
  fprintf(fptr, "  k2 = %12Lg", SC_PI_K2(C));
  fprintf(fptr, "  bias = %12Lg", SC_PI_BIAS(C));
  fprintf(fptr, "  safety = %12Lg", SC_PI_SAFETY(C));
  fprintf(fptr, "  ep = %12Lg", SC_PI_EP(C));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  fprintf(fptr, "  k1 = %12g", SC_PI_K1(C));
  fprintf(fptr, "  k2 = %12g", SC_PI_K2(C));
  fprintf(fptr, "  bias = %12g", SC_PI_BIAS(C));
  fprintf(fptr, "  safety = %12g", SC_PI_SAFETY(C));
  fprintf(fptr, "  ep = %12g", SC_PI_EP(C));
#else
  fprintf(fptr, "  k1 = %12g", SC_PI_K1(C));
  fprintf(fptr, "  k2 = %12g", SC_PI_K2(C));
  fprintf(fptr, "  bias = %12g", SC_PI_BIAS(C));
  fprintf(fptr, "  safety = %12g", SC_PI_SAFETY(C));
  fprintf(fptr, "  ep = %12g", SC_PI_EP(C));
#endif
  fprintf(fptr, "  p = %i", SC_PI_P(C));
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_PI(SUNControl C, int p)
{
  SC_PI_P(C) = p;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetSafetyFactor_PI(SUNControl C, realtype safety)
{
  SC_PI_SAFETY(C) = safety;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_PI(SUNControl C, realtype bias)
{
  SC_PI_BIAS(C) = bias;
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
  *leniw = 1;
  return SUNCONTROL_SUCCESS;
}
