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
#define SC_EXPGUS_SAFETY(C)      ( SC_EXPGUS_CONTENT(C)->safety )
#define SC_EXPGUS_EP(C)          ( SC_EXPGUS_CONTENT(C)->ep )
#define SC_EXPGUS_P(C)           ( SC_EXPGUS_CONTENT(C)->p )
#define SC_EXPGUS_FIRSTSTEP(C)   ( SC_EXPGUS_CONTENT(C)->firststep )

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
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid           = SUNControlGetID_ExpGus;
  C->ops->destroy         = SUNControlDestroy_ExpGus;
  C->ops->estimatestep    = SUNControlEstimateStep_ExpGus;
  C->ops->reset           = SUNControlReset_ExpGus;
  C->ops->setdefaults     = SUNControlSetDefaults_ExpGus;
  C->ops->write           = SUNControlWrite_ExpGus;
  C->ops->setmethodorder  = SUNControlSetMethodOrder_ExpGus;
  C->ops->setsafetyfactor = SUNControlSetSafetyFactor_ExpGus;
  C->ops->seterrorbias    = SUNControlSetErrorBias_ExpGus;
  C->ops->update          = SUNControlUpdate_ExpGus;
  C->ops->space           = SUNControlSpace_ExpGus;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_ExpGus)malloc(sizeof *content);
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
  SUNControlSetDefaults_ExpGus(C);
  SUNControlReset_ExpGus(C);

  /* Attach content and return */
  return (C);
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_ExpGus(SUNControl C) { return SUNDIALS_CONTROL_H; }

void SUNControlDestroy_ExpGus(SUNControl C)
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

int SUNControlEstimateStep_ExpGus(SUNControl C, realtype h,
                               realtype dsm, realtype* hnew)
{
  /* modified method for first step */
  if (SC_EXPGUS_FIRSTSTEP(C))
  {
    /* set usable time-step adaptivity parameters */
    const realtype k = -RCONST(1.0) / SC_EXPGUS_P(C);
    const realtype e = SUNMAX(SC_EXPGUS_BIAS(C) * dsm, RCONST(1.0e-10));

    /* compute estimated optimal time step size */
    *hnew = SC_EXPGUS_SAFETY(C) * h * SUNRpowerR(e,k);
  }
  else
  {
    /* set usable time-step adaptivity parameters */
    const realtype k1 = -SC_EXPGUS_K1(C) / SC_EXPGUS_P(C);
    const realtype k2 = -SC_EXPGUS_K2(C) / SC_EXPGUS_P(C);
    const realtype ecur = SC_EXPGUS_BIAS(C) * dsm;
    const realtype TINY = RCONST(1.0e-10);
    const realtype e1 = SUNMAX(ecur, TINY);
    const realtype e2 = e1 / SUNMAX(SC_EXPGUS_EP(C), TINY);

    /* compute estimated optimal time step size */
    *hnew = SC_EXPGUS_SAFETY(C) * h * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  }

  /* return with success */
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_ExpGus(SUNControl C)
{
  SC_EXPGUS_EP(C) = RCONST(1.0);
  SC_EXPGUS_FIRSTSTEP(C) = SUNTRUE;
}

int SUNControlSetDefaults_ExpGus(SUNControl C)
{
  SC_EXPGUS_K1(C)     = RCONST(0.367);
  SC_EXPGUS_K2(C)     = RCONST(0.268);
  SC_EXPGUS_BIAS(C)   = RCONST(1.5);
  SC_EXPGUS_SAFETY(C) = RCONST(0.96);
}

int SUNControlWrite_ExpGus(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "SUNControl_ExpGus module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %12Lg", SC_EXPGUS_K1(C));
  fprintf(fptr, "  k2 = %12Lg", SC_EXPGUS_K2(C));
  fprintf(fptr, "  bias = %12Lg", SC_EXPGUS_BIAS(C));
  fprintf(fptr, "  safety = %12Lg", SC_EXPGUS_SAFETY(C));
  fprintf(fptr, "  ep = %12Lg", SC_EXPGUS_EP(C));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  fprintf(fptr, "  k1 = %12g", SC_EXPGUS_K1(C));
  fprintf(fptr, "  k2 = %12g", SC_EXPGUS_K2(C));
  fprintf(fptr, "  bias = %12g", SC_EXPGUS_BIAS(C));
  fprintf(fptr, "  safety = %12g", SC_EXPGUS_SAFETY(C));
  fprintf(fptr, "  ep = %12g", SC_EXPGUS_EP(C));
#else
  fprintf(fptr, "  k1 = %12g", SC_EXPGUS_K1(C));
  fprintf(fptr, "  k2 = %12g", SC_EXPGUS_K2(C));
  fprintf(fptr, "  bias = %12g", SC_EXPGUS_BIAS(C));
  fprintf(fptr, "  safety = %12g", SC_EXPGUS_SAFETY(C));
  fprintf(fptr, "  ep = %12g", SC_EXPGUS_EP(C));
#endif
  fprintf(fptr, "  p = %i", SC_EXPGUS_P(C));
  fprintf(fptr, "  firststep = %i", SC_EXPGUS_FIRSTSTEP(C));
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_ExpGus(SUNControl C, int p)
{
  SC_EXPGUS_P(C) = p;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetSafetyFactor_ExpGus(SUNControl C, realtype safety)
{
  SC_EXPGUS_SAFETY(C) = safety;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_ExpGus(SUNControl C, realtype bias)
{
  SC_EXPGUS_BIAS(C) = bias;
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
  *leniw = 2;
  return SUNCONTROL_SUCCESS;
}
