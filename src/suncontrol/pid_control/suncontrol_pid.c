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
 * This is the implementation file for the SUNControl_PID module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_pid.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_PID_CONTENT(C)     ( (SUNControlContent_PID)(C->content) )
#define SC_PID_K1(C)          ( SC_PID_CONTENT(C)->k1 )
#define SC_PID_K2(C)          ( SC_PID_CONTENT(C)->k2 )
#define SC_PID_K3(C)          ( SC_PID_CONTENT(C)->k3 )
#define SC_PID_BIAS(C)        ( SC_PID_CONTENT(C)->bias )
#define SC_PID_SAFETY(C)      ( SC_PID_CONTENT(C)->safety )
#define SC_PID_EP(C)          ( SC_PID_CONTENT(C)->ep )
#define SC_PID_EPP(C)         ( SC_PID_CONTENT(C)->epp )
#define SC_PID_P(C)           ( SC_PID_CONTENT(C)->p )

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new PID controller
 */

SUNControl SUNControlPID(SUNContext sunctx)
{
  SUNControl C;
  SUNControlContent_PID content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid           = SUNControlGetID_PID;
  C->ops->destroy         = SUNControlDestroy_PID;
  C->ops->estimatestep    = SUNControlEstimateStep_PID;
  C->ops->reset           = SUNControlReset_PID;
  C->ops->setdefaults     = SUNControlSetDefaults_PID;
  C->ops->write           = SUNControlWrite_PID;
  C->ops->setmethodorder  = SUNControlSetMethodOrder_PID;
  C->ops->setsafetyfactor = SUNControlSetSafetyFactor_PID;
  C->ops->seterrorbias    = SUNControlSetErrorBias_PID;
  C->ops->update          = SUNControlUpdate_PID;
  C->ops->space           = SUNControlSpace_PID;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_PID)malloc(sizeof *content);
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
  SUNControlSetDefaults_PID(C);
  SUNControlReset_PID(C);

  /* Attach content and return */
  return (C);
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_PID(SUNControl C) { return SUNDIALS_CONTROL_H; }

void SUNControlDestroy_PID(SUNControl C)
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

int SUNControlEstimateStep_PID(SUNControl C, realtype h,
                               realtype dsm, realtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype k1 = -SC_PID_K1(C) / SC_PID_P(C);
  const realtype k2 =  SC_PID_K2(C) / SC_PID_P(C);
  const realtype k3 = -SC_PID_K3(C) / SC_PID_P(C);
  const realtype ecur = SC_PID_BIAS(C) * dsm;
  const realtype TINY = RCONST(1.0e-10);
  const realtype e1 = SUNMAX(ecur, TINY);
  const realtype e2 = SUNMAX(SC_PID_EP(C), TINY);
  const realtype e3 = SUNMAX(SC_PID_EPP(C), TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = SC_PID_SAFETY(C) * h * SUNRpowerR(e1,k1)
          * SUNRpowerR(e2,k2) * SUNRpowerR(e3,k3);
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_PI(SUNControl C)
{
  SC_PID_EP(C)  = RCONST(1.0);
  SC_PID_EPP(C) = RCONST(1.0);
}

int SUNControlSetDefaults_PID(SUNControl C)
{
  SC_PID_K1(C)     = RCONST(0.58);
  SC_PID_K2(C)     = RCONST(0.21);
  SC_PID_K3(C)     = RCONST(0.1);
  SC_PID_BIAS(C)   = RCONST(1.5);
  SC_PID_SAFETY(C) = RCONST(0.96);
}

int SUNControlWrite_PID(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "SUNControl_PID module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %12Lg", SC_PID_K1(C));
  fprintf(fptr, "  k2 = %12Lg", SC_PID_K2(C));
  fprintf(fptr, "  k3 = %12Lg", SC_PID_K3(C));
  fprintf(fptr, "  bias = %12Lg", SC_PID_BIAS(C));
  fprintf(fptr, "  safety = %12Lg", SC_PID_SAFETY(C));
  fprintf(fptr, "  ep = %12Lg", SC_PID_EP(C));
  fprintf(fptr, "  epp = %12Lg", SC_PID_EPP(C));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  fprintf(fptr, "  k1 = %12g", SC_PID_K1(C));
  fprintf(fptr, "  k2 = %12g", SC_PID_K2(C));
  fprintf(fptr, "  k3 = %12g", SC_PID_K3(C));
  fprintf(fptr, "  bias = %12g", SC_PID_BIAS(C));
  fprintf(fptr, "  safety = %12g", SC_PID_SAFETY(C));
  fprintf(fptr, "  ep = %12g", SC_PID_EP(C));
  fprintf(fptr, "  epp = %12g", SC_PID_EPP(C));
#else
  fprintf(fptr, "  k1 = %12g", SC_PID_K1(C));
  fprintf(fptr, "  k2 = %12g", SC_PID_K2(C));
  fprintf(fptr, "  k3 = %12g", SC_PID_K3(C));
  fprintf(fptr, "  bias = %12g", SC_PID_BIAS(C));
  fprintf(fptr, "  safety = %12g", SC_PID_SAFETY(C));
  fprintf(fptr, "  ep = %12g", SC_PID_EP(C));
  fprintf(fptr, "  epp = %12g", SC_PID_EPP(C));
#endif
  fprintf(fptr, "  p = %i", SC_PID_P(C));
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_PID(SUNControl C, int p)
{
  SC_PID_P(C) = p;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetSafetyFactor_PID(SUNControl C, realtype safety)
{
  SC_PID_SAFETY(C) = safety;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_PID(SUNControl C, realtype bias)
{
  SC_PID_BIAS(C) = bias;
  return SUNCONTROL_SUCCESS;
}

int SUNControlUpdate_PID(SUNControl C, realtype h, realtype dsm)
{
  SC_PID_EPP(C) = SC_PID_EP(C);
  SC_PID_EP(C) = SC_PID_BIAS(C) * dsm;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_PID(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 7;
  *leniw = 1;
  return SUNCONTROL_SUCCESS;
}
