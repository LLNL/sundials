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
#define SC_PID_EP(C)          ( SC_PID_CONTENT(C)->ep )
#define SC_PID_EPP(C)         ( SC_PID_CONTENT(C)->epp )
#define SC_PID_P(C)           ( SC_PID_CONTENT(C)->p )
#define SC_PID_PQ(C)          ( SC_PID_CONTENT(C)->pq )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1     RCONST(0.58)
#define DEFAULT_K2     RCONST(0.21)
#define DEFAULT_K3     RCONST(0.1)
#define DEFAULT_BIAS   RCONST(1.5)
#define DEFAULT_PQ     SUNFALSE
#define TINY           RCONST(1.0e-10)


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
  C = SUNControl_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype           = SUNControlGetType_PID;
  C->ops->estimatestep      = SUNControlEstimateStep_PID;
  C->ops->reset             = SUNControlReset_PID;
  C->ops->setdefaults       = SUNControlSetDefaults_PID;
  C->ops->write             = SUNControlWrite_PID;
  C->ops->setmethodorder    = SUNControlSetMethodOrder_PID;
  C->ops->setembeddingorder = SUNControlSetEmbeddingOrder_PID;
  C->ops->seterrorbias      = SUNControlSetErrorBias_PID;
  C->ops->update            = SUNControlUpdate_PID;
  C->ops->space             = SUNControlSpace_PID;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_PID)malloc(sizeof *content);
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
  SUNControlSetDefaults_PID(C);
  SUNControlReset_PID(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PID parameters
 */

int SUNControlPID_SetParams(SUNControl C, sunbooleantype pq,
                            realtype k1, realtype k2, realtype k3)
{
  /* store legal inputs, and return with success */
  SC_PID_PQ(C) = pq;
  if (k1 >= RCONST(0.0)) { SC_PID_K1(C) = k1; }
  if (k2 >= RCONST(0.0)) { SC_PID_K2(C) = k2; }
  if (k3 >= RCONST(0.0)) { SC_PID_K3(C) = k3; }
  return SUNCONTROL_SUCCESS;
}


/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_Type SUNControlGetType_PID(SUNControl C) { return SUNDIALS_CONTROL_H; }

int SUNControlEstimateStep_PID(SUNControl C, realtype h,
                               realtype dsm, realtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype k1 = -SC_PID_K1(C) / SC_PID_P(C);
  const realtype k2 =  SC_PID_K2(C) / SC_PID_P(C);
  const realtype k3 = -SC_PID_K3(C) / SC_PID_P(C);
  const realtype ecur = SC_PID_BIAS(C) * dsm;
  const realtype e1 = SUNMAX(ecur, TINY);
  const realtype e2 = SUNMAX(SC_PID_EP(C), TINY);
  const realtype e3 = SUNMAX(SC_PID_EPP(C), TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = h * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2) * SUNRpowerR(e3,k3);
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_PID(SUNControl C)
{
  SC_PID_EP(C)  = RCONST(1.0);
  SC_PID_EPP(C) = RCONST(1.0);
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_PID(SUNControl C)
{
  SC_PID_K1(C)     = DEFAULT_K1;
  SC_PID_K2(C)     = DEFAULT_K2;
  SC_PID_K3(C)     = DEFAULT_K3;
  SC_PID_BIAS(C)   = DEFAULT_BIAS;
  SC_PID_PQ(C)     = DEFAULT_PQ;
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_PID(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "PID-controller SUNControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SC_PID_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", SC_PID_K2(C));
  fprintf(fptr, "  k3 = %32Lg\n", SC_PID_K3(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SC_PID_BIAS(C));
  fprintf(fptr, "  previous errors = %32Lg  %32Lg\n", SC_PID_EP(C), SC_PID_EPP(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SC_PID_K1(C));
  fprintf(fptr, "  k2 = %16g\n", SC_PID_K2(C));
  fprintf(fptr, "  k3 = %16g\n", SC_PID_K3(C));
  fprintf(fptr, "  bias factor = %16g\n", SC_PID_BIAS(C));
  fprintf(fptr, "  previous errors = %16g  %16g\n", SC_PID_EP(C), SC_PID_EPP(C));
#endif
  if (SC_PID_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SC_PID_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SC_PID_P(C));
  }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_PID(SUNControl C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SC_PID_PQ(C)) { SC_PID_P(C) = q; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetEmbeddingOrder_PID(SUNControl C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SC_PID_PQ(C)) { SC_PID_P(C) = p; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_PID(SUNControl C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SC_PID_BIAS(C) = DEFAULT_BIAS;
  } else {
    SC_PID_BIAS(C) = bias;
  }

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
  *leniw = 2;
  return SUNCONTROL_SUCCESS;
}
