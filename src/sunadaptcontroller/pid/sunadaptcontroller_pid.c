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
 * This is the implementation file for the SUNAdaptController_PID
 * module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_pid.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SACPID_CONTENT(C)     ( (SUNAdaptControllerContent_PID)(C->content) )
#define SACPID_K1(C)          ( SACPID_CONTENT(C)->k1 )
#define SACPID_K2(C)          ( SACPID_CONTENT(C)->k2 )
#define SACPID_K3(C)          ( SACPID_CONTENT(C)->k3 )
#define SACPID_BIAS(C)        ( SACPID_CONTENT(C)->bias )
#define SACPID_EP(C)          ( SACPID_CONTENT(C)->ep )
#define SACPID_EPP(C)         ( SACPID_CONTENT(C)->epp )
#define SACPID_P(C)           ( SACPID_CONTENT(C)->p )
#define SACPID_PQ(C)          ( SACPID_CONTENT(C)->pq )

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

SUNAdaptController SUNAdaptControllerPID(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptControllerContent_PID content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype           = SUNAdaptControllerGetType_PID;
  C->ops->estimatestep      = SUNAdaptControllerEstimateStep_PID;
  C->ops->reset             = SUNAdaptControllerReset_PID;
  C->ops->setdefaults       = SUNAdaptControllerSetDefaults_PID;
  C->ops->write             = SUNAdaptControllerWrite_PID;
  C->ops->setmethodorder    = SUNAdaptControllerSetMethodOrder_PID;
  C->ops->setembeddingorder = SUNAdaptControllerSetEmbeddingOrder_PID;
  C->ops->seterrorbias      = SUNAdaptControllerSetErrorBias_PID;
  C->ops->update            = SUNAdaptControllerUpdate_PID;
  C->ops->space             = SUNAdaptControllerSpace_PID;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_PID)malloc(sizeof *content);
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
  SUNAdaptControllerSetDefaults_PID(C);
  SUNAdaptControllerReset_PID(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PID parameters
 */

int SUNAdaptControllerPID_SetParams(SUNAdaptController C, sunbooleantype pq,
                                    realtype k1, realtype k2, realtype k3)
{
  /* store legal inputs, and return with success */
  SACPID_PQ(C) = pq;
  if (k1 >= RCONST(0.0)) { SACPID_K1(C) = k1; }
  if (k2 >= RCONST(0.0)) { SACPID_K2(C) = k2; }
  if (k3 >= RCONST(0.0)) { SACPID_K3(C) = k3; }
  return SUNADAPTCONTROLLER_SUCCESS;
}


/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptControllerGetType_PID(SUNAdaptController C)
{ return SUN_ADAPTCONTROLLER_H; }

int SUNAdaptControllerEstimateStep_PID(SUNAdaptController C, realtype h,
                                       realtype dsm, realtype* hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype k1 = -SACPID_K1(C) / SACPID_P(C);
  const realtype k2 =  SACPID_K2(C) / SACPID_P(C);
  const realtype k3 = -SACPID_K3(C) / SACPID_P(C);
  const realtype ecur = SACPID_BIAS(C) * dsm;
  const realtype e1 = SUNMAX(ecur, TINY);
  const realtype e2 = SUNMAX(SACPID_EP(C), TINY);
  const realtype e3 = SUNMAX(SACPID_EPP(C), TINY);

  /* compute estimated optimal time step size and return with success */
  *hnew = h * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2) * SUNRpowerR(e3,k3);
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerReset_PID(SUNAdaptController C)
{
  SACPID_EP(C)  = RCONST(1.0);
  SACPID_EPP(C) = RCONST(1.0);
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetDefaults_PID(SUNAdaptController C)
{
  SACPID_K1(C)     = DEFAULT_K1;
  SACPID_K2(C)     = DEFAULT_K2;
  SACPID_K3(C)     = DEFAULT_K3;
  SACPID_BIAS(C)   = DEFAULT_BIAS;
  SACPID_PQ(C)     = DEFAULT_PQ;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerWrite_PID(SUNAdaptController C, FILE *fptr)
{
  fprintf(fptr, "PID SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SACPID_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", SACPID_K2(C));
  fprintf(fptr, "  k3 = %32Lg\n", SACPID_K3(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SACPID_BIAS(C));
  fprintf(fptr, "  previous errors = %32Lg  %32Lg\n", SACPID_EP(C), SACPID_EPP(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SACPID_K1(C));
  fprintf(fptr, "  k2 = %16g\n", SACPID_K2(C));
  fprintf(fptr, "  k3 = %16g\n", SACPID_K3(C));
  fprintf(fptr, "  bias factor = %16g\n", SACPID_BIAS(C));
  fprintf(fptr, "  previous errors = %16g  %16g\n", SACPID_EP(C), SACPID_EPP(C));
#endif
  if (SACPID_PQ(C))
  {
    fprintf(fptr, "  p = %i (method order)\n", SACPID_P(C));
  }
  else
  {
    fprintf(fptr, "  p = %i (embedding order)\n", SACPID_P(C));
  }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetMethodOrder_PID(SUNAdaptController C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (SACPID_PQ(C)) { SACPID_P(C) = q; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetEmbeddingOrder_PID(SUNAdaptController C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNADAPTCONTROLLER_ILL_INPUT; }

  /* store if "pq" specifies to use method order */
  if (!SACPID_PQ(C)) { SACPID_P(C) = p; }
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSetErrorBias_PID(SUNAdaptController C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SACPID_BIAS(C) = DEFAULT_BIAS;
  } else {
    SACPID_BIAS(C) = bias;
  }

  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerUpdate_PID(SUNAdaptController C, realtype h, realtype dsm)
{
  SACPID_EPP(C) = SACPID_EP(C);
  SACPID_EP(C) = SACPID_BIAS(C) * dsm;
  return SUNADAPTCONTROLLER_SUCCESS;
}

int SUNAdaptControllerSpace_PID(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 7;
  *leniw = 2;
  return SUNADAPTCONTROLLER_SUCCESS;
}
