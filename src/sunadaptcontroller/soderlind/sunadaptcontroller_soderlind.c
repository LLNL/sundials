/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the
 * SUNAdaptController_Soderlind module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>

#include "sundials_macros.h"

/* ---------------
 * Macro accessors
 * --------------- */

#define SODERLIND_CONTENT(C)    ((SUNAdaptControllerContent_Soderlind)(C->content))
#define SODERLIND_K1(C)         (SODERLIND_CONTENT(C)->k1)
#define SODERLIND_K2(C)         (SODERLIND_CONTENT(C)->k2)
#define SODERLIND_K3(C)         (SODERLIND_CONTENT(C)->k3)
#define SODERLIND_K4(C)         (SODERLIND_CONTENT(C)->k4)
#define SODERLIND_K5(C)         (SODERLIND_CONTENT(C)->k5)
#define SODERLIND_BIAS(C)       (SODERLIND_CONTENT(C)->bias)
#define SODERLIND_EP(C)         (SODERLIND_CONTENT(C)->ep)
#define SODERLIND_EPP(C)        (SODERLIND_CONTENT(C)->epp)
#define SODERLIND_HP(C)         (SODERLIND_CONTENT(C)->hp)
#define SODERLIND_HPP(C)        (SODERLIND_CONTENT(C)->hpp)
#define SODERLIND_FIRSTSTEPS(C) (SODERLIND_CONTENT(C)->firststeps)

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1        SUN_RCONST(1.25) /* H_{0}321 parameters */
#define DEFAULT_K2        SUN_RCONST(0.5)
#define DEFAULT_K3        SUN_RCONST(-0.75)
#define DEFAULT_K4        SUN_RCONST(0.25)
#define DEFAULT_K5        SUN_RCONST(0.75)
#define DEFAULT_PID_K1    SUN_RCONST(0.58) /* PID parameters */
#define DEFAULT_PID_K2    -SUN_RCONST(0.21)
#define DEFAULT_PID_K3    SUN_RCONST(0.1)
#define DEFAULT_PI_K1     SUN_RCONST(0.8) /* PI parameters */
#define DEFAULT_PI_K2     -SUN_RCONST(0.31)
#define DEFAULT_I_K1      SUN_RCONST(1.0)   /* I parameters */
#define DEFAULT_EXPGUS_K1 SUN_RCONST(0.367) /* Explicit Gustafsson parameters */
#define DEFAULT_EXPGUS_K2 SUN_RCONST(0.268)
#define DEFAULT_IMPGUS_K1 SUN_RCONST(0.98) /* Implicit Gustafsson parameters */
#define DEFAULT_IMPGUS_K2 SUN_RCONST(0.95)
#define DEFAULT_BIAS      SUN_RCONST(1.5)
#define TINY              SUN_RCONST(1.0e-10)

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new Soderlind controller
 */

SUNAdaptController SUNAdaptController_Soderlind(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  SUNAdaptControllerContent_Soderlind content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  C->ops->gettype      = SUNAdaptController_GetType_Soderlind;
  C->ops->estimatestep = SUNAdaptController_EstimateStep_Soderlind;
  C->ops->reset        = SUNAdaptController_Reset_Soderlind;
  C->ops->setdefaults  = SUNAdaptController_SetDefaults_Soderlind;
  C->ops->write        = SUNAdaptController_Write_Soderlind;
  C->ops->seterrorbias = SUNAdaptController_SetErrorBias_Soderlind;
  C->ops->updateh      = SUNAdaptController_UpdateH_Soderlind;
  C->ops->space        = SUNAdaptController_Space_Soderlind;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_Soderlind)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  C->content = content;

  /* Fill content with default/reset values */
  SUNCheckCallNull(SUNAdaptController_SetDefaults_Soderlind(C));
  SUNCheckCallNull(SUNAdaptController_Reset_Soderlind(C));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set Soderlind parameters
 */

SUNErrCode SUNAdaptController_SetParams_Soderlind(SUNAdaptController C,
                                                  sunrealtype k1,
                                                  sunrealtype k2, sunrealtype k3,
                                                  sunrealtype k4, sunrealtype k5)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_K1(C) = k1;
  SODERLIND_K2(C) = k2;
  SODERLIND_K3(C) = k3;
  SODERLIND_K4(C) = k4;
  SODERLIND_K5(C) = k5;
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to create a PID controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_PID(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(SUNAdaptController_SetParams_PID(C, DEFAULT_PID_K1,
                                                    DEFAULT_PID_K2,
                                                    DEFAULT_PID_K3));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PID parameters
 */

SUNErrCode SUNAdaptController_SetParams_PID(SUNAdaptController C, sunrealtype k1,
                                            sunrealtype k2, sunrealtype k3)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_K1(C) = k1;
  SODERLIND_K2(C) = k2;
  SODERLIND_K3(C) = k3;
  SODERLIND_K4(C) = SUN_RCONST(0.0);
  SODERLIND_K5(C) = SUN_RCONST(0.0);
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to create a PI controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_PI(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(
    SUNAdaptController_SetParams_PI(C, DEFAULT_PI_K1, DEFAULT_PI_K2));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PI parameters
 */

SUNErrCode SUNAdaptController_SetParams_PI(SUNAdaptController C, sunrealtype k1,
                                           sunrealtype k2)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_K1(C) = k1;
  SODERLIND_K2(C) = k2;
  SODERLIND_K3(C) = SUN_RCONST(0.0);
  SODERLIND_K4(C) = SUN_RCONST(0.0);
  SODERLIND_K5(C) = SUN_RCONST(0.0);
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to create an I controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_I(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(SUNAdaptController_SetParams_I(C, DEFAULT_I_K1));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PI parameters
 */

SUNErrCode SUNAdaptController_SetParams_I(SUNAdaptController C, sunrealtype k1)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_K1(C) = k1;
  SODERLIND_K2(C) = SUN_RCONST(0.0);
  SODERLIND_K3(C) = SUN_RCONST(0.0);
  SODERLIND_K4(C) = SUN_RCONST(0.0);
  SODERLIND_K5(C) = SUN_RCONST(0.0);
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to create an explicit Gustafsson controller
 */

SUNAdaptController SUNAdaptController_ExpGus(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(SUNAdaptController_SetParams_ExpGus(C, DEFAULT_EXPGUS_K1,
                                                       DEFAULT_EXPGUS_K2));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set explicit Gustafsson parameters
 */

SUNErrCode SUNAdaptController_SetParams_ExpGus(SUNAdaptController C,
                                               sunrealtype k1, sunrealtype k2)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_K1(C) = k1 + k2;
  SODERLIND_K2(C) = -k2;
  SODERLIND_K3(C) = SUN_RCONST(0.0);
  SODERLIND_K4(C) = SUN_RCONST(0.0);
  SODERLIND_K5(C) = SUN_RCONST(0.0);
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to create an implicit Gustafsson controller
 */

SUNAdaptController SUNAdaptController_ImpGus(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(SUNAdaptController_SetParams_ImpGus(C, DEFAULT_IMPGUS_K1,
                                                       DEFAULT_IMPGUS_K2));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set explicit Gustafsson parameters
 */

SUNErrCode SUNAdaptController_SetParams_ImpGus(SUNAdaptController C,
                                               sunrealtype k1, sunrealtype k2)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_K1(C) = k1 + k2;
  SODERLIND_K2(C) = -k2;
  SODERLIND_K3(C) = SUN_RCONST(0.0);
  SODERLIND_K4(C) = SUN_RCONST(1.0);
  SODERLIND_K5(C) = SUN_RCONST(0.0);
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_Soderlind(
  SUNDIALS_MAYBE_UNUSED SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_H;
}

SUNErrCode SUNAdaptController_EstimateStep_Soderlind(SUNAdaptController C,
                                                     sunrealtype h, int p,
                                                     sunrealtype dsm,
                                                     sunrealtype* hnew)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(hnew, SUN_ERR_ARG_CORRUPT);

  /* order parameter to use */
  const int ord = p + 1;

  /* set usable time-step adaptivity parameters */
  const sunrealtype k1    = -SODERLIND_K1(C) / ord;
  const sunrealtype k2    = -SODERLIND_K2(C) / ord;
  const sunrealtype k3    = -SODERLIND_K3(C) / ord;
  const sunrealtype k4    = SODERLIND_K4(C);
  const sunrealtype k5    = SODERLIND_K5(C);
  const sunrealtype e1    = SUNMAX(SODERLIND_BIAS(C) * dsm, TINY);
  const sunrealtype e2    = SUNMAX(SODERLIND_EP(C), TINY);
  const sunrealtype e3    = SUNMAX(SODERLIND_EPP(C), TINY);
  const sunrealtype hrat  = h / SODERLIND_HP(C);
  const sunrealtype hrat2 = SODERLIND_HP(C) / SODERLIND_HPP(C);

  /* compute estimated optimal time step size */
  if (SODERLIND_FIRSTSTEPS(C) < 1) { *hnew = h * SUNRpowerR(e1, k1); }
  else if (SODERLIND_FIRSTSTEPS(C) < 2)
  {
    *hnew = h * SUNRpowerR(e1, k1) * SUNRpowerR(e2, k2) * SUNRpowerR(hrat, k4);
  }
  else
  {
    *hnew = h * SUNRpowerR(e1, k1) * SUNRpowerR(e2, k2) * SUNRpowerR(e3, k3) *
            SUNRpowerR(hrat, k4) * SUNRpowerR(hrat2, k5);
  }

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Reset_Soderlind(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_EP(C)         = SUN_RCONST(1.0);
  SODERLIND_EPP(C)        = SUN_RCONST(1.0);
  SODERLIND_HP(C)         = SUN_RCONST(1.0);
  SODERLIND_HPP(C)        = SUN_RCONST(1.0);
  SODERLIND_FIRSTSTEPS(C) = 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetDefaults_Soderlind(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_K1(C)   = DEFAULT_K1;
  SODERLIND_K2(C)   = DEFAULT_K2;
  SODERLIND_K3(C)   = DEFAULT_K3;
  SODERLIND_K4(C)   = DEFAULT_K4;
  SODERLIND_K5(C)   = DEFAULT_K5;
  SODERLIND_BIAS(C) = DEFAULT_BIAS;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Write_Soderlind(SUNAdaptController C, FILE* fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "Soderlind SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SODERLIND_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", SODERLIND_K2(C));
  fprintf(fptr, "  k3 = %32Lg\n", SODERLIND_K3(C));
  fprintf(fptr, "  k4 = %32Lg\n", SODERLIND_K4(C));
  fprintf(fptr, "  k5 = %32Lg\n", SODERLIND_K5(C));
  fprintf(fptr, "  bias factor = %22Lg\n", SODERLIND_BIAS(C));
  fprintf(fptr, "  previous error = %22Lg\n", SODERLIND_EP(C));
  fprintf(fptr, "  previous-previous error = %22Lg\n", SODERLIND_EPP(C));
  fprintf(fptr, "  previous step = %22Lg\n", SODERLIND_HP(C));
  fprintf(fptr, "  previous-previous step = %22Lg\n", SODERLIND_HPP(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SODERLIND_K1(C));
  fprintf(fptr, "  k2 = %16g\n", SODERLIND_K2(C));
  fprintf(fptr, "  k3 = %16g\n", SODERLIND_K3(C));
  fprintf(fptr, "  k4 = %16g\n", SODERLIND_K4(C));
  fprintf(fptr, "  k5 = %16g\n", SODERLIND_K5(C));
  fprintf(fptr, "  bias factor = %16g\n", SODERLIND_BIAS(C));
  fprintf(fptr, "  previous error = %16g\n", SODERLIND_EP(C));
  fprintf(fptr, "  previous-previous error = %16g\n", SODERLIND_EPP(C));
  fprintf(fptr, "  previous step = %16g\n", SODERLIND_HP(C));
  fprintf(fptr, "  previous-previous step = %16g\n", SODERLIND_HPP(C));
#endif
  fprintf(fptr, "  firststeps = %i\n", SODERLIND_FIRSTSTEPS(C));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetErrorBias_Soderlind(SUNAdaptController C,
                                                     sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);

  /* set allowed value, otherwise set default */
  if (bias <= SUN_RCONST(0.0)) { SODERLIND_BIAS(C) = DEFAULT_BIAS; }
  else { SODERLIND_BIAS(C) = bias; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_UpdateH_Soderlind(SUNAdaptController C,
                                                sunrealtype h, sunrealtype dsm)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_EPP(C) = SODERLIND_EP(C);
  SODERLIND_EP(C)  = SODERLIND_BIAS(C) * dsm;
  SODERLIND_HPP(C) = SODERLIND_HP(C);
  SODERLIND_HP(C)  = h;
  if (SODERLIND_FIRSTSTEPS(C) < 2) { SODERLIND_FIRSTSTEPS(C) += 1; }
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Space_Soderlind(SUNAdaptController C,
                                              long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = 10;
  *leniw = 1;
  return SUN_SUCCESS;
}
