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
 * SUNAdaptController_MRIPID module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_mripid.h>
#include <sundials/sundials_core.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"

/* ---------------
 * Macro accessors
 * --------------- */

#define MRIPID_CONTENT(C) ((SUNAdaptControllerContent_MRIPID)(C->content))
#define MRIPID_K11(C)     (MRIPID_CONTENT(C)->k11)
#define MRIPID_K12(C)     (MRIPID_CONTENT(C)->k12)
#define MRIPID_K13(C)     (MRIPID_CONTENT(C)->k13)
#define MRIPID_K21(C)     (MRIPID_CONTENT(C)->k21)
#define MRIPID_K22(C)     (MRIPID_CONTENT(C)->k22)
#define MRIPID_K23(C)     (MRIPID_CONTENT(C)->k23)
#define MRIPID_BIAS(C)    (MRIPID_CONTENT(C)->bias)
#define MRIPID_ESP(C)     (MRIPID_CONTENT(C)->esp)
#define MRIPID_EFP(C)     (MRIPID_CONTENT(C)->efp)
#define MRIPID_ESPP(C)    (MRIPID_CONTENT(C)->espp)
#define MRIPID_EFPP(C)    (MRIPID_CONTENT(C)->efpp)
#define MRIPID_PFAST(C)   (MRIPID_CONTENT(C)->p)

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K11  SUN_RCONST(0.34)
#define DEFAULT_K12  SUN_RCONST(0.1)
#define DEFAULT_K13  SUN_RCONST(0.78)
#define DEFAULT_K21  SUN_RCONST(0.46)
#define DEFAULT_K22  SUN_RCONST(0.42)
#define DEFAULT_K23  SUN_RCONST(0.74)
#define DEFAULT_BIAS SUN_RCONST(1.5)
#define ONE          SUN_RCONST(1.0)
#define TINY         (SUN_RCONST(10.0)*SUN_UNIT_ROUNDOFF)

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRIPID controller
 */

SUNAdaptController SUNAdaptController_MRIPID(SUNContext sunctx, int p)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  SUNAdaptControllerContent_MRIPID content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  C->ops->gettype          = SUNAdaptController_GetType_MRIPID;
  C->ops->estimatemristeps = SUNAdaptController_EstimateMRISteps_MRIPID;
  C->ops->reset            = SUNAdaptController_Reset_MRIPID;
  C->ops->setdefaults      = SUNAdaptController_SetDefaults_MRIPID;
  C->ops->write            = SUNAdaptController_Write_MRIPID;
  C->ops->seterrorbias     = SUNAdaptController_SetErrorBias_MRIPID;
  C->ops->updatemrih       = SUNAdaptController_UpdateMRIH_MRIPID;
  C->ops->space            = SUNAdaptController_Space_MRIPID;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_MRIPID)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  C->content = content;

  /* Set fast method order */
  content->p = p;

  /* Fill content with default/reset values */
  SUNCheckCallNull(SUNAdaptController_SetDefaults_MRIPID(C));
  SUNCheckCallNull(SUNAdaptController_Reset_MRIPID(C));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set MRIPID parameters
 */

SUNErrCode SUNAdaptController_SetParams_MRIPID(SUNAdaptController C,
                                               sunrealtype k11, sunrealtype k12,
                                               sunrealtype k13, sunrealtype k21,
                                               sunrealtype k22, sunrealtype k23)
{
  SUNFunctionBegin(C->sunctx);
  MRIPID_K11(C) = k11;
  MRIPID_K12(C) = k12;
  MRIPID_K13(C) = k13;
  MRIPID_K21(C) = k21;
  MRIPID_K22(C) = k22;
  MRIPID_K23(C) = k23;
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_MRIPID(SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_MRI_H;
}

SUNErrCode SUNAdaptController_EstimateMRISteps_MRIPID(
  SUNAdaptController C, sunrealtype H, sunrealtype h, int P, sunrealtype DSM,
  sunrealtype dsm, sunrealtype* Hnew, sunrealtype* hnew)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Hnew, SUN_ERR_ARG_CORRUPT);
  SUNAssert(hnew, SUN_ERR_ARG_CORRUPT);

  /* set usable time-step adaptivity parameters */
  const int p           = MRIPID_PFAST(C);
  const sunrealtype k11 = MRIPID_K11(C);
  const sunrealtype k12 = MRIPID_K12(C);
  const sunrealtype k13 = MRIPID_K13(C);
  const sunrealtype k21 = MRIPID_K21(C);
  const sunrealtype k22 = MRIPID_K22(C);
  const sunrealtype k23 = MRIPID_K23(C);
  const sunrealtype a1  = (k11 + k12 + k13) / (3 * P);
  const sunrealtype a2  = -(k11 + k12) / (3 * P);
  const sunrealtype a3  = k11 / (3 * P);
  const sunrealtype b11 = (p + 1) * (k11 + k12 + k13) / (3 * P * p);
  const sunrealtype b12 = -(p + 1) * (k11 + k12) / (3 * P * p);
  const sunrealtype b13 = (p + 1) * k11 / (3 * P * p);
  const sunrealtype b21 = -(k21 + k22 + k23) / (3 * p);
  const sunrealtype b22 = (k21 + k22) / (3 * p);
  const sunrealtype b23 = -k21 / (3 * p);
  const sunrealtype es1 = ONE / SUNMAX(MRIPID_BIAS(C) * DSM, TINY);
  const sunrealtype es2 = ONE / MRIPID_ESP(C);
  const sunrealtype es3 = ONE / MRIPID_ESPP(C);
  const sunrealtype ef1 = ONE / SUNMAX(MRIPID_BIAS(C) * dsm, TINY);
  const sunrealtype ef2 = ONE / MRIPID_EFP(C);
  const sunrealtype ef3 = ONE / MRIPID_EFPP(C);
  const sunrealtype M   = SUNRceil(H / h);

  /* compute estimated optimal time step size */
  *Hnew = H * SUNRpowerR(es1, a1) * SUNRpowerR(es2, a2) * SUNRpowerR(es3, a3);
  const sunrealtype Mnew = M * SUNRpowerR(es1, b11) * SUNRpowerR(es2, b12) *
                           SUNRpowerR(es3, b13) * SUNRpowerR(ef1, b21) *
                           SUNRpowerR(ef2, b22) * SUNRpowerR(ef3, b23);
  *hnew = (*Hnew) / Mnew;

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Reset_MRIPID(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  MRIPID_ESP(C)  = SUN_RCONST(1.0);
  MRIPID_EFP(C)  = SUN_RCONST(1.0);
  MRIPID_ESPP(C) = SUN_RCONST(1.0);
  MRIPID_EFPP(C) = SUN_RCONST(1.0);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetDefaults_MRIPID(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  MRIPID_K11(C)  = DEFAULT_K11;
  MRIPID_K12(C)  = DEFAULT_K12;
  MRIPID_K13(C)  = DEFAULT_K13;
  MRIPID_K21(C)  = DEFAULT_K21;
  MRIPID_K22(C)  = DEFAULT_K22;
  MRIPID_K23(C)  = DEFAULT_K23;
  MRIPID_BIAS(C) = DEFAULT_BIAS;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Write_MRIPID(SUNAdaptController C, FILE* fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "Multirate PID SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k11 = %32Lg\n", MRIPID_K11(C));
  fprintf(fptr, "  k12 = %32Lg\n", MRIPID_K12(C));
  fprintf(fptr, "  k13 = %32Lg\n", MRIPID_K13(C));
  fprintf(fptr, "  k21 = %32Lg\n", MRIPID_K21(C));
  fprintf(fptr, "  k22 = %32Lg\n", MRIPID_K22(C));
  fprintf(fptr, "  k23 = %32Lg\n", MRIPID_K23(C));
  fprintf(fptr, "  bias factor = %32Lg\n", MRIPID_BIAS(C));
  fprintf(fptr, "  previous slow errors = %32Lg  %32Lg\n", MRIPID_ESP(C),
          MRIPID_ESPP(C));
  fprintf(fptr, "  previous fast errors = %32Lg  %32Lg\n", MRIPID_EFP(C),
          MRIPID_EFPP(C));
#else
  fprintf(fptr, "  k11 = %16g\n", MRIPID_K11(C));
  fprintf(fptr, "  k12 = %16g\n", MRIPID_K12(C));
  fprintf(fptr, "  k13 = %16g\n", MRIPID_K13(C));
  fprintf(fptr, "  k21 = %16g\n", MRIPID_K21(C));
  fprintf(fptr, "  k22 = %16g\n", MRIPID_K22(C));
  fprintf(fptr, "  k23 = %16g\n", MRIPID_K23(C));
  fprintf(fptr, "  bias factor = %16g\n", MRIPID_BIAS(C));
  fprintf(fptr, "  previous slow errors = %16g  %16g\n", MRIPID_ESP(C),
          MRIPID_ESPP(C));
  fprintf(fptr, "  previous fast errors = %16g  %16g\n", MRIPID_EFP(C),
          MRIPID_EFPP(C));
#endif
  fprintf(fptr, "  p = %i (fast method order)\n", MRIPID_PFAST(C));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetErrorBias_MRIPID(SUNAdaptController C,
                                                  sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);

  /* set allowed value, otherwise set default */
  if (bias <= SUN_RCONST(0.0)) { MRIPID_BIAS(C) = DEFAULT_BIAS; }
  else { MRIPID_BIAS(C) = bias; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_UpdateMRIH_MRIPID(SUNAdaptController C,
                                                sunrealtype H, sunrealtype h,
                                                sunrealtype DSM, sunrealtype dsm)
{
  SUNFunctionBegin(C->sunctx);
  MRIPID_ESPP(C) = MRIPID_ESP(C);
  MRIPID_EFPP(C) = MRIPID_EFP(C);
  MRIPID_ESP(C)  = SUNMAX(MRIPID_BIAS(C) * DSM, TINY);
  MRIPID_EFP(C)  = SUNMAX(MRIPID_BIAS(C) * dsm, TINY);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Space_MRIPID(SUNAdaptController C,
                                           long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = 11;
  *leniw = 1;
  return SUN_SUCCESS;
}
