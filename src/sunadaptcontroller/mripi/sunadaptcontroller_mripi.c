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
 * SUNAdaptController_MRIPI module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_mripi.h>
#include <sundials/sundials_core.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"


/* ---------------
 * Macro accessors
 * --------------- */

#define MRIPI_CONTENT(C) ( (SUNAdaptControllerContent_MRIPI)(C->content) )
#define MRIPI_K11(C)     ( MRIPI_CONTENT(C)->k11 )
#define MRIPI_K12(C)     ( MRIPI_CONTENT(C)->k12 )
#define MRIPI_K21(C)     ( MRIPI_CONTENT(C)->k21 )
#define MRIPI_K22(C)     ( MRIPI_CONTENT(C)->k22 )
#define MRIPI_BIAS(C)    ( MRIPI_CONTENT(C)->bias )
#define MRIPI_ESP(C)     ( MRIPI_CONTENT(C)->esp )
#define MRIPI_EFP(C)     ( MRIPI_CONTENT(C)->efp )
#define MRIPI_PFAST(C)   ( MRIPI_CONTENT(C)->p )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K11  SUN_RCONST(0.18)
#define DEFAULT_K12  SUN_RCONST(0.86)
#define DEFAULT_K21  SUN_RCONST(0.34)
#define DEFAULT_K22  SUN_RCONST(0.80)
#define DEFAULT_BIAS SUN_RCONST(1.5)
#define TINY         SUN_RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRIPI controller
 */

SUNAdaptController SUNAdaptController_MRIPI(SUNContext sunctx, int p)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  SUNAdaptControllerContent_MRIPI content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  C->ops->gettype          = SUNAdaptController_GetType_MRIPI;
  C->ops->estimatemristeps = SUNAdaptController_EstimateMRISteps_MRIPI;
  C->ops->reset            = SUNAdaptController_Reset_MRIPI;
  C->ops->setdefaults      = SUNAdaptController_SetDefaults_MRIPI;
  C->ops->write            = SUNAdaptController_Write_MRIPI;
  C->ops->seterrorbias     = SUNAdaptController_SetErrorBias_MRIPI;
  C->ops->updatemrih       = SUNAdaptController_UpdateMRIH_MRIPI;
  C->ops->space            = SUNAdaptController_Space_MRIPI;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_MRIPI)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  C->content = content;

  /* Set fast method order */
  content->p = p;

  /* Fill content with default/reset values */
  SUNCheckCallNull(SUNAdaptController_SetDefaults_MRIPI(C));
  SUNCheckCallNull(SUNAdaptController_Reset_MRIPI(C));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set MRIPI parameters
 */

SUNErrCode SUNAdaptController_SetParams_MRIPI(SUNAdaptController C,
                                              sunrealtype k11, sunrealtype k12,
                                              sunrealtype k21, sunrealtype k22)
{
  SUNFunctionBegin(C->sunctx);
  MRIPI_K11(C) = k11;
  MRIPI_K12(C) = k12;
  MRIPI_K21(C) = k21;
  MRIPI_K22(C) = k22;
  return SUN_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_MRIPI(SUNAdaptController C)
{
  return SUNDIALS_CONTROL_MRI_H;
}

SUNErrCode SUNAdaptController_EstimateMRISteps_MRIPI(SUNAdaptController C,
                                                     sunrealtype H, sunrealtype h,
                                                     int P, sunrealtype DSM,
                                                     sunrealtype dsm,
                                                     sunrealtype* Hnew,
                                                     sunrealtype* hnew)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Hnew, SUN_ERR_ARG_CORRUPT);
  SUNAssert(hnew, SUN_ERR_ARG_CORRUPT);

  /* set usable time-step adaptivity parameters */
  const int p = MRIPI_PFAST(C);
  const sunrealtype k11 = MRIPI_K11(C);
  const sunrealtype k12 = MRIPI_K12(C);
  const sunrealtype k21 = MRIPI_K21(C);
  const sunrealtype k22 = MRIPI_K22(C);
  const sunrealtype a1 = (k11 + k12) / (2 * P);
  const sunrealtype a2 = -k11 / (2 * P);
  const sunrealtype b11 = (p + 1) * (k11 + k12) / (2 * P * p);
  const sunrealtype b12 = -(p + 1) * k11 / (2 * P * p);
  const sunrealtype b21 = -(k21 + k22) / (2 * p);
  const sunrealtype b22 = k21 / (2 * p);
  const sunrealtype es1 = SUNMAX(MRIPI_BIAS(C) * DSM, TINY);
  const sunrealtype es2 = MRIPI_ESP(C);
  const sunrealtype ef1 = SUNMAX(MRIPI_BIAS(C) * dsm, TINY);
  const sunrealtype ef2 = MRIPI_EFP(C);
  const sunrealtype M = SUNRceil(H/h);

  /* compute estimated optimal time step size */
  *Hnew = H * SUNRpowerR(es1,a1) * SUNRpowerR(es2,a2);
  const sunrealtype Mnew = M * SUNRpowerR(es1,b11) * SUNRpowerR(es2,b12) *
                           SUNRpowerR(ef1,b21) * SUNRpowerR(ef2,b22);
  *hnew = (*Hnew) / Mnew;

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Reset_MRIPI(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  MRIPI_ESP(C) = SUN_RCONST(1.0);
  MRIPI_EFP(C) = SUN_RCONST(1.0);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetDefaults_MRIPI(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  MRIPI_K11(C)  = DEFAULT_K11;
  MRIPI_K12(C)  = DEFAULT_K12;
  MRIPI_K21(C)  = DEFAULT_K21;
  MRIPI_K22(C)  = DEFAULT_K22;
  MRIPI_BIAS(C) = DEFAULT_BIAS;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Write_MRIPI(SUNAdaptController C, FILE *fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "Multirate PI SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k11 = %32Lg\n", MRIPI_K11(C));
  fprintf(fptr, "  k12 = %32Lg\n", MRIPI_K12(C));
  fprintf(fptr, "  k21 = %32Lg\n", MRIPI_K21(C));
  fprintf(fptr, "  k22 = %32Lg\n", MRIPI_K22(C));
  fprintf(fptr, "  bias factor = %32Lg\n", MRIPI_BIAS(C));
  fprintf(fptr, "  previous slow error = %32Lg\n", MRIPI_ESP(C));
  fprintf(fptr, "  previous fast error = %32Lg\n", MRIPI_EFP(C));
#else
  fprintf(fptr, "  k11 = %16g\n", MRIPI_K11(C));
  fprintf(fptr, "  k12 = %16g\n", MRIPI_K12(C));
  fprintf(fptr, "  k21 = %16g\n", MRIPI_K21(C));
  fprintf(fptr, "  k22 = %16g\n", MRIPI_K22(C));
  fprintf(fptr, "  bias factor = %16g\n", MRIPI_BIAS(C));
  fprintf(fptr, "  previous slow error = %16g\n", MRIPI_ESP(C));
  fprintf(fptr, "  previous fast error = %16g\n", MRIPI_EFP(C));
#endif
  fprintf(fptr, "  p = %i (fast method order)\n", MRIPI_PFAST(C));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetErrorBias_MRIPI(SUNAdaptController C, sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);

  /* set allowed value, otherwise set default */
  if (bias <= SUN_RCONST(0.0)) { MRIPI_BIAS(C) = DEFAULT_BIAS; }
  else { MRIPI_BIAS(C) = bias; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_UpdateMRIH_MRIPI(SUNAdaptController C, sunrealtype H, sunrealtype h,
                                sunrealtype DSM, sunrealtype dsm)
{
  SUNFunctionBegin(C->sunctx);
  MRIPI_ESP(C) = SUNMAX(MRIPI_BIAS(C) * DSM, TINY);
  MRIPI_EFP(C) = SUNMAX(MRIPI_BIAS(C) * dsm, TINY);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Space_MRIPI(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = 7;
  *leniw = 1;
  return SUN_SUCCESS;
}
