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
 * SUNAdaptController_MRILL module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunadaptcontroller/sunadaptcontroller_mrill.h>
#include <sundials/sundials_core.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"


/* ---------------
 * Macro accessors
 * --------------- */

#define MRILL_CONTENT(C)   ( (SUNAdaptControllerContent_MRILL)(C->content) )
#define MRILL_K11(C)       ( MRILL_CONTENT(C)->k11 )
#define MRILL_K12(C)       ( MRILL_CONTENT(C)->k12 )
#define MRILL_K21(C)       ( MRILL_CONTENT(C)->k21 )
#define MRILL_K22(C)       ( MRILL_CONTENT(C)->k22 )
#define MRILL_BIAS(C)      ( MRILL_CONTENT(C)->bias )
#define MRILL_ESP(C)       ( MRILL_CONTENT(C)->esp )
#define MRILL_EFP(C)       ( MRILL_CONTENT(C)->efp )
#define MRILL_HSP(C)       ( MRILL_CONTENT(C)->hsp )
#define MRILL_HFP(C)       ( MRILL_CONTENT(C)->hfp )
#define MRILL_PFAST(C)     ( MRILL_CONTENT(C)->p )
#define MRILL_FIRSTSTEP(C) ( MRILL_CONTENT(C)->firststep )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K11  SUN_RCONST(0.82)
#define DEFAULT_K12  SUN_RCONST(0.54)
#define DEFAULT_K21  SUN_RCONST(0.94)
#define DEFAULT_K22  SUN_RCONST(0.9)
#define DEFAULT_BIAS SUN_RCONST(1.5)
#define TINY         SUN_RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRILL controller
 */

SUNAdaptController SUNAdaptController_MRILL(SUNContext sunctx, int p)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  SUNAdaptControllerContent_MRILL content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  C->ops->gettype          = SUNAdaptController_GetType_MRILL;
  C->ops->estimatemristeps = SUNAdaptController_EstimateMRISteps_MRILL;
  C->ops->reset            = SUNAdaptController_Reset_MRILL;
  C->ops->setdefaults      = SUNAdaptController_SetDefaults_MRILL;
  C->ops->write            = SUNAdaptController_Write_MRILL;
  C->ops->seterrorbias     = SUNAdaptController_SetErrorBias_MRILL;
  C->ops->updatemrih       = SUNAdaptController_UpdateMRIH_MRILL;
  C->ops->space            = SUNAdaptController_Space_MRILL;

  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_MRILL)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  C->content = content;

  /* Set fast method order */
  content->p = p;

  /* Fill content with default/reset values */
  SUNCheckCallNull(SUNAdaptController_SetDefaults_MRILL(C));
  SUNCheckCallNull(SUNAdaptController_Reset_MRILL(C));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set MRILL parameters
 */

SUNErrCode SUNAdaptController_SetParams_MRILL(SUNAdaptController C,
                                              sunrealtype k11, sunrealtype k12,
                                              sunrealtype k21, sunrealtype k22)
{
  SUNFunctionBegin(C->sunctx);
  MRILL_K11(C) = k11;
  MRILL_K12(C) = k12;
  MRILL_K21(C) = k21;
  MRILL_K22(C) = k22;
  return SUN_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_MRILL(SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_MRI_H;
}

SUNErrCode SUNAdaptController_EstimateMRISteps_MRILL(SUNAdaptController C,
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
  const int p = MRILL_PFAST(C);
  const sunrealtype k11 = MRILL_K11(C);
  const sunrealtype k12 = MRILL_K12(C);
  const sunrealtype k21 = MRILL_K21(C);
  const sunrealtype k22 = MRILL_K22(C);
  const sunrealtype a1  = (k11 + k12) / (2 * P);
  const sunrealtype a2  = -k11 / (2 * P);
  const sunrealtype b11 = (p + 1) * (k11 + k12) / (2 * P * p);
  const sunrealtype b12 = -(p + 1) * k11 / (2 * P * p);
  const sunrealtype b21 = -(k21 + k22) / (2 * p);
  const sunrealtype b22 = k21 / (2 * p);
  const sunrealtype es1 = SUNMAX(MRILL_BIAS(C) * DSM, TINY);
  const sunrealtype es2 = MRILL_ESP(C);
  const sunrealtype ef1 = SUNMAX(MRILL_BIAS(C) * dsm, TINY);
  const sunrealtype ef2 = MRILL_EFP(C);
  const sunrealtype M   = SUNRceil(H/h);

  /* handle first vs successive steps */
  sunrealtype Hfac, Mfac;
  if (MRILL_FIRSTSTEP(C))
  {
    Hfac = SUN_RCONST(1.0);
    Mfac = SUN_RCONST(1.0);
  }
  else
  {
    const sunrealtype Mp = SUNRceil(MRILL_HSP(C)/MRILL_HFP(C));
    Hfac = H / MRILL_HSP(C);
    Mfac = M / Mp;
  }

  /* compute estimated optimal time step size */
  *Hnew = H * Hfac * SUNRpowerR(es1,a1) * SUNRpowerR(es2,a2);
  const sunrealtype Mnew = M * Mfac * SUNRpowerR(es1,b11) *
                           SUNRpowerR(es2,b12) * SUNRpowerR(ef1,b21) *
                           SUNRpowerR(ef2,b22);
  *hnew = (*Hnew) / Mnew;

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Reset_MRILL(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  MRILL_ESP(C) = SUN_RCONST(1.0);
  MRILL_EFP(C) = SUN_RCONST(1.0);
  MRILL_HSP(C) = SUN_RCONST(1.0);
  MRILL_HFP(C) = SUN_RCONST(1.0);
  MRILL_FIRSTSTEP(C) = SUNTRUE;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetDefaults_MRILL(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  MRILL_K11(C)  = DEFAULT_K11;
  MRILL_K12(C)  = DEFAULT_K12;
  MRILL_K21(C)  = DEFAULT_K21;
  MRILL_K22(C)  = DEFAULT_K22;
  MRILL_BIAS(C) = DEFAULT_BIAS;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Write_MRILL(SUNAdaptController C, FILE *fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "Multirate LL SUNAdaptController module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k11 = %32Lg\n", MRILL_K11(C));
  fprintf(fptr, "  k12 = %32Lg\n", MRILL_K12(C));
  fprintf(fptr, "  k21 = %32Lg\n", MRILL_K21(C));
  fprintf(fptr, "  k22 = %32Lg\n", MRILL_K22(C));
  fprintf(fptr, "  bias factor = %32Lg\n", MRILL_BIAS(C));
  fprintf(fptr, "  previous slow error = %32Lg\n", MRILL_ESP(C));
  fprintf(fptr, "  previous fast error = %32Lg\n", MRILL_EFP(C));
  fprintf(fptr, "  previous slow step  = %32Lg\n", MRILL_HSP(C));
  fprintf(fptr, "  previous fast step  = %32Lg\n", MRILL_HFP(C));
#else
  fprintf(fptr, "  k11 = %16g\n", MRILL_K11(C));
  fprintf(fptr, "  k12 = %16g\n", MRILL_K12(C));
  fprintf(fptr, "  k21 = %16g\n", MRILL_K21(C));
  fprintf(fptr, "  k22 = %16g\n", MRILL_K22(C));
  fprintf(fptr, "  bias factor = %16g\n", MRILL_BIAS(C));
  fprintf(fptr, "  previous slow error = %16g\n", MRILL_ESP(C));
  fprintf(fptr, "  previous fast error = %16g\n", MRILL_EFP(C));
  fprintf(fptr, "  previous slow step  = %16g\n", MRILL_HSP(C));
  fprintf(fptr, "  previous fast step  = %16g\n", MRILL_HFP(C));
#endif
  fprintf(fptr, "  p = %i (fast method order)\n", MRILL_PFAST(C));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetErrorBias_MRILL(SUNAdaptController C, sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);

  /* set allowed value, otherwise set default */
  if (bias <= SUN_RCONST(0.0)) { MRILL_BIAS(C) = DEFAULT_BIAS; }
  else { MRILL_BIAS(C) = bias; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_UpdateMRIH_MRILL(SUNAdaptController C, sunrealtype H, sunrealtype h,
                                sunrealtype DSM, sunrealtype dsm)
{
  SUNFunctionBegin(C->sunctx);
  MRILL_ESP(C) = SUNMAX(MRILL_BIAS(C) * DSM, TINY);
  MRILL_EFP(C) = SUNMAX(MRILL_BIAS(C) * dsm, TINY);
  MRILL_HSP(C) = H;
  MRILL_HFP(C) = h;
  MRILL_FIRSTSTEP(C) = SUNFALSE;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Space_MRILL(SUNAdaptController C, long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = 9;
  *leniw = 2;
  return SUN_SUCCESS;
}
