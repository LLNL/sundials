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
 * This is the implementation file for the SUNControl_MRILL module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_mrill.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_MRILL_CONTENT(C)   ( (SUNControlContent_MRILL)(C->content) )
#define SC_MRILL_K11(C)       ( SC_MRILL_CONTENT(C)->k11 )
#define SC_MRILL_K12(C)       ( SC_MRILL_CONTENT(C)->k12 )
#define SC_MRILL_K21(C)       ( SC_MRILL_CONTENT(C)->k21 )
#define SC_MRILL_K22(C)       ( SC_MRILL_CONTENT(C)->k22 )
#define SC_MRILL_BIAS(C)      ( SC_MRILL_CONTENT(C)->bias )
#define SC_MRILL_ESP(C)       ( SC_MRILL_CONTENT(C)->esp )
#define SC_MRILL_EFP(C)       ( SC_MRILL_CONTENT(C)->efp )
#define SC_MRILL_HSP(C)       ( SC_MRILL_CONTENT(C)->hsp )
#define SC_MRILL_HFP(C)       ( SC_MRILL_CONTENT(C)->hfp )
#define SC_MRILL_PSLOW(C)     ( SC_MRILL_CONTENT(C)->P )
#define SC_MRILL_PFAST(C)     ( SC_MRILL_CONTENT(C)->p )
#define SC_MRILL_FIRSTSTEP(C) ( SC_MRILL_CONTENT(C)->firststep )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K11  RCONST(0.82)
#define DEFAULT_K12  RCONST(0.54)
#define DEFAULT_K21  RCONST(0.94)
#define DEFAULT_K22  RCONST(0.9)
#define DEFAULT_BIAS RCONST(1.5)
#define TINY         RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRILL controller
 */

SUNControl SUNControlMRILL(SUNContext sunctx, int P, int p)
{
  SUNControl C;
  SUNControlContent_MRILL content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid             = SUNControlGetID_MRILL;
  C->ops->estimatemristeps  = SUNControlEstimateMRISteps_MRILL;
  C->ops->reset             = SUNControlReset_MRILL;
  C->ops->setdefaults       = SUNControlSetDefaults_MRILL;
  C->ops->write             = SUNControlWrite_MRILL;
  C->ops->seterrorbias      = SUNControlSetErrorBias_MRILL;
  C->ops->updatemrih        = SUNControlUpdateMRIH_MRILL;
  C->ops->space             = SUNControlSpace_MRILL;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_MRILL)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNControlDestroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Initialize method orders */
  content->P = P;
  content->p = p;

  /* Fill content with default/reset values */
  SUNControlSetDefaults_MRILL(C);
  SUNControlReset_MRILL(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set MRILL parameters
 */

int SUNControlMRILL_SetParams(SUNControl C, realtype k11, realtype k12,
                              realtype k21, realtype k22)
{
  /* store legal inputs, and return with success */
  if (k11 >= RCONST(0.0)) { SC_MRILL_K11(C) = k11; }
  if (k12 >= RCONST(0.0)) { SC_MRILL_K12(C) = k12; }
  if (k21 >= RCONST(0.0)) { SC_MRILL_K21(C) = k21; }
  if (k22 >= RCONST(0.0)) { SC_MRILL_K22(C) = k22; }
  return SUNCONTROL_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_MRILL(SUNControl C) { return SUNDIALS_CONTROL_MRI_H; }

int SUNControlEstimateMRISteps_MRILL(SUNControl C, realtype H, realtype h,
                                     realtype DSM, realtype dsm,
                                     realtype* Hnew, realtype *hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype a1  = (SC_MRILL_K11(C) + SC_MRILL_K12(C))
    / (2 * SC_MRILL_PSLOW(C));
  const realtype a2  = -SC_MRILL_K11(C) / (2 * SC_MRILL_PSLOW(C));
  const realtype b11 = (SC_MRILL_PFAST(C) + 1) * (SC_MRILL_K11(C) + SC_MRILL_K12(C))
    / (2 * SC_MRILL_PSLOW(C) * SC_MRILL_PFAST(C));
  const realtype b12 = -(SC_MRILL_PFAST(C) + 1) * SC_MRILL_K11(C)
    / (2 * SC_MRILL_PSLOW(C) * SC_MRILL_PFAST(C));
  const realtype b21 = -(SC_MRILL_K21(C) + SC_MRILL_K22(C))
    / (2 * SC_MRILL_PFAST(C));
  const realtype b22 = SC_MRILL_K21(C) / (2 * SC_MRILL_PFAST(C));
  const realtype es1 = SUNMAX(SC_MRILL_BIAS(C) * DSM, TINY);
  const realtype es2 = SC_MRILL_ESP(C);
  const realtype ef1 = SUNMAX(SC_MRILL_BIAS(C) * dsm, TINY);
  const realtype ef2 = SC_MRILL_EFP(C);
  const realtype M   = SUNRceil(H/h);

  /* handle first vs successive steps */
  realtype Hfac, Mfac;
  if (SC_MRILL_FIRSTSTEP(C))
  {
    Hfac = RCONST(1.0);
    Mfac = RCONST(1.0);
  }
  else
  {
    const realtype Mp = SUNRceil(SC_MRILL_HSP(C)/SC_MRILL_HFP(C));
    Hfac = H / SC_MRILL_HSP(C);
    Mfac = M / Mp;
  }

  /* compute estimated optimal time step size */
  *Hnew = H * Hfac * SUNRpowerR(es1,a1) * SUNRpowerR(es2,a2);
  const realtype Mnew = M * Mfac * SUNRpowerR(es1,b11) * SUNRpowerR(es2,b12)
                                 * SUNRpowerR(ef1,b21) * SUNRpowerR(ef2,b22);
  *hnew = (*Hnew) / Mnew;

  /* return with success */
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_MRILL(SUNControl C)
{
  SC_MRILL_ESP(C)  = RCONST(1.0);
  SC_MRILL_EFP(C)  = RCONST(1.0);
  SC_MRILL_HSP(C)  = RCONST(1.0);
  SC_MRILL_HFP(C)  = RCONST(1.0);
  SC_MRILL_FIRSTSTEP(C) = SUNTRUE;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_MRILL(SUNControl C)
{
  SC_MRILL_K11(C)  = DEFAULT_K11;
  SC_MRILL_K12(C)  = DEFAULT_K12;
  SC_MRILL_K21(C)  = DEFAULT_K21;
  SC_MRILL_K22(C)  = DEFAULT_K22;
  SC_MRILL_BIAS(C) = DEFAULT_BIAS;
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_MRILL(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "Multirate LL SUNControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k11 = %32Lg\n", SC_MRILL_K11(C));
  fprintf(fptr, "  k12 = %32Lg\n", SC_MRILL_K12(C));
  fprintf(fptr, "  k21 = %32Lg\n", SC_MRILL_K21(C));
  fprintf(fptr, "  k22 = %32Lg\n", SC_MRILL_K22(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SC_MRILL_BIAS(C));
  fprintf(fptr, "  previous slow error = %32Lg\n", SC_MRILL_ESP(C));
  fprintf(fptr, "  previous fast error = %32Lg\n", SC_MRILL_EFP(C));
  fprintf(fptr, "  previous slow step  = %32Lg\n", SC_MRILL_HSP(C));
  fprintf(fptr, "  previous fast step  = %32Lg\n", SC_MRILL_HFP(C));
#else
  fprintf(fptr, "  k11 = %16g\n", SC_MRILL_K11(C));
  fprintf(fptr, "  k12 = %16g\n", SC_MRILL_K12(C));
  fprintf(fptr, "  k21 = %16g\n", SC_MRILL_K21(C));
  fprintf(fptr, "  k22 = %16g\n", SC_MRILL_K22(C));
  fprintf(fptr, "  bias factor = %16g\n", SC_MRILL_BIAS(C));
  fprintf(fptr, "  previous slow error = %16g\n", SC_MRILL_ESP(C));
  fprintf(fptr, "  previous fast error = %16g\n", SC_MRILL_EFP(C));
  fprintf(fptr, "  previous slow step  = %16g\n", SC_MRILL_HSP(C));
  fprintf(fptr, "  previous fast step  = %16g\n", SC_MRILL_HFP(C));
#endif
  fprintf(fptr, "  P = %i (slow method order)\n", SC_MRILL_PSLOW(C));
  fprintf(fptr, "  p = %i (fast method order)\n", SC_MRILL_PFAST(C));
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_MRILL(SUNControl C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SC_MRILL_BIAS(C) = DEFAULT_BIAS;
  } else {
    SC_MRILL_BIAS(C) = bias;
  }

  return SUNCONTROL_SUCCESS;
}

int SUNControlUpdateMRIH_MRILL(SUNControl C, realtype H, realtype h,
                                realtype DSM, realtype dsm)
{
  SC_MRILL_ESP(C) = SUNMAX(SC_MRILL_BIAS(C) * DSM, TINY);
  SC_MRILL_EFP(C) = SUNMAX(SC_MRILL_BIAS(C) * dsm, TINY);
  SC_MRILL_HSP(C) = H;
  SC_MRILL_HFP(C) = h;
  SC_MRILL_FIRSTSTEP(C) = SUNFALSE;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_MRILL(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 9;
  *leniw = 3;
  return SUNCONTROL_SUCCESS;
}
