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
 * This is the implementation file for the SUNControl_MRIPID module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_mripid.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_MRIPID_CONTENT(C) ( (SUNControlContent_MRIPID)(C->content) )
#define SC_MRIPID_K11(C)     ( SC_MRIPID_CONTENT(C)->k11 )
#define SC_MRIPID_K12(C)     ( SC_MRIPID_CONTENT(C)->k12 )
#define SC_MRIPID_K13(C)     ( SC_MRIPID_CONTENT(C)->k13 )
#define SC_MRIPID_K21(C)     ( SC_MRIPID_CONTENT(C)->k21 )
#define SC_MRIPID_K22(C)     ( SC_MRIPID_CONTENT(C)->k22 )
#define SC_MRIPID_K23(C)     ( SC_MRIPID_CONTENT(C)->k23 )
#define SC_MRIPID_BIAS(C)    ( SC_MRIPID_CONTENT(C)->bias )
#define SC_MRIPID_ESP(C)     ( SC_MRIPID_CONTENT(C)->esp )
#define SC_MRIPID_EFP(C)     ( SC_MRIPID_CONTENT(C)->efp )
#define SC_MRIPID_ESPP(C)    ( SC_MRIPID_CONTENT(C)->espp )
#define SC_MRIPID_EFPP(C)    ( SC_MRIPID_CONTENT(C)->efpp )
#define SC_MRIPID_PSLOW(C)   ( SC_MRIPID_CONTENT(C)->P )
#define SC_MRIPID_PFAST(C)   ( SC_MRIPID_CONTENT(C)->p )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K11  RCONST(0.34)
#define DEFAULT_K12  RCONST(0.1)
#define DEFAULT_K13  RCONST(0.78)
#define DEFAULT_K21  RCONST(0.46)
#define DEFAULT_K22  RCONST(0.42)
#define DEFAULT_K23  RCONST(0.74)
#define DEFAULT_BIAS RCONST(1.5)
#define TINY         RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRIPID controller
 */

SUNControl SUNControlMRIPID(SUNContext sunctx, int P, int p)
{
  SUNControl C;
  SUNControlContent_MRIPID content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid             = SUNControlGetID_MRIPID;
  C->ops->estimatemristeps  = SUNControlEstimateMRISteps_MRIPID;
  C->ops->reset             = SUNControlReset_MRIPID;
  C->ops->setdefaults       = SUNControlSetDefaults_MRIPID;
  C->ops->write             = SUNControlWrite_MRIPID;
  C->ops->seterrorbias      = SUNControlSetErrorBias_MRIPID;
  C->ops->updatemrih        = SUNControlUpdateMRIH_MRIPID;
  C->ops->space             = SUNControlSpace_MRIPID;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_MRIPID)malloc(sizeof *content);
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
  SUNControlSetDefaults_MRIPID(C);
  SUNControlReset_MRIPID(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set MRIPID parameters
 */

int SUNControlMRIPID_SetParams(SUNControl C, realtype k11,
                               realtype k12, realtype k13,
                               realtype k21, realtype k22,
                               realtype k23)
{
  /* store legal inputs, and return with success */
  if (k11 >= RCONST(0.0)) { SC_MRIPID_K11(C) = k11; }
  if (k12 >= RCONST(0.0)) { SC_MRIPID_K12(C) = k12; }
  if (k13 >= RCONST(0.0)) { SC_MRIPID_K13(C) = k13; }
  if (k21 >= RCONST(0.0)) { SC_MRIPID_K21(C) = k21; }
  if (k22 >= RCONST(0.0)) { SC_MRIPID_K22(C) = k22; }
  if (k23 >= RCONST(0.0)) { SC_MRIPID_K23(C) = k23; }
  return SUNCONTROL_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_MRIPID(SUNControl C) { return SUNDIALS_CONTROL_MRI_H; }

int SUNControlEstimateMRISteps_MRIPID(SUNControl C, realtype H, realtype h,
                                      realtype DSM, realtype dsm,
                                      realtype* Hnew, realtype *hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype a1 = (SC_MRIPID_K11(C) + SC_MRIPID_K12(C) + SC_MRIPID_K13(C))
                    / (3 * SC_MRIPID_PSLOW(C));
  const realtype a2 = -(SC_MRIPID_K11(C) + SC_MRIPID_K12(C))
                    / (3 * SC_MRIPID_PSLOW(C));
  const realtype a3 = SC_MRIPID_K11(C) / (3 * SC_MRIPID_PSLOW(C));
  const realtype b11 = (SC_MRIPID_PFAST(C) + 1)
                     * (SC_MRIPID_K11(C) + SC_MRIPID_K12(C) + SC_MRIPID_K13(C))
                     / (3 * SC_MRIPID_PSLOW(C) * SC_MRIPID_PFAST(C));
  const realtype b12 = -(SC_MRIPID_PFAST(C) + 1)
                     * (SC_MRIPID_K11(C) + SC_MRIPID_K12(C))
                     / (3 * SC_MRIPID_PSLOW(C) * SC_MRIPID_PFAST(C));
  const realtype b13 = (SC_MRIPID_PFAST(C) + 1) * SC_MRIPID_K11(C)
                     / (3 * SC_MRIPID_PSLOW(C) * SC_MRIPID_PFAST(C));
  const realtype b21 = -(SC_MRIPID_K21(C) + SC_MRIPID_K22(C) + SC_MRIPID_K23(C))
                     / (3 * SC_MRIPID_PFAST(C));
  const realtype b22 = (SC_MRIPID_K21(C) + SC_MRIPID_K22(C))
                     / (3 * SC_MRIPID_PFAST(C));
  const realtype b23 = -SC_MRIPID_K21(C) / (3 * SC_MRIPID_PFAST(C));
  const realtype es1 = SUNMAX(SC_MRIPID_BIAS(C) * DSM, TINY);
  const realtype es2 = SC_MRIPID_ESP(C);
  const realtype es3 = SC_MRIPID_ESPP(C);
  const realtype ef1 = SUNMAX(SC_MRIPID_BIAS(C) * dsm, TINY);
  const realtype ef2 = SC_MRIPID_EFP(C);
  const realtype ef3 = SC_MRIPID_EFPP(C);
  const realtype M = SUNRceil(H/h);

  /* compute estimated optimal time step size */
  *Hnew = H * SUNRpowerR(es1,a1) * SUNRpowerR(es2,a2) * SUNRpowerR(es3,a3);
  const realtype Mnew = M * SUNRpowerR(es1,b11) * SUNRpowerR(es2,b12)
                          * SUNRpowerR(es3,b13) * SUNRpowerR(ef1,b21)
                          * SUNRpowerR(ef2,b22) * SUNRpowerR(ef3,b23);
  *hnew = (*Hnew) / Mnew;

  /* return with success */
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_MRIPID(SUNControl C)
{
  SC_MRIPID_ESP(C)  = RCONST(1.0);
  SC_MRIPID_EFP(C)  = RCONST(1.0);
  SC_MRIPID_ESPP(C) = RCONST(1.0);
  SC_MRIPID_EFPP(C) = RCONST(1.0);
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_MRIPID(SUNControl C)
{
  SC_MRIPID_K11(C)  = DEFAULT_K11;
  SC_MRIPID_K12(C)  = DEFAULT_K12;
  SC_MRIPID_K13(C)  = DEFAULT_K13;
  SC_MRIPID_K21(C)  = DEFAULT_K21;
  SC_MRIPID_K22(C)  = DEFAULT_K22;
  SC_MRIPID_K23(C)  = DEFAULT_K23;
  SC_MRIPID_BIAS(C) = DEFAULT_BIAS;
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_MRIPID(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "Multirate PID SUNControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k11 = %32Lg\n", SC_MRIPID_K11(C));
  fprintf(fptr, "  k12 = %32Lg\n", SC_MRIPID_K12(C));
  fprintf(fptr, "  k13 = %32Lg\n", SC_MRIPID_K13(C));
  fprintf(fptr, "  k21 = %32Lg\n", SC_MRIPID_K21(C));
  fprintf(fptr, "  k22 = %32Lg\n", SC_MRIPID_K22(C));
  fprintf(fptr, "  k23 = %32Lg\n", SC_MRIPID_K23(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SC_MRIPID_BIAS(C));
  fprintf(fptr, "  previous slow errors = %32Lg  %32Lg\n",
          SC_MRIPID_ESP(C), SC_MRIPID_ESPP(C));
  fprintf(fptr, "  previous fast errors = %32Lg  %32Lg\n",
          SC_MRIPID_EFP(C), SC_MRIPID_EFPP(C));
#else
  fprintf(fptr, "  k11 = %16g\n", SC_MRIPID_K11(C));
  fprintf(fptr, "  k12 = %16g\n", SC_MRIPID_K12(C));
  fprintf(fptr, "  k13 = %16g\n", SC_MRIPID_K13(C));
  fprintf(fptr, "  k21 = %16g\n", SC_MRIPID_K21(C));
  fprintf(fptr, "  k22 = %16g\n", SC_MRIPID_K22(C));
  fprintf(fptr, "  k23 = %16g\n", SC_MRIPID_K23(C));
  fprintf(fptr, "  bias factor = %16g\n", SC_MRIPID_BIAS(C));
  fprintf(fptr, "  previous slow errors = %16g  %16g\n",
          SC_MRIPID_ESP(C), SC_MRIPID_ESPP(C));
  fprintf(fptr, "  previous fast errors = %16g  %16g\n",
          SC_MRIPID_EFP(C), SC_MRIPID_EFPP(C));
#endif
  fprintf(fptr, "  P = %i (slow method order)\n", SC_MRIPID_PSLOW(C));
  fprintf(fptr, "  p = %i (fast method order)\n", SC_MRIPID_PFAST(C));
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_MRIPID(SUNControl C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SC_MRIPID_BIAS(C) = DEFAULT_BIAS;
  } else {
    SC_MRIPID_BIAS(C) = bias;
  }

  return SUNCONTROL_SUCCESS;
}

int SUNControlUpdateMRIH_MRIPID(SUNControl C, realtype H, realtype h,
                                realtype DSM, realtype dsm)
{
  SC_MRIPID_ESPP(C) = SC_MRIPID_ESP(C);
  SC_MRIPID_EFPP(C) = SC_MRIPID_EFP(C);
  SC_MRIPID_ESP(C)  = SUNMAX(SC_MRIPID_BIAS(C) * DSM, TINY);
  SC_MRIPID_EFP(C)  = SUNMAX(SC_MRIPID_BIAS(C) * dsm, TINY);
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_MRIPID(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 11;
  *leniw = 2;
  return SUNCONTROL_SUCCESS;
}
