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
 * This is the implementation file for the SUNControl_MRICC module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_mricc.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_MRICC_CONTENT(C) ( (SUNControlContent_MRICC)(C->content) )
#define SC_MRICC_K1(C)      ( SC_MRICC_CONTENT(C)->k1 )
#define SC_MRICC_K2(C)      ( SC_MRICC_CONTENT(C)->k2 )
#define SC_MRICC_BIAS(C)    ( SC_MRICC_CONTENT(C)->bias )
#define SC_MRICC_PSLOW(C)   ( SC_MRICC_CONTENT(C)->P )
#define SC_MRICC_PFAST(C)   ( SC_MRICC_CONTENT(C)->p )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1   RCONST(0.42)
#define DEFAULT_K2   RCONST(0.44)
#define DEFAULT_BIAS RCONST(1.5)
#define TINY         RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRICC controller
 */

SUNControl SUNControlMRICC(SUNContext sunctx, int P, int p)
{
  SUNControl C;
  SUNControlContent_MRICC content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid            = SUNControlGetID_MRICC;
  C->ops->estimatemristeps = SUNControlEstimateMRISteps_MRICC;
  C->ops->setdefaults      = SUNControlSetDefaults_MRICC;
  C->ops->write            = SUNControlWrite_MRICC;
  C->ops->seterrorbias     = SUNControlSetErrorBias_MRICC;
  C->ops->space            = SUNControlSpace_MRICC;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_MRICC)malloc(sizeof *content);
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
  SUNControlSetDefaults_MRICC(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set MRICC parameters
 */

int SUNControlMRICC_SetParams(SUNControl C, realtype k1, realtype k2)
{
  /* store legal inputs, and return with success */
  if (k1 >= RCONST(0.0)) { SC_MRICC_K1(C) = k1; }
  if (k2 >= RCONST(0.0)) { SC_MRICC_K2(C) = k2; }
  return SUNCONTROL_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_MRICC(SUNControl C) { return SUNDIALS_CONTROL_MRI_H; }

int SUNControlEstimateMRISteps_MRICC(SUNControl C, realtype H, realtype h,
                                     realtype DSM, realtype dsm,
                                     realtype* Hnew, realtype *hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype al = SC_MRICC_K1(C) / SC_MRICC_PSLOW(C);
  const realtype b1 = (SC_MRICC_PFAST(C) + 1) * SC_MRICC_K1(C)
                    / (SC_MRICC_PSLOW(C) * SC_MRICC_PFAST(C));
  const realtype b2 = -SC_MRICC_K2(C) / SC_MRICC_PFAST(C);
  const realtype es = SUNMAX(SC_MRICC_BIAS(C) * DSM, TINY);
  const realtype ef = SUNMAX(SC_MRICC_BIAS(C) * dsm, TINY);
  const realtype M = SUNRceil(H/h);

  /* compute estimated optimal time step size */
  *Hnew = H * SUNRpowerR(es,al);
  const realtype Mnew = M * SUNRpowerR(es,b1) * SUNRpowerR(ef,b2);
  *hnew = (*Hnew) / Mnew;

  /* return with success */
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_MRICC(SUNControl C)
{
  SC_MRICC_K1(C)   = DEFAULT_K1;
  SC_MRICC_K2(C)   = DEFAULT_K2;
  SC_MRICC_BIAS(C) = DEFAULT_BIAS;
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_MRICC(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "Multirate constant-constant SUNControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k1 = %32Lg\n", SC_MRICC_K1(C));
  fprintf(fptr, "  k2 = %32Lg\n", SC_MRICC_K2(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SC_MRICC_BIAS(C));
#else
  fprintf(fptr, "  k1 = %16g\n", SC_MRICC_K1(C));
  fprintf(fptr, "  k2 = %16g\n", SC_MRICC_K2(C));
  fprintf(fptr, "  bias factor = %16g\n", SC_MRICC_BIAS(C));
#endif
  fprintf(fptr, "  P = %i (slow method order)\n", SC_MRICC_PSLOW(C));
  fprintf(fptr, "  p = %i (fast method order)\n", SC_MRICC_PFAST(C));
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_MRICC(SUNControl C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SC_MRICC_BIAS(C) = DEFAULT_BIAS;
  } else {
    SC_MRICC_BIAS(C) = bias;
  }

  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_MRICC(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 3;
  *leniw = 2;
  return SUNCONTROL_SUCCESS;
}
