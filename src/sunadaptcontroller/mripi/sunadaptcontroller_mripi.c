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
 * This is the implementation file for the SUNControl_MRIPI module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_mripi.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_MRIPI_CONTENT(C) ( (SUNControlContent_MRIPI)(C->content) )
#define SC_MRIPI_K11(C)     ( SC_MRIPI_CONTENT(C)->k11 )
#define SC_MRIPI_K12(C)     ( SC_MRIPI_CONTENT(C)->k12 )
#define SC_MRIPI_K21(C)     ( SC_MRIPI_CONTENT(C)->k21 )
#define SC_MRIPI_K22(C)     ( SC_MRIPI_CONTENT(C)->k22 )
#define SC_MRIPI_BIAS(C)    ( SC_MRIPI_CONTENT(C)->bias )
#define SC_MRIPI_ESP(C)     ( SC_MRIPI_CONTENT(C)->esp )
#define SC_MRIPI_EFP(C)     ( SC_MRIPI_CONTENT(C)->efp )
#define SC_MRIPI_PSLOW(C)   ( SC_MRIPI_CONTENT(C)->P )
#define SC_MRIPI_PFAST(C)   ( SC_MRIPI_CONTENT(C)->p )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K11  RCONST(0.18)
#define DEFAULT_K12  RCONST(0.86)
#define DEFAULT_K21  RCONST(0.34)
#define DEFAULT_K22  RCONST(0.80)
#define DEFAULT_BIAS RCONST(1.5)
#define TINY         RCONST(1.0e-10)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRIPI controller
 */

SUNControl SUNControlMRIPI(SUNContext sunctx, int P, int p)
{
  SUNControl C;
  SUNControlContent_MRIPI content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid             = SUNControlGetID_MRIPI;
  C->ops->estimatemristeps  = SUNControlEstimateMRISteps_MRIPI;
  C->ops->reset             = SUNControlReset_MRIPI;
  C->ops->setdefaults       = SUNControlSetDefaults_MRIPI;
  C->ops->write             = SUNControlWrite_MRIPI;
  C->ops->seterrorbias      = SUNControlSetErrorBias_MRIPI;
  C->ops->updatemrih        = SUNControlUpdateMRIH_MRIPI;
  C->ops->space             = SUNControlSpace_MRIPI;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_MRIPI)malloc(sizeof *content);
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
  SUNControlSetDefaults_MRIPI(C);
  SUNControlReset_MRIPI(C);

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set MRIPI parameters
 */

int SUNControlMRIPI_SetParams(SUNControl C, realtype k11, realtype k12,
                              realtype k21, realtype k22)
{
  /* store legal inputs, and return with success */
  if (k11 >= RCONST(0.0)) { SC_MRIPI_K11(C) = k11; }
  if (k12 >= RCONST(0.0)) { SC_MRIPI_K12(C) = k12; }
  if (k21 >= RCONST(0.0)) { SC_MRIPI_K21(C) = k21; }
  if (k22 >= RCONST(0.0)) { SC_MRIPI_K22(C) = k22; }
  return SUNCONTROL_SUCCESS;
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_MRIPI(SUNControl C) { return SUNDIALS_CONTROL_MRI_H; }

int SUNControlEstimateMRISteps_MRIPI(SUNControl C, realtype H, realtype h,
                                     realtype DSM, realtype dsm,
                                     realtype* Hnew, realtype *hnew)
{
  /* set usable time-step adaptivity parameters */
  const realtype a1 = (SC_MRIPI_K11(C) + SC_MRIPI_K12(C))
                    / (2 * SC_MRIPI_PSLOW(C));
  const realtype a2 = -SC_MRIPI_K11(C) / (2 * SC_MRIPI_PSLOW(C));
  const realtype b11 = (SC_MRIPI_PFAST(C) + 1) * (SC_MRIPI_K11(C) + SC_MRIPI_K12(C))
                     / (2 * SC_MRIPI_PSLOW(C) * SC_MRIPI_PFAST(C));
  const realtype b12 = -(SC_MRIPI_PFAST(C) + 1) * SC_MRIPI_K11(C)
                     / (2 * SC_MRIPI_PSLOW(C) * SC_MRIPI_PFAST(C));
  const realtype b21 = -(SC_MRIPI_K21(C) + SC_MRIPI_K22(C))
                     / (2 * SC_MRIPI_PFAST(C));
  const realtype b22 = SC_MRIPI_K21(C) / (2 * SC_MRIPI_PFAST(C));
  const realtype es1 = SUNMAX(SC_MRIPI_BIAS(C) * DSM, TINY);
  const realtype es2 = SC_MRIPI_ESP(C);
  const realtype ef1 = SUNMAX(SC_MRIPI_BIAS(C) * dsm, TINY);
  const realtype ef2 = SC_MRIPI_EFP(C);
  const realtype M = SUNRceil(H/h);

  /* compute estimated optimal time step size */
  *Hnew = H * SUNRpowerR(es1,a1) * SUNRpowerR(es2,a2);
  const realtype Mnew = M * SUNRpowerR(es1,b11) * SUNRpowerR(es2,b12)
                          * SUNRpowerR(ef1,b21) * SUNRpowerR(ef2,b22);
  *hnew = (*Hnew) / Mnew;

  /* return with success */
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_MRIPI(SUNControl C)
{
  SC_MRIPI_ESP(C) = RCONST(1.0);
  SC_MRIPI_EFP(C) = RCONST(1.0);
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_MRIPI(SUNControl C)
{
  SC_MRIPI_K11(C)  = DEFAULT_K11;
  SC_MRIPI_K12(C)  = DEFAULT_K12;
  SC_MRIPI_K21(C)  = DEFAULT_K21;
  SC_MRIPI_K22(C)  = DEFAULT_K22;
  SC_MRIPI_BIAS(C) = DEFAULT_BIAS;
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_MRIPI(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "Multirate PI SUNControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  k11 = %32Lg\n", SC_MRIPI_K11(C));
  fprintf(fptr, "  k12 = %32Lg\n", SC_MRIPI_K12(C));
  fprintf(fptr, "  k21 = %32Lg\n", SC_MRIPI_K21(C));
  fprintf(fptr, "  k22 = %32Lg\n", SC_MRIPI_K22(C));
  fprintf(fptr, "  bias factor = %32Lg\n", SC_MRIPI_BIAS(C));
  fprintf(fptr, "  previous slow error = %32Lg\n", SC_MRIPI_ESP(C));
  fprintf(fptr, "  previous fast error = %32Lg\n", SC_MRIPI_EFP(C));
#else
  fprintf(fptr, "  k11 = %16g\n", SC_MRIPI_K11(C));
  fprintf(fptr, "  k12 = %16g\n", SC_MRIPI_K12(C));
  fprintf(fptr, "  k21 = %16g\n", SC_MRIPI_K21(C));
  fprintf(fptr, "  k22 = %16g\n", SC_MRIPI_K22(C));
  fprintf(fptr, "  bias factor = %16g\n", SC_MRIPI_BIAS(C));
  fprintf(fptr, "  previous slow error = %16g\n", SC_MRIPI_ESP(C));
  fprintf(fptr, "  previous fast error = %16g\n", SC_MRIPI_EFP(C));
#endif
  fprintf(fptr, "  P = %i (slow method order)\n", SC_MRIPI_PSLOW(C));
  fprintf(fptr, "  p = %i (fast method order)\n", SC_MRIPI_PFAST(C));
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetErrorBias_MRIPI(SUNControl C, realtype bias)
{
  /* set allowed value, otherwise set default */
  if (bias <= RCONST(0.0)) {
    SC_MRIPI_BIAS(C) = DEFAULT_BIAS;
  } else {
    SC_MRIPI_BIAS(C) = bias;
  }

  return SUNCONTROL_SUCCESS;
}

int SUNControlUpdateMRIH_MRIPI(SUNControl C, realtype H, realtype h,
                                realtype DSM, realtype dsm)
{
  SC_MRIPI_ESP(C) = SUNMAX(SC_MRIPI_BIAS(C) * DSM, TINY);
  SC_MRIPI_EFP(C) = SUNMAX(SC_MRIPI_BIAS(C) * dsm, TINY);
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_MRIPI(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 7;
  *leniw = 2;
  return SUNCONTROL_SUCCESS;
}
