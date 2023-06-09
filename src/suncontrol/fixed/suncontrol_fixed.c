/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
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
 * This is the implementation file for the fixed step implementation
 * of the SUNControl module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_fixed.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_FIXED_CONTENT(C)     ( (SUNControlContent_Fixed)(C->content) )
#define SC_FIXED_H(C)           ( SC_FIXED_CONTENT(C)->hfixed )

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new fixed-step controller
 */

SUNControl SUNControlFixed(realtype hfixed, SUNContext sunctx)
{
  SUNControl C;
  SUNControlContent_Fixed content;

  /* return with NULL controller on illegal fixed step input */
  if (hfixed == RCONST(0.0))
    return (NULL);

  /* Create an empty controller object */
  C = NULL;
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid        = SUNControlGetID_Fixed;
  C->ops->destroy      = SUNControlDestroy_Fixed;
  C->ops->estimatestep = SUNControlEstimateStep_Fixed;
  C->ops->write        = SUNControlWrite_Fixed;
  C->ops->space        = SUNControlSpace_Fixed;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_Fixed)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNControlDestroy(C);
    return (NULL);
  }

  /* Fill content */
  content->hfixed = hfixed;

  /* Attach content and return */
  C->content = content;
  return (C);
}



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_Fixed(SUNControl C) { return SUNDIALS_CONTROL_H; }

void SUNControlDestroy_Fixed(SUNControl C)
{
  if (C == NULL) { return; }

  /* free content */
  if (C->content != NULL) {
    free(C->content);
    C->content = NULL;
  }

  /* free ops and controller */
  if (C->ops) {
    free(C->ops);
    C->ops = NULL;
  }
  free(C);
  C = NULL;

  return;
}

int SUNControlEstimateStep_Fixed(SUNControl C, realtype h,
                                 realtype dsm, realtype* hnew)
{
  /* set stored fixed step size and return with success */
  *hnew = SC_FIXED_H(C);
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_Fixed(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "SUNControl_Fixed module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  hfixed = %12Lg", SC_FIXED_H(C));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  fprintf(fptr, "  hfixed = %12g", SC_FIXED_H(C));
#else
  fprintf(fptr, "  hfixed = %12g", SC_FIXED_H(C));
#endif
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_Fixed(SUNControl C, long int* lenrw, long int* leniw)
{
  *lenrw = 1;
  *leniw = 0;
  return SUNCONTROL_SUCCESS;
}
