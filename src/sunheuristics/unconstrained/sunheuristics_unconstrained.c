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
 * This is the implementation file for the
 * SUNHeuristics_Unconstrained module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunheuristics/sunheuristics_unconstrained.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SH_CONTENT(H) ( (SUNHeuristicsContent_Unconstrained)(H->content) )
#define SH_NST_ACC(H) ( SH_CONTENT(H)->nst_acc )


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new unconstrained heuristics module
 */

SUNHeuristics SUNHeuristicsUnconstrained(SUNContext sunctx)
{
  SUNHeuristics H;
  SUNHeuristicsContent_Unconstrained content;

  /* Create an empty heuristics object */
  H = NULL;
  H = SUNHeuristics_NewEmpty(sunctx);
  if (H == NULL) { return (NULL); }

  /* Attach operations */
  H->ops->getid          = SUNHeuristicsGetID_Unconstrained;
  H->ops->constrainstep  = SUNHeuristicsConstrainStep_Unconstrained;
  H->ops->convfail       = SUNHeuristicsConvFail_Unconstrained;
  H->ops->reset          = SUNHeuristicsReset_Unconstrained;
  H->ops->write          = SUNHeuristicsWrite_Unconstrained;
  H->ops->getnumaccsteps = SUNHeuristicsGetNumAccSteps_Unconstrained;
  H->ops->space          = SUNHeuristicsSpace_Unconstrained;

  /* Create content */
  content = NULL;
  content = (SUNHeuristicsContent_Unconstrained)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNHeuristics_Destroy(H);
    return (NULL);
  }

  /* Attach content */
  H->content = content;

  /* Fill content with default/reset values */
  SUNHeuristicsReset_Unconstrained(H);

  return (H);
}


/* -----------------------------------------------------------------
 * implementation of heuristic operations
 * ----------------------------------------------------------------- */

SUNHeuristics_ID SUNHeuristicsGetID_Unconstrained(SUNHeuristics H)
{ return SUNDIALS_HEURISTICS_NULL; }

int SUNHeuristicsConstrainStep_Unconstrained(SUNHeuristics H, realtype hcur,
                                             realtype h_acc, realtype* hconstr)
{
  /* All steps are considered accuracy-limited */
  SH_NST_ACC(H)++;

  /* Pass recommended step size through and return */
  *hconstr = h_acc;
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsConvFail_Unconstrained(SUNHeuristics H, realtype hcur,
                                        realtype* hconstr)
{
  /* No recovery is possible upon an algebraic solver convergence failure */
  *hconstr = hcur;
  return SUNHEURISTICS_CANNOT_DECREASE;
}

int SUNHeuristicsReset_Unconstrained(SUNHeuristics H)
{
  SH_NST_ACC(H) = 0;
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsWrite_Unconstrained(SUNHeuristics H, FILE *fptr)
{
  fprintf(fptr, "Unconstrained SUNHeuristics module:\n");
  fprintf(fptr, "  Current step count = %li\n", SH_NST_ACC(H));
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsGetNumAccSteps_Unconstrained(SUNHeuristics H, long int* accsteps)
{
  *accsteps = SH_NST_ACC(H);
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSpace_Unconstrained(SUNHeuristics H, long int* lenrw, long int* leniw)
{
  *lenrw = 0;
  *leniw = 1;
  return SUNHEURISTICS_SUCCESS;
}
