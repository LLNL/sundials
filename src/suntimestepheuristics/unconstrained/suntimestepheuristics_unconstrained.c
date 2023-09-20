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
 * SUNTimestepHeuristics_Unconstrained module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suntimestepheuristics/suntimestepheuristics_unconstrained.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SH_CONTENT(H) ( (SUNTimestepHeuristicsContent_Unconstrained)(H->content) )
#define SH_NST_ACC(H) ( SH_CONTENT(H)->nst_acc )


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new unconstrained heuristics module
 */

SUNTimestepHeuristics SUNTimestepHeuristicsUnconstrained(SUNContext sunctx)
{
  SUNTimestepHeuristics H;
  SUNTimestepHeuristicsContent_Unconstrained content;

  /* Create an empty heuristics object */
  H = NULL;
  H = SUNTimestepHeuristics_NewEmpty(sunctx);
  if (H == NULL) { return (NULL); }

  /* Attach operations */
  H->ops->getid          = SUNTimestepHeuristicsGetID_Unconstrained;
  H->ops->constrainstep  = SUNTimestepHeuristicsConstrainStep_Unconstrained;
  H->ops->convfail       = SUNTimestepHeuristicsConvFail_Unconstrained;
  H->ops->reset          = SUNTimestepHeuristicsReset_Unconstrained;
  H->ops->write          = SUNTimestepHeuristicsWrite_Unconstrained;
  H->ops->getnumaccsteps = SUNTimestepHeuristicsGetNumAccSteps_Unconstrained;
  H->ops->space          = SUNTimestepHeuristicsSpace_Unconstrained;

  /* Create content */
  content = NULL;
  content = (SUNTimestepHeuristicsContent_Unconstrained)malloc(sizeof *content);
  if (content == NULL)
  {
    (void) SUNTimestepHeuristics_Destroy(H);
    return (NULL);
  }

  /* Attach content */
  H->content = content;

  /* Fill content with default/reset values */
  SUNTimestepHeuristicsReset_Unconstrained(H);

  return (H);
}


/* -----------------------------------------------------------------
 * implementation of heuristic operations
 * ----------------------------------------------------------------- */

SUNTimestepHeuristics_ID SUNTimestepHeuristicsGetID_Unconstrained(SUNTimestepHeuristics H)
{ return SUN_TIMESTEPHEURISTICS_NULL; }

int SUNTimestepHeuristicsConstrainStep_Unconstrained(SUNTimestepHeuristics H,
                                                     realtype hcur,
                                                     realtype h_acc,
                                                     realtype* hconstr)
{
  /* All steps are considered accuracy-limited */
  SH_NST_ACC(H)++;

  /* Pass recommended step size through and return */
  *hconstr = h_acc;
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristicsConvFail_Unconstrained(SUNTimestepHeuristics H,
                                                realtype hcur,
                                                realtype* hconstr)
{
  /* No recovery is possible upon an algebraic solver convergence failure */
  *hconstr = hcur;
  return SUNTIMESTEPHEURISTICS_CANNOT_DECREASE;
}

int SUNTimestepHeuristicsReset_Unconstrained(SUNTimestepHeuristics H)
{
  SH_NST_ACC(H) = 0;
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristicsWrite_Unconstrained(SUNTimestepHeuristics H,
                                             FILE *fptr)
{
  fprintf(fptr, "Unconstrained SUNTimestepHeuristics module:\n");
  fprintf(fptr, "  Current step count = %li\n", SH_NST_ACC(H));
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristicsGetNumAccSteps_Unconstrained(SUNTimestepHeuristics H,
                                                      long int* accsteps)
{
  *accsteps = SH_NST_ACC(H);
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristicsSpace_Unconstrained(SUNTimestepHeuristics H,
                                             long int* lenrw,
                                             long int* leniw)
{
  *lenrw = 0;
  *leniw = 1;
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}
