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
 * This is the header file for the
 * SUNTimestepHeuristics_Unconstrained module.
 * -----------------------------------------------------------------*/

#ifndef _SUNTIMESTEPHEURISTICS_UNCONSTRAINED_H
#define _SUNTIMESTEPHEURISTICS_UNCONSTRAINED_H

#include <stdio.h>
#include <sundials/sundials_timestepheuristics.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------
 * Unconstrained implementation of SUNTimestepHeuristics
 * ----------------------------------------------------- */

struct SUNTimestepHeuristicsContent_Unconstrained_ {
  long int nst_acc; /* num accuracy-limited internal steps */
};

typedef struct SUNTimestepHeuristicsContent_Unconstrained_ *SUNTimestepHeuristicsContent_Unconstrained;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNTimestepHeuristics SUNTimestepHeuristicsUnconstrained(SUNContext sunctx);
SUNDIALS_EXPORT
SUNTimestepHeuristics_ID SUNTimestepHeuristicsGetID_Unconstrained(SUNTimestepHeuristics H);
SUNDIALS_EXPORT
int SUNTimestepHeuristicsConstrainStep_Unconstrained(SUNTimestepHeuristics H,
                                                     realtype hcur,
                                                     realtype hnew,
                                                     realtype* hconstr);
SUNDIALS_EXPORT
int SUNTimestepHeuristicsConvFail_Unconstrained(SUNTimestepHeuristics H, realtype hcur,
                                                realtype* hconstr);
SUNDIALS_EXPORT
int SUNTimestepHeuristicsReset_Unconstrained(SUNTimestepHeuristics H);
SUNDIALS_EXPORT
int SUNTimestepHeuristicsWrite_Unconstrained(SUNTimestepHeuristics H, FILE* fptr);
int SUNTimestepHeuristicsGetNumAccSteps_Unconstrained(SUNTimestepHeuristics H,
                                                      long int* accsteps);
SUNDIALS_EXPORT
int SUNTimestepHeuristicsSpace_Unconstrained(SUNTimestepHeuristics H, long int *lenrw,
                                             long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNTIMESTEPHEURISTICS_UNCONSTRAINED_H */
