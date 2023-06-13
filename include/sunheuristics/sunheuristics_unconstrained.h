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
 * This is the header file for the SUNHeuristics_Unconstrained module.
 * -----------------------------------------------------------------*/

#ifndef _SUNHEURISTICS_UNCONSTRAINED_H
#define _SUNHEURISTICS_UNCONSTRAINED_H

#include <stdio.h>
#include <sundials/sundials_heuristics.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------
 * Unconstrained implementation of SUNHeuristics
 * --------------------------------------- */

struct _SUNHeuristicsContent_Unconstrained {
  long int nst_acc; /* num accuracy-limited internal steps */
};

typedef struct _SUNHeuristicsContent_Unconstrained *SUNHeuristicsContent_Unconstrained;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNHeuristics SUNHeuristicsUnconstrained(SUNContext sunctx);
SUNDIALS_EXPORT
SUNHeuristics_ID SUNHeuristicsGetID_Unconstrained(SUNHeuristics H);
SUNDIALS_EXPORT
int SUNHeuristicsConstrainStep_Unconstrained(SUNHeuristics H, realtype hcur,
                                       realtype hnew, realtype* hconstr);
SUNDIALS_EXPORT
int SUNHeuristicsConstrainCFail_Unconstrained(SUNHeuristics H, realtype hcur,
                                       realtype* hconstr);
SUNDIALS_EXPORT
int SUNHeuristicsReset_Unconstrained(SUNHeuristics H);
SUNDIALS_EXPORT
int SUNHeuristicsWrite_Unconstrained(SUNHeuristics H, FILE* fptr);
int SUNHeuristicsGetNumAccSteps_Unconstrained(SUNHeuristics H, long int* accsteps);
SUNDIALS_EXPORT
int SUNHeuristicsSpace_Unconstrained(SUNHeuristics H, long int *lenrw, long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNHEURISTICS_UNCONSTRAINED_H */
