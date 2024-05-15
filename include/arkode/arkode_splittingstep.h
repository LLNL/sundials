/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * TODO
 *--------------------------------------------------------------*/

#ifndef ARKODE_SPLITTINGSTEP_H_
#define ARKODE_SPLITTINGSTEP_H_

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <stdio.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
  Types : struct ARKodeSplittingCoeffsMem, ARKodeSplittingCoeffs
  ---------------------------------------------------------------*/
struct ARKodeSplittingCoeffsMem
{
  int sequential_methods; /* number of sequential splitting methods */
  int stages;     /* number of stages within each sequential splitting method */
  int partitions; /* number of RHS partitions */
  sunrealtype* alpha;  /* weights for sum over sequential splitting methods */
  sunrealtype*** beta; /* length of each flow, indexed by the sequential method, stage, and partition */
};

typedef _SUNDIALS_STRUCT_ ARKodeSplittingCoeffsMem* ARKodeSplittingCoeffs;

typedef enum
{
  ARKODE_SPLITTING_NONE        = -1, /* ensure enum is signed int */
  ARKODE_MIN_SPLITTING_NUM     = 300,
  ARKODE_SPLITTING_LIE_TROTTER = ARKODE_MIN_SPLITTING_NUM,
  ARKODE_SPLITTING_STRANG,
  ARKODE_SPLITTING_YOSHIDA,
  ARKODE_MAX_SPLITTING_NUM = ARKODE_SPLITTING_YOSHIDA
} ARKODE_SplittingCoeffsID;

/* TODO: remove once SUNStepper is ready */
typedef void* SUNStepper;

/* Coefficient memory management */
SUNDIALS_EXPORT ARKodeSplittingCoeffs
ARKodeSplittingCoeffs_Alloc(int sequential_methods, int stages, int partitions);
SUNDIALS_EXPORT void ARKodeSplittingCoeffs_Free(ARKodeSplittingCoeffs B);
SUNDIALS_EXPORT void ARKodeSplittingCoeffs_Space(ARKodeSplittingCoeffs B,
                                                 sunindextype* liw,
                                                 sunindextype* lrw);
SUNDIALS_EXPORT ARKodeSplittingCoeffs
ARKodeSplittingCoeffs_Copy(ARKodeSplittingCoeffs B);

/* Constructors for splitting coefficients */
SUNDIALS_EXPORT ARKodeSplittingCoeffs
ARKodeSplittingCoeffs_LieTrotter(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoeffs ARKodeSplittingCoeffs_Strang(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoeffs
ARKodeSplittingCoeffs_SymmetricParallel(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoeffs
ARKodeSplittingCoeffs_TripleJump(int partitions, int order);

/* Other coefficient functions */
SUNDIALS_EXPORT void ARKodeSplittingCoeffs_Write(ARKodeSplittingCoeffs B,
                                                 FILE* outfile);
SUNDIALS_EXPORT sunbooleantype ARKodeSplittingCoeffs_(ARKodeSplittingCoeffs B);

/* Parallelization Policy */
typedef int (*ARKParallelExecuteFn)(int i, N_Vector y, void *user_data);

struct ARKodeSplittingExecutionPolicyMem
{
  int (*setup)();

  int (*exectute)();

  int (*free)();

  void *data;
};

typedef _SUNDIALS_STRUCT_ ARKodeSplittingExecutionPolicyMem* ARKodeSplittingExecutionPolicy;


/* ARKODE functions */
void* ARKStepCreateSplitting(SUNStepper *steppers, int partitions, sunrealtype t0, N_Vector y0, SUNContext sunctx);

#ifdef __cplusplus
}
#endif

#endif
