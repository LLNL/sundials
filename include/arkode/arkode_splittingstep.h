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

#include <arkode/arkode_mristep.h>
#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
  Types : struct ARKodeSplittingCoefficientsMem, ARKodeSplittingCoefficients
  ---------------------------------------------------------------*/
struct ARKodeSplittingCoefficientsMem
{
  sunrealtype* alpha;  /* weights for sum over sequential splitting methods */
  sunrealtype*** beta; /* length of each flow, indexed by the sequential method, stage, and partition */
  int sequential_methods; /* number of sequential splitting methods */
  int stages;     /* number of stages within each sequential splitting method */
  int partitions; /* number of RHS partitions */
  int order;      /* order of convergence */
};

typedef _SUNDIALS_STRUCT_ ARKodeSplittingCoefficientsMem* ARKodeSplittingCoefficients;

typedef enum
{
  ARKODE_SPLITTING_NONE        = -1, /* ensure enum is signed int */
  ARKODE_MIN_SPLITTING_NUM     = 300,
  ARKODE_SPLITTING_LIE_TROTTER = ARKODE_MIN_SPLITTING_NUM,
  ARKODE_SPLITTING_STRANG,
  ARKODE_SPLITTING_YOSHIDA,
  ARKODE_MAX_SPLITTING_NUM = ARKODE_SPLITTING_YOSHIDA
} ARKODE_SplittingCoefficientsID;

/* TODO: remove once SUNStepper is ready */
typedef MRIStepInnerStepper SUNStepper;

/* Coefficient memory management */
SUNDIALS_EXPORT ARKodeSplittingCoefficients ARKodeSplittingCoefficients_Alloc(
  int sequential_methods, int stages, int partitions, int order);
SUNDIALS_EXPORT void ARKodeSplittingCoefficients_Free(ARKodeSplittingCoefficients B);
SUNDIALS_EXPORT void ARKodeSplittingCoefficients_Space(ARKodeSplittingCoefficients B,
                                                 sunindextype* liw,
                                                 sunindextype* lrw);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_Copy(ARKodeSplittingCoefficients B);

/* Constructors for splitting coefficients */
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_LieTrotter(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoefficients ARKodeSplittingCoefficients_Strang(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_SymmetricParallel(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_TripleJump(int partitions, int order);

/* Other coefficient functions */
SUNDIALS_EXPORT void ARKodeSplittingCoefficients_Write(ARKodeSplittingCoefficients B,
                                                 FILE* outfile);
SUNDIALS_EXPORT sunbooleantype ARKodeSplittingCoefficients_(ARKodeSplittingCoefficients B);

/* Parallelization Policy */
typedef int (*ARKParallelExecuteFn)(int i, N_Vector y, void* user_data);

struct ARKodeSplittingExecutionPolicyMem
{
  int (*setup)();

  int (*exectute)();

  int (*free)();

  void* data;
};

typedef _SUNDIALS_STRUCT_ ARKodeSplittingExecutionPolicyMem* ARKodeSplittingExecutionPolicy;

/* ARKODE functions */
SUNDIALS_EXPORT void* SplittingStepCreate(SUNStepper* steppers,
                                          const int partitions,
                                          const sunrealtype t0, N_Vector y0,
                                          SUNContext sunctx);

SUNDIALS_EXPORT int SplittingStep_SetCoefficients(
  void* arkode_mem, ARKodeSplittingCoefficients coefficients);

SUNDIALS_EXPORT int SplittingStep_SetExecutionPolicy(
  void* arkode_mem, ARKodeSplittingExecutionPolicy policy);

#ifdef __cplusplus
}
#endif

#endif
