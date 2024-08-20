// TODO: should filename be arkode_splittingcoefficients.h instead?

/* -----------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
* This is the header file for ARKode splitting coefficient structures.
 * ---------------------------------------------------------------------------*/

#ifndef ARKODE_SPLITTING_COEFFICIENTS_H_
#define ARKODE_SPLITTING_COEFFICIENTS_H_

#include <stdio.h>
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
  sunrealtype*** beta; /* subintegration nodes, indexed by the sequential method, stage, and partition */
  int sequential_methods; /* number of sequential splitting methods */
  int stages;     /* number of stages within each sequential splitting method */
  int partitions; /* number of RHS partitions */
  int order;      /* order of convergence */
};

typedef _SUNDIALS_STRUCT_ ARKodeSplittingCoefficientsMem* ARKodeSplittingCoefficients;

/* Splitting names use the convention ARKODE_SPLITTING_<name>_<partitions>_<order> */
typedef enum
{
  ARKODE_SPLITTING_NONE        = -1, /* ensure enum is signed int */
  ARKODE_MIN_SPLITTING_NUM     = 300,
  ARKODE_SPLITTING_LIE_TROTTER_2_1 = ARKODE_MIN_SPLITTING_NUM,
  ARKODE_SPLITTING_STRANG_2_2,
  ARKODE_SPLITTING_YOSHIDA_2_4,
  ARKODE_MAX_SPLITTING_NUM = ARKODE_SPLITTING_YOSHIDA_2_4
} ARKODE_SplittingCoefficientsID;

/* Coefficient memory management */
SUNDIALS_EXPORT ARKodeSplittingCoefficients ARKodeSplittingCoefficients_Alloc(
  int sequential_methods, int stages, int partitions);
/* TODO: Ideally, alpha and beta would be const, but that would be inconsistent
 * with other ARKODE function which accept arrays */
SUNDIALS_EXPORT ARKodeSplittingCoefficients ARKodeSplittingCoefficients_Create(
  int sequential_methods, int stages, int partitions, int order,
  sunrealtype* alpha, sunrealtype* beta);
SUNDIALS_EXPORT void ARKodeSplittingCoefficients_Free(ARKodeSplittingCoefficients B);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_Copy(ARKodeSplittingCoefficients B);
SUNDIALS_EXPORT void ARKodeSplittingCoefficients_Write(
  ARKodeSplittingCoefficients coefficients, FILE* outfile);

/* Load splitting coefficients */
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_LoadCoefficients(ARKODE_SplittingCoefficientsID id);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_LoadCoefficientsByName(const char* name);
SUNDIALS_EXPORT const char* ARKodeSplittingCoefficients_IDToName(
  ARKODE_SplittingCoefficientsID id);

/* Constructors for splitting coefficients */
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_LieTrotter(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_Strang(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_Parallel(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_SymmetricParallel(int partitions);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_TripleJump(int partitions, int order);
SUNDIALS_EXPORT ARKodeSplittingCoefficients
ARKodeSplittingCoefficients_SuzukiFractal(int partitions, int order);

#ifdef __cplusplus
}
#endif

#endif