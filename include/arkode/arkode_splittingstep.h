/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the header file for the ARKODE SplittingStep module.
 *--------------------------------------------------------------*/

#ifndef ARKODE_SPLITTINGSTEP_H_
#define ARKODE_SPLITTINGSTEP_H_

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_stepper.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
  Types : struct SplittingStepCoefficientsMem, SplittingStepCoefficients
  ---------------------------------------------------------------*/
struct SplittingStepCoefficientsMem
{
  sunrealtype* alpha;  /* weights for sum over sequential splitting methods */
  sunrealtype*** beta; /* subintegration nodes, indexed by the sequential method, stage, and partition */
  int sequential_methods; /* number of sequential splitting methods */
  int stages;     /* number of stages within each sequential splitting method */
  int partitions; /* number of RHS partitions */
  int order;      /* order of convergence */
};

typedef _SUNDIALS_STRUCT_ SplittingStepCoefficientsMem* SplittingStepCoefficients;

/* Splitting names use the convention
 * ARKODE_SPLITTING_<name>_<stages>_<order>_<partitions> */
typedef enum
{
  ARKODE_SPLITTING_NONE              = -1, /* ensure enum is signed int */
  ARKODE_MIN_SPLITTING_NUM           = 0,
  ARKODE_SPLITTING_LIE_TROTTER_1_1_2 = ARKODE_MIN_SPLITTING_NUM,
  ARKODE_SPLITTING_STRANG_2_2_2,
  ARKODE_SPLITTING_BEST_2_2_2,
  ARKODE_SPLITTING_SUZUKI_3_3_2,
  ARKODE_SPLITTING_RUTH_3_3_2,
  ARKODE_SPLITTING_YOSHIDA_4_4_2,
  ARKODE_SPLITTING_YOSHIDA_8_6_2,
  ARKODE_MAX_SPLITTING_NUM = ARKODE_SPLITTING_YOSHIDA_8_6_2
} ARKODE_SplittingCoefficientsID;

/* Coefficient memory management */
SUNDIALS_EXPORT SplittingStepCoefficients SplittingStepCoefficients_Alloc(
  int sequential_methods, int stages, int partitions);
SUNDIALS_EXPORT SplittingStepCoefficients SplittingStepCoefficients_Create(
  int sequential_methods, int stages, int partitions, int order,
  sunrealtype* alpha, sunrealtype* beta);
SUNDIALS_EXPORT void SplittingStepCoefficients_Destroy(
  SplittingStepCoefficients* coefficients);
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_Copy(SplittingStepCoefficients coefficients);
SUNDIALS_EXPORT void SplittingStepCoefficients_Write(
  SplittingStepCoefficients coefficients, FILE* outfile);

/* Load splitting coefficients */
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_LoadCoefficients(ARKODE_SplittingCoefficientsID id);
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_LoadCoefficientsByName(const char* name);
SUNDIALS_EXPORT const char* SplittingStepCoefficients_IDToName(
  ARKODE_SplittingCoefficientsID id);

/* Constructors for splitting coefficients */
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_LieTrotter(int partitions);
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_Strang(int partitions);
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_Parallel(int partitions);
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_SymmetricParallel(int partitions);
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_ThirdOrderSuzuki(int partitions);
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_TripleJump(int partitions, int order);
SUNDIALS_EXPORT SplittingStepCoefficients
SplittingStepCoefficients_SuzukiFractal(int partitions, int order);

/* Functions for SplittingStep integrator */
SUNDIALS_EXPORT void* SplittingStepCreate(SUNStepper* steppers, int partitions,
                                          sunrealtype t0, N_Vector y0,
                                          SUNContext sunctx);

SUNDIALS_EXPORT int SplittingStepReInit(void* arkode_mem, SUNStepper* steppers,
                                        int partitions, sunrealtype t0,
                                        N_Vector y0);

SUNDIALS_EXPORT int SplittingStepSetCoefficients(
  void* arkode_mem, SplittingStepCoefficients coefficients);

SUNDIALS_EXPORT int SplittingStepGetNumEvolves(void* arkode_mem, int partition,
                                               long int* evolves);

#ifdef __cplusplus
}
#endif

#endif
