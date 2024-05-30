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

#include <arkode/arkode_splittingstep.h>
#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_Alloc(
  const int sequential_methods, const int stages, const int partitions,
  const int order)
{
  if (sequential_methods < 1 || stages < 1 || partitions < 1) { return NULL; }

  const ARKodeSplittingCoefficients coefficients =
    (ARKodeSplittingCoefficients)malloc(sizeof(*coefficients));
  if (coefficients == NULL) { return NULL; }

  coefficients->sequential_methods = sequential_methods;
  coefficients->stages             = stages;
  coefficients->partitions         = partitions;
  coefficients->order              = order;

  coefficients->alpha = (sunrealtype*)calloc(sequential_methods,
                                             sizeof(*coefficients->alpha));
  if (coefficients->alpha == NULL)
  {
    ARKodeSplittingCoefficients_Free(coefficients);
    return NULL;
  }

  /* In order for beta to be indexed like beta[i][j][k], we can't use a single
     malloc. The j index requires allocating a matrix of pointers into beta. The
     i index requires allocating an array of pointers into that matrix. */

  /* Array of pointers for index i */
  coefficients->beta =
    (sunrealtype***)malloc(sequential_methods * sizeof(*coefficients->beta));
  if (coefficients->beta == NULL)
  {
    ARKodeSplittingCoefficients_Free(coefficients);
    return NULL;
  }

  /* Matrix of pointers for index j */
  sunrealtype** beta_cols =
    (sunrealtype**)malloc(sequential_methods * stages * sizeof(*beta_cols));
  if (beta_cols == NULL)
  {
    ARKodeSplittingCoefficients_Free(coefficients);
    return NULL;
  }

  /* Contiguous memory to store the beta tensor */
  sunrealtype* beta_mem =
    (sunrealtype*)calloc(sequential_methods * stages * partitions,
                         sizeof(*beta_mem));
  if (beta_mem == NULL)
  {
    ARKodeSplittingCoefficients_Free(coefficients);
    return NULL;
  }

  for (int i = 0; i < sequential_methods; i++)
  {
    coefficients->beta[i] = &beta_cols[i];
    for (int j = 0; j < stages; j++)
    {
      coefficients->beta[i][j] = &beta_mem[i * stages + j];
    }
  }

  return coefficients;
}

void ARKodeSplittingCoefficients_Free(ARKodeSplittingCoefficients coefficients)
{
  if (coefficients != NULL)
  {
    free(coefficients->alpha);
    if (coefficients->beta != NULL)
    {
      free(coefficients->beta[0][0]);
      free(coefficients->beta[0]);
      free(coefficients->beta);
    }
  }

  free(coefficients);
}

void ARKodeSplittingCoefficients_Space(ARKodeSplittingCoefficients coefficients,
                                       sunindextype* liw, sunindextype* lrw)
{
  /* initialize outputs and return if coefficients is not allocated */
  *liw = 0;
  *lrw = 0;
  if (coefficients == NULL) { return; }

  /* ints for # of sequential methods, stages, partitions, and order */
  *liw = 4;

  if (coefficients->beta)
  {
    /* pointers for index i of beta[i][j][k] */
    *liw += coefficients->sequential_methods;
  }

  if (coefficients->beta[0])
  {
    /* pointers for index j of beta[i][j][k] */
    *liw += coefficients->sequential_methods * coefficients->stages;
  }

  if (coefficients->alpha) { *lrw += coefficients->sequential_methods; }
  if (coefficients->beta[0][0])
  {
    *lrw += coefficients->sequential_methods * coefficients->stages *
            coefficients->partitions;
  }
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_Copy(
  ARKodeSplittingCoefficients coefficients)
{
  if (coefficients == NULL) { return NULL; }

  ARKodeSplittingCoefficients coefficientsCopy =
    ARKodeSplittingCoefficients_Alloc(coefficients->sequential_methods,
                                      coefficients->stages,
                                      coefficients->partitions,
                                      coefficients->order);
  if (coefficientsCopy == NULL) { return NULL; }

  memcpy(coefficientsCopy->alpha, coefficients->alpha,
         coefficients->sequential_methods * sizeof(*coefficients->alpha));

  /* beta[0][0] points to the contiguous memory allocation, so we can copy it
     with a single memcpy */
  memcpy(coefficientsCopy->beta[0][0], coefficients->beta[0][0],
         coefficients->sequential_methods * coefficients->stages *
           coefficients->partitions * sizeof(*coefficients->beta));

  return coefficientsCopy;
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_LoadCoefficients(
  ARKODE_SplittingCoefficientsID method)
{
  switch (method)
  {
#define ARK_SPLITTING_COEFFICIENTS(name, coeff) \
  case name: coeff break;
#include "arkode_splitting_coefficients.def"
#undef ARK_SPLITTING_COEFFICIENTS

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown splitting coefficients");
    return NULL;
  }
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_LoadCoefficientsByName(
  const char* method)
{
#define ARK_SPLITTING_COEFFICIENTS(name, coeff) \
  if (strcmp(#name, method) == 0) coeff
#include "arkode_splitting_coefficients.def"
#undef ARK_SPLITTING_COEFFICIENTS

  arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                  "Unknown splitting coefficients");

  return NULL;
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_LieTrotter(const int partitions)
{
  if (partitions < 1) { return NULL; }
  const ARKodeSplittingCoefficients coefficients =
    ARKodeSplittingCoefficients_Alloc(1, 1, partitions, 1);
  if (coefficients == NULL) { return NULL; }

  coefficients->alpha[0] = SUN_RCONST(1.0);
  for (int i = 0; i < partitions; i++)
  {
    coefficients->beta[0][0][i] = SUN_RCONST(1.0);
  }

  return coefficients;
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_Strang(const int partitions)
{
  return ARKodeSplittingCoefficients_TripleJump(partitions, 2);
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_Parallel(const int partitions)
{
  if (partitions < 1) { return NULL; }

  const ARKodeSplittingCoefficients coefficients =
    ARKodeSplittingCoefficients_Alloc(partitions + 1, 1, partitions, 1);
  if (coefficients == NULL) { return NULL; }

  for (int i = 0; i < partitions; i++)
  {
    coefficients->alpha[0]      = SUN_RCONST(1.0);
    coefficients->beta[i][0][i] = SUN_RCONST(1.0);
  }

  coefficients->alpha[partitions] = 1 - partitions;

  return coefficients;
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_SymmetricParallel(
  const int partitions)
{
  if (partitions < 1) { return NULL; }

  const ARKodeSplittingCoefficients coefficients =
    ARKodeSplittingCoefficients_Alloc(2, partitions, partitions, 2);
  if (coefficients == NULL) { return NULL; }

  coefficients->alpha[0] = SUN_RCONST(0.5);
  coefficients->alpha[1] = SUN_RCONST(0.5);

  for (int i = 0; i < partitions; i++)
  {
    coefficients->beta[0][0][i]                  = SUN_RCONST(1.0);
    coefficients->beta[1][i][partitions - i - 1] = SUN_RCONST(1.0);
  }

  return coefficients;
}

static sunrealtype** arkodeSplittingCoefficients_TripleJump(
  const int partitions, const int order, sunrealtype** beta,
  sunrealtype parent_length)
{
  // The base case is an order 2 Strang splitting
  if (order == 2)
  {
    for (int i = 0; i < partitions; i++)
    {
      // Use += here because the first integration may combine with the last
      // integration of the previous stage
      beta[0][i] += SUN_RCONST(0.5) * parent_length;
      beta[i][partitions - i - 1] += SUN_RCONST(0.5) * parent_length;
    }

    // Advance the beta pointer to the last stage of this Strang splitting
    return &beta[partitions - 1];
  }

  // Length of jump 1 and 3
  sunrealtype z1 = SUN_RCONST(1.0) /
                   (SUN_RCONST(2.0) -
                    SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0) / (order - 1)));
  // Length of jump 2
  sunrealtype z0 = SUN_RCONST(1.0) - SUN_RCONST(2.0) * z1;

  // Perform the triple jump
  beta = arkodeSplittingCoefficients_TripleJump(partitions, order - 2, beta, z1);
  beta = arkodeSplittingCoefficients_TripleJump(partitions, order - 2, beta, z0);
  beta = arkodeSplittingCoefficients_TripleJump(partitions, order - 2, beta, z1);

  return beta;
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_TripleJump(
  const int partitions, const int order)
{
  if (partitions < 1 || order < 2 || order % 2 != 0)
  {
    // Only even orders allowed
    return NULL;
  }

  const int stages = 1 + (partitions - 1) * SUNRpowerI(3, order / 2 - 1);
  const ARKodeSplittingCoefficients coefficients =
    ARKodeSplittingCoefficients_Alloc(1, stages, partitions, order);
  if (coefficients == NULL) { return NULL; }

  coefficients->alpha[0] = SUN_RCONST(1.0);
  arkodeSplittingCoefficients_TripleJump(partitions, order,
                                         coefficients->beta[0], SUN_RCONST(1.0));

  return coefficients;
}

/*---------------------------------------------------------------
  Routine to print a Butcher table structure
  ---------------------------------------------------------------*/
void ARKodeSplittingCoefficients_Write(ARKodeSplittingCoefficients coefficients,
                                       FILE* outfile)
{
  if (coefficients == NULL || coefficients->alpha == NULL ||
      coefficients->beta == NULL || coefficients->beta[0] == NULL ||
      coefficients->beta[0][0] == NULL)
  {
    return;
  }

  fprintf(outfile, "  sequential methods = %i\n",
          coefficients->sequential_methods);
  fprintf(outfile, "  stages = %i\n", coefficients->stages);
  fprintf(outfile, "  partitions = %i\n", coefficients->partitions);
  fprintf(outfile, "  order = %i\n", coefficients->order);
  fprintf(outfile, "  alpha = \n");
  for (int i = 0; i < coefficients->sequential_methods; i++)
  {
    fprintf(outfile, "%" RSYM "  ", coefficients->alpha[i]);
  }

  for (int i = 0; i < coefficients->sequential_methods; i++)
  {
    fprintf(outfile, "  beta[%i] = \n", i);
    for (int j = 0; j < coefficients->stages; j++)
    {
      fprintf(outfile, "      ");
      for (int k = 0; k < coefficients->partitions; k++)
      {
        fprintf(outfile, "%" RSYMW "  ", coefficients->beta[i][j][k]);
      }
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
  }
}
