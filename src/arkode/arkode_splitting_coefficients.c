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
  sunrealtype** beta_cols = (sunrealtype**)malloc(
    sequential_methods * (stages + 1) * sizeof(*beta_cols));
  if (beta_cols == NULL)
  {
    ARKodeSplittingCoefficients_Free(coefficients);
    return NULL;
  }

  /* Contiguous memory to store the beta tensor */
  sunrealtype* beta_mem =
    (sunrealtype*)calloc(sequential_methods * (stages + 1) * partitions,
                         sizeof(*beta_mem));
  if (beta_mem == NULL)
  {
    ARKodeSplittingCoefficients_Free(coefficients);
    return NULL;
  }

  /* Set pointers into the beta tensor */
  for (int i = 0; i < sequential_methods; i++)
  {
    coefficients->beta[i] = &beta_cols[i * (stages + 1)];
    for (int j = 0; j <= stages; j++)
    {
      coefficients->beta[i][j] = &beta_mem[(i * (stages + 1) + j) * partitions];
    }
  }

  return coefficients;
}

void ARKodeSplittingCoefficients_Free(const ARKodeSplittingCoefficients coefficients)
{
  // TODO: should argument be a pointer?
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
         coefficients->sequential_methods * (coefficients->stages + 1) *
           coefficients->partitions * sizeof(*coefficients->beta));

  return coefficientsCopy;
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_LoadCoefficients(
  const ARKODE_SplittingCoefficientsID method)
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
  const char* const method)
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
    coefficients->beta[0][1][i] = SUN_RCONST(1.0);
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
    coefficients->alpha[i] = SUN_RCONST(1.0);
    coefficients->beta[i][1][i] = SUN_RCONST(1.0);
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
    coefficients->beta[0][partitions][i] = SUN_RCONST(1.0);
    for (int j = partitions - i - 1; j < partitions; j++) {
      coefficients->beta[1][i + 1][j] = SUN_RCONST(1.0);
    }
  }

  return coefficients;
}

static sunrealtype* const* arkodeSplittingCoefficients_ComposeStrangHelper(
  const int partitions, const int order, const int composition_stages,
  const sunrealtype start, const sunrealtype end, sunrealtype* const* const beta)
{
  const sunrealtype diff = end - start;
  if (order == 2)
  {
    /* The base case is an order 2 Strang splitting */
    const sunrealtype mid = start + diff / SUN_RCONST(2.0);
    for (int j = 1; j <= partitions; j++)
    {
      for (int k = 0; k < partitions; k++)
      {
        beta[j][k] = (partitions - j > k) ? mid : end;
      }
    }

    return &beta[partitions - 1];
  }

  sunrealtype* const* beta_cur = beta;
  sunrealtype start_cur  = start;
  /* This is essentially the gamma coefficient from Geometric Numerical
   * Integration (https://doi.org/10.1007/3-540-30666-8) pg 44-45 scaled by the
   * current interval */
  const sunrealtype gamma = diff / (composition_stages - 1 - SUNRpowerR(composition_stages - 1, SUN_RCONST(1.0) / (order - 1)));
  for (int i = 1; i <= composition_stages; i++)
  {
    /* To avoid roundoff issues, this ensures end_cur=1 for the last value of i*/
    const sunrealtype end_cur = 2 * i < composition_stages ? (start + i * gamma) : (end + (i - composition_stages) * gamma);
    /* Recursively generate coefficients and shift beta_cur */
    beta_cur =
      arkodeSplittingCoefficients_ComposeStrangHelper(partitions, order - 2,
                                                      composition_stages,
                                                      start_cur, end_cur,
                                                      beta_cur);
    start_cur = end_cur;
  }

  return beta_cur;
}

static ARKodeSplittingCoefficients arkodeSplittingCoefficients_ComposeStrang(
  const int partitions, const int order, const int composition_stages)
{
  if (partitions < 1 || order < 2 || order % 2 != 0)
  {
    // Only even orders allowed
    return NULL;
  }

  const int stages = 1 + (partitions - 1) *
                           SUNRpowerI(composition_stages, order / 2 - 1);
  const ARKodeSplittingCoefficients coefficients =
    ARKodeSplittingCoefficients_Alloc(1, stages, partitions, order);
  if (coefficients == NULL) { return NULL; }

  arkodeSplittingCoefficients_ComposeStrangHelper(partitions, order,
                                                  composition_stages,
                                                  SUN_RCONST(0.0),
                                                  SUN_RCONST(1.0),
                                                  coefficients->beta[0]);

  return coefficients;
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_TripleJump(
  const int partitions, const int order)
{
  return arkodeSplittingCoefficients_ComposeStrang(partitions, order, 3);
}

ARKodeSplittingCoefficients ARKodeSplittingCoefficients_SuzukiFractal(
  const int partitions, const int order)
{
  return arkodeSplittingCoefficients_ComposeStrang(partitions, order, 5);
}

/*---------------------------------------------------------------
  Routine to print a splitting coefficients structure
  ---------------------------------------------------------------*/
void ARKodeSplittingCoefficients_Write(ARKodeSplittingCoefficients coefficients,
                                       FILE* outfile)
{
  // TODO: update when https://github.com/LLNL/sundials/pull/517 merged
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
  fprintf(outfile, "  alpha = ");
  for (int i = 0; i < coefficients->sequential_methods; i++)
  {
    fprintf(outfile, "%" RSYM "  ", coefficients->alpha[i]);
  }
  fprintf(outfile, "\n");

  for (int i = 0; i < coefficients->sequential_methods; i++)
  {
    fprintf(outfile, "  beta[%i] = \n", i);
    for (int j = 0; j <= coefficients->stages; j++)
    {
      fprintf(outfile, "      ");
      for (int k = 0; k < coefficients->partitions; k++)
      {
        fprintf(outfile, "%" RSYM "  ", coefficients->beta[i][j][k]);
      }
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
  }
}
