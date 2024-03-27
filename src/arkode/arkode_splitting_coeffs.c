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

#include <arkode/arkode_splitting.h>
#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"

ARKodeSplittingCoeffs ARKodeSplittingCoeffs_Alloc(const int sequential_methods,
                                                  const int stages,
                                                  const int partitions)
{
  if (sequential_methods < 1 || stages < 1 || partitions <= 1) { return NULL; }

  const ARKodeSplittingCoeffs coeffs =
    (ARKodeSplittingCoeffs)malloc(sizeof(*coeffs));
  if (coeffs == NULL) { return NULL; }

  coeffs->sequential_methods = sequential_methods;
  coeffs->stages             = stages;
  coeffs->partitions         = partitions;

  coeffs->alpha = (sunrealtype*)calloc(sequential_methods,
                                       sizeof(*coeffs->alpha));
  if (coeffs->alpha == NULL)
  {
    ARKodeSplittingCoeffs_Free(coeffs);
    return NULL;
  }

  /* In order for beta to be indexed like beta[i][j][k], we can't use a single
     malloc. The j index requires allocating a matrix of pointers into beta. The
     i index requires allocating an array of pointers into that matrix. */

  /* Array of pointers for index i */
  coeffs->beta =
    (sunrealtype***)malloc(sequential_methods * sizeof(*coeffs->beta));
  if (coeffs->beta == NULL)
  {
    ARKodeSplittingCoeffs_Free(coeffs);
    return NULL;
  }

  /* Matrix of pointers for index j */
  sunrealtype** beta_cols =
    (sunrealtype**)malloc(sequential_methods * stages * sizeof(*beta_cols));
  if (beta_cols == NULL)
  {
    ARKodeSplittingCoeffs_Free(coeffs);
    return NULL;
  }

  /* Contiguous memory to store the beta tensor */
  sunrealtype* beta_mem =
    (sunrealtype*)calloc(sequential_methods * stages * partitions,
                         sizeof(*beta_mem));
  if (beta_mem == NULL)
  {
    ARKodeSplittingCoeffs_Free(coeffs);
    return NULL;
  }

  for (int i = 0; i < sequential_methods; i++)
  {
    coeffs->beta[i] = &beta_cols[i];
    for (int j = 0; j < stages; j++)
    {
      coeffs->beta[i][j] = &beta_mem[i * stages + j];
    }
  }

  return coeffs;
}

void ARKodeSplittingCoeffs_Free(ARKodeSplittingCoeffs coeffs)
{
  if (coeffs != NULL)
  {
    free(coeffs->alpha);
    if (coeffs->beta != NULL)
    {
      free(coeffs->beta[0][0]);
      free(coeffs->beta[0]);
      free(coeffs->beta);
    }
  }

  free(coeffs);
}

void ARKodeSplittingCoeffs_Space(ARKodeSplittingCoeffs coeffs,
                                 sunindextype* liw, sunindextype* lrw)
{
  /* initialize outputs and return if coeffs is not allocated */
  *liw = 0;
  *lrw = 0;
  if (coeffs == NULL) { return; }

  /* ints for # of sequential methods, stages, and partitions */
  *liw = 3;

  if (coeffs->beta)
  {
    /* pointers for index i of beta[i][j][k] */
    *liw += coeffs->sequential_methods;
  }

  if (coeffs->beta[0])
  {
    /* pointers for index j of beta[i][j][k] */
    *liw += coeffs->sequential_methods * coeffs->stages;
  }

  if (coeffs->alpha) { *lrw += coeffs->sequential_methods; }
  if (coeffs->beta[0][0])
  {
    *lrw += coeffs->sequential_methods * coeffs->stages * coeffs->partitions;
  }
}

ARKodeSplittingCoeffs ARKodeSplittingCoeffs_Copy(ARKodeSplittingCoeffs coeffs)
{
  if (coeffs == NULL) { return NULL; }

  ARKodeSplittingCoeffs coeffsCopy =
    ARKodeSplittingCoeffs_Alloc(coeffs->sequential_methods, coeffs->stages,
                                coeffs->partitions);
  if (coeffsCopy == NULL) { return NULL; }

  memcpy(coeffsCopy->alpha, coeffs->alpha,
         coeffs->sequential_methods * sizeof(*coeffs->alpha));

  /* beta[0][0] points to the contiguous memory allocation, so we can copy it
     with a single memcpy */
  memcpy(coeffsCopy->beta[0][0], coeffs->beta[0][0],
         coeffs->sequential_methods * coeffs->stages * coeffs->partitions *
           sizeof(*coeffs->beta));

  return coeffsCopy;
}

ARKodeSplittingCoeffs ARKodeSplittingCoeffs_LoadCoeffs(ARKODE_SplittingCoeffsID method)
{
  switch (method)
  {
#define ARK_SPLITTING_COEFFS(name, coeff) \
  case name: coeff break;
#include "arkode_splitting_coeffs.def"
#undef ARK_SPLITTING_COEFFS

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown splitting coefficients");
    return NULL;
  }
}

ARKodeSplittingCoeffs ARKodeSplittingCoeffs_LoadCoeffsByName(const char* method)
{
#define ARK_SPLITTING_COEFFS(name, coeff) \
  if (strcmp(#name, method) == 0) coeff
#include "arkode_splitting_coeffs.def"
#undef ARK_SPLITTING_COEFFS

  arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                  "Unknown splitting coefficients");

  return NULL;
}

ARKodeSplittingCoeffs ARKodeSplittingCoeffs_LieTrotter(const int partitions)
{
  const ARKodeSplittingCoeffs coeffs = ARKodeSplittingCoeffs_Alloc(1, 1,
                                                                   partitions);
  if (coeffs == NULL) { return NULL; }

  coeffs->alpha[0] = SUN_RCONST(1.0);
  for (int i = 0; i < partitions; i++)
  {
    coeffs->beta[0][0][i] = SUN_RCONST(1.0);
  }

  return coeffs;
}

ARKodeSplittingCoeffs ARKodeSplittingCoeffs_Strang(const int partitions)
{
  return ARKodeSplittingCoeffs_TripleJump(partitions, 2);
}

ARKodeSplittingCoeffs ARKodeSplittingCoeffs_Parallel(const int partitions)
{
  const ARKodeSplittingCoeffs coeffs =
    ARKodeSplittingCoeffs_Alloc(partitions + 1, 1, partitions);
  if (coeffs == NULL) { return NULL; }

  for (int i = 0; i < partitions; i++)
  {
    coeffs->alpha[0]      = SUN_RCONST(1.0);
    coeffs->beta[i][0][i] = SUN_RCONST(1.0);
  }

  coeffs->alpha[partitions] = 1 - partitions;

  return coeffs;
}

ARKodeSplittingCoeffs ARKodeSplittingCoeffs_SymmetricParallel(const int partitions)
{
  const ARKodeSplittingCoeffs coeffs =
    ARKodeSplittingCoeffs_Alloc(2, partitions, partitions);
  if (coeffs == NULL) { return NULL; }

  coeffs->alpha[0] = SUN_RCONST(0.5);
  coeffs->alpha[1] = SUN_RCONST(0.5);

  for (int i = 0; i < partitions; i++)
  {
    coeffs->beta[0][0][i]                  = SUN_RCONST(1.0);
    coeffs->beta[1][i][partitions - i - 1] = SUN_RCONST(1.0);
  }

  return coeffs;
}

static sunrealtype** arkodeSplittingCoeffs_TripleJump(const int partitions,
                                                      const int order,
                                                      sunrealtype** beta,
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
  beta = arkodeSplittingCoeffs_TripleJump(partitions, order - 2, beta, z1);
  beta = arkodeSplittingCoeffs_TripleJump(partitions, order - 2, beta, z0);
  beta = arkodeSplittingCoeffs_TripleJump(partitions, order - 2, beta, z1);

  return beta;
}

ARKodeSplittingCoeffs ARKodeSplittingCoeffs_TripleJump(const int partitions,
                                                       const int order)
{
  if (order < 2 || order % 2 != 0)
  {
    // Only even orders allowed
    return NULL;
  }

  const int stages = 1 + (partitions - 1) * SUNRpowerI(3, order / 2 - 1);
  const ARKodeSplittingCoeffs coeffs = ARKodeSplittingCoeffs_Alloc(1, stages,
                                                                   partitions);
  if (coeffs == NULL) { return NULL; }

  coeffs->alpha[0] = SUN_RCONST(1.0);
  arkodeSplittingCoeffs_TripleJump(partitions, order, coeffs->beta[0],
                                   SUN_RCONST(1.0));

  return coeffs;
}

/*---------------------------------------------------------------
  Routine to print a Butcher table structure
  ---------------------------------------------------------------*/
void ARKodeSplittingCoeffs_Write(ARKodeSplittingCoeffs coeffs, FILE* outfile)
{
  if (coeffs == NULL || coeffs->alpha == NULL || coeffs->beta == NULL ||
      coeffs->beta[0] == NULL || coeffs->beta[0][0] == NULL)
  {
    return;
  }

  fprintf(outfile, "  sequential methods = %i\n", coeffs->sequential_methods);
  fprintf(outfile, "  stages = %i\n", coeffs->stages);
  fprintf(outfile, "  partitions = %i\n", coeffs->partitions);
  fprintf(outfile, "  alpha = \n");
  for (int i = 0; i < coeffs->sequential_methods; i++)
  {
    fprintf(outfile, "%" RSYM "  ", coeffs->alpha[i]);
  }

  for (int i = 0; i < coeffs->sequential_methods; i++)
  {
    fprintf(outfile, "  beta[%i] = \n", i);
    for (int j = 0; j < coeffs->stages; j++)
    {
      fprintf(outfile, "      ");
      for (int k = 0; k < coeffs->partitions; k++)
      {
        fprintf(outfile, "%" RSYMW "  ", coeffs->beta[i][j][k]);
      }
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
  }
}
