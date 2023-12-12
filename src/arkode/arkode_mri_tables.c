/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * This is the implementation file for ARKODE's MRIStepCoupling tables.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"

/* ===========================================================================
 * Exported Functions
 * ===========================================================================*/

/*---------------------------------------------------------------
  Returns MRIStepCoupling table structure for pre-set MRI methods.

  Input:  imeth -- integer key for the desired method
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_LoadTable(ARKODE_MRITableID method)
{
  switch (method)
  {
#define ARK_MRI_TABLE(name, coeff) \
  case name: coeff break;
#include "arkode_mri_tables.def"
#undef ARK_MRI_TABLE

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown coupling table");
    return NULL;
  }
}

/*---------------------------------------------------------------
  Returns MRIStepCoupling table structure for pre-set MRI methods.

  Input:  method -- string key for the desired method
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_LoadTableByName(const char* method)
{
#define ARK_MRI_TABLE(name, coeff) \
  if (strcmp(#name, method) == 0) coeff
#include "arkode_mri_tables.def"
#undef ARK_MRI_TABLE

  arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                  "Unknown coupling table");

  return NULL;
}

/*---------------------------------------------------------------
  Routine to allocate an empty MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages,
                                      MRISTEP_METHOD_TYPE type)
{
  int i, j;
  MRIStepCoupling MRIC = NULL;

  /* Check for legal input values */
  if (nmat < 1 || stages < 1) { return (NULL); }

  /* ------------------------------------------
   * Allocate and initialize coupling structure
   * ------------------------------------------ */

  MRIC = (MRIStepCoupling)malloc(sizeof(struct MRIStepCouplingMem));
  if (!MRIC) { return (NULL); }

  MRIC->nmat   = nmat;
  MRIC->stages = stages;
  MRIC->q      = 0;
  MRIC->p      = 0;
  MRIC->c      = NULL;
  MRIC->W      = NULL;
  MRIC->G      = NULL;

  /* --------------------------------------------
   * Allocate abscissae and coupling coefficients
   * -------------------------------------------- */

  MRIC->c = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
  if (!(MRIC->c))
  {
    MRIStepCoupling_Free(MRIC);
    return (NULL);
  }

  if (type == MRISTEP_EXPLICIT || type == MRISTEP_IMEX)
  {
    /* allocate W matrices */
    MRIC->W = (sunrealtype***)calloc(nmat, sizeof(sunrealtype**));
    if (!(MRIC->W))
    {
      MRIStepCoupling_Free(MRIC);
      return (NULL);
    }

    /* allocate rows of each matrix in W */
    for (i = 0; i < nmat; i++)
    {
      MRIC->W[i] = NULL;
      MRIC->W[i] = (sunrealtype**)calloc(stages, sizeof(sunrealtype*));
      if (!(MRIC->W[i]))
      {
        MRIStepCoupling_Free(MRIC);
        return (NULL);
      }
    }

    /* allocate columns of each matrix in W */
    for (i = 0; i < nmat; i++)
    {
      for (j = 0; j < stages; j++)
      {
        MRIC->W[i][j] = NULL;
        MRIC->W[i][j] = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
        if (!(MRIC->W[i][j]))
        {
          MRIStepCoupling_Free(MRIC);
          return (NULL);
        }
      }
    }
  }

  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
  {
    /* allocate G matrices */
    MRIC->G = (sunrealtype***)calloc(nmat, sizeof(sunrealtype**));
    if (!(MRIC->G))
    {
      MRIStepCoupling_Free(MRIC);
      return (NULL);
    }

    /* allocate rows of each matrix in G */
    for (i = 0; i < nmat; i++)
    {
      MRIC->G[i] = NULL;
      MRIC->G[i] = (sunrealtype**)calloc(stages, sizeof(sunrealtype*));
      if (!(MRIC->G[i]))
      {
        MRIStepCoupling_Free(MRIC);
        return (NULL);
      }
    }

    /* allocate columns of each matrix in G */
    for (i = 0; i < nmat; i++)
    {
      for (j = 0; j < stages; j++)
      {
        MRIC->G[i][j] = NULL;
        MRIC->G[i][j] = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
        if (!(MRIC->G[i][j]))
        {
          MRIStepCoupling_Free(MRIC);
          return (NULL);
        }
      }
    }
  }

  return (MRIC);
}

/*---------------------------------------------------------------
  Routine to allocate and fill a MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Create(int nmat, int stages, int q, int p,
                                       sunrealtype* W, sunrealtype* G,
                                       sunrealtype* c)
{
  int i, j, k;
  MRISTEP_METHOD_TYPE type;
  MRIStepCoupling MRIC = NULL;

  /* Check for legal inputs */
  if (nmat < 1 || stages < 1 || !c) { return (NULL); }

  /* Check for method coefficients and set method type */
  if (W && G) { type = MRISTEP_IMEX; }
  else if (W && !G) { type = MRISTEP_EXPLICIT; }
  else if (!W && G) { type = MRISTEP_IMPLICIT; }
  else { return (NULL); }

  /* Allocate MRIStepCoupling structure */
  MRIC = MRIStepCoupling_Alloc(nmat, stages, type);
  if (!MRIC) { return (NULL); }

  /* -------------------------
   * Copy the inputs into MRIC
   * ------------------------- */

  /* Method and embedding order */
  MRIC->q = q;
  MRIC->p = p;

  /* Abscissae */
  for (i = 0; i < stages; i++) { MRIC->c[i] = c[i]; }

  /* Coupling coefficients stored as 1D arrays of length nmat * stages * stages,
     with each stages * stages matrix stored in C (row-major) order */
  if (type == MRISTEP_EXPLICIT || type == MRISTEP_IMEX)
  {
    for (k = 0; k < nmat; k++)
    {
      for (i = 0; i < stages; i++)
      {
        for (j = 0; j < stages; j++)
        {
          MRIC->W[k][i][j] = W[stages * (stages * k + i) + j];
        }
      }
    }
  }
  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
  {
    for (k = 0; k < nmat; k++)
    {
      for (i = 0; i < stages; i++)
      {
        for (j = 0; j < stages; j++)
        {
          MRIC->G[k][i][j] = G[stages * (stages * k + i) + j];
        }
      }
    }
  }

  return (MRIC);
}

/*---------------------------------------------------------------
  Construct the MRI coupling matrix for an MIS method based on
  a given 'slow' Butcher table.
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B, int q, int p)
{
  int i, j, stages;
  sunbooleantype padding;
  sunrealtype Asum;
  sunrealtype*** C;
  MRISTEP_METHOD_TYPE type;
  MRIStepCoupling MRIC;

  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* Check that input table is non-NULL */
  if (!B) { return (NULL); }

  /* -----------------------------------
   * Check that the input table is valid
   * ----------------------------------- */

  /* First stage is just old solution */
  Asum = SUNRabs(B->c[0]);
  for (j = 0; j < B->stages; j++) { Asum += SUNRabs(B->A[0][j]); }
  if (Asum > tol) { return (NULL); }

  /* Last stage exceeds 1 */
  if (B->c[B->stages - 1] > ONE + tol) { return (NULL); }

  /* All stages are sorted */
  for (j = 1; j < B->stages; j++)
  {
    if ((B->c[j] - B->c[j - 1]) < -tol) { return (NULL); }
  }

  /* Each stage at most diagonally implicit */
  Asum = ZERO;
  for (i = 0; i < B->stages; i++)
  {
    for (j = i + 1; j < B->stages; j++) { Asum += SUNRabs(B->A[i][j]); }
  }
  if (Asum > tol) { return (NULL); }

  /* -----------------------------------------
   * determine whether the table needs padding
   * ----------------------------------------- */

  padding = SUNFALSE;

  /* Last stage time should equal 1 */
  if (SUNRabs(B->c[B->stages - 1] - ONE) > tol) { padding = SUNTRUE; }

  /* Last row of A should equal b */
  for (j = 0; j < B->stages; j++)
  {
    if (SUNRabs(B->A[B->stages - 1][j] - B->b[j]) > tol) { padding = SUNTRUE; }
  }
  stages = (padding) ? B->stages + 1 : B->stages;

  /* -------------------------
   * determine the method type
   * ------------------------- */

  /* Check if the table is strictly lower triangular (explicit) */
  type = MRISTEP_EXPLICIT;

  for (i = 0; i < B->stages; i++)
  {
    for (j = i; j < B->stages; j++)
    {
      if (SUNRabs(B->A[i][j]) > tol) { type = MRISTEP_IMPLICIT; }
    }
  }

  /* ----------------------------
   * construct coupling structure
   * ---------------------------- */

  MRIC = MRIStepCoupling_Alloc(1, stages, type);
  if (!MRIC) { return (NULL); }

  /* Copy method/embedding orders */
  MRIC->q = q;
  MRIC->p = p;

  /* Copy abscissae, padding if needed */
  for (i = 0; i < B->stages; i++) { MRIC->c[i] = B->c[i]; }

  if (padding) { MRIC->c[stages - 1] = ONE; }

  /* Construct the coupling table */
  if (type == MRISTEP_EXPLICIT) { C = MRIC->W; }
  else { C = MRIC->G; }

  /* First row is identically zero */
  for (i = 0; i < stages; i++)
  {
    for (j = 0; j < stages; j++) { C[0][i][j] = ZERO; }
  }

  /* Remaining rows = A(2:end,:) - A(1:end-1,:) */
  for (i = 1; i < B->stages; i++)
  {
    for (j = 0; j < B->stages; j++)
    {
      C[0][i][j] = B->A[i][j] - B->A[i - 1][j];
    }
  }

  /* Padded row = b(:) - A(end,:) */
  if (padding)
  {
    for (j = 0; j < B->stages; j++)
    {
      C[0][stages - 1][j] = B->b[j] - B->A[B->stages - 1][j];
    }
  }

  return (MRIC);
}

/*---------------------------------------------------------------
  Routine to copy a MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling MRIC)
{
  int i, j, k, nmat, stages;
  MRISTEP_METHOD_TYPE type;
  MRIStepCoupling MRICcopy;

  /* Check for legal input */
  if (!MRIC) { return (NULL); }

  /* Check for method coefficients and set method type */
  if (MRIC->W && MRIC->G) { type = MRISTEP_IMEX; }
  else if (MRIC->W && !(MRIC->G)) { type = MRISTEP_EXPLICIT; }
  else if (!(MRIC->W) && MRIC->G) { type = MRISTEP_IMPLICIT; }
  else { return (NULL); }

  /* Check for stage times */
  if (!(MRIC->c)) { return (NULL); }

  /* Get the number of coupling matrices and stages */
  nmat   = MRIC->nmat;
  stages = MRIC->stages;

  /* Allocate coupling structure */
  MRICcopy = MRIStepCoupling_Alloc(nmat, stages, type);
  if (!MRICcopy) { return (NULL); }

  /* Copy method and embedding orders */
  MRICcopy->q = MRIC->q;
  MRICcopy->p = MRIC->p;

  /* Copy abscissae */
  for (i = 0; i < stages; i++) { MRICcopy->c[i] = MRIC->c[i]; }

  /* Copy explicit coupling matrices W */
  if (MRIC->W)
  {
    for (k = 0; k < nmat; k++)
    {
      for (i = 0; i < stages; i++)
      {
        for (j = 0; j < stages; j++)
        {
          MRICcopy->W[k][i][j] = MRIC->W[k][i][j];
        }
      }
    }
  }

  /* Copy implicit coupling matrices G */
  if (MRIC->G)
  {
    for (k = 0; k < nmat; k++)
    {
      for (i = 0; i < stages; i++)
      {
        for (j = 0; j < stages; j++)
        {
          MRICcopy->G[k][i][j] = MRIC->G[k][i][j];
        }
      }
    }
  }

  return (MRICcopy);
}

/*---------------------------------------------------------------
  Routine to query the MRIStepCoupling structure workspace size
  ---------------------------------------------------------------*/
void MRIStepCoupling_Space(MRIStepCoupling MRIC, sunindextype* liw,
                           sunindextype* lrw)
{
  /* initialize outputs and return if MRIC is not allocated */
  *liw = 0;
  *lrw = 0;
  if (!MRIC) { return; }

  /* fill outputs based on MRIC */
  *liw = 4;
  if (MRIC->c) { *lrw += MRIC->stages; }
  if (MRIC->W) { *lrw += MRIC->nmat * MRIC->stages * MRIC->stages; }
  if (MRIC->G) { *lrw += MRIC->nmat * MRIC->stages * MRIC->stages; }
}

/*---------------------------------------------------------------
  Routine to free a MRIStepCoupling structure
  ---------------------------------------------------------------*/
void MRIStepCoupling_Free(MRIStepCoupling MRIC)
{
  int k, i;

  /* Free each field within MRIStepCoupling structure, and then
     free structure itself */
  if (MRIC)
  {
    if (MRIC->c) { free(MRIC->c); }

    if (MRIC->W)
    {
      for (k = 0; k < MRIC->nmat; k++)
      {
        if (MRIC->W[k])
        {
          for (i = 0; i < MRIC->stages; i++)
          {
            if (MRIC->W[k][i])
            {
              free(MRIC->W[k][i]);
              MRIC->W[k][i] = NULL;
            }
          }
          free(MRIC->W[k]);
          MRIC->W[k] = NULL;
        }
      }
      free(MRIC->W);
    }

    if (MRIC->G)
    {
      for (k = 0; k < MRIC->nmat; k++)
      {
        if (MRIC->G[k])
        {
          for (i = 0; i < MRIC->stages; i++)
          {
            if (MRIC->G[k][i])
            {
              free(MRIC->G[k][i]);
              MRIC->G[k][i] = NULL;
            }
          }
          free(MRIC->G[k]);
          MRIC->G[k] = NULL;
        }
      }
      free(MRIC->G);
    }

    free(MRIC);
  }
}

/*---------------------------------------------------------------
  Routine to print a MRIStepCoupling structure
  ---------------------------------------------------------------*/
void MRIStepCoupling_Write(MRIStepCoupling MRIC, FILE* outfile)
{
  int i, j, k;

  /* check for vaild coupling structure */
  if (!MRIC) { return; }
  if (!(MRIC->W) && !(MRIC->G)) { return; }
  if (!(MRIC->c)) { return; }

  if (MRIC->W)
  {
    for (i = 0; i < MRIC->nmat; i++)
    {
      if (!(MRIC->W[i])) { return; }
      for (j = 0; j < MRIC->stages; j++)
      {
        if (!(MRIC->W[i][j])) { return; }
      }
    }
  }

  if (MRIC->G)
  {
    for (i = 0; i < MRIC->nmat; i++)
    {
      if (!(MRIC->G[i])) { return; }
      for (j = 0; j < MRIC->stages; j++)
      {
        if (!(MRIC->G[i][j])) { return; }
      }
    }
  }

  fprintf(outfile, "  nmat = %i\n", MRIC->nmat);
  fprintf(outfile, "  stages = %i\n", MRIC->stages);
  fprintf(outfile, "  method order (q) = %i\n", MRIC->q);
  fprintf(outfile, "  embedding order (p) = %i\n", MRIC->p);
  fprintf(outfile, "  c = ");
  for (i = 0; i < MRIC->stages; i++)
  {
    fprintf(outfile, "%" RSYM "  ", MRIC->c[i]);
  }
  fprintf(outfile, "\n");

  if (MRIC->W)
  {
    for (k = 0; k < MRIC->nmat; k++)
    {
      fprintf(outfile, "  W[%i] = \n", k);
      for (i = 0; i < MRIC->stages; i++)
      {
        fprintf(outfile, "      ");
        for (j = 0; j < MRIC->stages; j++)
        {
          fprintf(outfile, "%" RSYMW "  ", MRIC->W[k][i][j]);
        }
        fprintf(outfile, "\n");
      }
      fprintf(outfile, "\n");
    }
  }

  if (MRIC->G)
  {
    for (k = 0; k < MRIC->nmat; k++)
    {
      fprintf(outfile, "  G[%i] = \n", k);
      for (i = 0; i < MRIC->stages; i++)
      {
        fprintf(outfile, "      ");
        for (j = 0; j < MRIC->stages; j++)
        {
          fprintf(outfile, "%" RSYMW "  ", MRIC->G[k][i][j]);
        }
        fprintf(outfile, "\n");
      }
      fprintf(outfile, "\n");
    }
  }
}

/* ===========================================================================
 * Private Functions
 * ===========================================================================*/

/* ---------------------------------------------------------------------------
 * Stage type identifier: returns one of the constants
 *
 * MRISTAGE_ERK_FAST    -- standard MIS-like stage
 * MRISTAGE_ERK_NOFAST  -- standard ERK stage
 * MRISTAGE_DIRK_NOFAST -- standard DIRK stage
 * MRISTAGE_DIRK_FAST   -- coupled DIRK + MIS-like stage
 *
 * for each nontrivial stage in an MRI-like method. Otherwise (i.e., stage is
 * not in [1,MRIC->stages-1]), returns ARK_INVALID_TABLE (<0).
 *
 * The stage type is determined by 2 factors:
 * (a) Sum |MRIC->G[:][is][is]| (nonzero => DIRK)
 * (b) MRIC->c[is] - MRIC->c[is-1]  (nonzero => fast)
 * ---------------------------------------------------------------------------*/

int mriStepCoupling_GetStageType(MRIStepCoupling MRIC, int is)
{
  int i;
  sunrealtype Gabs, cdiff;
  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  if ((is < 1) || (is >= MRIC->stages)) { return ARK_INVALID_TABLE; }

  /* sum of stage diagonal entries across implicit tables */
  Gabs = ZERO;
  if (MRIC->G)
  {
    for (i = 0; i < MRIC->nmat; i++) { Gabs += SUNRabs(MRIC->G[i][is][is]); }
  }

  /* abscissae difference */
  cdiff = MRIC->c[is] - MRIC->c[is - 1];

  if (Gabs > tol)
  { /* DIRK */
    if (cdiff > tol)
    { /* Fast */
      return (MRISTAGE_DIRK_FAST);
    }
    else { return (MRISTAGE_DIRK_NOFAST); }
  }
  else
  { /* ERK */
    if (cdiff > tol)
    { /* Fast */
      return (MRISTAGE_ERK_FAST);
    }
    else { return (MRISTAGE_ERK_NOFAST); }
  }
}

/* ---------------------------------------------------------------------------
 * Computes the stage RHS vector storage maps. With repeated abscissae the
 * first stage of the pair generally corresponds to a column of zeros and so
 * does not need to be computed and stored. The stage_map indicates if the RHS
 * needs to be computed and where to store it i.e., stage_map[i] > -1.
 * ---------------------------------------------------------------------------*/

int mriStepCoupling_GetStageMap(MRIStepCoupling MRIC, int* stage_map,
                                int* nstages_active)
{
  int i, j, k, idx;
  sunrealtype Wsum, Gsum;
  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* ----------------------
   * Check for valid inputs
   * ---------------------- */

  if (!MRIC) { return (ARK_ILL_INPUT); }
  if (!(MRIC->W) && !(MRIC->G)) { return (ARK_ILL_INPUT); }
  if (!stage_map || !nstages_active) { return (ARK_ILL_INPUT); }

  /* -------------------
   * Compute storage map
   * ------------------- */

  /* Number of stage RHS vectors active */
  *nstages_active = 0;

  /* Initial storage index */
  idx = 0;

  /* Check if a stage corresponds to a column of zeros for all coupling
   * matrices by computing the column sums */
  for (j = 0; j < MRIC->stages; j++)
  {
    Wsum = ZERO;
    Gsum = ZERO;

    if (MRIC->W)
    {
      for (k = 0; k < MRIC->nmat; k++)
      {
        for (i = 0; i < MRIC->stages; i++)
        {
          Wsum += SUNRabs(MRIC->W[k][i][j]);
        }
      }
    }

    if (MRIC->G)
    {
      for (k = 0; k < MRIC->nmat; k++)
      {
        for (i = 0; i < MRIC->stages; i++)
        {
          Gsum += SUNRabs(MRIC->G[k][i][j]);
        }
      }
    }

    if (Wsum > tol || Gsum > tol)
    {
      stage_map[j] = idx;
      idx++;
    }
    else { stage_map[j] = -1; }
  }

  /* Check and set number of stage RHS vectors active */
  if (idx < 1) { return (ARK_ILL_INPUT); }

  *nstages_active = idx;

  return (ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
