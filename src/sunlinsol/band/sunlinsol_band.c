/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
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
 * This is the implementation file for the band implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_band.h>
#include <sundials/sundials_math.h>
#include "sundials/sundials_context.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_linearsolver.h"

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define ROW(i,j,smu) (i-j+smu)

/*
 * -----------------------------------------------------------------
 * Band solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define BAND_CONTENT(S)   ( (SUNLinearSolverContent_Band)(S->content) )
#define PIVOTS(S)         ( BAND_CONTENT(S)->pivots )
#define LASTFLAG(S)       ( BAND_CONTENT(S)->last_flag )

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new band linear solver
 */

SUNLinearSolver SUNLinSol_Band(N_Vector y, SUNMatrix A, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_Band content;
  sunindextype MatrixRows;

  SUNAssertContext(sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE, sunctx);
  SUNAssert(SUNBandMatrix_Rows(A) == SUNBandMatrix_Columns(A), SUN_ERR_ARG_DIMSMISMATCH, sunctx);
  SUNAssert(y->ops->nvgetarraypointer, SUN_ERR_ARG_INCOMPATIBLE, sunctx);

  /* Check that A has appropriate storage upper bandwidth for factorization */
  MatrixRows = SUNBandMatrix_Rows(A);
  SUNAssert(SUNBandMatrix_StoredUpperBandwidth(A) >=
              SUNMIN(MatrixRows - 1, SUNBandMatrix_LowerBandwidth(A) +
                                       SUNBandMatrix_UpperBandwidth(A)),
            SUN_ERR_ARG_INCOMPATIBLE, sunctx);
  SUNAssert(MatrixRows == N_VGetLength(y), SUN_ERR_ARG_DIMSMISMATCH, sunctx);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNCheckCallLastErrReturnNull(SUNLinSolNewEmpty(sunctx), sunctx);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_Band;
  S->ops->getid      = SUNLinSolGetID_Band;
  S->ops->initialize = SUNLinSolInitialize_Band;
  S->ops->setup      = SUNLinSolSetup_Band;
  S->ops->solve      = SUNLinSolSolve_Band;
  S->ops->lastflag   = SUNLinSolLastFlag_Band;
  S->ops->space      = SUNLinSolSpace_Band;
  S->ops->free       = SUNLinSolFree_Band;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_Band) malloc(sizeof *content);
  SUNAssert(content, SUN_ERR_MALLOC_FAIL, sunctx);

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->N         = MatrixRows;
  content->last_flag = 0;
  content->pivots    = NULL;

  /* Allocate content */
  content->pivots = (sunindextype *) malloc(MatrixRows * sizeof(sunindextype));
  SUNAssert(content->pivots, SUN_ERR_MALLOC_FAIL, sunctx);

  return(S);
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_Band(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_Band(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_BAND);
}

SUNErrCode SUNLinSolInitialize_Band(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  LASTFLAG(S) = SUNLS_SUCCESS;
  return SUN_SUCCESS;
}

SUNLsStatus SUNLinSolSetup_Band(SUNLinearSolver S, SUNMatrix A)
{
  realtype **A_cols;
  sunindextype *pivots;
  SUNContext sunctx = S->sunctx;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE, sunctx);

  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  pivots = NULL;
  A_cols = SM_COLS_B(A);
  pivots = PIVOTS(S);
  if ( (A_cols == NULL) || (pivots == NULL) ) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(SUNLS_MEM_FAIL);
  }

  /* ensure that storage upper bandwidth is sufficient for fill-in */
  if (SM_SUBAND_B(A) < SUNMIN(SM_COLUMNS_B(A)-1, SM_UBAND_B(A) + SM_LBAND_B(A))) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(SUNLS_MEM_FAIL);
  }

  /* perform LU factorization of input matrix */
  LASTFLAG(S) = SUNDlsMat_bandGBTRF(A_cols, SM_COLUMNS_B(A), SM_UBAND_B(A),
                                    SM_LBAND_B(A), SM_SUBAND_B(A), pivots);

  /* store error flag (if nonzero, that row encountered zero-valued pivod) */
  if (LASTFLAG(S) > 0)
    return(SUNLS_LUFACT_FAIL);
  return(SUNLS_SUCCESS);
}

SUNLsStatus SUNLinSolSolve_Band(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                N_Vector b, realtype tol)
{
  realtype **A_cols, *xdata;
  sunindextype *pivots;

  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  xdata = NULL;
  pivots = NULL;
  A_cols = SUNBandMatrix_Cols(A);
  xdata = N_VGetArrayPointer(x);
  pivots = PIVOTS(S);
  if ( (A_cols == NULL) || (xdata == NULL)  || (pivots == NULL) ) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(SUNLS_MEM_FAIL);
  }

  /* solve using LU factors */
  SUNDlsMat_bandGBTRS(A_cols, SM_COLUMNS_B(A), SM_SUBAND_B(A),
                      SM_LBAND_B(A), pivots, xdata);
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(SUNLS_SUCCESS);
}

sunindextype SUNLinSolLastFlag_Band(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return(LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_Band(SUNLinearSolver S, long int* lenrwLS,
                               long int* leniwLS)
{
  SUNAssert(SUNLinSolGetID(S) == SUNLINEARSOLVER_BAND, SUN_ERR_ARG_WRONGTYPE, S->sunctx);
  *leniwLS = 2 + BAND_CONTENT(S)->N;
  *lenrwLS = 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_Band(SUNLinearSolver S)
{
  /* return if S is already free */
  if (S == NULL) return SUN_SUCCESS;

  /* delete items from contents, then delete generic structure */
  if (S->content) {
    if (PIVOTS(S)) {
      free(PIVOTS(S));
      PIVOTS(S) = NULL;
    }
    free(S->content);
    S->content = NULL;
  }
  if (S->ops) {
    free(S->ops);
    S->ops = NULL;
  }
  free(S); S = NULL;
  return SUN_SUCCESS;
}
