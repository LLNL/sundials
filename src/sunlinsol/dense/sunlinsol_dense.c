/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ SMU
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
 * This is the implementation file for the dense implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_math.h>
#include "sundials/sundials_errors.h"
#include "sundials/sundials_types.h"

#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Dense solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define DENSE_CONTENT(S)  ( (SUNLinearSolverContent_Dense)(S->content) )
#define PIVOTS(S)         ( DENSE_CONTENT(S)->pivots )
#define LASTFLAG(S)       ( DENSE_CONTENT(S)->last_flag )

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new dense linear solver
 */

SUNLinearSolver SUNLinSol_Dense(N_Vector y, SUNMatrix A, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_Dense content;
  sunindextype MatrixRows;

  SUNAssertContext(sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, sunctx);
  SUNAssert(SUNDenseMatrix_Rows(A) == SUNDenseMatrix_Columns(A), SUN_ERR_ARG_DIMSMISMATCH, sunctx);
  SUNAssert(y->ops->nvgetarraypointer, SUN_ERR_ARG_INCOMPATIBLE, sunctx);

  MatrixRows = SUNDenseMatrix_Rows(A);
  SUNAssert(MatrixRows == N_VGetLength(y), SUN_ERR_ARG_DIMSMISMATCH, sunctx);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNCheckCallLastErrReturnNull(SUNLinSolNewEmpty(sunctx), sunctx);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_Dense;
  S->ops->getid      = SUNLinSolGetID_Dense;
  S->ops->initialize = SUNLinSolInitialize_Dense;
  S->ops->setup      = SUNLinSolSetup_Dense;
  S->ops->solve      = SUNLinSolSolve_Dense;
  S->ops->lastflag   = SUNLinSolLastFlag_Dense;
  S->ops->space      = SUNLinSolSpace_Dense;
  S->ops->free       = SUNLinSolFree_Dense;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_Dense) malloc(sizeof *content);
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

SUNLinearSolver_Type SUNLinSolGetType_Dense(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_Dense(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DENSE);
}

SUNErrCode SUNLinSolInitialize_Dense(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  LASTFLAG(S) = SUNLS_SUCCESS;
  return SUN_SUCCESS;
}

SUNLsStatus SUNLinSolSetup_Dense(SUNLinearSolver S, SUNMatrix A)
{
  realtype **A_cols;
  sunindextype *pivots;
  SUNContext sunctx = S->sunctx;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, sunctx);

  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  pivots = NULL;
  A_cols = SUNDenseMatrix_Cols(A);
  pivots = PIVOTS(S);
  if ( (A_cols == NULL) || (pivots == NULL) ) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(SUNLS_MEM_FAIL);
  }

  /* perform LU factorization of input matrix */
  LASTFLAG(S) = SUNDlsMat_denseGETRF(A_cols, SUNDenseMatrix_Rows(A),
                                     SUNDenseMatrix_Columns(A), pivots);

  /* store error flag (if nonzero, this row encountered zero-valued pivod) */
  if (LASTFLAG(S) > 0)
    return(SUNLS_LUFACT_FAIL);
  return(SUNLS_SUCCESS);
}

SUNLsStatus SUNLinSolSolve_Dense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                 N_Vector b, realtype tol)
{
  realtype **A_cols, *xdata;
  sunindextype *pivots;
  SUNContext sunctx = S->sunctx;

  /* copy b into x */
  SUNCheckCallLastErr(N_VScale(ONE, b, x), sunctx);

  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  xdata = NULL;
  pivots = NULL;
  A_cols = SUNCheckCallLastErr(SUNDenseMatrix_Cols(A), sunctx);
  xdata = SUNCheckCallLastErr(N_VGetArrayPointer(x), sunctx);
  pivots = PIVOTS(S);
  if ( (A_cols == NULL) || (xdata == NULL)  || (pivots == NULL) ) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(SUNLS_MEM_FAIL);
  }

  /* solve using LU factors */
  SUNDlsMat_denseGETRS(A_cols, SUNDenseMatrix_Rows(A), pivots, xdata);
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(SUNLS_SUCCESS);
}

sunindextype SUNLinSolLastFlag_Dense(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return(LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_Dense(SUNLinearSolver S, long int* lenrwLS,
                                long int* leniwLS)
{
  SUNAssert(SUNLinSolGetID(S) == SUNLINEARSOLVER_DENSE, SUN_ERR_ARG_WRONGTYPE, S->sunctx);
  *leniwLS = 2 + DENSE_CONTENT(S)->N;
  *lenrwLS = 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_Dense(SUNLinearSolver S)
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
