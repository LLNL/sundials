/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_dense.c by: Scott D. Cohen,
 *     Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * the SUNMATRIX package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "sundials/sundials_context.h"
#include "sundials/sundials_errors.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Private function prototypes */
static booleantype compatibleMatrices(SUNMatrix A, SUNMatrix B);
static booleantype compatibleMatrixAndVectors(SUNMatrix A, N_Vector x, N_Vector y);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new dense matrix
 */

SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N, SUNContext sunctx)
{
  SUNMatrix A;
  SUNMatrixContent_Dense content;
  sunindextype j;

  /* return with NULL matrix on illegal dimension input */
  SUNAssertContext(sunctx);
  SUNAssert(N > 0 && M > 0, SUN_ERR_ARG_OUTOFRANGE, sunctx);

  /* Create an empty matrix object */
  A = NULL;
  A = SUNCheckCallLastErrNull(SUNMatNewEmpty(sunctx), sunctx);

  /* Attach operations */
  A->ops->getid     = SUNMatGetID_Dense;
  A->ops->clone     = SUNMatClone_Dense;
  A->ops->destroy   = SUNMatDestroy_Dense;
  A->ops->zero      = SUNMatZero_Dense;
  A->ops->copy      = SUNMatCopy_Dense;
  A->ops->scaleadd  = SUNMatScaleAdd_Dense;
  A->ops->scaleaddi = SUNMatScaleAddI_Dense;
  A->ops->matvec    = SUNMatMatvec_Dense;
  A->ops->space     = SUNMatSpace_Dense;

  /* Create content */
  content = NULL;
  content = (SUNMatrixContent_Dense)malloc(sizeof *content);
  SUNAssert(content, SUN_ERR_MALLOC_FAIL, sunctx);

  /* Attach content */
  A->content = content;

  /* Fill content */
  content->M     = M;
  content->N     = N;
  content->ldata = M * N;
  content->data  = NULL;
  content->cols  = NULL;

  /* Allocate content */
  content->data = (realtype*)calloc(M * N, sizeof(realtype));
  SUNAssert(content->data, SUN_ERR_MALLOC_FAIL, sunctx);

  content->cols = (realtype**)malloc(N * sizeof(realtype*));
  SUNAssert(content->cols, SUN_ERR_MALLOC_FAIL, sunctx);
  for (j = 0; j < N; j++)
    content->cols[j] = content->data + j * M;

  return (A);
}

/* ----------------------------------------------------------------------------
 * Function to print the dense matrix
 */

void SUNDenseMatrix_Print(SUNMatrix A, FILE* outfile)
{
  sunindextype i, j;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);

  /* perform operation */
  fprintf(outfile, "\n");
  for (i = 0; i < SM_ROWS_D(A); i++) {
    for (j = 0; j < SM_COLUMNS_D(A); j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(outfile, "%12Lg  ", SM_ELEMENT_D(A, i, j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      fprintf(outfile, "%12g  ", SM_ELEMENT_D(A, i, j));
#else
      fprintf(outfile, "%12g  ", SM_ELEMENT_D(A, i, j));
#endif
    }
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "\n");
  return;
}

/* ----------------------------------------------------------------------------
 * Functions to access the contents of the dense matrix structure
 */

sunindextype SUNDenseMatrix_Rows(SUNMatrix A)
{
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);
  return SM_ROWS_D(A);
}

sunindextype SUNDenseMatrix_Columns(SUNMatrix A)
{
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);
  return SM_COLUMNS_D(A);
}

sunindextype SUNDenseMatrix_LData(SUNMatrix A)
{
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);
  return SM_LDATA_D(A);
}

realtype* SUNDenseMatrix_Data(SUNMatrix A)
{
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);
  return SM_DATA_D(A);
}

realtype** SUNDenseMatrix_Cols(SUNMatrix A)
{
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);
  return SM_COLS_D(A);
}

realtype* SUNDenseMatrix_Column(SUNMatrix A, sunindextype j)
{
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);
  return SM_COLUMN_D(A, j);
}

/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_Dense(SUNMatrix A) { return SUNMATRIX_DENSE; }

SUNMatrix SUNMatClone_Dense(SUNMatrix A)
{
  SUNContext sunctx = A->sunctx;
  SUNMatrix B =
    SUNCheckCallLastErrNull(SUNDenseMatrix(SM_ROWS_D(A), SM_COLUMNS_D(A),
                                                 sunctx),
                                  sunctx);
  return (B);
}

void SUNMatDestroy_Dense(SUNMatrix A) SUNDIALS_NOEXCEPT
{
  if (A == NULL)
    return;

  /* free content */
  if (A->content != NULL) {
    /* free data array */
    if (SM_DATA_D(A) != NULL) {
      free(SM_DATA_D(A));
      SM_DATA_D(A) = NULL;
    }
    /* free column pointers */
    if (SM_CONTENT_D(A)->cols != NULL) {
      free(SM_CONTENT_D(A)->cols);
      SM_CONTENT_D(A)->cols = NULL;
    }
    /* free content struct */
    free(A->content);
    A->content = NULL;
  }

  /* free ops and matrix */
  if (A->ops) {
    free(A->ops);
    A->ops = NULL;
  }
  free(A);
  A = NULL;

  return;
}

SUNErrCode SUNMatZero_Dense(SUNMatrix A)
{
  sunindextype i;
  realtype* Adata;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);

  /* Perform operation A_ij = 0 */
  Adata = SM_DATA_D(A);
  for (i = 0; i < SM_LDATA_D(A); i++)
    Adata[i] = ZERO;

  return SUN_SUCCESS;
}

SUNErrCode SUNMatCopy_Dense(SUNMatrix A, SUNMatrix B)
{
  sunindextype i, j;
  SUNContext sunctx = A->sunctx;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, sunctx);
  SUNAssert(SUNMatGetID(B) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, sunctx);
  SUNCheck(compatibleMatrices(A, B), SUN_ERR_ARG_DIMSMISMATCH, sunctx);

  /* Perform operation B_ij = A_ij */
  for (j = 0; j < SM_COLUMNS_D(A); j++)
    for (i = 0; i < SM_ROWS_D(A); i++)
      SM_ELEMENT_D(B, i, j) = SM_ELEMENT_D(A, i, j);

  return SUN_SUCCESS;
}

SUNErrCode SUNMatScaleAddI_Dense(realtype c, SUNMatrix A)
{
  sunindextype i, j;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);

  /* Perform operation A = c*A + I */
  for (j = 0; j < SM_COLUMNS_D(A); j++) {
    for (i = 0; i < SM_ROWS_D(A); i++) {
      SM_ELEMENT_D(A, i, j) *= c;
      if (i == j)
        SM_ELEMENT_D(A, i, j) += ONE;
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatScaleAdd_Dense(realtype c, SUNMatrix A, SUNMatrix B)
{
  sunindextype i, j;
  SUNContext sunctx = A->sunctx;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, sunctx);
  SUNCheck(compatibleMatrices(A, B), SUN_ERR_ARG_DIMSMISMATCH, sunctx);

  /* Perform operation A = c*A + B */
  for (j = 0; j < SM_COLUMNS_D(A); j++)
    for (i = 0; i < SM_ROWS_D(A); i++)
      SM_ELEMENT_D(A, i, j) = c * SM_ELEMENT_D(A, i, j) + SM_ELEMENT_D(B, i, j);

  return SUN_SUCCESS;
}

SUNErrCode SUNMatMatvec_Dense(SUNMatrix A, N_Vector x, N_Vector y)
{
  sunindextype i, j;
  realtype *col_j, *xd, *yd;
  SUNContext sunctx = A->sunctx;

  SUNCheck(compatibleMatrixAndVectors(A, x, y), SUN_ERR_ARG_DIMSMISMATCH, sunctx);

  /* access vector data (return if NULL data pointers) */
  xd = SUNCheckCallLastErr(N_VGetArrayPointer(x), sunctx);
  yd = SUNCheckCallLastErr(N_VGetArrayPointer(y), sunctx);

  /* Perform operation y = Ax */
  for (i = 0; i < SM_ROWS_D(A); i++)
    yd[i] = ZERO;
  for (j = 0; j < SM_COLUMNS_D(A); j++) {
    col_j = SM_COLUMN_D(A, j);
    for (i = 0; i < SM_ROWS_D(A); i++)
      yd[i] += col_j[i] * xd[j];
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMatSpace_Dense(SUNMatrix A, long int* lenrw, long int* leniw)
{
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE, A->sunctx);
  *lenrw = SM_LDATA_D(A);
  *leniw = 3 + SM_COLUMNS_D(A);
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static booleantype compatibleMatrices(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must have the same shape */
  if ((SM_ROWS_D(A) != SM_ROWS_D(B)) || (SM_COLUMNS_D(A) != SM_COLUMNS_D(B))) {
    return SUNFALSE;
  }

  return SUNTRUE;
}

static booleantype compatibleMatrixAndVectors(SUNMatrix A, N_Vector x, N_Vector y)
{
  /* Vectors must provide nvgetarraypointer and cannot be a parallel vector */
  if (!x->ops->nvgetarraypointer || !y->ops->nvgetarraypointer) {
    return SUNFALSE;
  }

  /* Check that the dimensions agree */
  if ((N_VGetLength(x) != SM_COLUMNS_D(A)) || (N_VGetLength(y) != SM_ROWS_D(A))) {
    return SUNFALSE;
  }

  return SUNTRUE;
}
