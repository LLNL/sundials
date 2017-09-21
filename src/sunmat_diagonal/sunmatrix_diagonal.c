/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the diagonal implementation of 
 * the SUNMATRIX package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunmatrix/sunmatrix_diagonal.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)


/* Private function prototypes */
static booleantype SMCompatible_Diagonal(SUNMatrix A, SUNMatrix B);
static booleantype SMCompatible2_Diagonal(SUNMatrix A, N_Vector x, N_Vector y);


/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new diagonal matrix
 */

SUNMatrix SUNDiagonalMatrix(N_Vector tmpl)
{
  SUNMatrix A;
  SUNMatrix_Ops ops;
  SUNMatrixContent_Diagonal content;

  /* return with NULL matrix on illegal input */
  if ( tmpl==NULL ) return(NULL);

  /* return with NULL if template vector is missing any require operations */
  if (tmpl->ops->nvclone == NULL)       return(NULL);
  if (tmpl->ops->nvgetvectorid == NULL) return(NULL);
  if (tmpl->ops->nvconst == NULL)       return(NULL);
  if (tmpl->ops->nvscale == NULL)       return(NULL);
  if (tmpl->ops->nvaddconst == NULL)    return(NULL);
  if (tmpl->ops->nvlinearsum == NULL)   return(NULL);
  if (tmpl->ops->nvprod == NULL)        return(NULL);
  
  /* Create matrix */
  A = NULL;
  A = (SUNMatrix) malloc(sizeof *A);
  if (A == NULL) return(NULL);
  
  /* Create matrix operation structure */
  ops = NULL;
  ops = (SUNMatrix_Ops) malloc(sizeof(struct _generic_SUNMatrix_Ops));
  if (ops == NULL) { free(A); return(NULL); }

  /* Attach operations */
  ops->getid     = SUNMatGetID_Diagonal;
  ops->clone     = SUNMatClone_Diagonal;
  ops->destroy   = SUNMatDestroy_Diagonal;
  ops->zero      = SUNMatZero_Diagonal;
  ops->copy      = SUNMatCopy_Diagonal;
  ops->scaleadd  = SUNMatScaleAdd_Diagonal;
  ops->scaleaddi = SUNMatScaleAddI_Diagonal;
  ops->matvec    = SUNMatMatvec_Diagonal;
  ops->space     = SUNMatSpace_Diagonal;

  /* Create content */
  content = NULL;
  content = (SUNMatrixContent_Diagonal)
    malloc(sizeof(struct _SUNMatrixContent_Diagonal));
  if (content == NULL) { free(ops); free(A); return(NULL); }

  /* Fill content */
  content->d = NULL;
  content->d = N_VClone(tmpl);
  if (content->d == NULL) {
    free(content); free(ops); free(A); return(NULL);
  }
  
  /* Attach content and ops */
  A->content = content;
  A->ops     = ops;

  return(A);
}


/* ----------------------------------------------------------------------------
 * Functions to access the contents of the diagonal matrix structure
 */

N_Vector SUNDiagonalMatrix_Diag(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_DIAGONAL)
    return SM_DATA_DIAG(A);
  else
    return NULL;
}


/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_Diagonal(SUNMatrix A)
{
  return SUNMATRIX_DIAGONAL;
}

SUNMatrix SUNMatClone_Diagonal(SUNMatrix A)
{
  SUNMatrix B = SUNDiagonalMatrix(SM_DATA_DIAG(A));
  return(B);
}

void SUNMatDestroy_Diagonal(SUNMatrix A)
{
  /* perform operation */
  if (SM_DATA_DIAG(A))
    N_VDestroy(SM_DATA_DIAG(A));
  SM_DATA_DIAG(A) = NULL;
  free(A->content);  A->content = NULL;
  free(A->ops);  A->ops = NULL;
  free(A); A = NULL;
  return;
}

int SUNMatZero_Diagonal(SUNMatrix A)
{
  /* Perform operation */
  N_VConst(ZERO, SM_DATA_DIAG(A));
  return 0;
}

int SUNMatCopy_Diagonal(SUNMatrix A, SUNMatrix B)
{
  /* Verify that A and B are compatible */
  if (!SMCompatible_Diagonal(A, B))
    return 1;

  /* Perform operation */
  N_VScale(ONE, SM_DATA_DIAG(A), SM_DATA_DIAG(B));
  return 0;
}

int SUNMatScaleAddI_Diagonal(realtype c, SUNMatrix A)
{
  /* Perform operation */
  N_VScale(c, SM_DATA_DIAG(A), SM_DATA_DIAG(A));
  N_VAddConst(SM_DATA_DIAG(A), ONE, SM_DATA_DIAG(A));
  return 0;
}

int SUNMatScaleAdd_Diagonal(realtype c, SUNMatrix A, SUNMatrix B)
{
  /* Verify that A and B are compatible */
  if (!SMCompatible_Diagonal(A, B))
    return 1;

  /* Perform operation */
  N_VLinearSum(c, SM_DATA_DIAG(A), ONE, SM_DATA_DIAG(B), SM_DATA_DIAG(A));
  return 0;
}

int SUNMatMatvec_Diagonal(SUNMatrix A, N_Vector x, N_Vector y)
{
  /* Verify that A, x and y are compatible */
  if (!SMCompatible2_Diagonal(A, x, y))
    return 1;

  /* Perform operation */
  N_VProd(SM_DATA_DIAG(A), x, y);
  return 0;
}

int SUNMatSpace_Diagonal(SUNMatrix A, long int *lenrw, long int *leniw)
{
  sunindextype lrw1, liw1;
  N_VSpace(SM_DATA_DIAG(A), &lrw1, &liw1);
  *lenrw = lrw1;
  *leniw = liw1;
  return 0;
}


/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static booleantype SMCompatible_Diagonal(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must be SUNMATRIX_DIAGONAL */
  if (SUNMatGetID(A) != SUNMATRIX_DIAGONAL)
    return FALSE;
  if (SUNMatGetID(B) != SUNMATRIX_DIAGONAL)
    return FALSE;

  /* both matrix diagonals must have the same N_Vector ID */
  if (N_VGetVectorID(SM_DATA_DIAG(A)) != N_VGetVectorID(SM_DATA_DIAG(B)))
    return FALSE;

  /* it would be nice to verify that the vectors all have the same dimensions, 
     but there is no generic 'length' N_Vector routine */
  
  return TRUE;
}


static booleantype SMCompatible2_Diagonal(SUNMatrix A, N_Vector x, N_Vector y)
{
  /* the matrix diagonal and vectors must have the same N_Vector ID */ 
  if ( (N_VGetVectorID(SM_DATA_DIAG(A)) != N_VGetVectorID(x)) ||
       (N_VGetVectorID(SM_DATA_DIAG(A)) != N_VGetVectorID(y)) )
    return FALSE;

  /* it would be nice to verify that the vectors all have the same dimensions, 
     but there is no generic 'length' N_Vector routine */
  
  return TRUE;
}

