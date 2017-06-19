/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
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
 * This is the implementation file for the dense implementation of 
 * the SUNMATRIX package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)


/* Private function prototypes */
static booleantype SMCompatible_Dense(SUNMatrix A, SUNMatrix B);
static booleantype SMCompatible2_Dense(SUNMatrix A, N_Vector x, N_Vector y);


/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns the SUNMatrix type ID. Used to identify matrix 
 * implementations from the abstract SUNMatrix interface.
 */

SUNMatrix_ID SUNMatrixGetID_Dense(SUNMatrix A)
{
  return SUNMATRIX_DENSE;
}

/* ----------------------------------------------------------------------------
 * Function to create a new dense matrix
 */

SUNMatrix SUNMatrixNew_Dense(long int M, long int N)
{
  SUNMatrix A;
  SUNMatrix_Ops ops;
  SUNMatrixContent_Dense content;
  long int j;

  /* return with NULL matrix on illegal dimension input */
  if ( (M <= 0) || (N <= 0) ) return(NULL);

  /* Create matrix */
  A = NULL;
  A = (SUNMatrix) malloc(sizeof *A);
  if (A == NULL) return(NULL);
  
  /* Create matrix operation structure */
  ops = NULL;
  ops = (SUNMatrix_Ops) malloc(sizeof(struct _generic_SUNMatrix_Ops));
  if (ops == NULL) { free(A); return(NULL); }

  /* Attach operations */
  ops->getid       = SUNMatrixGetID_Dense;
  ops->clone       = SUNMatrixClone_Dense;
  ops->destroy     = SUNMatrixDestroy_Dense;
  ops->zero        = SUNMatrixZero_Dense;
  ops->scale       = SUNMatrixScale_Dense;
  ops->copy        = SUNMatrixCopy_Dense;
  ops->addidentity = SUNMatrixAddIdentity_Dense;
  ops->add         = SUNMatrixAdd_Dense;
  ops->matvec      = SUNMatrixMatvec_Dense;

  /* Create content */
  content = NULL;
  content = (SUNMatrixContent_Dense) malloc(sizeof(struct _SUNMatrixContent_Dense));
  if (content == NULL) { free(ops); free(A); return(NULL); }

  /* Fill content */
  content->M = M;
  content->N = N;
  content->ldata = M*N;
  content->data = NULL;
  content->data = (realtype *) malloc(M * N * sizeof(realtype));
  if (content->data == NULL) {
    free(content); free(ops); free(A); return(NULL);
  }
  content->cols = NULL;
  content->cols = (realtype **) malloc(N * sizeof(realtype *));
  if (content->cols == NULL) {
    free(content->data); free(content); free(ops); free(A); return(NULL);
  }
  for (j=0; j<N; j++) content->cols[j] = content->data + j * M;
  
  /* Attach content and ops */
  A->content = content;
  A->ops     = ops;

  return(A);
}


/* ----------------------------------------------------------------------------
 * Function to free a matrix created with SUNMatrixNew_Dense
 */

void SUNMatrixDestroy_Dense(SUNMatrix A)
{
  /* should not be called unless A is a dense matrix; 
     otherwise return immediately */
  if (SUNMatrixGetID(A) != SUNMATRIX_DENSE)
    return;

  /* perform operation */
  free(SM_DATA_D(A));  SM_DATA_D(A) = NULL;
  free(SM_CONTENT_D(A)->cols);  SM_CONTENT_D(A)->cols = NULL;
  free(A->content);  A->content = NULL;
  free(A->ops);  A->ops = NULL;
  free(A); A = NULL;
  return;
}

/* ----------------------------------------------------------------------------
 * Function to print the dense matrix 
 */
 
void SUNMatrixPrint_Dense(SUNMatrix A, FILE* outfile)
{
  long int i, j;
  
  /* should not be called unless A is a dense matrix; 
     otherwise return immediately */
  if (SUNMatrixGetID(A) != SUNMATRIX_DENSE)
    return;

  /* perform operation */
  fprintf(outfile,"\n");
  for (i=0; i<SM_ROWS_D(A); i++) {
    for (j=0; j<SM_COLUMNS_D(A); j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(outfile,"%12Lg  ", SM_ELEMENT_D(A,i,j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      fprintf(outfile,"%12g  ", SM_ELEMENT_D(A,i,j));
#else
      fprintf(outfile,"%12g  ", SM_ELEMENT_D(A,i,j));
#endif
    }
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
  return;
}


/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix SUNMatrixClone_Dense(SUNMatrix A)
{
  SUNMatrix B = SUNMatrixNew_Dense(SM_ROWS_D(A), SM_COLUMNS_D(A));
  return(B);
}

int SUNMatrixZero_Dense(SUNMatrix A)
{
  long int i;
  realtype *Adata;

  /* Verify that A is a dense matrix */
  if (SUNMatrixGetID(A) != SUNMATRIX_DENSE)
    return 1;

  /* Perform operation */
  Adata = SM_DATA_D(A);
  for (i=0; i<SM_LDATA_D(A); i++)
    Adata[i] = ZERO;
  return 0;
}

int SUNMatrixScale_Dense(realtype c, SUNMatrix A)
{
  long int i;
  realtype *Adata;

  /* Verify that A is a dense matrix */
  if (SUNMatrixGetID(A) != SUNMATRIX_DENSE)
    return 1;

  /* Perform operation */
  Adata = SM_DATA_D(A);
  for (i=0; i<SM_LDATA_D(A); i++)
    Adata[i] *= c;
  return 0;
}

int SUNMatrixCopy_Dense(SUNMatrix A, SUNMatrix B)
{
  long int i, j;

  /* Verify that A and B are compatible */
  if (!SMCompatible_Dense(A, B))
    return 1;

  /* Perform operation */
  for (j=0; j<SM_COLUMNS_D(A); j++)
    for (i=0; i<SM_ROWS_D(A); i++)
      SM_ELEMENT_D(A,i,j) = SM_ELEMENT_D(B,i,j);
  return 0;
}

int SUNMatrixAddIdentity_Dense(SUNMatrix A)
{
  long int i;

  /* Verify that A is a dense matrix */
  if (SUNMatrixGetID(A) != SUNMATRIX_DENSE)
    return 1;

  /* Perform operation */
  for (i=0; i<SM_COLUMNS_D(A); i++)
    SM_ELEMENT_D(A,i,i) += ONE;
  return 0;
}

int SUNMatrixAdd_Dense(SUNMatrix A, SUNMatrix B)
{
  long int i, j;

  /* Verify that A and B are compatible */
  if (!SMCompatible_Dense(A, B))
    return 1;

  /* Perform operation */
  for (j=0; j<SM_COLUMNS_D(A); j++)
    for (i=0; i<SM_ROWS_D(A); i++)
      SM_ELEMENT_D(A,i,j) += SM_ELEMENT_D(B,i,j);
  return 0;
}

int SUNMatrixMatvec_Dense(SUNMatrix A, N_Vector x, N_Vector y)
{
  long int i, j;
  realtype *col_j, *xd, *yd;
  
  /* Verify that A, x and y are compatible */
  if (!SMCompatible2_Dense(A, x, y))
    return 1;

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  yd = N_VGetArrayPointer(y);
  if ((xd == NULL) || (yd == NULL))
    return 1;

  /* Perform operation */
  for (i=0; i<SM_ROWS_D(A); i++)
    yd[i] = ZERO;
  for(j=0; j<SM_COLUMNS_D(A); j++) {
    col_j = SM_COLUMN_D(A,j);
    for (i=0; i<SM_ROWS_D(A); i++)
      yd[i] += col_j[i]*xd[j];
  }
  return 0;
}


/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static booleantype SMCompatible_Dense(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must be SUNMATRIX_DENSE */
  if (SUNMatrixGetID(A) != SUNMATRIX_DENSE)
    return FALSE;
  if (SUNMatrixGetID(B) != SUNMATRIX_DENSE)
    return FALSE;

  /* both matrices must have the same shape */
  if (SM_ROWS_D(A) != SM_ROWS_D(B))
    return FALSE;
  if (SM_COLUMNS_D(A) != SM_COLUMNS_D(B))
    return FALSE;

  return TRUE;
}


static booleantype SMCompatible2_Dense(SUNMatrix A, N_Vector x, N_Vector y)
{
  /* matrix must be SUNMATRIX_DENSE */
  if (SUNMatrixGetID(A) != SUNMATRIX_DENSE)
    return FALSE;

  /* vectors must be one of {SERIAL, OPENMP, PTHREADS}, and 
     have compatible dimensions */ 
  if (N_VGetVectorID(x) == SUNDIALS_NVEC_SERIAL) {
    if (NV_LENGTH_S(x) != SM_COLUMNS_D(A))
      return FALSE;
  } else if (N_VGetVectorID(x) == SUNDIALS_NVEC_OPENMP) {
    if (NV_LENGTH_OMP(x) != SM_COLUMNS_D(A))
      return FALSE;
  } else if (N_VGetVectorID(x) == SUNDIALS_NVEC_PTHREADS) {
    if (NV_LENGTH_PT(x) != SM_COLUMNS_D(A))
      return FALSE;
  } else {   /* incompatible type */
    return FALSE;
  }
  if (N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL) {
    if (NV_LENGTH_S(y) != SM_ROWS_D(A))
      return FALSE;
  } else if (N_VGetVectorID(y) == SUNDIALS_NVEC_OPENMP) {
    if (NV_LENGTH_OMP(y) != SM_ROWS_D(A))
      return FALSE;
  } else if (N_VGetVectorID(y) == SUNDIALS_NVEC_PTHREADS) {
    if (NV_LENGTH_PT(y) != SM_ROWS_D(A))
      return FALSE;
  } else {   /* incompatible type */
    return FALSE;
  }

  return TRUE;
}

