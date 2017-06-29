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
 * This is the implementation file for the sparse implementation of 
 * the SUNMATRIX package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunmatrix/sunmatrix_sparse.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Private function prototypes */
static booleantype SMCompatible_Sparse(SUNMatrix A, SUNMatrix B);
static booleantype SMCompatible2_Sparse(SUNMatrix A, N_Vector x, N_Vector y);
int Matvec_SparseCSC(SUNMatrix A, N_Vector x, N_Vector y);
int Matvec_SparseCSR(SUNMatrix A, N_Vector x, N_Vector y);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/*
 * ==================================================================
 * Private function prototypes (functions working on SlsMat)
 * ==================================================================
 */

/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix
 */

SUNMatrix SUNSparseMatrix(long int M, long int N,
                          long int NNZ, int sparsetype)
{
  SUNMatrix A;
  SUNMatrix_Ops ops;
  SUNMatrixContent_Sparse content;

  /* return with NULL matrix on illegal input */
  if ( (M <= 0) || (N <= 0) || (NNZ < 0) ) return(NULL);
  if ( (sparsetype != CSC_MAT) && (sparsetype != CSR_MAT) ) return(NULL);

  /* Create matrix */
  A = NULL;
  A = (SUNMatrix) malloc(sizeof *A);
  if (A == NULL) return(NULL);
  
  /* Create matrix operation structure */
  ops = NULL;
  ops = (SUNMatrix_Ops) malloc(sizeof(struct _generic_SUNMatrix_Ops));
  if (ops == NULL) { free(A); return(NULL); }

  /* Attach operations */
  ops->getid       = SUNMatGetID_Sparse;
  ops->clone       = SUNMatClone_Sparse;
  ops->destroy     = SUNMatDestroy_Sparse;
  ops->zero        = SUNMatZero_Sparse;
  ops->copy        = SUNMatCopy_Sparse;
  ops->scaleadd    = SUNMatScaleAdd_Sparse;
  ops->scaleaddi   = SUNMatScaleAddI_Sparse;
  ops->matvec      = SUNMatMatvec_Sparse;

  /* Create content */
  content = NULL;
  content = (SUNMatrixContent_Sparse) malloc(sizeof(struct _SUNMatrixContent_Sparse));
  if (content == NULL) { free(ops); free(A); return(NULL); }

  /* Fill content */
  content->sparsetype = sparsetype;
  content->M = M;
  content->N = N;
  content->NNZ = NNZ;
  switch(sparsetype){
  case CSC_MAT:
    content->NP = N;
    content->rowvals = &(content->indexvals);
    content->colptrs = &(content->indexptrs);
    /* CSR indices */
    content->colvals = NULL;
    content->rowptrs = NULL;
    break;
  case CSR_MAT:
    content->NP = M;
    content->colvals = &(content->indexvals);
    content->rowptrs = &(content->indexptrs);
    /* CSC indices */
    content->rowvals = NULL;
    content->colptrs = NULL;
  }
  content->data = (realtype *) calloc(NNZ, sizeof(realtype));
  if (content->data == NULL) {
    free(content); free(ops); free(A); return(NULL);
  }
  content->indexvals = (long int *) calloc(NNZ, sizeof(long int));
  if (content->indexvals == NULL) {
    free(content->data); free(content); free(ops); free(A); return(NULL);
  }
  content->indexptrs = (long int *) calloc((content->NP + 1), sizeof(long int));
  if (content->indexptrs == NULL) {
    free(content->indexvals);
    free(content->data);
    free(content);
    free(ops);
    free(A);
    return(NULL);
  }
  content->indexptrs[content->NP] = 0;

  /* Attach content and ops */
  A->content = content;
  A->ops     = ops;

  return(A);
}



/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix from an existing dense matrix 
 * by copying all nonzero values into the sparse matrix structure.  Returns NULL 
 * if the request for matrix storage cannot be satisfied.
 */

SUNMatrix SUNSparseFromDenseMatrix(SUNMatrix Ad, realtype droptol, int sparsetype)
{
  long int i, j, nnz;
  long int M, N;
  SUNMatrix As;

  /* check for legal sparsetype, droptol and input matrix type */
  if ( (sparsetype != CSR_MAT) && (sparsetype != CSC_MAT) )
    return NULL;
  if ( droptol < ZERO )
    return NULL;
  if (SUNMatGetID(Ad) != SUNMATRIX_DENSE)
    return NULL;
  
  /* set size of new matrix */
  M = SM_ROWS_D(Ad);
  N = SM_COLUMNS_D(Ad);

  /* determine total number of nonzeros */
  nnz = 0;
  for (j=0; j<N; j++)
    for (i=0; i<M; i++)
      nnz += (SUNRabs(SM_ELEMENT_D(Ad,i,j)) > droptol);
    
  /* allocate sparse matrix */
  As = SUNSparseMatrix(M, N, nnz, sparsetype);
  if (As == NULL)  return NULL;
  
  /* copy nonzeros from Ad into As, based on CSR/CSC type */
  nnz = 0;
  if (sparsetype == CSC_MAT) {
    for (j=0; j<N; j++) {
      (SM_INDEXPTRS_S(As))[j] = nnz;
      for (i=0; i<M; i++) {
        if ( SUNRabs(SM_ELEMENT_D(Ad,i,j)) > droptol ) { 
          (SM_INDEXVALS_S(As))[nnz] = i;
          (SM_DATA_S(As))[nnz++] = SM_ELEMENT_D(Ad,i,j);
        }
      }
    }
    (SM_INDEXPTRS_S(As))[N] = nnz;
  } else {       /* CSR_MAT */
    for (i=0; i<M; i++) {
      (SM_INDEXPTRS_S(As))[i] = nnz;
      for (j=0; j<N; j++) {
        if ( SUNRabs(SM_ELEMENT_D(Ad,i,j)) > droptol ) { 
          (SM_INDEXVALS_S(As))[nnz] = j;
          (SM_DATA_S(As))[nnz++] = SM_ELEMENT_D(Ad,i,j);
        }
      }
    }
    (SM_INDEXPTRS_S(As))[M] = nnz;
  }
    
  return(As);
}


/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix from an existing band matrix 
 * by copying all nonzero values into the sparse matrix structure.  Returns NULL 
 * if the request for matrix storage cannot be satisfied.
 */

SUNMatrix SUNSparseFromBandMatrix(SUNMatrix Ad, realtype droptol, int sparsetype)
{
  long int i, j, nnz;
  long int M, N;
  SUNMatrix As;

  /* check for legal sparsetype, droptol and input matrix type */
  if ( (sparsetype != CSR_MAT) && (sparsetype != CSC_MAT) )
    return NULL;
  if ( droptol < ZERO )
    return NULL;
  if (SUNMatGetID(Ad) != SUNMATRIX_BAND)
    return NULL;
  
  /* set size of new matrix */
  M = SM_ROWS_B(Ad);
  N = SM_COLUMNS_B(Ad);

  /* determine total number of nonzeros */
  nnz = 0;
  for (j=0; j<N; j++)
    for (i=SUNMAX(0,j-SM_UBAND_B(Ad)); i<=SUNMIN(M-1,j+SM_LBAND_B(Ad)); i++)
      nnz += (SUNRabs(SM_ELEMENT_B(Ad,i,j)) > droptol);

  /* allocate sparse matrix */
  As = SUNSparseMatrix(M, N, nnz, sparsetype);
  if (As == NULL)  return NULL;

  /* copy nonzeros from Ad into As, based on CSR/CSC type */
  nnz = 0;
  if (sparsetype == CSC_MAT) {
    for (j=0; j<N; j++) {
      (SM_INDEXPTRS_S(As))[j] = nnz;
      for (i=SUNMAX(0,j-SM_UBAND_B(Ad)); i<=SUNMIN(M-1,j+SM_LBAND_B(Ad)); i++) {
        if ( SUNRabs(SM_ELEMENT_B(Ad,i,j)) > droptol ) { 
          (SM_INDEXVALS_S(As))[nnz] = i;
          (SM_DATA_S(As))[nnz++] = SM_ELEMENT_B(Ad,i,j);
        }
      }
    }
    (SM_INDEXPTRS_S(As))[N] = nnz;
  } else {       /* CSR_MAT */
    for (i=0; i<M; i++) {
      (SM_INDEXPTRS_S(As))[i] = nnz;
      for (j=SUNMAX(0,i-SM_LBAND_B(Ad)); j<=SUNMIN(N-1,i+SM_UBAND_B(Ad)); j++) {
        if ( SUNRabs(SM_ELEMENT_B(Ad,i,j)) > droptol ) { 
          (SM_INDEXVALS_S(As))[nnz] = j;
          (SM_DATA_S(As))[nnz++] = SM_ELEMENT_B(Ad,i,j);
        }
      }
    }
    (SM_INDEXPTRS_S(As))[M] = nnz;
  }

  return(As);
}


/* ----------------------------------------------------------------------------
 * Function to reallocate internal sparse matrix storage arrays so that the
 * resulting sparse matrix holds indexptrs[NP] nonzeros.  Returns 0 on success
 * and 1 on failure (e.g. if A does not have sparse type, or if nnz is negative)
 */

int SUNSparseMatrix_Realloc(SUNMatrix A)
{
  /* check for valid matrix type */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE)
    return 1;

  /* get total number of nonzeros (return with failure if illegal) */
  long int nzmax; 
  nzmax = (SM_INDEXPTRS_S(A))[SM_NP_S(A)];
  if (nzmax < 0)
    return 1;

  /* perform reallocation */
  SM_INDEXVALS_S(A) = realloc(SM_INDEXVALS_S(A), nzmax*sizeof(long int));
  SM_DATA_S(A) = realloc(SM_DATA_S(A), nzmax*sizeof(realtype));
  SM_NNZ_S(A) = nzmax;

  return 0;
}


/* ----------------------------------------------------------------------------
 * Function to print the sparse matrix 
 */
 
void SUNSparseMatrix_Print(SUNMatrix A, FILE* outfile)
{
  long int i, j;
  char *matrixtype;
  char *indexname;
  
  /* should not be called unless A is a sparse matrix; 
     otherwise return immediately */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE)
    return;

  /* perform operation */
  switch(SM_SPARSETYPE_S(A)) {
  case CSC_MAT:
    indexname = (char*) "col";
    matrixtype = (char*) "CSC";
    break;
  case CSR_MAT:
    indexname = (char*) "row";
    matrixtype = (char*) "CSR";
    break;
  }
  fprintf(outfile, "\n");
  fprintf(outfile, "%ld by %ld %s matrix, NNZ: %ld \n", SM_ROWS_S(A),
          SM_COLUMNS_S(A), matrixtype, SM_NNZ_S(A));
  for (j=0; j<SM_NP_S(A); j++) {
    fprintf(outfile, "%s %ld : locations %ld to %ld\n", indexname,
            j, (SM_INDEXPTRS_S(A))[j], (SM_INDEXPTRS_S(A))[j+1]-1);
    fprintf(outfile, "  ");
    for (i=(SM_INDEXPTRS_S(A))[j]; i<(SM_INDEXPTRS_S(A))[j+1]; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(outfile, "%ld: %Lg   ", (SM_INDEXVALS_S(A))[i], (SM_DATA_S(A))[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      fprintf(outfile, "%ld: %g   ", (SM_INDEXVALS_S(A))[i], (SM_DATA_S(A))[i]);
#else
      fprintf(outfile, "%ld: %g   ", (SM_INDEXVALS_S(A))[i], (SM_DATA_S(A))[i]);
#endif
    }
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "\n");
  return;
}


/* ----------------------------------------------------------------------------
 * Functions to access the contents of the sparse matrix structure
 */

long int SUNSparseMatrix_Rows(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
    return SM_ROWS_S(A);
  else
    return -1;
}

long int SUNSparseMatrix_Columns(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
    return SM_COLUMNS_S(A);
  else
    return -1;
}

long int SUNSparseMatrix_NNZ(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
    return SM_NNZ_S(A);
  else
    return -1;
}

long int SUNSparseMatrix_NP(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
    return SM_NP_S(A);
  else
    return -1;
}

int SUNSparseMatrix_SparseType(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
    return SM_SPARSETYPE_S(A);
  else
    return -1;
}

realtype* SUNSparseMatrix_Data(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
    return SM_DATA_S(A);
  else
    return NULL;
}

long int* SUNSparseMatrix_IndexValues(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
    return SM_INDEXVALS_S(A);
  else
    return NULL;
}

long int* SUNSparseMatrix_IndexPointers(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
    return SM_INDEXPTRS_S(A);
  else
    return NULL;
}


/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_Sparse(SUNMatrix A)
{
  return SUNMATRIX_SPARSE;
}

SUNMatrix SUNMatClone_Sparse(SUNMatrix A)
{
  SUNMatrix B = SUNSparseMatrix(SM_ROWS_S(A), SM_COLUMNS_S(A),
                                SM_NNZ_S(A), SM_SPARSETYPE_S(A));
  return(B);
}

void SUNMatDestroy_Sparse(SUNMatrix A)
{
  /* perform operation */
  if (SM_DATA_S(A)) {
    free(SM_DATA_S(A));  SM_DATA_S(A) = NULL;
  }
  if (SM_INDEXVALS_S(A)) {
    free(SM_INDEXVALS_S(A));
    SM_INDEXVALS_S(A) = NULL;
    SM_CONTENT_S(A)->rowvals = NULL;
    SM_CONTENT_S(A)->colvals = NULL;
  }
  if (SM_INDEXPTRS_S(A)) {
    free(SM_INDEXPTRS_S(A));
    SM_INDEXPTRS_S(A) = NULL;
    SM_CONTENT_S(A)->colptrs = NULL;
    SM_CONTENT_S(A)->rowptrs = NULL;
  }
  free(A->content); A->content = NULL;
  free(A->ops);  A->ops = NULL;
  free(A); A = NULL;
  return;
}
  
int SUNMatZero_Sparse(SUNMatrix A)
{
  long int i;

  /* Perform operation */
  for (i=0; i<SM_NNZ_S(A); i++) {
    (SM_DATA_S(A))[i] = ZERO;
    (SM_INDEXVALS_S(A))[i] = 0;
  }
  for (i=0; i<SM_NP_S(A); i++) 
    (SM_INDEXPTRS_S(A))[i] = 0;
  (SM_INDEXPTRS_S(A))[SM_NP_S(A)] = 0;
  return 0;
}

int SUNMatCopy_Sparse(SUNMatrix B, SUNMatrix A)
{
  long int i, A_nz;

  /* Verify that A and B are compatible */
  if (!SMCompatible_Sparse(A, B))
    return 1;

  /* Perform operation */
  A_nz = (SM_INDEXPTRS_S(A))[SM_NP_S(A)];
  
  /* ensure that B is allocated with at least as 
     much memory as we have nonzeros in A */
  if (SM_NNZ_S(B) < A_nz) {
    SM_INDEXVALS_S(B) = realloc(SM_INDEXVALS_S(B), A_nz*sizeof(int));
    SM_DATA_S(B) = realloc(SM_DATA_S(B), A_nz*sizeof(realtype));
    SM_NNZ_S(B) = A_nz;
  }

  /* zero out B so that copy works correctly */
  SUNMatZero_Sparse(B);

  /* copy the data and row indices over */
  for (i=0; i<A_nz; i++){
    (SM_DATA_S(B))[i] = (SM_DATA_S(A))[i];
    (SM_INDEXVALS_S(B))[i] = (SM_INDEXVALS_S(A))[i];
  }

  /* copy the column pointers over */
  for (i=0; i<SM_NP_S(A); i++) {
    (SM_INDEXPTRS_S(B))[i] = (SM_INDEXPTRS_S(A))[i];
  }
  (SM_INDEXPTRS_S(B))[SM_NP_S(A)] = A_nz;
  
  return 0;
}

int SUNMatScaleAddI_Sparse(realtype c, SUNMatrix A)
{
  long int j, i, p, nz;
  booleantype newmat, found;
  long int *w, *Ap, *Ai, *Cp, *Ci;
  realtype *x, *Ax, *Cx;
  SUNMatrix C;
  long int M;
  long int N;

  /* Perform operation */

  /* determine if A already contains values on the diagonal (hence 
     no memory allocation necessary) */
  newmat=FALSE;
  for (j=0; j < SUNMIN(SM_COLUMNS_S(A),SM_ROWS_S(A)); j++) {
    /* scan column (row if CSR) of A, searching for diagonal value */
    found = FALSE;
    for (i=(SM_INDEXPTRS_S(A))[j]; i<(SM_INDEXPTRS_S(A))[j+1]; i++) {
      if ((SM_INDEXVALS_S(A))[i] == j) {
        found = TRUE;
        break;
      }
    }
    /* if no diagonal found, signal new matrix */
    if (!found) {
      newmat=TRUE;
      break;
    }
  }

  /* perform operation */

  /*   case 1: A already contains a diagonal */
  if (!newmat) {

    /* iterate through columns, adding 1.0 to diagonal */
    for (j=0; j < SUNMIN(SM_COLUMNS_S(A),SM_ROWS_S(A)); j++)
      for (i=(SM_INDEXPTRS_S(A))[j]; i<(SM_INDEXPTRS_S(A))[j+1]; i++)
        if ((SM_INDEXVALS_S(A))[i] == j) {
          (SM_DATA_S(A))[i] = ONE + c*(SM_DATA_S(A))[i];
        } else {
          (SM_DATA_S(A))[i] = c*(SM_DATA_S(A))[i];
        }

  /*   case 2: A does not already contain a diagonal */
  } else {
    
    if (SM_SPARSETYPE_S(A) == CSC_MAT) {
      M = SM_ROWS_S(A);
      N = SM_COLUMNS_S(A);
    }
    else {
      M = SM_COLUMNS_S(A);
      N = SM_ROWS_S(A);
    }
  
    /* create work arrays for row indices and nonzero column values */
    w = (long int *) malloc(SM_ROWS_S(A) * sizeof(long int));
    x = (realtype *) malloc(SM_ROWS_S(A) * sizeof(realtype));

    /* create new matrix for sum (overestimate nnz as sum of each) */
    C = SUNSparseMatrix(SM_ROWS_S(A), SM_COLUMNS_S(A),
                        (SM_INDEXPTRS_S(A))[SM_NP_S(A)]
                        + SUNMIN(SM_ROWS_S(A), SM_COLUMNS_S(A)),
                        SM_SPARSETYPE_S(A));

    /* access data from CSR structures (return if failure) */
    Cp = Ci = Ap = Ai = NULL;
    Cx = Ax = NULL;
    if (SM_INDEXPTRS_S(C))  Cp = SM_INDEXPTRS_S(C);
    else  return (-1);
    if (SM_INDEXVALS_S(C))  Ci = SM_INDEXVALS_S(C);
    else  return (-1);
    if (SM_DATA_S(C))       Cx = SM_DATA_S(C);
    else  return (-1);
    if (SM_INDEXPTRS_S(A))  Ap = SM_INDEXPTRS_S(A);
    else  return (-1);
    if (SM_INDEXVALS_S(A))  Ai = SM_INDEXVALS_S(A);
    else  return (-1);
    if (SM_DATA_S(A))       Ax = SM_DATA_S(A);
    else  return (-1);

    /* initialize total nonzero count */
    nz = 0;

    /* iterate through columns (rows for CSR) */
    for (j=0; j<N; j++) {

      /* set current column (row) pointer to current # nonzeros */
      Cp[j] = nz;

      /* clear out temporary arrays for this column (row) */
      for (i=0; i<M; i++) {
        w[i] = 0;
        x[i] = 0.0;
      }

      /* iterate down column (along row) of A, collecting nonzeros */
      for (p=Ap[j]; p<Ap[j+1]; p++) {
        w[Ai[p]] += 1;         /* indicate that row is filled */
        x[Ai[p]] = c*Ax[p];    /* collect/scale value */
      }

      /* add identity to this column (row) */
      if (j < M) {
        w[j] += 1;     /* indicate that row is filled */
        x[j] += ONE;   /* update value */
      }

      /* fill entries of C with this column's (row's) data */
      for (i=0; i<M; i++) {
        if ( w[i] > 0 ) { 
          Ci[nz] = i;  
          Cx[nz++] = x[i];
        }
      }
    }

    /* indicate end of data */
    Cp[N] = nz;

    /* update A's structure with C's values; nullify C's pointers */
    SM_NNZ_S(A) = SM_NNZ_S(C);

    if (SM_DATA_S(A))
      free(SM_DATA_S(A));  
    SM_DATA_S(A) = SM_DATA_S(C);
    SM_DATA_S(C) = NULL;

    if (SM_INDEXVALS_S(A))
      free(SM_INDEXVALS_S(A));
    SM_INDEXVALS_S(A) = SM_INDEXVALS_S(C);
    SM_INDEXVALS_S(C) = NULL;

    if (SM_INDEXPTRS_S(A))
      free(SM_INDEXPTRS_S(A));
    SM_INDEXPTRS_S(A) = SM_INDEXPTRS_S(C);
    SM_INDEXPTRS_S(C) = NULL;

    /* clean up */
    SUNMatDestroy_Sparse(C); 
    free(w);
    free(x);

    /* reallocate the new matrix to remove extra space */
    SUNSparseMatrix_Realloc(A);
  }
  return 0;

}

int SUNMatScaleAdd_Sparse(realtype c, SUNMatrix A, SUNMatrix B)
{
  long int j, i, p, nz;
  booleantype newmat;
  long int *w, *Ap, *Ai, *Bp, *Bi, *Cp, *Ci;
  realtype *x, *Ax, *Bx, *Cx;
  SUNMatrix C;
  long int M;
  long int N;

  /* Verify that A and B are compatible */
  if (!SMCompatible_Sparse(A, B))
    return 1;

  /* Perform operation */

  /* if A is CSR matrix, transpose M and N */
  if (SM_SPARSETYPE_S(A) == CSC_MAT) {
    M = SM_ROWS_S(A);
    N = SM_COLUMNS_S(A);
  } else {
    M = SM_COLUMNS_S(A);
    N = SM_ROWS_S(A);
  }
  
  /* create work arrays for row indices and nonzero column values */
  w = (long int *) malloc(M * sizeof(long int));
  x = (realtype *) malloc(M * sizeof(realtype));

  /* determine if A already contains the sparsity pattern of B */
  newmat=FALSE;
  for (j=0; j<N; j++) {

    /* clear work array */
    for (i=0; i<M; i++)  w[i] = 0;

    /* scan column of A, incrementing w by one */
    for (i=(SM_INDEXPTRS_S(A))[j]; i<(SM_INDEXPTRS_S(A))[j+1]; i++)
      w[(SM_INDEXVALS_S(A))[i]] += 1;

    /* scan column of B, decrementing w by one */
    for (i=(SM_INDEXPTRS_S(B))[j]; i<(SM_INDEXPTRS_S(B))[j+1]; i++)
      w[(SM_INDEXVALS_S(B))[i]] -= 1;

    /* if any entry of w is negative, A doesn't contain B's sparsity */
    for (i=0; i<M; i++)
      if (w[i] < 0) {
        newmat = TRUE;
        break;
      }
    if (newmat) break;

  }

  /* perform operation */

  /*   case 1: A already contains sparsity pattern of B */
  if (!newmat) {

    /* iterate through columns, adding matrices */
    for (j=0; j<N; j++) {

      /* clear work array */
      for (i=0; i<M; i++)
        x[i] = ZERO;

      /* scan column of B, updating work array */
      for (i = (SM_INDEXPTRS_S(B))[j]; i < (SM_INDEXPTRS_S(B))[j+1]; i++)
        x[(SM_INDEXVALS_S(B))[i]] = (SM_DATA_S(B))[i];

      /* scan column of A, updating entries appropriately array */
      for (i = (SM_INDEXPTRS_S(A))[j]; i < (SM_INDEXPTRS_S(A))[j+1]; i++)
        (SM_DATA_S(A))[i] = c*(SM_DATA_S(A))[i] + x[(SM_INDEXVALS_S(A))[i]];

    }

  /*   case 2: A does not already contain B's sparsity */
  } else {

    /* create new matrix for sum (overestimate nnz as sum of each) */
    C = SUNSparseMatrix(M, N, (SM_INDEXPTRS_S(A))[N]+(SM_INDEXPTRS_S(B))[N],
                        SM_SPARSETYPE_S(A));

    /* access data from CSR structures (return if failure) */
    Cp = Ci = Ap = Ai = Bp = Bi = NULL;
    Cx = Ax = Bx = NULL;
    if (SM_INDEXPTRS_S(C))  Cp = SM_INDEXPTRS_S(C);
    else  return(-1);
    if (SM_INDEXVALS_S(C))  Ci = SM_INDEXVALS_S(C);
    else  return(-1);
    if (SM_DATA_S(C))       Cx = SM_DATA_S(C);
    else  return(-1);
    if (SM_INDEXPTRS_S(A))  Ap = SM_INDEXPTRS_S(A);
    else  return(-1);
    if (SM_INDEXVALS_S(A))  Ai = SM_INDEXVALS_S(A);
    else  return(-1);
    if (SM_DATA_S(A))       Ax = SM_DATA_S(A);
    else  return(-1);
    if (SM_INDEXPTRS_S(B))  Bp = SM_INDEXPTRS_S(B);
    else  return(-1);
    if (SM_INDEXVALS_S(B))  Bi = SM_INDEXVALS_S(B);
    else  return(-1);
    if (SM_DATA_S(B))       Bx = SM_DATA_S(B);
    else  return(-1);

    /* initialize total nonzero count */
    nz = 0;

    /* iterate through columns */
    for (j=0; j<N; j++) {

      /* set current column pointer to current # nonzeros */
      Cp[j] = nz;

      /* clear out temporary arrays for this column */
      for (i=0; i<M; i++) {
        w[i] = 0;
        x[i] = RCONST(0.0);
      }

      /* iterate down column of A, collecting nonzeros */
      for (p=Ap[j]; p<Ap[j+1]; p++) {
        w[Ai[p]] += 1;         /* indicate that row is filled */
        x[Ai[p]] = c*Ax[p];    /* collect/scale value */
      }

      /* iterate down column of B, collecting nonzeros */
      for (p=Bp[j]; p<Bp[j+1]; p++) {
        w[Bi[p]] += 1;       /* indicate that row is filled */
        x[Bi[p]] += Bx[p];   /* collect value */
      }

      /* fill entries of C with this column's data */
      for (i=0; i<M; i++) {
        if ( w[i] > 0 ) { 
          Ci[nz] = i;  
          Cx[nz++] = x[i];
        }
      }
    }

    /* indicate end of data */
    Cp[N] = nz;

    /* update A's structure with C's values; nullify C's pointers */
    SM_NNZ_S(A) = SM_NNZ_S(C);

    free(SM_DATA_S(A));  
    SM_DATA_S(A) = SM_DATA_S(C);
    SM_DATA_S(C) = NULL;

    free(SM_INDEXVALS_S(A));
    SM_INDEXVALS_S(A) = SM_INDEXVALS_S(C);
    SM_INDEXVALS_S(C) = NULL;

    free(SM_INDEXPTRS_S(A));
    SM_INDEXPTRS_S(A) = SM_INDEXPTRS_S(C);
    SM_INDEXPTRS_S(C) = NULL;

    /* clean up */
    SUNMatDestroy_Sparse(C); 

    /* reallocate the new matrix to remove extra space */
    SUNSparseMatrix_Realloc(A);

  }

  /* clean up */
  free(w);
  free(x);

  /* return success */
  return(0);

}

int SUNMatMatvec_Sparse(SUNMatrix A, N_Vector x, N_Vector y)
{
  /* Verify that A, x and y are compatible */
  if (!SMCompatible2_Sparse(A, x, y))
    return 1;

  /* Perform operation */
  if(SM_SPARSETYPE_S(A) == CSC_MAT)
    return Matvec_SparseCSC(A, x, y);
  else 
    return Matvec_SparseCSR(A, x, y);
}


/*
 * =================================================================
 * private functions
 * =================================================================
 */

/* -----------------------------------------------------------------
 * Function to check compatibility of two sparse SUNMatrix objects
 */

static booleantype SMCompatible_Sparse(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must be sparse */
  if ( (SUNMatGetID(A) != SUNMATRIX_SPARSE) ||
       (SUNMatGetID(B) != SUNMATRIX_SPARSE) )
    return FALSE;

  /* both matrices must have the same shape and sparsity type */
  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Rows(B))
    return FALSE;
  if (SUNSparseMatrix_Columns(A) != SUNSparseMatrix_Columns(B))
    return FALSE;
  if (SM_SPARSETYPE_S(A) != SM_SPARSETYPE_S(B))
    return FALSE;

  return TRUE;
}


/* -----------------------------------------------------------------
 * Function to check compatibility of a SUNMatrix object with two
 * N_Vectors (A*x = b)
 */

static booleantype SMCompatible2_Sparse(SUNMatrix A, N_Vector x, N_Vector y)
{

  /*   vectors must be one of {SERIAL, OPENMP, PTHREADS} */ 
  if ( (N_VGetVectorID(x) != SUNDIALS_NVEC_SERIAL) &&
       (N_VGetVectorID(x) != SUNDIALS_NVEC_OPENMP) &&
       (N_VGetVectorID(x) != SUNDIALS_NVEC_PTHREADS) )
    return FALSE;

/*   /\* vectors must be one of {SERIAL, OPENMP, PTHREADS}, and  */
/*      have compatible dimensions *\/  */
/*   if (N_VGetVectorID(x) == SUNDIALS_NVEC_SERIAL) { */
/*     if (N_VGetLength_Serial(x) != SUNSparseMatrix_Columns(A)) */
/*       return FALSE; */
/*   } */
/* #ifdef SUNDIALS_OPENMP_ENABLED */
/*   else if (N_VGetVectorID(x) == SUNDIALS_NVEC_OPENMP) { */
/*     if (N_VGetLength_OpenMP(x) != SUNSparseMatrix_Columns(A)) */
/*       return FALSE; */
/*   } */
/* #endif */
/* #ifdef SUNDIALS_PTHREADS_ENABLED */
/*   else if (N_VGetVectorID(x) == SUNDIALS_NVEC_PTHREADS) { */
/*     if (N_VGetLength_Pthreads(x) != SUNSparseMatrix_Columns(A)) */
/*       return FALSE; */
/*   } */
/* #endif */
/*   else {   /\* incompatible type *\/ */
/*     return FALSE; */
/*   } */
  
/*   if (N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL) { */
/*     if (N_VGetLength_Serial(y) != SUNSparseMatrix_Rows(A)) */
/*       return FALSE; */
/*   } */
/* #ifdef SUNDIALS_OPENMP_ENABLED */
/*   else if (N_VGetVectorID(y) == SUNDIALS_NVEC_OPENMP) { */
/*     if (N_VGetLength_OpenMP(y) != SUNSparseMatrix_Rows(A)) */
/*       return FALSE; */
/*   } */
/* #endif */
/* #ifdef SUNDIALS_PTHREADS_ENABLED */
/*   else if (N_VGetVectorID(y) == SUNDIALS_NVEC_PTHREADS) { */
/*     if (N_VGetLength_Pthreads(y) != SUNSparseMatrix_Rows(A)) */
/*       return FALSE; */
/*   } */
/* #endif */
/*   else {   /\* incompatible type *\/ */
/*     return FALSE; */
/*   } */

  return TRUE;
}


/* -----------------------------------------------------------------
 * Computes y=A*x, where A is a CSC SUNMatrix_Sparse of dimension MxN, x is a 
 * compatible N_Vector object of length N, and y is a compatible 
 * N_Vector object of length M.
 * 
 * Returns 0 if successful, 1 if unsuccessful (failed memory access, or both 
 * x and y are the same vector).
 */
int Matvec_SparseCSC(SUNMatrix A, N_Vector x, N_Vector y)
{
  long int i, j;
  long int *Ap, *Ai;
  realtype *Ax, *xd, *yd;

  /* access data from CSC structure (return if failure) */
  Ap = SM_INDEXPTRS_S(A);
  Ai = SM_INDEXVALS_S(A);
  Ax = SM_DATA_S(A);
  if ((Ap == NULL) || (Ai == NULL) || (Ax == NULL))
    return 1;

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  yd = N_VGetArrayPointer(y);
  if ((xd == NULL) || (yd == NULL) || (xd == yd) )
    return 1;

  /* initialize result */
  for (i=0; i<SM_ROWS_S(A); i++)
    yd[i] = 0.0;

  /* iterate through matrix columns */
  for (j=0; j<SM_COLUMNS_S(A); j++) {

    /* iterate down column of A, performing product */
    for (i=Ap[j]; i<Ap[j+1]; i++)
      yd[Ai[i]] += Ax[i]*xd[j];

  }

  return 0;
}


/* -----------------------------------------------------------------
 * Computes y=A*x, where A is a CSR SUNMatrix_Sparse of dimension MxN, x is a 
 * compatible N_Vector object of length N, and y is a compatible 
 * N_Vector object of length M.
 * 
 * Returns 0 if successful, 1 if unsuccessful (failed memory access).
 */
int Matvec_SparseCSR(SUNMatrix A, N_Vector x, N_Vector y)
{
  long int i, j;
  long int *Ap, *Aj;
  realtype *Ax, *xd, *yd;

  /* access data from CSR structure (return if failure) */
  Ap = SM_INDEXPTRS_S(A);
  Aj = SM_INDEXVALS_S(A);
  Ax = SM_DATA_S(A);
  if ((Ap == NULL) || (Aj == NULL) || (Ax == NULL))
    return 1;

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  yd = N_VGetArrayPointer(y);
  if ((xd == NULL) || (yd == NULL) || (xd == yd))
    return 1;

  /* initialize result */
  for (i=0; i<SM_ROWS_S(A); i++)
    yd[i] = 0.0;

  /* iterate through matrix rows */
  for (i=0; i<SM_ROWS_S(A); i++) {

    /* iterate along row of A, performing product */
    for (j=Ap[i]; j<Ap[i+1]; j++)
      yd[i] += Ax[j]*xd[Aj[j]];

  }

  return(0);
}


