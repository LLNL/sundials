/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:34 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a generic package of dense
 * matrix operations based on BLAS/LAPACK.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_lapack.h>
#include <sundials/sundials_math.c>

LapackMat LapackAllocDenseMat(int M, int N)
{
  LapackMat A;
  int j;

  if ( (M <= 0) || (N <= 0) ) return(NULL);

  A = NULL;
  A = (LapackMat) malloc(sizeof *A);
  if (A==NULL) return (NULL);
  
  A->data = (realtype *) malloc(M * N * sizeof(realtype));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }
  A->cols = (realtype **) malloc(N * sizeof(realtype *));
  if (A->cols == NULL) {
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }

  for (j=0; j < N; j++) A->cols[j] = A->data + j * M;

  A->M = M;
  A->N = N;
  A->ldim = M;
  A->ldata = M*N;

  A->type = LAPACK_DENSE;

  return(A);
}

LapackMat LapackAllocBandMat(int N, int mu, int ml, int smu)
{
  LapackMat A;
  int j, colSize;

  if (N <= 0) return(NULL);
  
  A = NULL;
  A = (LapackMat) malloc(sizeof *A);
  if (A == NULL) return (NULL);

  colSize = smu + ml + 1;
  A->data = NULL;
  A->data = (realtype *) malloc(N * colSize * sizeof(realtype));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }

  A->cols = NULL;
  A->cols = (realtype **) malloc(N * sizeof(realtype *));
  if (A->cols == NULL) {
    free(A->data);
    free(A); A = NULL;
    return(NULL);
  }

  for (j=0; j < N; j++) A->cols[j] = A->data + j * colSize;

  A->M = N;
  A->N = N;
  A->mu = mu;
  A->ml = ml;
  A->storage_mu = smu;
  A->ldim =  colSize;
  A->ldata = N * colSize;

  A->type = LAPACK_BAND;

  return(A);
}

void LapackFreeMat(LapackMat A)
{
  free(A->data);  A->data = NULL;
  free(A->cols);
  free(A); A = NULL;
}

int *LapackAllocIntArray(int N)
{
  int *vec;

  if (N <= 0) return(NULL);

  vec = NULL;
  vec = (int *) malloc(N * sizeof(int));

  return(vec);
}

realtype *LapackAllocRealArray(int N)
{
  realtype *vec;

  if (N <= 0) return(NULL);

  vec = NULL;
  vec = (realtype *) malloc(N * sizeof(realtype));

  return(vec);
}

void LapackFreeArray(void *vec)
{ 
  free(vec); 
  vec = NULL;
}

void LapackPrintMat(LapackMat A)
{
  int i, j, start, finish;;
  realtype **a;

  switch (A->type) {

  case LAPACK_DENSE:

    printf("\n");
    for (i=0; i < A->M; i++) {
      for (j=0; j < A->N; j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%12Lg  ", LAPACK_DENSE_ELEM(A,i,j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%12lg  ", LAPACK_DENSE_ELEM(A,i,j));
#else
        printf("%12g  ", LAPACK_DENSE_ELEM(A,i,j));
#endif
      }
      printf("\n");
    }
    printf("\n");

    break;

  case LAPACK_BAND:

    a = A->cols;
    printf("\n");
    for (i=0; i < A->N; i++) {
      start = MAX(0,i-A->ml);
      finish = MIN(A->N-1,i+A->mu);
      for (j=0; j < start; j++) printf("%12s  ","");
      for (j=start; j <= finish; j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%12Lg  ", a[j][i-j+A->storage_mu]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%12lg  ", a[j][i-j+A->storage_mu]);
#else
        printf("%12g  ", a[j][i-j+A->storage_mu]);
#endif
      }
      printf("\n");
    }
    printf("\n");
    
    
    break;

  }
}
