/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Programmer: Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for operations on the SUNDIALS
 * sparse matrix structure.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_sparse.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

SlsMat NewSparseMat(int M, int N, int NNZ)
{
  SlsMat A;

  if ( (M <= 0) || (N <= 0) ) return(NULL);

  A = NULL;
  A = (SlsMat) malloc(sizeof *A);
  if (A==NULL) return (NULL);
  
  A->data = (realtype *) malloc(NNZ * sizeof(realtype));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }
  A->rowvals = (int *) malloc(NNZ * sizeof(int));
  if (A->rowvals == NULL) {
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }
  A->colptrs = (int *) malloc((N+1) * sizeof(int));
  if (A->colptrs == NULL) {
    free(A->rowvals);
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }

  A->M = M;
  A->N = N;
  A->NNZ = NNZ;
  A->colptrs[N] = NNZ;

  return(A);
}


void DestroySparseMat(SlsMat A)
{
  if (A->data) {
    free(A->data);  
    A->data = NULL;
  }
  if (A->rowvals) {
    free(A->rowvals);
    A->rowvals = NULL;
  }
  if (A->colptrs) {
    free(A->colptrs);
    A->colptrs = NULL;
  }
  free(A); A = NULL;
}


void SlsSetToZero(SlsMat A)
{
  int i;

  for (i=0; i<A->NNZ; i++) {
    A->data[i] = ZERO;
    A->rowvals[i] = 0;
  }

  for (i=0; i<A->N; i++) {
    A->colptrs[i] = 0;
  }
  A->colptrs[A->N] = A->NNZ;

}

/* 
Copy sparse matrix A into sparse matrix B.  
Assumed that A and B are the same size. 
*/
void CopySparseMat(SlsMat A, SlsMat B)
{
  int i;

  for (i=0; i<A->NNZ; i++){
    B->data[i] = A->data[i];
    B->rowvals[i] = A->rowvals[i];
  }

  for (i=0; i<A->N; i++) {
    B->colptrs[i] = A->colptrs[i];
  }
  B->colptrs[A->N] = A->NNZ;

}

/* 
Scale sparse matrix A by coefficient b.  
*/
void ScaleSparseMat(realtype b, SlsMat A)
{
  int i;

  for (i=0; i<A->NNZ; i++){
    A->data[i] = b * (A->data[i]);
  }
}

/* 
Add the identity to a sparse matrix.  May need to resize if df/dy has a 
0-valued diagonal entry.
*/
void AddIdentitySparseMat(SlsMat A)
{
  int j, i, p, M, N, nz;
  int *w, *Cp, *Ap, *Ai, *Ci;
  realtype *x, *Ax, *Cx;
  SlsMat C;

  M = A->M;
  N = A->N;

  w = (int *)malloc(M * sizeof(int));
  x = (realtype *)malloc(M * sizeof(realtype));
  C = NewSparseMat(A->M, A->N, (A->NNZ)+M);

  Cp = C->colptrs;
  Ci = C->rowvals;
  Cx = C->data;
  Ap = A->colptrs;
  Ai = A->rowvals;
  Ax = A->data;

  /* Initialize values */
  nz = 0;
  for (j=0; j<M; j++) {
    w[j] = 0;
    x[j] = 0.0;
  }

  for (j=0; j<N; j++) {
    Cp[j] = nz;
    for (p=Ap[j]; p<Ap[j+1]; p++) {
      i = Ai[p];
      w[i] = j+1;
      Ci[nz] = i;
      nz++;
      x[i] = Ax[p];
    }
    if (w[j] < j+1) {
      Ci[nz] = j;
      nz++;
      w[j] = 1.0;
    } else {
      x[j] += 1.0;
    }
    for (p=Cp[j]; p<nz; p++) {
      Cx[p] = x[Ci[p]];
    }

  }
  Cp[N] = nz;
  
  if (A->data) {
    free(A->data);  
    A->data = C->data;
    C->data = NULL;
  }
  if (A->rowvals) {
    free(A->rowvals);
    A->rowvals = C->rowvals;
    C->rowvals = NULL;
  }
  if (A->colptrs) {
    free(A->colptrs);
    A->colptrs = C->colptrs;
    C->data = NULL;
  }
  DestroySparseMat(C); 

}

void PrintSparseMat(SlsMat A)
{
  int i,j, M, N, NNZ;
  int *colptrs;

  colptrs = A->colptrs;

  M = A->M;
  N = A->N;
  NNZ = A->NNZ;

  printf("\n");
  
  printf("%d by %d NNZ: %d \n", M, N, NNZ);
  for (j=0; j < A->N; j++) {
    printf("  col %d : locations %d to %d\n", j, colptrs[j], colptrs[j+1]-1);
    for (i = colptrs[j]; i < colptrs[j+1]; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("%d  %12Lg  ", A->rowvals[i], A->data[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("%d  %12g  ", A->rowvals[i], A->data[i]);
#else
      printf("%d  %12g  ", A->rowvals[i], A->data[i]);
#endif
    }
    printf("\n");
  }
  printf("\n");
    
}


