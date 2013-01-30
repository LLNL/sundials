/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Programmer: Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2013, 
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for operations on teh SUNDIALS
 * sparse matrix structure.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_sparse.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

SlsMat NewSparseMat(long int M, long int N, long int NNZ)
{
  SlsMat A;
  long int j;

  if ( (M <= 0) || (N <= 0) ) return(NULL);

  A = NULL;
  A = (SlsMat) malloc(sizeof *A);
  if (A==NULL) return (NULL);
  
  A->data = (realtype *) malloc(NNZ * sizeof(realtype));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }
  A->rowvals = (long int *) malloc(NNZ * sizeof(long int));
  if (A->rowvals == NULL) {
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }
  A->colptrs = (long int *) malloc((N+1) * sizeof(long int));
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
  free(A->data);  A->data = NULL;
  free(A->rowvals);
  free(A->colptrs);
  free(A); A = NULL;
}


void SlsSetToZero(SlsMat A)
{
  long int i;

  for (i=0; i<A->NNZ; i++) {
    A->data[i] = ZERO;
    A->rowvals[i] = 0;
  }

  for (i=0; i<A->N; i++) {
    A->colptrs[i] = 0;
  }
  A->colptrs[A->N] = A->NNZ;

}


void PrintSparseMat(SlsMat A)
{
  long int i,j, M, N, NNZ;
  long int *colptrs;

  colptrs = A->colptrs;
  M = A->M;
  N = A->N;
  NNZ = A->NNZ;

  printf("\n");
  
  printf("%ld by %ld NNZ: %ld \n", M, N, NNZ);
  for (j=0; j < A->N; j++) {
    printf("  col %ld : locations %ld to %ld\n", j, colptrs[j], colptrs[j+1]-1);
    for (i = colptrs[j]; i < colptrs[j+1]; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("%ld  %12Lg  ", A->rowvals[i], A->data[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("%ld  %12lg  ", A->rowvals[i], A->data[i]);
#else
      printf("%ld  %12g  ", A->rowvals[i], A->data[i]);
#endif
    }
  }
  printf("\n");
    
}


