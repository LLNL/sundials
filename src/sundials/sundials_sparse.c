/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Programmer: Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * begincopyright(llns)
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * endcopyright(llns)
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


