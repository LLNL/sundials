/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-06-15 15:39:53 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a generic DENSE linear
 * solver package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include "sundials_dense.h"
#include "sundials_math.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Implementation */

DenseMat DenseAllocMat(long int N)
{
  DenseMat A;

  if (N <= 0) return(NULL);

  A = NULL;
  A = (DenseMat) malloc(sizeof *A);
  if (A==NULL) return (NULL);
  
  A->data = NULL;
  A->data = denalloc(N);
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }

  A->size = N;

  return(A);
}

long int *DenseAllocPiv(long int N)
{
  return(denallocpiv(N));
}

long int DenseFactor(DenseMat A, long int *p)
{
  return(gefa(A->data, A->size, p));
}

void DenseBacksolve(DenseMat A, long int *p, realtype *b)
{
  gesl(A->data, A->size, p, b);
}

void DenseZero(DenseMat A)
{
  denzero(A->data, A->size);
}

void DenseCopy(DenseMat A, DenseMat B)
{
  dencopy(A->data, B->data, A->size);
}

void DenseScale(realtype c, DenseMat A)
{
  denscale(c, A->data, A->size);
}

void DenseAddI(DenseMat A)
{
  denaddI(A->data, A->size);
}

void DenseFreeMat(DenseMat A)
{
  denfree(A->data);
  free(A); A = NULL;
}

void DenseFreePiv(long int *p)
{  
  denfreepiv(p);
}

void DensePrint(DenseMat A)
{
  denprint(A->data, A->size);
}
