/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-10-11 16:34:20 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a generic package of dense
 * matrix operations.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_dense.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Implementation */

DenseMat DenseAllocMat(long int M, long int N)
{
  DenseMat A;

  if ( (M <=0) || (N<=0) ) return(NULL);

  A = NULL;
  A = (DenseMat) malloc(sizeof *A);
  if (A==NULL) return (NULL);
  
  A->data = NULL;
  A->data = denalloc(M, N);
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }

  A->rows = M;
  A->cols = N;

  return(A);
}

long int *DenseAllocPiv(long int M)
{
  return(denallocpiv(M));
}

long int DenseGETRF(DenseMat A, long int *p)
{
  return(denGETRF(A->data, A->rows, A->cols, p));
}

void DenseGETRS(DenseMat A, long int *p, realtype *b)
{
  denGETRS(A->data, A->rows, p, b);
}

void DenseZero(DenseMat A)
{
  denzero(A->data, A->rows, A->cols);
}

void DenseCopy(DenseMat A, DenseMat B)
{
  dencopy(A->data, B->data, A->rows, A->cols);
}

void DenseScale(realtype c, DenseMat A)
{
  denscale(c, A->data, A->rows, A->cols);
}

void DenseAddI(DenseMat A)
{
  denaddI(A->data, A->rows);
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
  denprint(A->data, A->rows, A->cols);
}
