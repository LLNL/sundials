/******************************************************************
 *                                                                *
 * File          : dense.c                                        *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                             *
 * Version of    : 4 March 2002                                   *
 *----------------------------------------------------------------*
 * This is the implementation file for a generic DENSE linear     *
 * solver package.                                                *
 *                                                                *
 ******************************************************************/ 

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h"
#include "llnlmath.h"
#include "dense.h"
#include "smalldense.h"


#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)


/* Implementation */


DenseMat DenseAllocMat(integer N)
{
  DenseMat A;

  if (N <= 0) return(NULL);

  A = (DenseMat) malloc(sizeof *A);
  if (A==NULL) return (NULL);
  
  A->data = denalloc(N);
  if (A->data == NULL) {
    free(A);
    return(NULL);
  }

  A->size = N;

  return(A);
}


integer *DenseAllocPiv(integer N)
{
  if (N <= 0) return(NULL);

  return((integer *) malloc(N * sizeof(integer)));
}


integer DenseFactor(DenseMat A, integer *p)
{
  return(gefa(A->data, A->size, p));
}


void DenseBacksolve(DenseMat A, integer *p, real *b)
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

void DenseScale(real c, DenseMat A)
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
  free(A);
}

void DenseFreePiv(integer *p)
{  
  free(p);
}

void DensePrint(DenseMat A)
{
  denprint(A->data, A->size);
}

