/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-11-08 01:01:34 $
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

  /* Note that M and N are tested in denalloc */

  A = NULL;
  A = (DenseMat) malloc(sizeof *A);
  if (A==NULL) return (NULL);
  
  A->data = NULL;
  A->data = denalloc(M, N);
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }

  A->M = M;
  A->N = N;

  return(A);
}

void DenseFreeMat(DenseMat A)
{
  denfree(A->data);
  free(A); A = NULL;
}

long int *DenseAllocPiv(long int N)
{
  return(denallocpiv(N));
}

void DenseFreePiv(long int *p)
{  
  denfreepiv(p);
}

realtype *DenseAllocBeta(long int M)
{
  return(denallocbeta(M));
}

void DenseFreeBeta(realtype *beta)
{
  denfreebeta(beta);
}

long int DenseGETRF(DenseMat A, long int *p)
{
  return(denGETRF(A->data, A->M, A->N, p));
}

void DenseGETRS(DenseMat A, long int *p, realtype *b)
{
  denGETRS(A->data, A->N, p, b);
}

int DensePOTRF(DenseMat A)
{
  return(denPOTRF(A->data, A->M));
}

void DensePOTRS(DenseMat A, realtype *b)
{
  denPOTRS(A->data, A->M, b);
}

int DenseGEQRF(DenseMat A, realtype *beta, realtype *wrk)
{
  return(denGEQRF(A->data, A->M, A->N, beta, wrk));
}

int DenseORMQR(DenseMat A, realtype *beta, realtype *vn, realtype *vm, realtype *wrk)
{
  return(denORMQR(A->data, A->M, A->N, beta, vn, vm, wrk));
}

void DenseZero(DenseMat A)
{
  denzero(A->data, A->M, A->N);
}

void DenseCopy(DenseMat A, DenseMat B)
{
  dencopy(A->data, B->data, A->M, A->N);
}

void DenseScale(realtype c, DenseMat A)
{
  denscale(c, A->data, A->M, A->N);
}

void DenseAddI(DenseMat A)
{
  denaddI(A->data, A->N);
}

void DensePrint(DenseMat A)
{
  denprint(A->data, A->M, A->N);
}
