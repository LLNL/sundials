/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol Dense module
 * implementation.
 * ----------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_magmadense.h>
#include <sunmatrix/sunmatrix_magmadense.h>
#include <sundials/sundials_math.h>
#include "test_sunlinsol.h"

#if defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#define HIP_OR_CUDA(a,b) a
#elif defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#define HIP_OR_CUDA(a,b) b
#else
#define HIP_OR_CUDA(a,b) ((void)0);
#endif

#if defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#include <nvector/nvector_cuda.h>
#include <sunmemory/sunmemory_cuda.h>
#elif defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#include <nvector/nvector_hip.h>
#include <sunmemory/sunmemory_hip.h>
#endif

/* ----------------------------------------------------------------------
 * SUNLinSol_MagmaDense Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int             fails = 0;          /* counter for test failures  */
  sunindextype    cols, rows;         /* matrix columns, rows       */
  sunindextype    nblocks;            /* number of matrix blocks    */
  SUNLinearSolver LS;                 /* solver object              */
  SUNMatrix       A, I;               /* test matrices              */
  N_Vector        x, b;               /* test vectors               */
  int             print_timing;
  sunindextype    i, j, k;
  realtype        *Adata, *Idata, *xdata;
  SUNContext      sunctx;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  SUNMemoryHelper memhelper = HIP_OR_CUDA( SUNMemoryHelper_Hip(sunctx);,
                                           SUNMemoryHelper_Cuda(sunctx); )

  /* check input and set matrix dimensions */
  if (argc < 4){
    printf("ERROR: THREE (3) Inputs required: matrix cols, number of blocks, print timing \n");
    return(-1);
  }

  cols = (sunindextype) atol(argv[1]);
  if (cols <= 0) {
    printf("ERROR: number of matrix columns must be a positive integer \n");
    return(-1);
  }
  rows = cols;

  nblocks = (sunindextype) atol(argv[2]);
  if (nblocks <= 0) {
    printf("ERROR: number of blocks must be a positive integer \n");
    return(-1);
  }

  print_timing = atoi(argv[3]);
  SetTiming(print_timing);

  printf("\n MAGMA dense linear solver test: size %ld, blocks %ld\n\n",
         (long int) cols, (long int) nblocks);

  /* Create matrices and vectors */
  if (nblocks > 1)
    A = SUNMatrix_MagmaDenseBlock(nblocks, rows, cols, SUNMEMTYPE_DEVICE, memhelper, NULL, sunctx);
  else
    A = SUNMatrix_MagmaDense(rows, cols, SUNMEMTYPE_DEVICE, memhelper, NULL, sunctx);
  I = SUNMatClone(A);
  x = HIP_OR_CUDA( N_VNew_Hip(cols*nblocks, sunctx);,
                   N_VNew_Cuda(cols*nblocks, sunctx); )
  b = N_VClone(x);

  /* Allocate host data */
  Adata = (realtype*) malloc(sizeof(realtype)*SUNMatrix_MagmaDense_LData(A));
  Idata = (realtype*) malloc(sizeof(realtype)*SUNMatrix_MagmaDense_LData(I));

  /* Fill A matrix with uniform random data in [0,1/cols] */
  for (k=0; k<nblocks; k++)
    for (j=0; j<cols; j++)
      for (i=0; i<rows; i++)
        Adata[k*cols*rows + j*rows + i] =
            (realtype) rand() / (realtype) RAND_MAX / cols;

  /* Create anti-identity matrix */
  for (k=0; k<nblocks; k++)
    for(j=0; j<cols; j++)
      for (i=0; i<rows; i++)
        Idata[k*cols*rows + j*rows + i] =
            ((rows-1-i) == j) ? RCONST(1.0) : RCONST(0.0);

  /* Add anti-identity to ensure the solver needs to do row-swapping */
  for (k=0; k<nblocks; k++)
    for (i=0; i<rows; i++)
      for(j=0; j<cols; j++)
        Adata[k*cols*rows + j*rows + i] +=
            Idata[k*cols*rows + j*rows + i];

  SUNMatrix_MagmaDense_CopyToDevice(A, Adata);
  SUNMatrix_MagmaDense_CopyToDevice(I, Idata);

  /* Fill x vector with uniform random data in [0,1] */
  xdata = N_VGetArrayPointer(x);
  for (j=0; j<cols*nblocks; j++)
    xdata[j] = (realtype) rand() / (realtype) RAND_MAX;
  HIP_OR_CUDA( N_VCopyToDevice_Hip(x);,
               N_VCopyToDevice_Cuda(x); )

  /* create right-hand side vector for linear solve */
  fails += SUNMatMatvecSetup(A);
  fails += SUNMatMatvec(A, x, b);
  if (fails) {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");

    /* Free matrices and vectors */
    SUNMatDestroy(A);
    SUNMatDestroy(I);
    N_VDestroy(x);
    N_VDestroy(b);

    free(Adata);
    free(Idata);

    return(1);
  }

  /* Create dense linear solver */
  LS = SUNLinSol_MagmaDense(x, A, sunctx);
  if (LS == NULL) {
    printf("FAIL: SUNLinSol_MagmaDense failure\n");

    /* Free matrices and vectors */
    SUNMatDestroy(A);
    SUNMatDestroy(I);
    N_VDestroy(x);
    N_VDestroy(b);

    free(Adata);
    free(Idata);

    return(1);
  }

  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSetup(LS, A, 0);
  fails += Test_SUNLinSolSolve(LS, A, x, b, RCONST(1e-10), SUNTRUE, 0);
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
  fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_MAGMADENSE, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
  } else {
    printf("SUCCESS: SUNLinSol module passed all tests \n \n");
  }

  /* Free solver, matrix and vectors */
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNMatDestroy(I);
  N_VDestroy(x);
  N_VDestroy(b);

  free(Adata);
  free(Idata);

  SUNContext_Free(&sunctx);

  return(fails);
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_vector(N_Vector actual, N_Vector expected, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* copy vectors to host */
  HIP_OR_CUDA( N_VCopyFromDevice_Hip(actual);,
               N_VCopyFromDevice_Cuda(actual); )
  HIP_OR_CUDA( N_VCopyFromDevice_Hip(expected);,
               N_VCopyFromDevice_Cuda(expected); )

  /* get vector data */
  xdata = N_VGetArrayPointer(actual);
  ydata = N_VGetArrayPointer(expected);

  /* check data lengths */
  xldata = N_VGetLength(actual);
  yldata = N_VGetLength(expected);


  if (xldata != yldata) {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return(1);
  }

  /* check vector data */
  for(i=0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO) {
    printf("Check_vector failures:\n");
    for(i=0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        printf("  actual[%ld] = %g != %e (err = %g)\n", (long int) i,
               xdata[i], ydata[i], SUNRabs(xdata[i]-ydata[i]));
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

void sync_device()
{
  HIP_OR_CUDA( hipDeviceSynchronize();,
               cudaDeviceSynchronize(); )
}
