/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to evaluate the performance of the
 * PETSc NVECTOR module implementation.
 * -----------------------------------------------------------------*/

#include <mpi.h>
#include <nvector/nvector_petsc.h>
#include <petscvec.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "test_nvector_performance.h"

/* private functions */
static int InitializeClearCache(int cachesize);
static int FinalizeClearCache();

/* private data for clearing cache */
static sunindextype N;    /* data length */
static sunrealtype* data; /* host data   */

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  SUNContext sunctx;
  N_Vector X;          /* test vector        */
  sunindextype veclen; /* vector length      */

  Vec Xpetsc;          /* PETSc vector     */
  PetscErrorCode ierr; /* PETSc error code */

  int print_timing; /* output timings     */
  int ntests;       /* number of tests    */
  int nvecs;        /* number of tests    */
  int nsums;        /* number of sums     */
  int cachesize;    /* size of cache (MB) */
  int flag;         /* return flag        */

  MPI_Comm comm;    /* MPI Communicator   */
  int nprocs, myid; /* Num procs, proc id */

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);
  ierr = PetscInitializeNoArguments();
  CHKERRQ(ierr);

  if (SUNContext_Create(comm, &sunctx))
  {
    printf("ERROR: SUNContext_Create returned nonzero\n");
    return (-1);
  }

  if (myid == 0)
  {
    printf("\nStart Tests\n");
    printf("Vector Name: PETSc\n");
  }

  /* check input and set vector length */
  if (argc < 7)
  {
    printf("ERROR: SIX (6) arguments required: ");
    printf("<vector length> <number of vectors> <number of sums> <number of "
           "tests> ");
    printf("<cache size (MB)> <print timing>\n");
    return (-1);
  }

  veclen = (sunindextype)atol(argv[1]);
  if (veclen <= 0)
  {
    printf("ERROR: local vector length must be a positive integer \n");
    return (-1);
  }

  nvecs = (int)atol(argv[2]);
  if (nvecs < 1) { printf("WARNING: Fused operation tests disabled\n"); }

  nsums = (int)atol(argv[3]);
  if (nsums < 1) { printf("WARNING: Some fused operation tests disabled\n"); }

  ntests = (int)atol(argv[4]);
  if (ntests <= 0)
  {
    printf("ERROR: number of tests must be a positive integer \n");
    return (-1);
  }

  cachesize = (int)atol(argv[5]);
  if (cachesize < 0)
  {
    printf("ERROR: cache size (MB) must be a non-negative integer \n");
    return (-1);
  }
  InitializeClearCache(cachesize);

  print_timing = atoi(argv[6]);
  SetTiming(print_timing, myid);

  if (myid == 0)
  {
    printf("\nRunning with: \n");
    printf("  local vector length   %ld \n", (long int)veclen);
    printf("  max number of vectors %d  \n", nvecs);
    printf("  max number of sums    %d  \n", nsums);
    printf("  number of tests       %d  \n", ntests);
    printf("  timing on/off         %d  \n", print_timing);
    printf("  number of MPI procs   %d  \n", nprocs);
  }

  /* Create vectors */
  VecCreate(comm, &Xpetsc);
  VecSetSizes(Xpetsc, veclen, nprocs * veclen);
  VecSetFromOptions(Xpetsc);
  X = N_VMake_Petsc(Xpetsc, sunctx);

  /* run tests */
  if (myid == 0 && print_timing)
  {
    printf("\n\n standard operations:\n");
    PrintTableHeader(1);
  }
  flag = Test_N_VLinearSum(X, veclen, ntests);
  flag = Test_N_VConst(X, veclen, ntests);
  flag = Test_N_VProd(X, veclen, ntests);
  flag = Test_N_VDiv(X, veclen, ntests);
  flag = Test_N_VScale(X, veclen, ntests);
  flag = Test_N_VAbs(X, veclen, ntests);
  flag = Test_N_VInv(X, veclen, ntests);
  flag = Test_N_VAddConst(X, veclen, ntests);
  flag = Test_N_VDotProd(X, veclen, ntests);
  flag = Test_N_VMaxNorm(X, veclen, ntests);
  flag = Test_N_VWrmsNorm(X, veclen, ntests);
  flag = Test_N_VWrmsNormMask(X, veclen, ntests);
  flag = Test_N_VMin(X, veclen, ntests);
  flag = Test_N_VWL2Norm(X, veclen, ntests);
  flag = Test_N_VL1Norm(X, veclen, ntests);
  flag = Test_N_VCompare(X, veclen, ntests);
  flag = Test_N_VInvTest(X, veclen, ntests);
  flag = Test_N_VConstrMask(X, veclen, ntests);
  flag = Test_N_VMinQuotient(X, veclen, ntests);

  if (nvecs > 0)
  {
    if (myid == 0 && print_timing)
    {
      printf("\n\n fused operations 1: nvecs= %d\n", nvecs);
      PrintTableHeader(2);
    }
    flag = Test_N_VLinearCombination(X, veclen, nvecs, ntests);
    flag = Test_N_VScaleAddMulti(X, veclen, nvecs, ntests);
    flag = Test_N_VDotProdMulti(X, veclen, nvecs, ntests);
    flag = Test_N_VLinearSumVectorArray(X, veclen, nvecs, ntests);
    flag = Test_N_VScaleVectorArray(X, veclen, nvecs, ntests);
    flag = Test_N_VConstVectorArray(X, veclen, nvecs, ntests);
    flag = Test_N_VWrmsNormVectorArray(X, veclen, nvecs, ntests);
    flag = Test_N_VWrmsNormMaskVectorArray(X, veclen, nvecs, ntests);

    if (nsums > 0)
    {
      if (myid == 0 && print_timing)
      {
        printf("\n\n fused operations 2: nvecs= %d nsums= %d\n", nvecs, nsums);
        PrintTableHeader(2);
      }
      flag = Test_N_VScaleAddMultiVectorArray(X, veclen, nvecs, nsums, ntests);
      flag = Test_N_VLinearCombinationVectorArray(X, veclen, nvecs, nsums,
                                                  ntests);
    }
  }

  /* Free vectors */
  N_VDestroy(X);

  FinalizeClearCache();

  if (myid == 0) { printf("\nFinished Tests\n"); }

  ierr = PetscFinalize();
  CHKERRQ(ierr);
  SUNContext_Free(&sunctx);
  MPI_Finalize();

  return (flag);
}

/* ----------------------------------------------------------------------
 * Functions required by testing routines to fill vector data
 * --------------------------------------------------------------------*/

/* random data between lower and upper */
void N_VRand(N_Vector Xvec, sunindextype Xlen, sunrealtype lower,
             sunrealtype upper)
{
  Vec Xpetsc;
  PetscScalar* Xdata;

  Xpetsc = N_VGetVector_Petsc(Xvec);
  VecGetArray(Xpetsc, &Xdata);
  rand_realtype(Xdata, Xlen, lower, upper);
  VecRestoreArray(Xpetsc, &Xdata);
}

/* series of 0 and 1 */
void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen)
{
  Vec Xpetsc;
  PetscScalar* Xdata;

  Xpetsc = N_VGetVector_Petsc(Xvec);
  VecGetArray(Xpetsc, &Xdata);
  rand_realtype_zero_one(Xdata, Xlen);
  VecRestoreArray(Xpetsc, &Xdata);
}

/* random values for constraint array */
void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen)
{
  Vec Xpetsc;
  PetscScalar* Xdata;

  Xpetsc = N_VGetVector_Petsc(Xvec);
  VecGetArray(Xpetsc, &Xdata);
  rand_realtype_constraints(Xdata, Xlen);
  VecRestoreArray(Xpetsc, &Xdata);
}

/* ----------------------------------------------------------------------
 * Functions required for MPI or GPU testing
 * --------------------------------------------------------------------*/

void collect_times(N_Vector X, double* times, int ntimes)
{
  MPI_Comm comm;
  int myid;

  comm = ((N_VectorContent_Petsc)(X->content))->comm;
  MPI_Comm_rank(comm, &myid);

  if (myid == 0)
  {
    MPI_Reduce(MPI_IN_PLACE, times, ntimes, MPI_DOUBLE, MPI_MAX, 0, comm);
  }
  else { MPI_Reduce(times, times, ntimes, MPI_DOUBLE, MPI_MAX, 0, comm); }

  return;
}

void sync_device(N_Vector x)
{
  /* not running on GPU, just return */
  return;
}

/* ----------------------------------------------------------------------
 * Functions required for clearing cache
 * --------------------------------------------------------------------*/

static int InitializeClearCache(int cachesize)
{
  size_t nbytes; /* cache size in bytes */

  /* determine size of vector to clear cache, N = ceil(2 * nbytes/sunrealtype) */
  nbytes = (size_t)(2 * cachesize * 1024 * 1024);
  N = (sunindextype)((nbytes + sizeof(sunrealtype) - 1) / sizeof(sunrealtype));

  /* allocate data and fill random values */
  data = (sunrealtype*)malloc(N * sizeof(sunrealtype));
  rand_realtype(data, N, SUN_RCONST(-1.0), SUN_RCONST(1.0));

  return (0);
}

static int FinalizeClearCache()
{
  free(data);
  return (0);
}

void ClearCache()
{
  sunrealtype sum;
  sunindextype i;

  sum = SUN_RCONST(0.0);
  for (i = 0; i < N; i++) { sum += data[i]; }

  return;
}
