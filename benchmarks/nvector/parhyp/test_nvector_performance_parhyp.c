/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to evaluate the performance of the
 * MPI parallel NVECTOR module implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_math.h>
#include "test_nvector_performance.h"

#include <mpi.h>
static realtype* Tdata;
/* private functions */
static int InitializeClearCache(int cachesize);
static int FinalizeClearCache();

/* private data for clearing cache */
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
static sunindextype N;  /* data length */
static realtype* data;  /* host data   */
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
//static realtype* Tdata;
static sunindextype N;    /* data length */
static realtype* h_data;  /* host data   */
static realtype* h_sum;   /* host sum    */
static realtype* d_data;  /* device data */
static realtype* d_sum;   /* device sum  */
static int blocksPerGrid;

/* cuda reduction kernel to clearing cache between tests */
__global__
void ClearCacheKernel(sunindextype N, realtype* data, realtype* out)
{
  __shared__ realtype shared[256];

  int sharedidx = blockIdx.x;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  realtype tmp = 0;
  while (tid < N) {
    tmp += data[tid];
    tid += blockDim.x * gridDim.x;
  }
  shared[sharedidx] = tmp;
  __syncthreads();

  /* assues blockDim is a power of 2 */
  int i = blockDim.x/2;
  while (i != 0) {
    if (sharedidx < i)
      shared[sharedidx] += shared[sharedidx + i];
    __syncthreads();
    i /= 2;
  }

  if (sharedidx == 0)
    out[sharedidx] = shared[0];
}
#endif

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  SUNContext   ctx = NULL;  /* SUNDIALS context */
  N_Vector     X   = NULL;  /* test vector      */
  sunindextype veclen;      /* vector length    */

  HYPRE_ParVector Xhyp;    /* hypre parallel vector */
  HYPRE_Int *partitioning; /* vector partitioning   */

  int print_timing;    /* output timings     */
  int ntests;          /* number of tests    */
  int nvecs;           /* number of tests    */
  int nsums;           /* number of sums     */
  int cachesize;       /* size of cache (MB) */
  int flag;            /* return flag        */

  MPI_Comm comm;          /* MPI Communicator   */
  int      nprocs, myid;  /* Num procs, proc id */

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /* check input and set vector length */
  if (argc < 7){
    printf("ERROR: SIX (6) arguments required: ");
    printf("<vector length> <number of vectors> <number of sums> <number of tests> ");
    printf("<cache size (MB)> <print timing>\n");
    return(-1);
  }

  veclen = (sunindextype) atol(argv[1]);
  if (veclen <= 0) {
    printf("ERROR: local vector length must be a positive integer \n");
    return(-1);
  }

  nvecs = (int) atol(argv[2]);
  if (nvecs < 1) {
    printf("WARNING: Fused operation tests disabled\n");
  }

  nsums = (int) atol(argv[3]);
  if (nsums < 1) {
    printf("WARNING: Some fused operation tests disabled\n");
  }

  ntests = (int) atol(argv[4]);
  if (ntests <= 0) {
    printf("ERROR: number of tests must be a positive integer \n");
    return(-1);
  }

  cachesize = (int) atol(argv[5]);
  if (cachesize < 0) {
    printf("ERROR: cache size (MB) must be a non-negative integer \n");
    return(-1);
  }
  InitializeClearCache(cachesize);

  print_timing = atoi(argv[6]);
  SetTiming(print_timing, myid);

  Tdata = NULL;
  Tdata = (realtype*) malloc(veclen * sizeof(realtype));
  if (!Tdata) {
    printf("ERROR: malloc failed\n");
    return(-1);
  }

  if (myid == 0)
  {
    printf("\nStart Tests\n");
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
    printf("Vector Name: ParHyp\n");
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
    printf("Vector Name: ParHyp+CUDA_Unmanaged\n");
#endif
    printf("\nRunning with: \n");
    printf("  local vector length   %ld \n", (long int) veclen);
    printf("  max number of vectors %d  \n", nvecs);
    printf("  max number of sums    %d  \n", nsums);
    printf("  number of tests       %d  \n", ntests);
    printf("  timing on/off         %d  \n", print_timing);
    printf("  number of MPI procs   %d  \n", nprocs);
  }

  flag = SUNContext_Create(&comm, &ctx);
  if (flag) return flag;

  /* set partitioning */
  if(HYPRE_AssumedPartitionCheck()) {
    partitioning = (HYPRE_Int*) malloc(2*sizeof(HYPRE_Int));
    partitioning[0] = myid*veclen;
    partitioning[1] = (myid+1)*veclen;
  } else {
    partitioning = (HYPRE_Int*) malloc((nprocs+1)*sizeof(HYPRE_Int));
    if (veclen <= 0) {
      printf("Error: Using global partition not supported.\n");
      return -1;
    }
  }

  /* create template hypre vector */
  HYPRE_ParVectorCreate(comm, nprocs*veclen, partitioning, &Xhyp);
  HYPRE_ParVectorInitialize(Xhyp);

  /* Create vectors */
  X = N_VMake_ParHyp(Xhyp, ctx);

  /* run tests */
  if (myid == 0 && print_timing) {
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
      flag = Test_N_VLinearCombinationVectorArray(X, veclen, nvecs, nsums, ntests);
    }
  }

  /* Free vectors */
  N_VDestroy(X);
  HYPRE_ParVectorDestroy(Xhyp);

  FinalizeClearCache();

  flag = SUNContext_Free(&ctx);
  if (flag) return flag;

  if (myid == 0)
    printf("\nFinished Tests\n");

  free(Tdata);

  MPI_Finalize();

  return(flag);
}


/* ----------------------------------------------------------------------
 * Functions required by testing routines to fill vector data
 * --------------------------------------------------------------------*/

/* random data between lower and upper */
void N_VRand(N_Vector Xvec, sunindextype Xlen, realtype lower, realtype upper)
{
  HYPRE_ParVector Xhyp;
  realtype *Xdata;

  Xhyp  = N_VGetVector_ParHyp(Xvec);
  Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xhyp));
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  rand_realtype(Xdata, Xlen, lower, upper);
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
  rand_realtype(Tdata, Xlen, lower, upper);
  cudaMemcpy(Xdata, Tdata,
             Xlen * sizeof(realtype),
             cudaMemcpyHostToDevice);
#endif
}

/* series of 0 and 1 */
void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen)
{
  HYPRE_ParVector Xhyp;
  realtype *Xdata;

  Xhyp  = N_VGetVector_ParHyp(Xvec);
  Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xhyp));
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  rand_realtype_zero_one(Xdata, Xlen);
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
  rand_realtype_zero_one(Tdata, Xlen);
  cudaMemcpy(Xdata, Tdata,
             Xlen * sizeof(realtype),
             cudaMemcpyHostToDevice);
#endif
}

/* random values for constraint array */
void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen)
{
  HYPRE_ParVector Xhyp;
  realtype *Xdata;

  Xhyp  = N_VGetVector_ParHyp(Xvec);
  Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xhyp));
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  rand_realtype_constraints(Xdata, Xlen);
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
  rand_realtype_constraints(Tdata, Xlen);
  cudaMemcpy(Xdata, Tdata,
             Xlen * sizeof(realtype),
             cudaMemcpyHostToDevice);
#endif
}


/* ----------------------------------------------------------------------
 * Functions required for MPI or GPU testing
 * --------------------------------------------------------------------*/

void collect_times(N_Vector X, double *times, int ntimes)
{
  MPI_Comm comm;
  int myid;

  comm = ((N_VectorContent_ParHyp)(X->content))->comm;
  MPI_Comm_rank(comm, &myid);

  if (myid == 0)
    MPI_Reduce(MPI_IN_PLACE, times, ntimes, MPI_DOUBLE, MPI_MAX, 0, comm);
  else
    MPI_Reduce(times, times, ntimes, MPI_DOUBLE, MPI_MAX, 0, comm);

  return;
}

void sync_device(N_Vector x)
{
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
  cudaDeviceSynchronize();
#endif
  return;
}


/* ----------------------------------------------------------------------
 * Functions required for clearing cache
 * --------------------------------------------------------------------*/

static int InitializeClearCache(int cachesize)
{
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)

  size_t nbytes;  /* cache size in bytes */

  /* determine size of vector to clear cache, N = ceil(2 * nbytes/realtype) */
  nbytes = (size_t) (2 * cachesize * 1024 * 1024);
  N = (sunindextype) ((nbytes + sizeof(realtype) - 1)/sizeof(realtype));

  /* allocate data and fill random values */
  data = (realtype*) malloc(N*sizeof(realtype));
  rand_realtype(data, N, RCONST(-1.0), RCONST(1.0));

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)

  cudaError_t err;     /* cuda error flag     */
  size_t      nbytes;  /* cache size in bytes */

  /* determine size of vector to clear cache, N = ceil(2 * nbytes/realtype) */
  nbytes = (size_t) (2 * cachesize * 1024 * 1024);
  N = (sunindextype) ((nbytes + sizeof(realtype) - 1)/sizeof(realtype));

  /* allocate host data */
  blocksPerGrid = SUNMIN(32,(N+255)/256);

  h_data = (realtype*) malloc(N*sizeof(realtype));
  h_sum  = (realtype*) malloc(blocksPerGrid*sizeof(realtype));

  /* allocate device data */
  err = cudaMalloc((void**) &d_data, N*sizeof(realtype));
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to allocate device vector (error code %d )!\n",err);
    return(-1);
  }

  err = cudaMalloc((void**) &d_sum, blocksPerGrid*sizeof(realtype));
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to allocate device vector (error code %d )!\n",err);
    return(-1);
  }

  /* fill host vector with random data and copy to device */
  rand_realtype(h_data, N, RCONST(-1.0), RCONST(1.0));

  err = cudaMemcpy(d_data, h_data, N*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to copy data from host to device (error code %d )!\n",err);
    return(-1);
  }

#endif

  return(0);
}

static int FinalizeClearCache()
{
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)

  free(data);

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)

  cudaError_t err;  /* cuda error flag */

  free(h_data);
  free(h_sum);

  err = cudaFree(d_data);
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to free device data (error code %d )!\n",err);
    return(-1);
  }

  err = cudaFree(d_sum);
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to free device data (error code %d )!\n",err);
    return(-1);
  }

#endif

  return(0);
}

void ClearCache()
{
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)

  realtype     sum;
  sunindextype i;

  sum = RCONST(0.0);
  for (i=0; i<N; i++)
    sum += data[i];

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)

  /* call cuda kernel to clear the cache */
  ClearCacheKernel<<<SUNMIN(32,(N+255)/256), 256>>>(N, d_data, d_sum);
  cudaMemcpy(h_sum, d_sum, blocksPerGrid*sizeof(realtype), cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

#endif

  return;
}
