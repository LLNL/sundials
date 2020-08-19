/* -----------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the timing routine to time various data transfer amounts
 * for the NVECTOR Parallel Multiple Dot Product function. 
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <nvector/nvector_parallel.h>

/* define constants */
#define ZERO     RCONST(0.0)
#define ONE      RCONST(1.0)
#define TWO      RCONST(2.0)

/* macro for printing timings */
#define FMT "%s Time: %22.15e\n\n"
#define PRINT_TIME(test, time) printf(FMT, test, time)

/* private functions */
double max_time(double time);

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          retval;                /* function return value     */
  sunindextype local_length;          /* local vector length       */
  sunindextype global_length;         /* global vector length      */
  N_Vector     X;                     /* test vector               */
  N_Vector     *V;                    /* nvector array containing clones of X      */
  realtype     *dotprods;             /* array containing dotprods of vectors in V */
  realtype     *xd, *yd;              /* arrays used in dotprod calculations       */
  int          num_vecs;              /* number of vectors in V    */ 
  MPI_Comm     comm;                  /* MPI Communicator          */
  int          nprocs, myid;          /* Number of procs, proc id  */
  double       start_time, stop_time; /* start and stop times for processes   */ 
  double       maxt;                  /* maximum time required by any process */
  int i, j;

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /* check inputs */
  if (argc < 3) {
    if (myid == 0)
      printf("ERROR: TWO (2) Inputs required: vector length   number of vectors \n");
    MPI_Abort(comm, -1);
  }

  local_length = (sunindextype) atol(argv[1]);
  if (local_length < 1) {
    if (myid == 0)
      printf("ERROR: local vector length must be a positive integer \n");
    MPI_Abort(comm, -1);
  }

  num_vecs = atoi(argv[2]);
  if (num_vecs < 1) {
    if (myid == 0)
      printf("ERROR: number of vectors must be a positive integer \n");
    MPI_Abort(comm, -1);
  }

  /* global length */
  global_length = nprocs*local_length;

  if (myid == 0) {
    printf("Timing the parallel (MPI) N_Vector All Reduce\n");
    printf("Vector global length %ld \n", (long int) global_length);
    printf("MPI processes %d \n", nprocs);
  }

  /* Create new vector */
  X = N_VNew_Parallel(comm, local_length, global_length);
  if (X == NULL) {
    if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
    MPI_Abort(comm, 1);
  }
  
  /* Enable fused vector operations */ 
  retval = N_VEnableFusedOps_Parallel(X, SUNTRUE);

  /* Perform All Reduce Timing Test */ 

  /* Create vectors for testing */
  V = N_VCloneVectorArray(num_vecs, X);
  /* Create array to store dot products */
  dotprods = NULL;
  dotprods = (realtype *) malloc(num_vecs * sizeof(realtype));

  /* fill vector data */
  N_VConst(TWO, X);
  for (i = 0; i < num_vecs; i++)
    N_VConst(ONE, V[i]);
  
  /* get data array */
  xd   = NV_DATA_P(X);

  MPI_Barrier(MPI_COMM_WORLD);
  /* compute multiple dot products */
  for (i=0; i<num_vecs; i++) {
    yd = NV_DATA_P(V[i]);
    dotprods[i] = ZERO;
    for (j=0; j<local_length; j++) {
      dotprods[i] += xd[j] * yd[j];
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  start_time = MPI_Wtime();
  retval = MPI_Allreduce(MPI_IN_PLACE, dotprods, num_vecs, MPI_SUNREALTYPE, MPI_SUM, MPI_COMM_WORLD);
  stop_time = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);

  /* find max time across all processes */
  maxt = max_time(stop_time - start_time);

  /* Print result */
  if (myid == 0) PRINT_TIME("N_VDotProdMulti", maxt);

  /* Free memory */
  N_VDestroyVectorArray(V, num_vecs);
  N_VDestroy(X);
  free(dotprods);

  MPI_Finalize();

  return(0);
}

/* ----------------------------------------------------------------------
 * Grab max time of all mpi processes 
 * --------------------------------------------------------------------*/
double max_time(double time)
{
  double maxt;

  /* get max time across all MPI ranks */
  (void) MPI_Reduce(&time, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  return(maxt);
}
