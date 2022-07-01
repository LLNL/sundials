#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <mpi.h>

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails = 0;         /* counter for test failures */
  int          globfails = 0;     /* counter for test failures */
  int          retval;            /* function return value     */
  sunindextype local_length;      /* local vector length       */
  sunindextype global_length;     /* global vector length      */
  N_Vector     U, V, X, Y, Z;     /* test vectors              */
  int          print_timing;      /* turn timing on/off        */
  MPI_Comm     comm;              /* MPI Communicator          */
  int          nprocs, myid;      /* Number of procs, proc id  */

  HYPRE_Int       *partitioning;  /* Vector Partitioning               */
  HYPRE_ParVector Xhyp;           /* Instantiate hypre parallel vector */

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);

  comm = MPI_COMM_WORLD;
  Test_Init(&comm);

  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /* check inputs */
  if (argc < 3) {
    if (myid == 0)
      printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    Test_AbortMPI(&comm, -1);
  }

  local_length = (sunindextype) atol(argv[1]);
  if (local_length < 1) {
    if (myid == 0)
      printf("ERROR: local vector length must be a positive integer \n");
    Test_AbortMPI(&comm, -1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, myid);

  /* global length */
  global_length = nprocs*local_length;

  if (myid == 0) {
    printf("Testing the hypre ParVector N_Vector wrapper \n");
    printf("Vector global length %ld \n", (long int) global_length);
    printf("MPI processes %d \n\n", nprocs);
  }

  /* set partitioning */
  if(HYPRE_AssumedPartitionCheck()) {
    partitioning    = (HYPRE_Int*) malloc(2*sizeof(HYPRE_Int));
    partitioning[0] = myid*local_length;
    partitioning[1] = (myid+1)*local_length;
  } else {
    partitioning = (HYPRE_Int*) malloc((nprocs+1)*sizeof(HYPRE_Int));
    if (local_length < 1) {
      printf("Using global partition.\n");
      printf("I don't do this stuff. Now exiting...\n");
      Test_AbortMPI(&comm, 1);
    }
  }

  /* Create hypre vector */
  HYPRE_ParVectorCreate(comm, global_length, partitioning, &Xhyp);
  HYPRE_ParVectorInitialize(Xhyp);

  /* Create hypre ParVector N_Vector wrapper and test */
  X = N_VMake_ParHyp(Xhyp, sunctx);
  fails += Test_N_VMake(X, local_length, myid);
  if (fails != 0) {
    N_VDestroy(X);
    HYPRE_ParVectorDestroy(Xhyp);
    if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
    Test_AbortMPI(&comm, 1);
  }

  /* Check vector ID */
  fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_PARHYP, myid);

 /* Check vector length */
  fails += Test_N_VGetLength(X, myid);


  /* Free vectors */
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);
  N_VDestroy(U);
  N_VDestroy(V);
  HYPRE_ParVectorDestroy(Xhyp);
#ifndef hypre_ParVectorOwnsPartitioning
  free(partitioning);
#endif

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests, Proc %d \n\n", fails, myid);
  } else {
    if (myid == 0)
      printf("SUCCESS: NVector module passed all tests \n\n");
  }

  /* check if any other process failed */
  (void) MPI_Allreduce(&fails, &globfails, 1, MPI_INT, MPI_MAX, comm);

  Test_Finalize();
  MPI_Finalize();
  return(globfails);
}

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  int             failure = 0;
  sunindextype    i;
  HYPRE_ParVector Xvec;
  realtype        *Xdata;

  Xvec  = N_VGetVector_ParHyp(X);
  Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xvec));

  /* check vector data */
  for (i = 0; i < local_length; i++) {
    failure += SUNRCompare(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  /* check if wrapped hypre ParVector is non-null */
  return (N_VGetVector_ParHyp(X) == NULL) ? SUNFALSE : SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  /* set i-th element of data array */
  set_element_range(X, i, i, val);
}

void set_element_range(N_Vector X, sunindextype is, sunindextype ie,
                       realtype val)
{
  HYPRE_ParVector  Xvec;
  realtype        *Xdata;
  sunindextype     i;

  /* set elements [is,ie] of the data array */
  Xvec  = N_VGetVector_ParHyp(X);
  Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xvec));

  for(i = is; i <= ie; i++) Xdata[i] = val;
}

realtype get_element(N_Vector X, sunindextype i)
{
  HYPRE_ParVector Xvec;
  realtype        *Xdata;

  /* get i-th element of data array */
  Xvec  = N_VGetVector_ParHyp(X);
  Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xvec));

  return Xdata[i];
}

double max_time(N_Vector X, double time)
{
  MPI_Comm comm;
  double maxt;

  /* get max time across all MPI ranks */
  comm = hypre_ParVectorComm(N_VGetVector_ParHyp(X));
  (void) MPI_Reduce(&time, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  return(maxt);
}

void sync_device(N_Vector x)
{
  /* not running on GPU, just return */
  return;
}
