/* Header files */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <mpi.h>

/* Main Test */

int main(int argc, char *argv[])
{
  int fails = 0;
  int globfails = 0;
  int retval;
  sunindextype local_length;
  sunindextype global_length;
  N_Vector U, V, X, Y, Z;
  int print_timing;
  MPI_Comm comm;
  int nprocs, myid;

  HYPRE_Int *partitioning;
  HYPRE_ParVector Xhyp;

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);

  comm = MPI_COMM_WORLD;
  Test_Init(&comm);
  
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /* Check inputs */

  if (argc < 3) {
    if (myid == 0)
      printf("ERROR: local vector length must be a positive intiger \n");
    Test_AbortMPI(&comm, -1);
  }

  print_timing = atoi(argv[2]);
  setTiming(print_timing, myid);

  /* global length */
  global_length = nprocs*local_length;

  if (myid == 0) {
    printf("Testing the hypre ParVector N_Vector wrapper \n");
    printf("Vector global length &ld \n", (long int) global_length);
    printf("MPI processes %d \n\n", nprocs);
  }

  /* Set partitioning */

  if(HYPRE_AssumedPartitionCheck()) {
    partitioning = (HYPRE_Int*) malloc(2*sizeof(HYPRE_Int));
    partitioning[0] = myid*local_length;
    partitioning[1] = (myid+1)*local_length;
  } else {
    partitioning = (HYPRE_Int*) malloc((nprocs+1)*sizeof(HYPRE_Int));
    if (local_length < 1) {
      printf("Using global partitioning. \n");
      printf("I don't do this stuff. Now exiting...\n");
      Test_AbortMPI(&comm, 1);
    }
  }

  /* Create a hypre vector */

  HYPRE_ParVectorCreate(comm, global_length, partitioning, &Xhyp);
  HYPRE_ParVectorInitialize(Xhyp);

  /* Create a hypre ParVector N_Vector wrapper and test */

  X = N_VMake_ParHyp(Xhyp, sunctx);
  fails += Test_N_VMake(X, local_length, myid);
  if (fails != 0) {
    N_VDestroy(X);
    HYPRE_ParVectorDestroy(Xhyp);
    if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
    Test_AbortMPI(&comm, 1);
  }


  /* check vector ID */

  fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC, PARHYP, myid);

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
  
  /* Print results */

  if (fails) {
    printf("FAIL: NVector module failed %i tests, Proc %d \n\n", fails, myid);
  } else {
    if (myid == 0)
      printf("SUCCESS: NVector module passed all tests \n\n");
  }

  /*Check if any other processes failed*/

  (void) MPI_Allreduce(&fails, &globfails, 1, MPI_INT, MPI_MAX, comm);
  
  Test_Finalize();
  MPI_Finalize();
  return(globfails);
}


