/*
 * -----------------------------------------------------------------
 * $Revision: 4137 $
 * $Date: 2014-06-15 12:26:15 -0700 (Sun, 15 Jun 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the NVECTOR Parallel module 
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <petscvec.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_petsc.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <mpi.h>

/* local vector length */
#define VECLEN 10000

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int      fails = 0;                   /* counter for test failures */
  long int local_length, global_length; /* vector lengths            */
  N_Vector W, X, Y, Z;                  /* test vectors              */
  MPI_Comm comm;                        /* MPI Communicator          */
  int      nprocs, myid;                /* Number of procs, proc id  */
  PetscErrorCode ierr;                  /* PETSc error code          */

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);
  ierr = PetscInitializeNoArguments();
  CHKERRQ(ierr);
  
  /* set local and global lengths */
  local_length = VECLEN;
  global_length = nprocs*local_length;

  /* Create vectors */
  //printf("Creating W...\n");
  W = N_VNewEmpty_petsc(comm, local_length, global_length);
  //printf("Creating X...\n");
  X = N_VNew_petsc(comm, local_length, global_length);
  Y = N_VNew_petsc(comm, local_length, global_length);
  Z = N_VNew_petsc(comm, local_length, global_length);

  /* NVector Test */
  fails += Test_N_VSetArrayPointer(W, local_length, myid);
  fails += Test_N_VGetArrayPointer(X, local_length, myid);
  fails += Test_N_VLinearSum(X, Y, Z, local_length, myid);
  fails += Test_N_VConst(X, local_length, myid);
  fails += Test_N_VProd(X, Y, Z, local_length, myid);
  fails += Test_N_VDiv(X, Y, Z, local_length, myid);
  fails += Test_N_VScale(X, Z, local_length, myid);
  fails += Test_N_VAbs(X, Z, local_length, myid);
  fails += Test_N_VInv(X, Z, local_length, myid);
  fails += Test_N_VAddConst(X, Z, local_length, myid);
  fails += Test_N_VDotProd(X, Y, local_length, global_length, myid);
  fails += Test_N_VMaxNorm(X, local_length, myid);
  fails += Test_N_VWrmsNorm(X, Y, local_length, myid);
  fails += Test_N_VWrmsNormMask(X, Y, Z, local_length, global_length, myid);
  fails += Test_N_VMin(X, local_length, myid);
  fails += Test_N_VWL2Norm(X, Y, local_length, global_length, myid);
  fails += Test_N_VL1Norm(X, local_length, global_length, myid);
  fails += Test_N_VCompare(X, Z, local_length, myid);
  fails += Test_N_VInvTest(X, Z, local_length, myid);
  fails += Test_N_VConstrMask(X, Y, Z, local_length, myid);
  fails += Test_N_VMinQuotient(X, Y, local_length, myid);
  fails += Test_N_VCloneVectorArray(5, X, local_length, myid);
  fails += Test_N_VCloneEmptyVectorArray(5, X, myid);
  fails += Test_N_VCloneEmpty(X, myid);
  fails += Test_N_VClone(X, local_length, myid);

  /* Free vectors */
  //printf("Destroying W...\n");
  N_VDestroy_petsc(W);
  //printf("Destroying X...\n");
  N_VDestroy_petsc(X);
  //printf("Destroying Y...\n");
  N_VDestroy_petsc(Y);
  //printf("Destroying Z...\n");
  N_VDestroy_petsc(Z);

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests, Proc %d \n \n", fails, myid);
  } else {
     if(myid == 0) {
	printf("SUCCESS: NVector module passed all tests, Proc %d \n \n",myid);
     }
  }
  
  //printf("Finalizing PETSc...\n");
  ierr = PetscFinalize();
  CHKERRQ(ierr);
  MPI_Finalize();
  return(0);
}
