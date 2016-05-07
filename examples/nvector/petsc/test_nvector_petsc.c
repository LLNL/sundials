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

static int Test_N_VMake(Vec* W, int myid);

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
  Vec xvec;                             /* PETSc vector              */
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
  
  /* Allocate and initialize PETSc vector */
  VecCreate(comm, &xvec);
  VecSetSizes(xvec, PETSC_DECIDE, global_length);
  VecSetFromOptions(xvec);
  
  /* Create vectors */
  W = N_VNewEmpty_petsc(comm, local_length, global_length);

  /* NVector Test */

  /* PETSc specific tests */
  //fails += Test_N_VMake(&xvec, myid);
  
  X = N_VMake_petsc(&xvec);

  /* Memory allocation tests */
  fails += Test_N_VCloneVectorArray(5, X, local_length, myid);
  fails += Test_N_VCloneEmptyVectorArray(5, X, myid);
  fails += Test_N_VCloneEmpty(X, myid);
  fails += Test_N_VClone(X, local_length, myid);

  Y = N_VClone_petsc(X);
  Z = N_VClone_petsc(X);

  /* Skipped tests */
  /* fails += Test_N_VSetArrayPointer(W, local_length, myid); */
  /* fails += Test_N_VGetArrayPointer(X, local_length, myid); */
  
  /* Vector operations tests */
  fails += Test_N_VConst(X, local_length, myid);
  fails += Test_N_VLinearSum(X, Y, Z, local_length, myid);
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

  /* Free vectors */
  N_VDestroy_petsc(W);
  N_VDestroy_petsc(X);
  N_VDestroy_petsc(Y);
  N_VDestroy_petsc(Z);

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests, Proc %d \n \n", fails, myid);
  } else {
     if(myid == 0) {
       printf("SUCCESS: NVector module passed all tests, Proc %d \n \n",myid);
     }
  }
  
  ierr = PetscFinalize();
  CHKERRQ(ierr);
  MPI_Finalize();
  return(0);
}

/* ----------------------------------------------------------------------
 * Check vector
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, long int local_length)
{
  int failure = 0;
  long int i;
  Vec *xv = N_VGetVector_petsc(X);
  PetscScalar *a;

  failure = 0;
  /* check PETSc vector */
  VecGetArray(*xv, &a);
  for (i = 0; i < local_length; ++i){
    failure += FNEQ(a[i], ans);
  }
  VecRestoreArray(*xv, &a);

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(N_Vector X)
{
  if(N_VGetVector_petsc(X) == NULL)
    return FALSE;
  else
    return TRUE;
}

void set_element(N_Vector X, long int i, realtype val)
{
  PetscScalar *a;
  Vec *xv = N_VGetVector_petsc(X);
  
  VecGetArray(*xv, &a);
  a[i] = val;
  VecRestoreArray(*xv, &a);
}

realtype get_element(N_Vector X, long int i)
{
  PetscScalar *a;
  Vec *xv = N_VGetVector_petsc(X);
  realtype val;
  
  VecGetArray(*xv, &a);
  val = a[i];
  VecRestoreArray(*xv, &a);
  
  return val;    
}

/* ----------------------------------------------------------------------
 * N_VMake Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
static int Test_N_VMake(Vec* W, int myid)
{
  int failure;
  /* double   start_time, stop_time; */
  N_Vector X;
  int local_length = 0;

  VecGetLocalSize(*W, &local_length);

  /* clone vector */
  /* start_time = get_time(); */  
  X = N_VMake_petsc(W);
  /* stop_time = get_time();  */

  /* check cloned vector */
  if (X == NULL) {
    printf(">>> FAILED test -- N_VMake, Proc %d \n", myid);
    printf("    After N_VMakeEmpty, X == NULL \n \n");
    return(1);
  } 

  /* check cloned vector data */
  if (!has_data(X)) {
    printf(">>> FAILED test -- N_VMake, Proc %d \n", myid);
    printf("    Vector data == NULL \n \n");
    N_VDestroy(X);
    return(1);
  }    

  N_VConst(ONE,X);
  failure = check_ans(ONE, X, local_length);
  if (failure) {
    printf(">>> FAILED test -- N_VMake, Proc %d \n", myid);
    printf("    Failed N_VConst check \n \n");
    N_VDestroy(X);
    return(1);
  }    

  N_VDestroy(X); 

  if (myid == 0) {
    printf("    PASSED test -- N_VMake \n");
    /* PRINT_TIME("    N_VMake Time: %22.15e \n \n", stop_time - start_time); */
  }

  return(0);
}

