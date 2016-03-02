/*
 * -----------------------------------------------------------------
 * $Revision: 4137 $
 * $Date: 2014-06-15 12:26:15 -0700 (Sun, 15 Jun 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner, Slaven Peles @ LLNL
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

#include <sundials/sundials_types.h>
#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <mpi.h>
#include <omp.h>


/* default local vector length */
#define VECLEN 5000

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
  long int veclen;                      /* vector length             */
  int      print_timing;

  /* check input and set vector length */
  if (argc < 3) {
    //printf("ERROR: TWO (2) arguments required: <vector length> <print timing>\n");
    //return(-1);
    SetTiming(0);
  } else {
   print_timing = atoi(argv[2]);
   SetTiming(print_timing);
  }
  
  if (argc < 2) {
    veclen = VECLEN;
  } else {
    veclen = atol(argv[1]); 
    if (veclen <= 0) {
      printf("ERROR: length of vector must be a positive integer \n");
      return(-1); 
    }
  }

//   print_timing = atoi(argv[2]);
//   SetTiming(print_timing);
// 
//   printf("\nRunning with vector length %ld \n \n", veclen);
  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  /* omp_set_num_threads(4); */
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /* set local and global lengths */
  local_length = veclen;
  global_length = nprocs*local_length;

  /* Create vectors */
  W = N_VNewEmpty_ParHyp(comm, local_length, global_length);
  X = N_VNew_ParHyp(comm, local_length, global_length);
  Y = N_VNew_ParHyp(comm, local_length, global_length);
  Z = N_VNew_ParHyp(comm, local_length, global_length);

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
  fails += Test_N_VCloneEmpty(X, myid);
  fails += Test_N_VClone(X, local_length, myid);
  fails += Test_N_VCloneVectorArray(5, X, local_length, myid);
  fails += Test_N_VCloneEmptyVectorArray(5, X, myid);

  /* Free vectors (Can't free W due to a bug in the code) */
  //N_VDestroy_ParHyp(W);
  N_VDestroy_ParHyp(X);
  N_VDestroy_ParHyp(Y);
  N_VDestroy_ParHyp(Z);

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests, Proc %d \n \n", fails, myid);
  } else {
     if(myid == 0) {
       printf("SUCCESS: NVector module passed all tests, Proc %d \n \n",myid);
     }
  }
  
  MPI_Finalize();
  return(0);
}

/* ----------------------------------------------------------------------
 * Check vector
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, long int local_length)
{
  int      failure = 0;
  long int i;
  realtype *Xdata;
  
  Xdata = N_VGetArrayPointer(X);

  /* check vector data */
  for(i=0; i < local_length; i++){
    failure += FNEQ(Xdata[i], ans);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(N_Vector X)
{
  realtype *Xdata = N_VGetArrayPointer(X);
  if (Xdata == NULL)
    return FALSE;
  else
    return TRUE;
}


void set_element(N_Vector X, long int i, realtype val)
{
  NV_Ith_PH(X,i) = val;
}
 
realtype get_element(N_Vector X, long int i)
{
  return NV_Ith_PH(X,i);
}

