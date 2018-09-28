/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the testing routine to check the NVECTOR Raja module
 * implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_mpiraja.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int           fails = 0;    /* counter for test failures  */
  N_Vector      W, X, Y, Z;   /* test vectors               */
  MPI_Comm      comm;
  int           nprocs, myid; /* Number of procs, proc id  */
  sunindextype  global_length, local_length;

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /* set local and global lengths */
  local_length  = 1000;
  global_length = nprocs*local_length;

  if (myid == 0)
    printf("\nRunning with vector length %ld \n\n", (long) global_length);

  /* Create vectors */
  W = N_VNewEmpty_Raja(local_length);
  X = N_VNew_Raja(comm, local_length, global_length);

  /* NVector Tests */

  /* Raja specific tests */

  /* Memory allocation tests */
  fails += Test_N_VCloneEmpty(X, myid);
  fails += Test_N_VClone(X, local_length, myid);
  fails += Test_N_VCloneEmptyVectorArray(5, X, myid);
  fails += Test_N_VCloneVectorArray(5, X, local_length, myid);

  Y = N_VClone_Raja(X);
  Z = N_VClone_Raja(X);

  /* Skipped tests */
  /*   fails += Test_N_VSetArrayPointer(W, local_length, myid); */
  /*   fails += Test_N_VGetArrayPointer(X, local_length, myid); */

  /* Vector operation tests */
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
  N_VDestroy(W);
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);

  /* Print result */
  if (myid == 0) {
    if (fails) {
      printf("FAIL: NVector module failed %i tests \n \n", fails);
    } else {
      printf("SUCCESS: NVector module passed all tests \n \n");
    }
  }

  MPI_Finalize();

  return(fails);
}

/* ----------------------------------------------------------------------
 * Check vector
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  int failure = 0;
  sunindextype i;
  realtype *xdata;

  N_VCopyFromDevice_Raja(X);

  xdata = N_VGetHostArrayPointer_Raja(X);
  /* check vector data */
  for(i=0; i < local_length; i++){
    failure += FNEQ(xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  return (X->content == NULL ? SUNFALSE : SUNTRUE);
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  N_VCopyFromDevice_Raja(X);
  (N_VGetHostArrayPointer_Raja(X))[i] = val;
  N_VCopyToDevice_Raja(X);
}

realtype get_element(N_Vector X, sunindextype i)
{
  N_VCopyFromDevice_Raja(X);
  return (N_VGetHostArrayPointer_Raja(X))[i];
}
