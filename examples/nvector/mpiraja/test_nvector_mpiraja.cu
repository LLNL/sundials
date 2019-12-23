/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the NVECTOR Raja module
 * implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_mpiplusx.h>
#include <nvector/nvector_raja.h>
#include <nvector/raja/Vector.hpp>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <mpi.h>

using namespace sunrajavec;

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails = 0;           /* counter for test failures */
  int          globfails = 0;       /* counter for test failures */
  int          retval;              /* function return value     */
  sunindextype local_length;        /* local vector length       */
  sunindextype global_length;       /* global vector length      */
  N_Vector     U, V, X;             /* local test vectors        */
  N_Vector     plusU, plusV, plusX; /* MPIPlusX test vectors     */
  N_Vector     plusY, plusZ;        /* MPIPlusX test vectors     */
  int          print_timing;        /* turn timing on/off        */
  MPI_Comm     comm;                /* MPI Communicator          */
  int          nprocs, myid;        /* Number of procs, proc id  */

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /* check inputs */
  if (argc < 3) {
    if (myid == 0)
      printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    MPI_Abort(comm, -1);
  }

  local_length = (sunindextype) atol(argv[1]);
  if (local_length < 1) {
    if (myid == 0)
      printf("ERROR: local vector length must be a positive integer \n");
    MPI_Abort(comm, -1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, myid);

  /* global length */
  global_length = nprocs*local_length;

  if (myid == 0) {
    printf("Testing the RAJA N_Vector \n");
    printf("Vector global length %ld \n", (long int) global_length);
    printf("MPI processes %d \n", nprocs);
  }

  /* Create new vectors */
  X = N_VNew_Raja(local_length);
  if (X == NULL) {
    if (myid == 0) printf("FAIL: Unable to create a new RAJA vector \n\n");
    MPI_Abort(comm, 1);
  }

  /* Create the MPI+X vector */
  plusX = N_VMake_MPIPlusX(comm, X);
  if (plusX == NULL) {
    N_VDestroy(X);
    if (myid == 0) printf("FAIL: Unable to create a new MPIPlusX vector \n\n");
    MPI_Abort(comm, 1);
  }

  /* Check vector ID */
  fails += Test_N_VGetVectorID(plusX, SUNDIALS_NVEC_MPIPLUSX, myid);

  /* Check vector length */
  fails += Test_N_VGetLength(plusX, myid);

  /* Check vector communicator */
  fails += Test_N_VGetCommunicator(plusX, &comm, myid);

  /* Test clone functions */
  fails += Test_N_VCloneEmpty(plusX, myid);
  fails += Test_N_VClone(plusX, local_length, myid);
  fails += Test_N_VCloneEmptyVectorArray(5, plusX, myid);
  fails += Test_N_VCloneVectorArray(5, plusX, local_length, myid);

  /* Clone additional vectors for testing */
  plusY = N_VClone(plusX);
  if (plusY == NULL) {
    N_VDestroy(X);
    N_VDestroy(plusX);
    if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
    MPI_Abort(comm, 1);
  }

  plusZ = N_VClone(plusX);
  if (plusZ == NULL) {
    N_VDestroy(X);
    N_VDestroy(plusX);
    N_VDestroy(plusY);
    if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
    MPI_Abort(comm, 1);
  }

  /* Standard vector operation tests */
  if (myid == 0) printf("\nTesting standard vector operations:\n\n");

  fails += Test_N_VConst(plusX, local_length, myid);
  fails += Test_N_VLinearSum(plusX, plusY, plusZ, local_length, myid);
  fails += Test_N_VProd(plusX, plusY, plusZ, local_length, myid);
  fails += Test_N_VDiv(plusX, plusY, plusZ, local_length, myid);
  fails += Test_N_VScale(plusX, plusZ, local_length, myid);
  fails += Test_N_VAbs(plusX, plusZ, local_length, myid);
  fails += Test_N_VInv(plusX, plusZ, local_length, myid);
  fails += Test_N_VAddConst(plusX, plusZ, local_length, myid);
  fails += Test_N_VDotProd(plusX, plusY, local_length, myid);
  fails += Test_N_VMaxNorm(plusX, local_length, myid);
  fails += Test_N_VWrmsNorm(plusX, plusY, local_length, myid);
  fails += Test_N_VWrmsNormMask(plusX, plusY, plusZ, local_length, myid);
  fails += Test_N_VMin(plusX, local_length, myid);
  fails += Test_N_VWL2Norm(plusX, plusY, local_length, myid);
  fails += Test_N_VL1Norm(plusX, local_length, myid);
  fails += Test_N_VCompare(plusX, plusZ, local_length, myid);
  fails += Test_N_VInvTest(plusX, plusZ, local_length, myid);
  fails += Test_N_VConstrMask(plusX, plusY, plusZ, local_length, myid);
  fails += Test_N_VMinQuotient(plusX, plusY, local_length, myid);

  /* Fused and vector array operations tests (disabled) */
  if (myid == 0) printf("\nTesting fused and vector array operations (disabled):\n\n");

  /* create vector and disable all fused and vector array operations */
  U = N_VNew_Raja(local_length);
  retval = N_VEnableFusedOps_Raja(U, SUNFALSE);
  if (U == NULL || retval != 0) {
    N_VDestroy(X);
    N_VDestroy(plusX);
    N_VDestroy(plusY);
    N_VDestroy(plusZ);
    if (myid == 0) printf("FAIL: Unable to create a new RAJA vector \n\n");
    MPI_Abort(comm, 1);
  }

  plusU = N_VMake_MPIPlusX(comm, U);
  if (U == NULL || retval != 0) {
    N_VDestroy(X);
    N_VDestroy(U);
    N_VDestroy(plusX);
    N_VDestroy(plusY);
    N_VDestroy(plusZ);
    if (myid == 0) printf("FAIL: Unable to create a new MPIPlusX vector \n\n");
    MPI_Abort(comm, 1);
  }

  /* fused operations */
  fails += Test_N_VLinearCombination(plusU, local_length, myid);
  fails += Test_N_VScaleAddMulti(plusU, local_length, myid);
  fails += Test_N_VDotProdMulti(plusU, local_length, myid);

  /* vector array operations */
  fails += Test_N_VLinearSumVectorArray(plusU, local_length, myid);
  fails += Test_N_VScaleVectorArray(plusU, local_length, myid);
  fails += Test_N_VConstVectorArray(plusU, local_length, myid);
  fails += Test_N_VWrmsNormVectorArray(plusU, local_length, myid);
  fails += Test_N_VWrmsNormMaskVectorArray(plusU, local_length, myid);
  fails += Test_N_VScaleAddMultiVectorArray(plusU, local_length, myid);
  fails += Test_N_VLinearCombinationVectorArray(plusU, local_length, myid);

  /* Fused and vector array operations tests (enabled) */
  if (myid == 0) printf("\nTesting fused and vector array operations (enabled):\n\n");

  /* create vector and enable all fused and vector array operations */
  V = N_VNew_Raja(local_length);
  retval = N_VEnableFusedOps_Raja(V, SUNTRUE);
  if (V == NULL || retval != 0) {
    N_VDestroy(X);
    N_VDestroy(U);
    N_VDestroy(plusX);
    N_VDestroy(plusY);
    N_VDestroy(plusZ);
    N_VDestroy(plusU);
    if (myid == 0) printf("FAIL: Unable to create a new RAJA vector \n\n");
    MPI_Abort(comm, 1);
  }

  /* create the MPIPlusX vector */
  plusV = N_VMake_MPIPlusX(comm, V);
  if (V == NULL || retval != 0) {
    N_VDestroy(X);
    N_VDestroy(U);
    N_VDestroy(V);
    N_VDestroy(plusU);
    N_VDestroy(plusX);
    N_VDestroy(plusY);
    N_VDestroy(plusZ);
    if (myid == 0) printf("FAIL: Unable to create a new MPIPlusX vector \n\n");
    MPI_Abort(comm, 1);
  }

  /* fused operations */
  fails += Test_N_VLinearCombination(plusV, local_length, myid);
  fails += Test_N_VScaleAddMulti(plusV, local_length, myid);
  fails += Test_N_VDotProdMulti(plusV, local_length, myid);

  /* vector array operations */
  fails += Test_N_VLinearSumVectorArray(plusV, local_length, myid);
  fails += Test_N_VScaleVectorArray(plusV, local_length, myid);
  fails += Test_N_VConstVectorArray(plusV, local_length, myid);
  fails += Test_N_VWrmsNormVectorArray(plusV, local_length, myid);
  fails += Test_N_VWrmsNormMaskVectorArray(plusV, local_length, myid);
  fails += Test_N_VScaleAddMultiVectorArray(plusV, local_length, myid);
  fails += Test_N_VLinearCombinationVectorArray(plusV, local_length, myid);

  /* local reduction operations */
  printf("\nTesting local reduction operations:\n\n");

  fails += Test_N_VDotProdLocal(plusX, plusY, local_length, myid);
  fails += Test_N_VMaxNormLocal(plusX, local_length, myid);
  fails += Test_N_VMinLocal(plusX, local_length, myid);
  fails += Test_N_VL1NormLocal(plusX, local_length, myid);
  fails += Test_N_VWSqrSumLocal(plusX, plusY, local_length, myid);
  fails += Test_N_VWSqrSumMaskLocal(plusX, plusY, plusZ, local_length, myid);
  fails += Test_N_VInvTestLocal(plusX, plusZ, local_length, myid);
  fails += Test_N_VConstrMaskLocal(plusX, plusY, plusZ, local_length, myid);
  fails += Test_N_VMinQuotientLocal(plusX, plusY, local_length, myid);

  /* Free vectors */
  N_VDestroy(X);
  N_VDestroy(U);
  N_VDestroy(V);
  N_VDestroy(plusX);
  N_VDestroy(plusY);
  N_VDestroy(plusZ);
  N_VDestroy(plusU);
  N_VDestroy(plusV);

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests, Proc %d \n\n", fails, myid);
  } else {
    if (myid == 0)
      printf("SUCCESS: NVector module passed all tests \n\n");
  }

  /* check if any other process failed */
  (void) MPI_Allreduce(&fails, &globfails, 1, MPI_INT, MPI_MAX, comm);

  MPI_Finalize();

  return(globfails);
}

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector plusX, sunindextype local_length)
{
  int          failure = 0;
  sunindextype i;
  realtype     *Xdata;
  N_Vector     X;

  X = N_VGetLocalVector_MPIPlusX(plusX);
  N_VCopyFromDevice_Raja(X);
  Xdata = N_VGetHostArrayPointer_Raja(X);

  /* check vector data */
  for (i = 0; i < local_length; i++) {
    failure += FNEQ(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector plusX)
{
  return (N_VGetLocalVector_MPIPlusX(plusX)->content == NULL) ? SUNFALSE : SUNTRUE;
}

void set_element(N_Vector plusX, sunindextype i, realtype val)
{
  /* set i-th element of data array */
  set_element_range(plusX, i, i, val);
}

void set_element_range(N_Vector plusX, sunindextype is, sunindextype ie,
                       realtype val)
{
  sunindextype i;
  realtype*    xd;
  N_Vector     X;

  X = N_VGetLocalVector_MPIPlusX(plusX);

  /* set elements [is,ie] of the data array */
  N_VCopyFromDevice_Raja(X);
  xd = N_VGetHostArrayPointer_Raja(X);
  for(i = is; i <= ie; i++) xd[i] = val;
  N_VCopyToDevice_Raja(X);
}

realtype get_element(N_Vector plusX, sunindextype i)
{
  N_Vector X = N_VGetLocalVector_MPIPlusX(plusX);

  /* get i-th element of data array */
  N_VCopyFromDevice_Raja(X);
  return (N_VGetHostArrayPointer_Raja(X))[i];
}

double max_time(N_Vector plusX, double time)
{
  MPI_Comm *comm;
  double maxt;

  comm = (MPI_Comm*) N_VGetCommunicator(plusX);

  /* get max time across all MPI ranks */
  (void) MPI_Reduce(&time, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, *comm);
  return(maxt);
}

void sync_device()
{
  /* sync with GPU */
  cudaDeviceSynchronize();
  return;
}
