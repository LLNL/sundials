/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, Daniel McGreer @ LLNL
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
 * This is the testing routine to check the NVECTOR Kokkos module
 * implementation.
 * -----------------------------------------------------------------*/

#include <nvector/nvector_kokkos.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "test_nvector.h"

/* KOKKOS vector specific tests */
static int Test_N_VMake_Kokkos(N_Vector X, sunindextype length, int myid);

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails = 0;          /* counter for test failures */
  int retval;             /* function return value     */
  sunindextype length;    /* vector length             */
  N_Vector U, V, X, Y, Z; /* test vectors              */
  int print_timing;       /* turn timing on/off        */
  int memtype;

  Test_Init(NULL);

  /* check input and set vector length */
  if (argc < 3) {
    printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    return (-1);
  }

  length = (sunindextype)atol(argv[1]);
  if (length <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    return (-1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, 0);

  printf("Testing KOKKOS N_Vector \n");

  printf("Vector length %ld \n\n", (long int)length);

  Kokkos::initialize(argc, argv);
  {
    /* Create new vectors */
    X = N_VNew_Kokkos(length);

    if (X == NULL) {
      printf("FAIL: Unable to create a new vector \n\n");
      return (1);
    }

    /* Check vector ID */
    fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_KOKKOS, 0);

    /* Check vector length */
    fails += Test_N_VGetLength(X, 0);

    /* Check vector communicator */
    fails += Test_N_VGetCommunicator(X, NULL, 0);

    /* Test clone functions */
    fails += Test_N_VCloneEmpty(X, 0);
    fails += Test_N_VClone(X, length, 0);
    fails += Test_N_VCloneEmptyVectorArray(5, X, 0);
    fails += Test_N_VCloneVectorArray(5, X, length, 0);

    /* Clone additional vectors for testing */
    Y = N_VClone(X);
    if (Y == NULL) {
      N_VDestroy(X);
      printf("FAIL: Unable to create a new vector \n\n");
      return (1);
    }

    Z = N_VClone(X);
    if (Z == NULL) {
      N_VDestroy(X);
      N_VDestroy(Y);
      printf("FAIL: Unable to create a new vector \n\n");
      return (1);
    }

    /* Standard vector operation tests */
    printf("\nTesting standard vector operations:\n\n");

    fails += Test_N_VConst(X, length, 0);
    fails += Test_N_VLinearSum(X, Y, Z, length, 0);
    fails += Test_N_VProd(X, Y, Z, length, 0);
    fails += Test_N_VDiv(X, Y, Z, length, 0);
    fails += Test_N_VScale(X, Z, length, 0);
    fails += Test_N_VAbs(X, Z, length, 0);
    fails += Test_N_VInv(X, Z, length, 0);
    fails += Test_N_VAddConst(X, Z, length, 0);
    fails += Test_N_VDotProd(X, Y, length, 0);
    fails += Test_N_VMaxNorm(X, length, 0);
    fails += Test_N_VWrmsNorm(X, Y, length, 0);
    fails += Test_N_VWrmsNormMask(X, Y, Z, length, 0);
    fails += Test_N_VMin(X, length, 0);
    fails += Test_N_VWL2Norm(X, Y, length, 0);
    fails += Test_N_VL1Norm(X, length, 0);
    fails += Test_N_VCompare(X, Z, length, 0);
    fails += Test_N_VInvTest(X, Z, length, 0);
    fails += Test_N_VConstrMask(X, Y, Z, length, 0);
    fails += Test_N_VMinQuotient(X, Y, length, 0);

    /* Fused and vector array operations tests (disabled) */
    printf("\nTesting fused and vector array operations (disabled):\n\n");

    /* create vector and disable all fused and vector array operations */
    U = N_VNew_Kokkos(length);

    retval = N_VEnableFusedOps_Kokkos(U, SUNFALSE);
    if (U == NULL || retval != 0) {
      N_VDestroy(X);
      N_VDestroy(Y);
      N_VDestroy(Z);
      printf("FAIL: Unable to create a new vector \n\n");
      return (1);
    }

    /* fused operations */
    fails += Test_N_VLinearCombination(U, length, 0);
    fails += Test_N_VScaleAddMulti(U, length, 0);
    fails += Test_N_VDotProdMulti(U, length, 0);

    /* vector array operations */
    fails += Test_N_VLinearSumVectorArray(U, length, 0);
    fails += Test_N_VScaleVectorArray(U, length, 0);
    fails += Test_N_VConstVectorArray(U, length, 0);
    fails += Test_N_VWrmsNormVectorArray(U, length, 0);
    fails += Test_N_VWrmsNormMaskVectorArray(U, length, 0);
    fails += Test_N_VScaleAddMultiVectorArray(U, length, 0);
    fails += Test_N_VLinearCombinationVectorArray(U, length, 0);

    /* Fused and vector array operations tests (enabled) */
    printf("\nTesting fused and vector array operations (enabled):\n\n");

    /* create vector and enable all fused and vector array operations */
    V = N_VNew_Kokkos(length);

    retval = N_VEnableFusedOps_Kokkos(V, SUNTRUE);
    if (V == NULL || retval != 0) {
      N_VDestroy(X);
      N_VDestroy(Y);
      N_VDestroy(Z);
      N_VDestroy(U);
      printf("FAIL: Unable to create a new vector \n\n");
      return (1);
    }

    /* fused operations */
    fails += Test_N_VLinearCombination(V, length, 0);
    fails += Test_N_VScaleAddMulti(V, length, 0);
    fails += Test_N_VDotProdMulti(V, length, 0);

    /* vector array operations */
    fails += Test_N_VLinearSumVectorArray(V, length, 0);
    fails += Test_N_VScaleVectorArray(V, length, 0);
    fails += Test_N_VConstVectorArray(V, length, 0);
    fails += Test_N_VWrmsNormVectorArray(V, length, 0);
    fails += Test_N_VWrmsNormMaskVectorArray(V, length, 0);
    fails += Test_N_VScaleAddMultiVectorArray(V, length, 0);
    fails += Test_N_VLinearCombinationVectorArray(V, length, 0);

    /* local reduction operations */
    printf("\nTesting local reduction operations:\n\n");

    fails += Test_N_VDotProdLocal(X, Y, length, 0);
    fails += Test_N_VMaxNormLocal(X, length, 0);
    fails += Test_N_VMinLocal(X, length, 0);
    fails += Test_N_VL1NormLocal(X, length, 0);
    fails += Test_N_VWSqrSumLocal(X, Y, length, 0);
    fails += Test_N_VWSqrSumMaskLocal(X, Y, Z, length, 0);
    fails += Test_N_VInvTestLocal(X, Z, length, 0);
    fails += Test_N_VConstrMaskLocal(X, Y, Z, length, 0);
    fails += Test_N_VMinQuotientLocal(X, Y, length, 0);

    /* KOKKOS specific tests */
    printf("\nTesting kokkos vector specific operations:\n\n");
    fails += Test_N_VMake_Kokkos(X, length, 0);

    /* Free vectors */
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(U);
    N_VDestroy(V);
  }
  Kokkos::finalize();

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests \n\n", fails);
  }
  else {
    printf("SUCCESS: NVector module passed all tests \n\n");
  }

  Test_Finalize();

  return (fails);
}

/* ----------------------------------------------------------------------
 * KOKKOS specific tests
 * --------------------------------------------------------------------*/

/* --------------------------------------------------------------------
 * Test for the KOKKOS N_Vector N_VMake_Kokkos function. Requires N_VConst
 * to check data.
 */

int Test_N_VMake_Kokkos(N_Vector X, sunindextype length, int myid)
{
  int failure = 0;
  realtype *h_data, *d_data;
  N_Vector Y;

  N_VConst(NEG_HALF, X);
  N_VCopyFromDevice_Kokkos(X);

  h_data = N_VGetHostArrayPointer_Kokkos(X);
  d_data = N_VGetDeviceArrayPointer_Kokkos(X);

  /* Case 1: h_data and d_data are not null */
  Y = N_VMake_Kokkos(length, h_data, d_data);
  if (Y == NULL) {
    printf(">>> FAILED test -- N_VMake_Kokkos, Proc %d \n", myid);
    printf("    Vector is NULL \n \n");
    return (1);
  }

  if (N_VGetHostArrayPointer_Kokkos(Y) == NULL) {
    printf(">>> FAILED test -- N_VMake_Kokkos, Proc %d \n", myid);
    printf("    Vector host data == NULL \n \n");
    N_VDestroy(Y);
    return (1);
  }

  if (N_VGetDeviceArrayPointer_Kokkos(Y) == NULL) {
    printf(">>> FAILED test -- N_VMake_Kokkos, Proc %d \n", myid);
    printf("    Vector device data -= NULL \n \n");
    N_VDestroy(Y);
    return (1);
  }

  failure += check_ans(NEG_HALF, Y, length);

  if (failure) {
    printf(">>> FAILED test -- N_VMake_Kokkos Case 1, Proc %d \n", myid);
    printf("    Failed N_VConst check \n \n");
    N_VDestroy(Y);
    return (1);
  }

  if (myid == 0) {
    printf("PASSED test -- N_VMake_Kokkos Case 1 \n");
  }

  N_VDestroy(Y);

  /* Case 2: data is null */
  Y = N_VMake_Kokkos(length, NULL, NULL);
  if (Y != NULL) {
    printf(">>> FAILED test -- N_VMake_Kokkos Case 2, Proc %d \n", myid);
    printf("    Vector is not NULL \n \n");
    return (1);
  }

  if (myid == 0) {
    printf("PASSED test -- N_VMake_Kokkos Case 2 \n");
  }

  N_VDestroy(Y);

  return (failure);
}

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/
extern "C" int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  int failure = 0;
  sunindextype i;
  realtype* Xdata;

  N_VCopyFromDevice_Kokkos(X);
  Xdata = N_VGetHostArrayPointer_Kokkos(X);

  /* check vector data */
  for (i = 0; i < local_length; i++) {
    failure += FNEQ(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

extern "C" booleantype has_data(N_Vector X)
{
  /* check if vector data is non-null */
  if ((N_VGetHostArrayPointer_Kokkos(X) == NULL) && (N_VGetDeviceArrayPointer_Kokkos(X) == NULL)) return SUNFALSE;
  return SUNTRUE;
}

extern "C" void set_element(N_Vector X, sunindextype i, realtype val)
{
  /* set i-th element of data array */
  set_element_range(X, i, i, val);
}

extern "C" void set_element_range(N_Vector X, sunindextype is, sunindextype ie, realtype val)
{
  sunindextype i;
  realtype* xd;

  /* set elements [is,ie] of the data array */
  N_VCopyFromDevice_Kokkos(X);
  xd = N_VGetHostArrayPointer_Kokkos(X);
  for (i = is; i <= ie; i++) xd[i] = val;
  N_VCopyToDevice_Kokkos(X);
}

extern "C" realtype get_element(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  N_VCopyFromDevice_Kokkos(X);
  return (N_VGetHostArrayPointer_Kokkos(X))[i];
}

extern "C" double max_time(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return (time);
}

extern "C" void sync_device()
{
  /* sync with GPU */
  Kokkos::fence();
  return;
}
