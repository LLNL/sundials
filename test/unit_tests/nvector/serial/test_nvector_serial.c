/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the NVECTOR Serial module
 * implementation.
 * -----------------------------------------------------------------*/

#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "test_nvector.h"
#include "test_nvector_complex.h"

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails = 0;             /* counter for test failures */
  int retval;                /* function return value     */
  sunindextype length;       /* vector length             */
  N_Vector U, V, W, X, Y, Z; /* test vectors              */
  int print_timing;          /* turn timing on/off        */

  Test_Init(SUN_COMM_NULL);

  /* check input and set vector length */
  if (argc < 3)
  {
    printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    Test_Finalize();
    return (-1);
  }

  length = (sunindextype)atol(argv[1]);
  if (length <= 0)
  {
    printf("ERROR: length of vector must be a positive integer \n");
    Test_Finalize();
    return (-1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, 0);

  printf("Testing serial N_Vector \n");
  printf("Vector length %ld \n", (long int)length);

  /* Create new vectors */
  W = N_VNewEmpty_Serial(length, sunctx);
  if (W == NULL)
  {
    printf("FAIL: Unable to create a new empty vector \n\n");
    Test_Finalize();
    return (1);
  }

  X = N_VNew_Serial(length, sunctx);
  if (X == NULL)
  {
    N_VDestroy(W);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Finalize();
    return (1);
  }

  /* Check vector ID */
  fails += Test_N_VGetVectorID_Z(X, SUNDIALS_NVEC_SERIAL, 0);

  /* Check vector length */
  fails += Test_N_VGetLength_Z(X, 0);

  /* Check vector communicator */
  fails += Test_N_VGetCommunicator_Z(X, SUN_COMM_NULL, 0);

  /* Test clone functions */
  fails += Test_N_VCloneEmpty_Z(X, 0);
  fails += Test_N_VClone_Z(X, length, 0);
  fails += Test_N_VCloneEmptyVectorArray_Z(5, X, 0);
  fails += Test_N_VCloneVectorArray_Z(5, X, length, 0);

  /* Test setting/getting array data */
  fails += Test_N_VSetArrayPointer_Z(W, length, 0);
  fails += Test_N_VGetArrayPointer_Z(X, length, 0);

  /* Clone additional vectors for testing */
  Y = N_VClone(X);
  if (Y == NULL)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Finalize();
    return (1);
  }

  Z = N_VClone(X);
  if (Z == NULL)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Finalize();
    return (1);
  }

  /* Standard vector operation tests */
  printf("\nTesting standard vector operations:\n\n");

  fails += Test_N_VConst_Z(X, length, 0);
  fails += Test_N_VLinearSum_Z(X, Y, Z, length, 0);
  fails += Test_N_VProd_Z(X, Y, Z, length, 0);
  fails += Test_N_VDiv_Z(X, Y, Z, length, 0);
  fails += Test_N_VScale_Z(X, Z, length, 0);
  fails += Test_N_VAbs_Z(X, Z, length, 0);
  fails += Test_N_VInv_Z(X, Z, length, 0);
  fails += Test_N_VAddConst_Z(X, Z, length, 0);
  fails += Test_N_VDotProd_Z(X, Y, length, 0);
  fails += Test_N_VMaxNorm_Z(X, length, 0);
  fails += Test_N_VWrmsNorm_Z(X, Y, length, 0);
  fails += Test_N_VWrmsNormMask_Z(X, Y, Z, length, 0);
  fails += Test_N_VMin_Z(X, length, 0);
  fails += Test_N_VWL2Norm_Z(X, Y, length, 0);
  fails += Test_N_VL1Norm_Z(X, length, 0);
  fails += Test_N_VCompare_Z(X, Z, length, 0);
  fails += Test_N_VInvTest_Z(X, Z, length, 0);
  fails += Test_N_VConstrMask_Z(X, Y, Z, length, 0);
  fails += Test_N_VMinQuotient_Z(X, Y, length, 0);

  /* Fused and vector array operations tests (disabled) */
  printf("\nTesting fused and vector array operations (disabled):\n\n");

  /* create vector and disable all fused and vector array operations */
  U      = N_VNew_Serial(length, sunctx);
  retval = N_VEnableFusedOps_Serial(U, SUNFALSE);
  if (U == NULL || retval != 0)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Finalize();
    return (1);
  }

  /* fused operations */
  fails += Test_N_VLinearCombination_Z(U, length, 0);
  fails += Test_N_VScaleAddMulti_Z(U, length, 0);
  fails += Test_N_VDotProdMulti_Z(U, length, 0);

  /* vector array operations */
  fails += Test_N_VLinearSumVectorArray_Z(U, length, 0);
  fails += Test_N_VScaleVectorArray_Z(U, length, 0);
  fails += Test_N_VConstVectorArray_Z(U, length, 0);
  fails += Test_N_VWrmsNormVectorArray_Z(U, length, 0);
  fails += Test_N_VWrmsNormMaskVectorArray_Z(U, length, 0);
  fails += Test_N_VScaleAddMultiVectorArray_Z(U, length, 0);
  fails += Test_N_VLinearCombinationVectorArray_Z(U, length, 0);

  /* Fused and vector array operations tests (enabled) */
  printf("\nTesting fused and vector array operations (enabled):\n\n");

  /* create vector and enable all fused and vector array operations */
  V      = N_VNew_Serial(length, sunctx);
  retval = N_VEnableFusedOps_Serial(V, SUNTRUE);
  if (V == NULL || retval != 0)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(U);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Finalize();
    return (1);
  }

  /* fused operations */
  fails += Test_N_VLinearCombination_Z(V, length, 0);
  fails += Test_N_VScaleAddMulti_Z(V, length, 0);
  fails += Test_N_VDotProdMulti_Z(V, length, 0);

  /* vector array operations */
  fails += Test_N_VLinearSumVectorArray_Z(V, length, 0);
  fails += Test_N_VScaleVectorArray_Z(V, length, 0);
  fails += Test_N_VConstVectorArray_Z(V, length, 0);
  fails += Test_N_VWrmsNormVectorArray_Z(V, length, 0);
  fails += Test_N_VWrmsNormMaskVectorArray_Z(V, length, 0);
  fails += Test_N_VScaleAddMultiVectorArray_Z(V, length, 0);
  fails += Test_N_VLinearCombinationVectorArray_Z(V, length, 0);

  /* local reduction operations */
  printf("\nTesting local reduction operations:\n\n");

  fails += Test_N_VDotProdLocal_Z(X, Y, length, 0);
  fails += Test_N_VMaxNormLocal_Z(X, length, 0);
  fails += Test_N_VMinLocal_Z(X, length, 0);
  fails += Test_N_VL1NormLocal_Z(X, length, 0);
  fails += Test_N_VWSqrSumLocal_Z(X, Y, length, 0);
  fails += Test_N_VWSqrSumMaskLocal_Z(X, Y, Z, length, 0);
  fails += Test_N_VInvTestLocal_Z(X, Z, length, 0);
  fails += Test_N_VConstrMaskLocal_Z(X, Y, Z, length, 0);
  fails += Test_N_VMinQuotientLocal_Z(X, Y, length, 0);

  /* local fused reduction operations */
  printf("\nTesting local fused reduction operations:\n\n");
  fails += Test_N_VDotProdMultiLocal_Z(V, length, 0);

  /* XBraid interface operations */
  printf("\nTesting XBraid interface operations:\n\n");

  fails += Test_N_VBufSize_Z(X, length, 0);
  fails += Test_N_VBufPack_Z(X, length, 0);
  fails += Test_N_VBufUnpack_Z(X, length, 0);

  /* Free vectors */
  N_VDestroy(W);
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);
  N_VDestroy(U);
  N_VDestroy(V);

  /* Print result */
  if (fails) { printf("FAIL: NVector module failed %i tests \n\n", fails); }
  else { printf("SUCCESS: NVector module passed all tests \n\n"); }

  Test_Finalize();
  return (fails);
}

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/
int check_ans(sunrealtype ans, N_Vector X, sunindextype local_length)
{
  return check_ans_Z(ans, X, local_length);
}

int check_ans_Z(sunscalartype ans, N_Vector X, sunindextype local_length)
{
  int failure = 0;
  sunindextype i;
  sunscalartype* Xdata;

  Xdata = N_VGetArrayPointer(X);

  /* check vector data */
  for (i = 0; i < local_length; i++) { failure += SUNCompare(Xdata[i], ans); }

  return (failure > ZERO) ? (1) : (0);
}

sunbooleantype has_data(N_Vector X)
{
  /* check if data array is non-null */
  return (N_VGetArrayPointer(X) == NULL) ? SUNFALSE : SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, sunrealtype val)
{
  /* set i-th element of data array */
  set_element_range(X, i, i, val);
}

void set_element_Z(N_Vector X, sunindextype i, sunscalartype val)
{
  /* set i-th element of data array */
  set_element_range_Z(X, i, i, val);
}

void set_element_range(N_Vector X, sunindextype is, sunindextype ie,
                       sunrealtype val)
{
  return set_element_range_Z(X, is, ie, val);
}

void set_element_range_Z(N_Vector X, sunindextype is, sunindextype ie,
                         sunscalartype val)
{
  sunindextype i;

  /* set elements [is,ie] of the data array */
  sunscalartype* xd = N_VGetArrayPointer(X);
  for (i = is; i <= ie; i++) { xd[i] = val; }
}

sunrealtype get_element(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  return NV_Ith_S(X, i);
}

double max_time(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return (time);
}

void sync_device(N_Vector x)
{
  /* not running on GPU, just return */
  return;
}
