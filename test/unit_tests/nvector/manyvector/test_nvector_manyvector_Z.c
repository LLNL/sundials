/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * This is the testing routine to check the NVECTOR ManyVector
 * (serial) module implementation.
 * -----------------------------------------------------------------*/

#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "test_nvector_complex.h"

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails = 0;             /* counter for test failures */
  int retval;                /* function return value     */
  sunindextype len1, len2;   /* vector lengths            */
  sunindextype length;       /* overall vector length     */
  N_Vector Xsub[2];          /* subvector pointer array   */
  N_Vector U, V, W, X, Y, Z; /* test vectors              */
  int print_timing;          /* turn timing on/off        */

  Test_Init_Z(SUN_COMM_NULL);

  /* check input and set vector length */
  if (argc < 3)
  {
    printf("ERROR: TWO (2) Inputs required: subvector 1 length, subvector 2 "
           "length\n");
    Test_Abort_Z(1);
  }

  len1 = (sunindextype)atol(argv[1]);
  if (len1 <= 0)
  {
    printf("ERROR: length of subvector 1 must be a positive integer \n");
    Test_Abort_Z(1);
  }

  len2 = (sunindextype)atol(argv[2]);
  if (len2 <= 0)
  {
    printf("ERROR: length of subvector 2 must be a positive integer \n");
    Test_Abort_Z(1);
  }

  /* overall length */
  length = len1 + len2;

  printf("Testing ManyVector (serial) N_Vector \n");
  printf("Vector lengths: %ld %ld \n", (long int)len1, (long int)len2);

  /* Create subvectors */
  Xsub[0] = N_VNew_Serial(len1, sunctx);
  if (Xsub[0] == NULL)
  {
    printf("FAIL: Unable to create a new serial subvector \n\n");
    Test_Abort_Z(1);
  }
  Xsub[1] = N_VNew_Serial(len2, sunctx);
  if (Xsub[1] == NULL)
  {
    N_VDestroy(Xsub[0]);
    printf("FAIL: Unable to create a new serial subvector \n\n");
    Test_Abort_Z(1);
  }

  /* Create a new ManyVector */
  X = N_VNew_ManyVector(2, Xsub, sunctx);

  /* Check vector ID */
  fails += Test_N_VGetVectorID_Z(X, SUNDIALS_NVEC_MANYVECTOR, 0);

  /* Check vector length */
  fails += Test_N_VGetLength_Z(X, 0);

  /* Check vector communicator */
  fails += Test_N_VGetCommunicator_Z(X, SUN_COMM_NULL, 0);

  /* Test subvector accessors */
  if (N_VGetNumSubvectors_ManyVector(X) != 2)
  {
    printf(">>> FAILED test -- N_VGetNumSubvectors_ManyVector\n");
    fails += 1;
  }
  U = N_VGetSubvector_ManyVector(X, 0);
  if (N_VGetLength(U) != len1)
  {
    printf(">>> FAILED test -- N_VGetSubvector_ManyVector\n");
    fails += 1;
  }
  U = N_VGetSubvector_ManyVector(X, 1);
  if (N_VGetLength(U) != len2)
  {
    printf(">>> FAILED test -- N_VGetSubvector_ManyVector\n");
    fails += 1;
  }

  /* Clone additional vectors for testing */
  W = N_VClone(X);
  if (W == NULL)
  {
    N_VDestroy(X);
    printf("FAIL: Unable to create a new vector\n\n");
    Test_Abort_Z(1);
  }

  Y = N_VClone(X);
  if (Y == NULL)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    printf("FAIL: Unable to create a new vector\n\n");
    Test_Abort_Z(1);
  }

  Z = N_VClone(X);
  if (Z == NULL)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    printf("FAIL: Unable to create a new vector\n\n");
    Test_Abort_Z(1);
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
  U      = N_VClone(X);
  retval = N_VEnableFusedOps_ManyVector(U, SUNFALSE);
  if (U == NULL || retval != 0)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Abort_Z(1);
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
  V      = N_VClone(X);
  retval = N_VEnableFusedOps_ManyVector(V, SUNTRUE);
  if (V == NULL || retval != 0)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(U);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Abort_Z(1);
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
  N_VDestroy(Xsub[0]);
  N_VDestroy(Xsub[1]);

  /* Print result */
  if (fails) { printf("FAIL: NVector module failed %i tests \n\n", fails); }
  else { printf("SUCCESS: NVector module passed all tests \n\n"); }

  Test_Finalize_Z();
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
  N_Vector Xsub[2];
  sunscalartype *x0, *x1;
  sunindextype x0len, x1len;

  Xsub[0] = N_VGetSubvector_ManyVector(X, 0);
  Xsub[1] = N_VGetSubvector_ManyVector(X, 1);
  x0len   = N_VGetLength(Xsub[0]);
  x1len   = N_VGetLength(Xsub[1]);
  x0      = N_VGetSubvectorArrayPointer_ManyVector(X, 0);
  x1      = N_VGetSubvectorArrayPointer_ManyVector(X, 1);

  /* ensure that local_length = x0len + x1len */
  if (local_length != x0len + x1len) { return (1); }

  /* check vector data */
  for (i = 0; i < x0len; i++) { failure += SUNCompare(x0[i], ans); }
  for (i = 0; i < x1len; i++) { failure += SUNCompare(x1[i], ans); }

  return (failure > ZERO) ? (1) : (0);
}

sunbooleantype has_data_Z(N_Vector X)
{
  /* should not be called in these tests */
  return SUNTRUE;
}

void set_element_Z(N_Vector X, sunindextype i, sunscalartype val)
{
  N_Vector Xsub[2];
  sunindextype x0len;

  Xsub[0] = N_VGetSubvector_ManyVector(X, 0);
  Xsub[1] = N_VGetSubvector_ManyVector(X, 1);
  x0len   = N_VGetLength(Xsub[0]);

  /* set i-th element of data array (in appropriate subvector) */
  if (i < x0len) { NV_Ith_S(Xsub[0], i) = val; }
  else { NV_Ith_S(Xsub[1], i - x0len) = val; }
}

void set_element_range_Z(N_Vector X, sunindextype is, sunindextype ie,
                         sunscalartype val)
{
  N_Vector Xsub[2];
  sunindextype x0len, i;

  Xsub[0] = N_VGetSubvector_ManyVector(X, 0);
  Xsub[1] = N_VGetSubvector_ManyVector(X, 1);
  x0len   = N_VGetLength(Xsub[0]);

  /* set i-th element of data array (in appropriate subvector) */
  for (i = is; i < x0len; i++) { NV_Ith_S(Xsub[0], i) = val; }
  for (i = x0len; i <= ie; i++) { NV_Ith_S(Xsub[1], i - x0len) = val; }
}

sunscalartype get_element_Z(N_Vector X, sunindextype i)
{
  N_Vector Xsub[2];
  sunindextype x0len;

  Xsub[0] = N_VGetSubvector_ManyVector(X, 0);
  Xsub[1] = N_VGetSubvector_ManyVector(X, 1);
  x0len   = N_VGetLength(Xsub[0]);

  /* get i-th element of data array (from appropriate subvector) */
  if (i < x0len) { return NV_Ith_S(Xsub[0], i); }
  else { return NV_Ith_S(Xsub[1], i - x0len); }
}

double max_time_Z(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return (time);
}

void sync_device_Z(N_Vector x)
{
  /* not running on GPU, just return */
  return;
}
