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
 * This is the testing routine to check the NVECTOR CUDA module
 * implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_cuda.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int           fails = 0;    /* counter for test failures  */
  N_Vector      W, X, Y, Z;   /* test vectors               */
  int           print_timing;
  sunindextype  length;

  /* check input and set vector length */
  if (argc < 3){
    printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    return(-1);
  }

  length = atol(argv[1]);
  if (length <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    return(-1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing);

  printf("\nRunning with vector length %ld \n\n", (long) length);

  /* Create vectors */
  W = N_VNewEmpty_Cuda(length);
  X = N_VNew_Cuda(length);

  /* NVector Tests */

  /* CUDA specific tests */

  /* Memory allocation tests */
  fails += Test_N_VCloneEmpty(X, 0);
  fails += Test_N_VClone(X, length, 0);
  fails += Test_N_VCloneEmptyVectorArray(5, X, 0);
  fails += Test_N_VCloneVectorArray(5, X, length, 0);

  Y = N_VClone_Cuda(X);
  Z = N_VClone_Cuda(X);

  /* Skipped tests */
  /*   fails += Test_N_VSetArrayPointer(W, length, 0); */
  /*   fails += Test_N_VGetArrayPointer(X, length, 0); */

  /* Vector operation tests */
  fails += Test_N_VConst(X, length, 0);
  fails += Test_N_VLinearSum(X, Y, Z, length, 0);
  fails += Test_N_VProd(X, Y, Z, length, 0);
  fails += Test_N_VDiv(X, Y, Z, length, 0);
  fails += Test_N_VScale(X, Z, length, 0);
  fails += Test_N_VAbs(X, Z, length, 0);
  fails += Test_N_VInv(X, Z, length, 0);
  fails += Test_N_VAddConst(X, Z, length, 0);
  fails += Test_N_VDotProd(X, Y, length, length, 0);
  fails += Test_N_VMaxNorm(X, length, 0);
  fails += Test_N_VWrmsNorm(X, Y, length, 0);
  fails += Test_N_VWrmsNormMask(X, Y, Z, length, length, 0);
  fails += Test_N_VMin(X, length, 0);
  fails += Test_N_VWL2Norm(X, Y, length, length, 0);
  fails += Test_N_VL1Norm(X, length, length, 0);
  fails += Test_N_VCompare(X, Z, length, 0);
  fails += Test_N_VInvTest(X, Z, length, 0);
  fails += Test_N_VConstrMask(X, Y, Z, length, 0);
  fails += Test_N_VMinQuotient(X, Y, length, 0);

  /* Free vectors */
  N_VDestroy(W);
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests \n \n", fails);
  } else {
    printf("SUCCESS: NVector module passed all tests \n \n");
  }

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

  N_VCopyFromDevice_Cuda(X);

  xdata = N_VGetHostArrayPointer_Cuda(X);
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
  N_VCopyFromDevice_Cuda(X);
  (N_VGetHostArrayPointer_Cuda(X))[i] = val;
  N_VCopyToDevice_Cuda(X);
}

realtype get_element(N_Vector X, sunindextype i)
{
  N_VCopyFromDevice_Cuda(X);
  return (N_VGetHostArrayPointer_Cuda(X))[i];
}
