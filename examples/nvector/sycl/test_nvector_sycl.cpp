/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the NVECTOR SYCL module
 * implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_sycl.h>

#include "custom_memory_helper_sycl.h"
#include "test_nvector.h"

/* SYCL vector variants */
enum mem_type { UNMANAGED, MANAGED, SUNMEMORY };
enum pol_type { DEFAULT_POL, GRID_STRIDE };

/* --------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int             fails = 0;         /* counter for test failures  */
  int             retval;            /* function return value      */
  sunindextype    length;            /* vector length              */
  N_Vector        U, V, X, Y, Z;     /* test vectors               */
  int             print_timing;      /* turn timing on/off         */
  int             threadsPerBlock;   /* sycl block size            */
  int             memtype, policy;

  Test_Init(NULL);

  /* check input and set vector length */
  if (argc < 4)
  {
    printf("ERROR: THREE (3) Inputs required: vector length, SYCL threads per block (0 for default), print timing \n");
    Test_Abort(-1);
  }

  length = (sunindextype) atol(argv[1]);
  if (length <= 0)
  {
    printf("ERROR: length of vector must be a positive integer\n");
    Test_Abort(-1);
  }

  threadsPerBlock = (int) atoi(argv[2]);
  if (threadsPerBlock < 0)
  {
    printf("ERROR: SYCL threads per block must be a positive value or 0 to use the default\n");
    Test_Abort(-1);
  }

  print_timing = atoi(argv[3]);
  SetTiming(print_timing, 0);

  /* Create an in-order GPU queue */
#if SYCL_LANGUAGE_VERSION >= 2020
  sycl::queue myQueue(sycl::gpu_selector_v,
                      sycl::property_list{sycl::property::queue::in_order{}});
#else
  sycl::gpu_selector selector;
  sycl::queue myQueue(selector,
                      sycl::property_list{sycl::property::queue::in_order{}});
#endif

  sycl::device dev = myQueue.get_device();
  std::cout << "Running on "
            << (dev.get_info<sycl::info::device::name>())
            << std::endl;
  std::cout << " is cpu? "
            << (dev.is_cpu() ? "Yes" : "No")
            << std::endl;
  std::cout << " is gpu? "
            << (dev.is_gpu() ? "Yes" : "No")
            << std::endl;
  std::cout << " is accelerator? "
            << (dev.is_accelerator() ? "Yes" : "No")
            << std::endl;
  std::cout << " is the queue in order? "
            << (myQueue.is_in_order() ? "Yes" : "No")
            << std::endl;
  std::cout << " supports usm host allocations? "
            << (dev.get_info<sycl::info::device::usm_host_allocations>() ?
                "Yes" : "No")
            << std::endl;
  std::cout << " supports usm device allocations? "
            << (dev.get_info<sycl::info::device::usm_device_allocations>() ?
                "Yes" : "No")
            << std::endl;
  std::cout << " suports usm shared allocations? "
            << (dev.get_info<sycl::info::device::usm_shared_allocations>() ?
                "Yes" : "No")
            << std::endl;
  std::cout << " max work group size: "
            << dev.get_info<sycl::info::device::max_work_group_size>()
            << std::endl;
  std::cout << " max global memory size (bytes): "
            << dev.get_info<sycl::info::device::global_mem_size>()
            << std::endl;
  std::cout << " max local memory size (bytes): "
            << dev.get_info<sycl::info::device::local_mem_size>()
            << std::endl;
  std::cout << std::endl;

  /* initialize vectors to NULL */
  U = V = X = Y = Z = NULL;

  /* test with all policy variants */
  for (policy = DEFAULT_POL; policy <=GRID_STRIDE; ++policy)
  {
    int actualThreadsPerBlock = threadsPerBlock ? threadsPerBlock :
      myQueue.get_device().get_info<sycl::info::device::max_work_group_size>();

    SUNSyclExecPolicy* stream_exec_policy = NULL;
    SUNSyclExecPolicy* reduce_exec_policy = NULL;

    if (policy == GRID_STRIDE)
    {
      stream_exec_policy = new SUNSyclGridStrideExecPolicy(actualThreadsPerBlock, 1);
      reduce_exec_policy = new SUNSyclBlockReduceExecPolicy(actualThreadsPerBlock, 1);
    }

    /* test with all memory variants */
    for (memtype = UNMANAGED; memtype <= SUNMEMORY; ++memtype)
    {
      SUNMemoryHelper mem_helper = NULL;

      printf("=====> Beginning setup\n\n");

      /* Create new vectors */
      if (memtype == UNMANAGED)
      {
        printf("Testing SYCL N_Vector, policy %d\n", policy);
        X = N_VNew_Sycl(length, &myQueue, sunctx);
      }
      else if (memtype == MANAGED)
      {
        printf("Testing SYCL N_Vector with managed memory, policy %d\n", policy);
        X = N_VNewManaged_Sycl(length, &myQueue, sunctx);
      }
      else if (memtype == SUNMEMORY)
      {
        printf("Testing SYCL N_Vector with SUNMemoryHelper, policy %d\n", policy);
        mem_helper = MyMemoryHelper(sunctx);
        X = N_VNewWithMemHelp_Sycl(length, SUNFALSE, mem_helper, &myQueue, sunctx);
      }
      printf("Vector length: %ld \n", (long int) length);

      if (X == NULL) {
        delete stream_exec_policy;
        delete reduce_exec_policy;
        if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
        printf("FAIL: Unable to create a new vector \n\n");
        Test_Abort(1);
      }

      if (stream_exec_policy != NULL && reduce_exec_policy != NULL) {
        if (N_VSetKernelExecPolicy_Sycl(X, stream_exec_policy, reduce_exec_policy)) {
          N_VDestroy(X);
          delete stream_exec_policy;
          delete reduce_exec_policy;
          if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
          printf("FAIL: Unable to set kernel execution policy \n\n");
          Test_Abort(1);
        }
        printf("Using non-default kernel execution policy\n");
        printf("Threads per block: %d\n\n", actualThreadsPerBlock);
      }

      /* Fill vector with uniform random data in [-1,1] */
      realtype* xdata = N_VGetHostArrayPointer_Sycl(X);
      for (sunindextype j=0; j<length; j++)
        xdata[j] = ((realtype) rand() / (realtype) RAND_MAX)*2-1;
      N_VCopyToDevice_Sycl(X);

      /* Clone additional vectors for testing */
      Y = N_VClone(X);
      if (Y == NULL) {
        N_VDestroy(X);
        printf("FAIL: Unable to create a new vector \n\n");
        delete stream_exec_policy;
        delete reduce_exec_policy;
        if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
        Test_Abort(1);
      }

      Z = N_VClone(X);
      if (Z == NULL) {
        N_VDestroy(X);
        N_VDestroy(Y);
        delete stream_exec_policy;
        delete reduce_exec_policy;
        if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
        printf("FAIL: Unable to create a new vector \n\n");
        Test_Abort(1);
      }

      /* Fill vectors with uniform random data in [-1,1] */
      realtype* ydata = N_VGetHostArrayPointer_Sycl(Y);
      realtype* zdata = N_VGetHostArrayPointer_Sycl(Z);
      for (sunindextype j=0; j<length; j++) {
        ydata[j] = ((realtype) rand() / (realtype) RAND_MAX)*2-1;
        zdata[j] = ((realtype) rand() / (realtype) RAND_MAX)*2-1;
      }
      N_VCopyToDevice_Sycl(Y);
      N_VCopyToDevice_Sycl(Z);

      printf("=====> Setup complete\n");
      printf("=====> Beginning tests\n\n");

      /* Standard vector operation tests */
      printf("\nTesting standard vector operations:\n\n");

      /* Check vector ID */
      fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_SYCL, 0);

      /* Check vector length */
      fails += Test_N_VGetLength(X, 0);

      /* Check vector communicator */
      fails += Test_N_VGetCommunicator(X, NULL, 0);

      /* Test clone functions */
      fails += Test_N_VCloneEmpty(X, 0);
      fails += Test_N_VClone(X, length, 0);
      fails += Test_N_VCloneEmptyVectorArray(5, X, 0);
      fails += Test_N_VCloneVectorArray(5, X, length, 0);

      /* Test vector math kernels */
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
      if (length >= 3) fails += Test_N_VCompare(X, Z, length, 0);
      fails += Test_N_VInvTest(X, Z, length, 0);
      if (length >= 7) fails += Test_N_VConstrMask(X, Y, Z, length, 0);
      fails += Test_N_VMinQuotient(X, Y, length, 0);

      /* Fused and vector array operations tests (disabled) */
      printf("\nTesting fused and vector array operations (disabled):\n\n");

      /* create vector and disable all fused and vector array operations */
      U = N_VClone(X);
      if (U == NULL) {
        N_VDestroy(X);
        N_VDestroy(Y);
        delete stream_exec_policy;
        delete reduce_exec_policy;
        if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
        printf("FAIL: Unable to create a new vector \n\n");
        Test_Abort(1);
      }
      retval = N_VEnableFusedOps_Sycl(U, SUNFALSE);
      if (retval != 0) {
        N_VDestroy(X);
        N_VDestroy(Y);
        N_VDestroy(Z);
        N_VDestroy(U);
        delete stream_exec_policy;
        delete reduce_exec_policy;
        if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
        printf("FAIL: Unable to create a new vector \n\n");
        Test_Abort(1);
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
      V = N_VClone(X);
      retval = N_VEnableFusedOps_Sycl(V, SUNTRUE);
      if (V == NULL) {
        N_VDestroy(X);
        N_VDestroy(Y);
        N_VDestroy(Z);
        N_VDestroy(U);
        delete stream_exec_policy;
        delete reduce_exec_policy;
        if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
        printf("FAIL: Unable to create a new vector \n\n");
        Test_Abort(1);
      }
      if (retval != 0) {
        N_VDestroy(X);
        N_VDestroy(Y);
        N_VDestroy(Z);
        N_VDestroy(U);
        N_VDestroy(V);
        delete stream_exec_policy;
        delete reduce_exec_policy;
        if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
        printf("FAIL: Unable to create a new vector \n\n");
        Test_Abort(1);
      }

      /* fused operations */
      fails += Test_N_VLinearCombination(V, length, 0);
      fails += Test_N_VScaleAddMulti(V, length, 0);
      // fails += Test_N_VDotProdMulti(V, length, 0);

      /* vector array operations */
      fails += Test_N_VLinearSumVectorArray(V, length, 0);
      fails += Test_N_VScaleVectorArray(V, length, 0);
      fails += Test_N_VConstVectorArray(V, length, 0);
      // fails += Test_N_VWrmsNormVectorArray(V, length, 0);
      // fails += Test_N_VWrmsNormMaskVectorArray(V, length, 0);
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
      if (length >= 7) fails += Test_N_VConstrMaskLocal(X, Y, Z, length, 0);
      fails += Test_N_VMinQuotientLocal(X, Y, length, 0);

      /* XBraid interface operations */
      printf("\nTesting XBraid interface operations:\n\n");

      fails += Test_N_VBufSize(X, length, 0);
      fails += Test_N_VBufPack(X, length, 0);
      fails += Test_N_VBufUnpack(X, length, 0);

      printf("\n=====> Beginning teardown\n");

      /* Synchronize */
      myQueue.wait();

      /* Free vectors */
      N_VDestroy(X);
      N_VDestroy(Y);
      N_VDestroy(Z);
      N_VDestroy(U);
      N_VDestroy(V);

      if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);

      printf("=====> Teardown complete\n\n");
    }

    /* Print result */
    if (fails) {
      printf("\n\nFAIL: NVector module failed %i tests \n\n", fails);
    } else {
      printf("\n\nSUCCESS: NVector module passed all tests \n\n");
    }

    delete stream_exec_policy;
    delete reduce_exec_policy;
  }

  /* Synchronize */
  myQueue.wait();

  Test_Finalize();
  return(fails);
}

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, sunindextype length)
{
  int          failure = 0;
  sunindextype i;
  realtype     *Xdata;

  N_VCopyFromDevice_Sycl(X);
  Xdata = N_VGetHostArrayPointer_Sycl(X);

  /* check vector data */
  for (i = 0; i < length; i++) {
    if (failure += SUNRCompare(Xdata[i], ans)) {
      printf("check_ans fail: Xdata[%ld] = %f, expected Xdata[%ld] = %f\n",
             (long int)i, Xdata[i], (long int)i, ans);
    }
  }

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  /* check if vector data is non-null */
  if ((N_VGetHostArrayPointer_Sycl(X) == NULL) &&
      (N_VGetDeviceArrayPointer_Sycl(X) == NULL))
    return SUNFALSE;
  return SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  /* set i-th element of data array */
  set_element_range(X, i, i, val);
}

void set_element_range(N_Vector X, sunindextype is, sunindextype ie,
                       realtype val)
{
  sunindextype i;
  realtype*    xd;

  /* set elements [is,ie] of the data array */
  N_VCopyFromDevice_Sycl(X);
  xd = N_VGetHostArrayPointer_Sycl(X);
  for(i = is; i <= ie; i++) xd[i] = val;
  N_VCopyToDevice_Sycl(X);
}

realtype get_element(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  N_VCopyFromDevice_Sycl(X);
  return (N_VGetHostArrayPointer_Sycl(X))[i];
}

double max_time(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return(time);
}

void sync_device(N_Vector x)
{
  /* sync with host and device */
  ((N_VectorContent_Sycl)(x->content))->queue->wait();
  return;
}
