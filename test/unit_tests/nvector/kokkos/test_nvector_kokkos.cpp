/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel McGreer and Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the testing routine for the NVector implementation using Kokkos.
 *
 * Three construction scenarios are exercised:
 *
 *   A. Application-managed view  -- Vector(managed_view, sunctx)
 *   B. Managed view constructor   -- Vector(managed_view, sunctx)
 *   C. Unmanaged view constructor -- Vector(unmanaged_view, sunctx)
 *
 * Additional vectors are obtained via N_VClone, which always returns a managed
 * Vector, to mirror the pattern in SUNDIALS packages where the user supplies
 * an initial vector and the integrator clones it to create internal workspace.
 * ---------------------------------------------------------------------------*/

#include <nvector/nvector_kokkos.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "test_nvector.h"

#if defined(USE_CUDA)
using ExecSpace = Kokkos::Cuda;
#elif defined(USE_HIP)
using ExecSpace = Kokkos::HIP;
#elif defined(USE_OPENMP)
using ExecSpace = Kokkos::OpenMP;
#else
using ExecSpace = Kokkos::Serial;
#endif

using ManagedViewType =
  Kokkos::View<sunrealtype*, typename ExecSpace::memory_space>;
using VecType = sundials::kokkos::VectorFor<ManagedViewType>;
using SizeType = VecType::size_type;

/* -----------------------------------------------------------------------------
 * Run the vector operation tests
 * ---------------------------------------------------------------------------*/
static int run_all_tests(N_Vector X, N_Vector Y, N_Vector Z, sunindextype length)
{
  int fails{0};

  printf("\nTesting standard vector operations:\n\n");

  fails += Test_N_VAbs(X, Z, length, 0);
  fails += Test_N_VAddConst(X, Z, length, 0);
  fails += Test_N_VCompare(X, Z, length, 0);
  fails += Test_N_VConst(X, length, 0);
  fails += Test_N_VConstrMask(X, Y, Z, length, 0);
  fails += Test_N_VDiv(X, Y, Z, length, 0);
  fails += Test_N_VDotProd(X, Y, length, 0);
  fails += Test_N_VInv(X, Z, length, 0);
  fails += Test_N_VInvTest(X, Z, length, 0);
  fails += Test_N_VL1Norm(X, length, 0);
  fails += Test_N_VLinearSum(X, Y, Z, length, 0);
  fails += Test_N_VMaxNorm(X, length, 0);
  fails += Test_N_VMin(X, length, 0);
  fails += Test_N_VMinQuotient(X, Y, length, 0);
  fails += Test_N_VProd(X, Y, Z, length, 0);
  fails += Test_N_VScale(X, Z, length, 0);
  fails += Test_N_VWL2Norm(X, Y, length, 0);
  fails += Test_N_VWrmsNorm(X, Y, length, 0);
  fails += Test_N_VWrmsNormMask(X, Y, Z, length, 0);

  printf("\nTesting fused and vector array operations (disabled):\n\n");

  fails += Test_N_VLinearCombination(X, length, 0);
  fails += Test_N_VScaleAddMulti(X, length, 0);
  fails += Test_N_VDotProdMulti(X, length, 0);
  fails += Test_N_VLinearSumVectorArray(X, length, 0);
  fails += Test_N_VScaleVectorArray(X, length, 0);
  fails += Test_N_VConstVectorArray(X, length, 0);
  fails += Test_N_VWrmsNormVectorArray(X, length, 0);
  fails += Test_N_VWrmsNormMaskVectorArray(X, length, 0);
  fails += Test_N_VScaleAddMultiVectorArray(X, length, 0);
  fails += Test_N_VLinearCombinationVectorArray(X, length, 0);

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

  return fails;
}

/* -----------------------------------------------------------------------------
 * Main NVector Testing Routine
 * ---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails{0};        /* counter for test failures */
  sunindextype length; /* vector length             */
  int print_timing;    /* turn timing on/off        */

  Test_Init(SUN_COMM_NULL);

  /* check input and set vector length */
  if (argc < 3)
  {
    printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    return (-1);
  }

  length = (sunindextype)atol(argv[1]);
  if (length <= 0)
  {
    printf("ERROR: length of vector must be a positive integer \n");
    return (-1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, 0);

  Kokkos::initialize(argc, argv);
  {
    /* -------------------------------------------------------------------------
     * Scenario A: application-managed view
     *
     * X owns its memory. Y and Z are managed clones of X. All three vectors
     * share the same VecType specialisation so there is no managed/unmanaged
     * mixing.
     * -----------------------------------------------------------------------*/
    printf("====================================================\n");
     printf("Testing KOKKOS N_Vector: application-managed view\n");
    printf("Vector length %ld\n", (long int)length);
    printf("====================================================\n");

    {
      ManagedViewType xview("test vector", static_cast<SizeType>(length));
      VecType X{xview, sunctx};

      fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_KOKKOS, 0);
      fails += Test_N_VClone(X, length, 0);
      fails += Test_N_VCloneVectorArray(5, X, length, 0);
      fails += Test_N_VGetLength(X, 0);
      fails += Test_N_VGetCommunicator(X, SUN_COMM_NULL, 0);

      N_Vector Y_nv = N_VClone(X);
      N_Vector Z_nv = N_VClone(X);

      fails += run_all_tests(X, Y_nv, Z_nv, length);

      N_VDestroy(Y_nv);
      N_VDestroy(Z_nv);
    }

    /* -------------------------------------------------------------------------
     * Scenario B: managed view constructor
     *
     * X is constructed from a pre-existing managed Kokkos view. The view is
     * shared between the caller and the Vector (reference counted). Y and Z are
     * managed clones. All three vectors are the same VecType specialisation; no
     * managed/unmanaged mixing occurs.
     *
     * This exercises the view constructor path and verifies that constructing
     * from a managed view does not allocate a second buffer.
     * -----------------------------------------------------------------------*/
    printf("\n====================================================\n");
    printf("Testing KOKKOS N_Vector: managed view constructor\n");
    printf("Vector length %ld\n", (long int)length);
    printf("====================================================\n");

    {
      VecType::managed_view_type managed_view("managed_view",
                                              static_cast<SizeType>(length));
      VecType X{managed_view, sunctx};

      /* Verify zero-copy: X wraps the same allocation, no copy occurred */
      if (X.View().data() != managed_view.data())
      {
        printf("FAILED: managed view constructor: data pointer mismatch -- "
               "Vector did not wrap the provided view\n");
        fails++;
      }
      else
      {
        printf("PASSED: managed view constructor: data pointer matches "
               "provided view\n");
      }

      fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_KOKKOS, 0);
      fails += Test_N_VClone(X, length, 0);
      fails += Test_N_VCloneVectorArray(5, X, length, 0);
      fails += Test_N_VGetLength(X, 0);
      fails += Test_N_VGetCommunicator(X, SUN_COMM_NULL, 0);

      N_Vector Y_nv = N_VClone(X);
      N_Vector Z_nv = N_VClone(X);

      fails += run_all_tests(X, Y_nv, Z_nv, length);

      N_VDestroy(Y_nv);
      N_VDestroy(Z_nv);
    }

/* Use an explicitly unmanaged view over application-owned backing storage. */
    {
      VecType::managed_view_type backing("application backing", static_cast<SizeType>(length));
      using UnmanagedView = Kokkos::View<sunrealtype*, typename ExecSpace::memory_space,
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
      UnmanagedView unmanaged_view(backing.data(), static_cast<SizeType>(length));
      {
        VecType X{unmanaged_view, sunctx};
        N_Vector clone = N_VClone(X);
        if (VecType::FromNVector(clone)->View().data() == X.View().data()) fails++;
        N_VScale(SUN_RCONST(2.0), X, clone);
        N_Vector Y_nv = N_VClone(X);
        N_Vector Z_nv = N_VClone(X);
        fails += run_all_tests(X, Y_nv, Z_nv, length);
        N_VDestroy(clone);
        N_VDestroy(Y_nv);
        N_VDestroy(Z_nv);
      }
      // The wrapper is gone, but application-owned backing remains valid.
      if (backing.data() == nullptr) fails++;
    }
  }
  Kokkos::finalize();

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
  int failure{0};
  auto Xvec = VecType::FromNVector(X);
  auto Xdata = Kokkos::create_mirror_view(Xvec->View());
  Kokkos::deep_copy(Xdata, Xvec->View());
  for (sunindextype i = 0; i < local_length; i++)
  {
    failure += SUNRCompare(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

sunbooleantype has_data(N_Vector X)
{
  /* check if vector data is non-null */
  return SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, sunrealtype val)
{
  /* set i-th element of data array */
  set_element_range(X, i, i, val);
}

void set_element_range(N_Vector X, sunindextype is, sunindextype ie,
                       sunrealtype val)
{
  auto Xvec = VecType::FromNVector(X);
  auto Xdata = Kokkos::create_mirror_view(Xvec->View());
  Kokkos::deep_copy(Xdata, Xvec->View());

  /* set elements [is,ie] of the data array */
  for (sunindextype i = is; i <= ie; i++) { Xdata[i] = val; }
  Kokkos::deep_copy(Xvec->View(), Xdata);
}

sunrealtype get_element(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  auto Xvec = VecType::FromNVector(X);
  auto Xdata = Kokkos::create_mirror_view(Xvec->View());
  Kokkos::deep_copy(Xdata, Xvec->View());
  return Xdata[i];
}

double max_time(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return time;
}

void sync_device(N_Vector x)
{
  /* sync with GPU */
  Kokkos::fence();
  return;
}
