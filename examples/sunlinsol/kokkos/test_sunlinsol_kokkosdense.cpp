/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the testing routine for the dense SUNLinearSolver using Kokkos.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include <nvector/nvector_kokkos.h>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>
#include <sunlinsol/sunlinsol_kokkosdense.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "test_sunlinsol.h"

#define VEC_CONTENT(x) ( static_cast<N_VectorContent_Kokkos>(x->content) )
#define VEC_VIEW(x)    ( *VEC_CONTENT(x)->device_data )

#define MAT_CONTENT(A) ( static_cast<SUNMatrixContent_KokkosDense>(A->content) )
#define MAT_VIEW(A)    ( MAT_CONTENT(A)->data_view )

/* -----------------------------------------------------------------------------
 * SUNLinearSolver Testing
 * ---------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  // counter for test failures
  int fails = 0;

  // check input and set matrix dimensions
  if (argc < 4) {
    printf("ERROR: THREE (3) Inputs required: matrix cols, number of blocks, print timing \n");
    return -1;
  }

  sunindextype cols = (sunindextype) atol(argv[1]);
  if (cols <= 0) {
    printf("ERROR: number of matrix columns must be a positive integer \n");
    return -1;
  }
  sunindextype rows = cols;

  sunindextype nblocks = (sunindextype) atol(argv[2]);
  if (nblocks <= 0) {
    printf("ERROR: number of blocks must be a positive integer \n");
    return -1;
  }

  int print_timing = atoi(argv[3]);
  SetTiming(print_timing);

  SUNContext sunctx = nullptr;
  if (SUNContext_Create(nullptr, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return -1;
  }

  printf("\nKokkos dense linear solver test: size %ld, blocks %ld\n\n",
         (long int) cols, (long int) nblocks);

  Kokkos::initialize( argc, argv );
  {
    // Create matrices and vectors
    SUNMatrix A;
    if (nblocks > 1)
      A = SUNMatrix_KokkosDenseBlock(nblocks, rows, cols, sunctx);
    else
      A = SUNMatrix_KokkosDense(rows, cols, sunctx);
    N_Vector x = N_VNew_Kokkos(cols * nblocks, sunctx);
    N_Vector b = N_VClone(x);

    using RandPoolType = Kokkos::Random_XorShift64_Pool<>;
    using GenType = RandPoolType::generator_type;

    RandPoolType rand_pool(5374857);

    // Fill A matrix with uniform random data 1 + rand[0,1]
    auto A_data = MAT_VIEW(A);

    Kokkos::parallel_for("fill_A",
                         Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nblocks, rows, cols}),
                         KOKKOS_LAMBDA(const sunindextype i,
                                       const sunindextype j,
                                       const sunindextype k)
                         {
                           auto rgen = rand_pool.get_state();
                           auto rval = Kokkos::rand<GenType, sunrealtype>::draw(rgen, ONE);
                           A_data(i, j, k) = ONE + rval;
                           rand_pool.free_state(rgen);
                         });

    // Fill x vector with uniform random data in 1 + rand[0,1]
    auto x_data = VEC_VIEW(x);

    Kokkos::parallel_for("fill_x",
                         Kokkos::RangePolicy<>(0, cols * nblocks),
                         KOKKOS_LAMBDA(const sunindextype j)
                         {
                           auto rgen = rand_pool.get_state();
                           auto rval = Kokkos::rand<GenType, sunrealtype>::draw(rgen, ONE);
                           x_data(j) = ONE + rval;
                           rand_pool.free_state(rgen);
                         });

    // create right-hand side vector for linear solve
    fails += SUNMatMatvecSetup(A);
    fails += SUNMatMatvec(A, x, b);
    if (fails)
    {
      printf("FAIL: SUNLinSol SUNMatMatvec failure\n");

      // Free matrices and vectors
      SUNMatDestroy(A);
      N_VDestroy(x);
      N_VDestroy(b);

      return 1;
    }

    // Create dense linear solver
    SUNLinearSolver LS = SUNLinSol_KokkosDense(x, A, sunctx);
    if (!LS) {
      printf("FAIL: SUNLinSol_KokkosDense failure\n");

      // Free matrices and vectors
      SUNMatDestroy(A);
      N_VDestroy(x);
      N_VDestroy(b);

      return 1;
    }

    // Run Tests
    fails += Test_SUNLinSolInitialize(LS, 0);
    fails += Test_SUNLinSolSetup(LS, A, 0);
    fails += Test_SUNLinSolSolve(LS, A, x, b, SUN_RCONST(1.0e-10), SUNTRUE, 0);
    fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
    fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_KOKKOSDENSE, 0);
    fails += Test_SUNLinSolLastFlag(LS, 0);
    fails += Test_SUNLinSolSpace(LS, 0);

    // Print result
    if (fails) {
      printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
    } else {
      printf("SUCCESS: SUNLinSol module passed all tests \n \n");
    }

    // Free solver, matrix and vectors
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(x);
    N_VDestroy(b);
  }
  Kokkos::finalize();

  SUNContext_Free(&sunctx);

  return fails;
}

/* -----------------------------------------------------------------------------
 * Check functions
 * ---------------------------------------------------------------------------*/

KOKKOS_FUNCTION
int CompareTol(sunrealtype a, sunrealtype b, sunrealtype tol)
{
  if (a == b) return 0;
  if (std::isnan(a) || std::isnan(b)) return 1;
  if (std::isinf(a) || std::isinf(b)) return 1;
  sunrealtype diff = std::abs(a - b);
  sunrealtype norm = std::min(std::abs(a + b),
                              std::numeric_limits<sunrealtype>::max());
  return diff >= std::max(10 * std::numeric_limits<sunrealtype>::epsilon(),
                          tol * norm);
}

int check_vector(N_Vector expected, N_Vector computed, realtype tol)
{
  int failure = 0;

  auto e_data = VEC_VIEW(expected);
  auto c_data = VEC_VIEW(computed);
  auto length = VEC_CONTENT(expected)->length;

  // check vector data
  Kokkos::parallel_reduce("check_vector",
                          Kokkos::RangePolicy<>(0, length),
                          KOKKOS_LAMBDA(const sunindextype i, int &l_fail)
                          {
                            l_fail += CompareTol(c_data(i), e_data(i), tol);
                          }, failure);

  if (failure > ZERO)
  {
    constexpr auto max_precision {std::numeric_limits<long double>::digits10};
    std::cout << std::setprecision(max_precision);
    for (int i = 0; i < length; i++)
      if (CompareTol(c_data(i), e_data(i), tol))
        std::cout << c_data(i) << " vs " << e_data(i) << " diff " << std::abs(c_data(i) - e_data(i)) << std::endl;
    return 1;
  }
  else
    return 0;
}

void sync_device()
{
  Kokkos::fence();
}
