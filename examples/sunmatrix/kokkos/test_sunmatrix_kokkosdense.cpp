/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the testing routine for the dense SUNMatrix using Kokkos.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <nvector/nvector_kokkos.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>

#include "test_sunmatrix.h"

#define VEC_CONTENT(x) ( static_cast<N_VectorContent_Kokkos>(x->content) )
#define VEC_VIEW(x)    ( *VEC_CONTENT(x)->device_data )

#define MAT_CONTENT(A) ( static_cast<SUNMatrixContent_KokkosDense>(A->content) )
#define MAT_VIEW(A)    ( MAT_CONTENT(A)->data_view )

/* -----------------------------------------------------------------------------
 * SUNMatrix Testing
 * ---------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  // counter for test failures
  int fails = 0;

  // check input and set vector length
  if (argc < 5) {
    printf("ERROR: FOUR (4) Input required: matrix rows, matrix cols, number of matrix blocks, print timing \n");
    return -1;
  }

  sunindextype matrows = (sunindextype) atol(argv[1]);
  if (matrows <= 0) {
    printf("ERROR: number of rows must be a positive integer \n");
    return -1;
  }

  sunindextype matcols = (sunindextype) atol(argv[2]);
  if (matcols <= 0) {
    printf("ERROR: number of cols must be a positive integer \n");
    return -1;
  }

  sunindextype nblocks = (sunindextype) atol(argv[3]);
  if (nblocks <= 0) {
    printf("ERROR: number of blocks must be a positive integer \n");
    return -1;
  }

  int print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  SUNContext sunctx = nullptr;
  if (SUNContext_Create(nullptr, &sunctx))
  {
    printf("ERROR: SUNContext_Create failed\n");
    return -1;
  }

  bool square = (matrows == matcols);
  printf("\nKokkos dense matrix test: size %ld by %ld\n\n",
         (long int) matrows, (long int) matcols);

  Kokkos::initialize( argc, argv );
  {
    // Create vectors and matrices
    N_Vector  x = N_VNew_Kokkos(matcols * nblocks, sunctx);
    N_Vector  y = N_VNew_Kokkos(matrows * nblocks, sunctx);
    SUNMatrix A = SUNMatrix_KokkosDenseBlock(nblocks, matrows, matcols, sunctx);
    SUNMatrix I = square ? SUNMatClone(A) : nullptr;

    // Fill matrices and vectors
    auto A_data = MAT_VIEW(A);

    Kokkos::parallel_for("fill_A",
                         Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nblocks, matrows, matcols}),
                         KOKKOS_LAMBDA(const sunindextype i,
                                       const sunindextype j,
                                       const sunindextype k)
                         {
                           A_data(i, j, k) = (k + 1) * (j + k);
                         });

    if (square)
    {
      auto I_data = MAT_VIEW(I);
      Kokkos::parallel_for("fill_I",
                           Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nblocks, matrows, matcols}),
                           KOKKOS_LAMBDA(const sunindextype i,
                                         const sunindextype j,
                                         const sunindextype k)
                           {
                             I_data(i, j, k) = (j == k) ? ONE : ZERO;
                           });
    }

    auto x_data = VEC_VIEW(x);

    Kokkos::parallel_for("fill_x",
                         Kokkos::RangePolicy<>(0, matcols * nblocks),
                         KOKKOS_LAMBDA(const sunindextype j)
                         {
                           x_data(j) = ONE / ((j % matcols) + 1);
                         });

    auto y_data = VEC_VIEW(y);

    Kokkos::parallel_for("fill_y",
                         Kokkos::RangePolicy<>(0, matrows * nblocks),
                         KOKKOS_LAMBDA(const sunindextype j)
                         {
                           auto m = j % matrows;
                           auto n = m + matcols - 1;
                           y_data(j) = HALF * (n + 1 - m) * (n + m);
                         });

    // SUNMatrix Tests
    fails += Test_SUNMatGetID(A, SUNMATRIX_KOKKOSDENSE, 0);
    fails += Test_SUNMatClone(A, 0);
    fails += Test_SUNMatCopy(A, 0);
    fails += Test_SUNMatZero(A, 0);
    if (square)
    {
      fails += Test_SUNMatScaleAdd(A, I, 0);
      fails += Test_SUNMatScaleAddI(A, I, 0);
    }
    fails += Test_SUNMatMatvec(A, x, y, 0);

    // Print result
    if (fails)
      printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
    else
      printf("SUCCESS: SUNMatrix module passed all tests \n \n");

    // Free vectors and matrices
    N_VDestroy(x);
    N_VDestroy(y);
    SUNMatDestroy(A);
    if (square)
      SUNMatDestroy(I);
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

int check_matrix(SUNMatrix A, SUNMatrix B, sunrealtype tol)
{
  int failure = 0;

  const auto nblocks = MAT_CONTENT(A)->nblocks;
  const auto matrows = MAT_CONTENT(A)->M;
  const auto matcols = MAT_CONTENT(A)->N;

  auto A_data = MAT_VIEW(A);
  auto B_data = MAT_VIEW(B);

  // compare data
  Kokkos::parallel_reduce("check_matrix",
                          Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nblocks, matrows, matcols}),
                          KOKKOS_LAMBDA(const sunindextype i, const sunindextype j, const sunindextype k, int &l_fail)
                          {
                            l_fail += CompareTol(A_data(i, j, k), B_data(i , j, k), tol);
                          }, failure);

  if (failure > ZERO)
    return 1;
  else
    return 0;
}

int check_matrix_entry(SUNMatrix A, sunrealtype val, sunrealtype tol)
{
  int failure = 0;

  const auto nblocks = MAT_CONTENT(A)->nblocks;
  const auto matrows = MAT_CONTENT(A)->M;
  const auto matcols = MAT_CONTENT(A)->N;

  auto A_data = MAT_VIEW(A);

  // compare data
  Kokkos::parallel_reduce("check_matrix_entry",
                          Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nblocks, matrows, matcols}),
                          KOKKOS_LAMBDA(const sunindextype i, const sunindextype j, const sunindextype k, int &l_fail)
                          {
                            l_fail += CompareTol(A_data(i, j, k), val, tol);
                          }, failure);

  if (failure > ZERO)
    return 1;
  else
    return 0;
}

int check_vector(N_Vector actual, N_Vector expected, sunrealtype tol)
{
  int failure = 0;

  auto a_data = VEC_VIEW(actual);
  auto e_data = VEC_VIEW(expected);
  auto length = VEC_CONTENT(actual)->length;

  // check vector data
  Kokkos::parallel_reduce("check_vector",
                          Kokkos::RangePolicy<>(0, length),
                          KOKKOS_LAMBDA(const sunindextype i, int &l_fail)
                          {
                            l_fail += CompareTol(a_data(i), e_data(i), tol);
                          }, failure);

  if (failure > ZERO)
    return 1;
  else
    return 0;
}

booleantype has_data(SUNMatrix A)
{
  return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  if (SUNMatrix_KokkosDense_Rows(A) == SUNMatrix_KokkosDense_Columns(A))
    return SUNTRUE;
  else
    return SUNFALSE;
}

void sync_device(SUNMatrix A)
{
  Kokkos::fence();
}
