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

//#include <nvector/nvector_kokkos.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>

#include "test_sunmatrix.h"

#if defined(USE_CUDA)
using ExecSpace = Kokkos::Cuda;
#elif defined(USE_HIP)
using ExecSpace = Kokkos::HIP;
#elif defined(USE_OPENMP)
using ExecSpace = Kokkos::OpenMP;
#else
using ExecSpace = Kokkos::Serial;
#endif

/* -----------------------------------------------------------------------------
 * SUNMatrix Testing
 * ---------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  // Create SUNDIALS context
  sundials::Context sunctx;

  // counter for test failures
  int fails = 0;

  // check input and set vector length
  if (argc < 5)
  {
    std::cerr << "ERROR: FOUR (4) Input required:\n"
              << " (1) matrix rows\n"
              << " (2) matrix cols\n"
              << " (3) number of matrix blocks\n"
              << " (4) print timing\n";
    return -1;
  }

  sunindextype matrows = (sunindextype) atol(argv[1]);
  if (matrows <= 0)
  {
    std::cerr << "ERROR: number of rows must be a positive integer\n";
    return -1;
  }

  sunindextype matcols = (sunindextype) atol(argv[2]);
  if (matcols <= 0)
  {
    std::cerr << "ERROR: number of cols must be a positive integer\n";
    return -1;
  }

  sunindextype nblocks = (sunindextype) atol(argv[3]);
  if (nblocks <= 0)
  {
    std::cerr << "ERROR: number of blocks must be a positive integer\n";
    return -1;
  }

  int print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  bool square = (matrows == matcols);
  std::cout << "\nKokkos dense matrix test:\n"
            << "  Size:   " << matrows << " by " << matcols << "\n"
            << "  Blocks: " << nblocks << "\n";

  Kokkos::initialize( argc, argv );
  {
    auto exec_instance = ExecSpace();

    // Create vectors and matrices
    // N_Vector x = N_VNew_Kokkos(matcols * nblocks, sunctx);
    // N_Vector y = N_VNew_Kokkos(matrows * nblocks, sunctx);
    sundials::kokkos::DenseMatrix<ExecSpace> A{nblocks, matrows, matcols,
                                               exec_instance, sunctx};
    sundials::kokkos::DenseMatrix<ExecSpace> I{nblocks, matrows, matcols,
                                               exec_instance, sunctx};

    // Fill matrices and vectors
    auto A_data = A.view();

    Kokkos::parallel_for("fill_A",
                         Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>
                         ({0, 0, 0}, {nblocks, matrows, matcols}),
                         KOKKOS_LAMBDA(const sunindextype i,
                                       const sunindextype j,
                                       const sunindextype k)
                         {
                           A_data(i, j, k) = (k + 1) * (j + k);
                         });

    if (square)
    {
      auto I_data = I.view();
      Kokkos::parallel_for("fill_I",
                           Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>
                           ({0, 0, 0}, {nblocks, matrows, matcols}),
                           KOKKOS_LAMBDA(const sunindextype i,
                                         const sunindextype j,
                                         const sunindextype k)
                           {
                             I_data(i, j, k) = (j == k) ? ONE : ZERO;
                           });
    }

    // auto x_data = VEC_VIEW(x);

    // Kokkos::parallel_for("fill_x",
    //                      Kokkos::RangePolicy<ExecSpace>
    //                      (0, matcols * nblocks),
    //                      KOKKOS_LAMBDA(const sunindextype j)
    //                      {
    //                        x_data(j) = ONE / ((j % matcols) + 1);
    //                      });

    // auto y_data = VEC_VIEW(y);

    // Kokkos::parallel_for("fill_y",
    //                      Kokkos::RangePolicy<ExecSpace>
    //                      (0, matrows * nblocks),
    //                      KOKKOS_LAMBDA(const sunindextype j)
    //                      {
    //                        auto m = j % matrows;
    //                        auto n = m + matcols - 1;
    //                        y_data(j) = HALF * (n + 1 - m) * (n + m);
    //                      });

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
    // fails += Test_SUNMatMatvec(A, x, y, 0);

    // Print result
    if (fails)
      std::cout << "FAIL: SUNMatrix module failed " << fails << " tests\n\n";
    else
      std::cout << "SUCCESS: SUNMatrix module passed all tests\n\n";
  }
  Kokkos::finalize();

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

extern "C" int check_matrix(SUNMatrix A, SUNMatrix B, sunrealtype tol)
{
  int failure = 0;

  auto A_mat{sundials::kokkos::GetDenseMat<ExecSpace>(A)};
  auto B_mat{sundials::kokkos::GetDenseMat<ExecSpace>(B)};

  const auto nblocks = A_mat->blocks();
  const auto matrows = A_mat->block_rows();
  const auto matcols = A_mat->block_cols();

  auto A_data{A_mat->view()};
  auto B_data{B_mat->view()};

  // compare data
  Kokkos::parallel_reduce("check_matrix",
                          Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>
                          ({0, 0, 0}, {nblocks, matrows, matcols}),
                          KOKKOS_LAMBDA(const sunindextype i,
                                        const sunindextype j,
                                        const sunindextype k, int &l_fail)
                          {
                            l_fail += CompareTol(A_data(i, j, k), B_data(i , j, k), tol);
                          }, failure);

  if (failure > ZERO)
    return 1;
  else
    return 0;
}

extern "C" int check_matrix_entry(SUNMatrix A, sunrealtype val, sunrealtype tol)
{
  int failure = 0;

  auto A_mat{sundials::kokkos::GetDenseMat<ExecSpace>(A)};

  const auto nblocks = A_mat->blocks();
  const auto matrows = A_mat->block_rows();
  const auto matcols = A_mat->block_cols();

  auto A_data{A_mat->view()};

  // compare data
  Kokkos::parallel_reduce("check_matrix_entry",
                          Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>
                          ({0, 0, 0}, {nblocks, matrows, matcols}),
                          KOKKOS_LAMBDA(const sunindextype i,
                                        const sunindextype j,
                                        const sunindextype k, int &l_fail)
                          {
                            l_fail += CompareTol(A_data(i, j, k), val, tol);
                          }, failure);

  if (failure > ZERO)
    return 1;
  else
    return 0;
}

extern "C" int check_vector(N_Vector actual, N_Vector expected, sunrealtype tol)
{
  int failure = 0;

  // auto a_data = VEC_VIEW(actual);
  // auto e_data = VEC_VIEW(expected);
  // auto length = VEC_CONTENT(actual)->length;

  // // check vector data
  // Kokkos::parallel_reduce("check_vector",
  //                         Kokkos::RangePolicy<ExecSpace>(0, length),
  //                         KOKKOS_LAMBDA(const sunindextype i, int &l_fail)
  //                         {
  //                           l_fail += CompareTol(a_data(i), e_data(i), tol);
  //                         }, failure);

  if (failure > ZERO)
    return 1;
  else
    return 0;
}

extern "C" booleantype has_data(SUNMatrix A)
{
  return SUNTRUE;
}

extern "C" booleantype is_square(SUNMatrix A)
{
  auto A_mat{sundials::kokkos::GetDenseMat<ExecSpace>(A)};

  const auto matrows = A_mat->rows();
  const auto matcols = A_mat->cols();

  if (matrows == matcols)
    return SUNTRUE;
  else
    return SUNFALSE;
}

extern "C" void sync_device(SUNMatrix A)
{
  Kokkos::fence();
}
