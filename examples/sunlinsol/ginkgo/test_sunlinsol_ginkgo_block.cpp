/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details->
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol Dense module
 * implementation.
 * ----------------------------------------------------------------- */

#include <algorithm>
#include <cctype>
#include <map>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_ginkgoblock.hpp>

#include "test_sunlinsol.h"

#if defined(USE_HIP)
#define OMP_OR_HIP_OR_CUDA(a, b, c) b
#elif defined(USE_CUDA)
#define OMP_OR_HIP_OR_CUDA(a, b, c) c
#else
#define OMP_OR_HIP_OR_CUDA(a, b, c) a
#endif

#if defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#elif defined(USE_HIP)
#include <nvector/nvector_hip.h>
#else
#include <nvector/nvector_serial.h>
#endif

using namespace sundials::ginkgo;

#if defined(USE_CSR)
using GkoMatrixType = gko::matrix::Csr<sunrealtype>;
using GkoBatchMatrixType =
  gko::batch::matrix::Csr<sunrealtype, gkoblock_indextype>;
using SUNMatrixType = BlockMatrix<GkoBatchMatrixType>;
#else
using GkoMatrixType      = gko::matrix::Dense<sunrealtype>;
using GkoBatchMatrixType = gko::batch::matrix::Dense<sunrealtype>;
using SUNMatrixType      = BlockMatrix<GkoBatchMatrixType>;
#endif

const std::unordered_map<std::string, int> methods{{"bicgstab", 0},
                                                   {"cg", 2},
                                                   {"direct", 3},
                                                   {"gmres", 4},
                                                   {"idr", 5},
                                                   {"richardson", 6}};

constexpr sunrealtype solve_tolerance =
  1000 * std::numeric_limits<sunrealtype>::epsilon();

void fill_matrix_data(gko::matrix_data<sunrealtype>& data)
{
  auto num_rows = data.size[0];
  for (int row = 0; row < num_rows; ++row)
  {
    if (row > 0)
    {
      data.nonzeros.emplace_back(row - 1, row, sunrealtype{-1.0});
    }
    data.nonzeros.emplace_back(row, row, sunrealtype{5.0});
    if (row < num_rows - 1)
    {
      data.nonzeros.emplace_back(row, row + 1, sunrealtype{-1.0});
    }
  }
}

/* ----------------------------------------------------------------------
 * SUNLinSol_Ginkgo Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int argi  = 0;
  int fails = 0; /* counter for test failures    */
  std::string method;
  sunrealtype matcond;           /* matrix condition number      */
  sunindextype matcols, matrows; /* matrix columns, matrows      */
  sunindextype num_blocks;       /* number of matrix blocks      */
  sunindextype max_iters;        /* maximum number of iterations */
  N_Vector x, y, b, s1, s2;      /* test vectors                 */
  int print_timing;
  realtype* xdata;
  sundials::Context sunctx;

  auto gko_exec =
    OMP_OR_HIP_OR_CUDA(gko::OmpExecutor::create(),
                       gko::HipExecutor::create(0, gko::OmpExecutor::create(),
                                                true),
                       gko::CudaExecutor::create(0, gko::OmpExecutor::create(),
                                                 true));

  /* check input and set matrix dimensions */
  if (argc < 7)
  {
    printf("ERROR: SIX (6) Inputs required: method, number of matrix columns, "
           "condition number, number of blocks, "
           "max iterations, print timing \n");
    return (-1);
  }

  method = argv[++argi];
  std::transform(method.begin(), method.end(), method.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (methods.count(method) == 0)
  {
    printf("ERROR: method must be one of ");
    for (const auto& m : methods) { std::cout << m.first << ", "; }
    printf("\n");
    return (-1);
  }

  matcols = static_cast<sunindextype>(atol(argv[++argi]));
  if (matcols <= 0)
  {
    printf("ERROR: number of matrix columns must be a positive integer \n");
    return (-1);
  }
  matrows = matcols;

  matcond = static_cast<sunrealtype>(atof(argv[++argi]));
  if (matcond <= 0)
  {
    printf("ERROR: matrix condition number must be positive\n");
    return (-1);
  }

  num_blocks = static_cast<sunindextype>(atol(argv[++argi]));
  if (num_blocks <= 0)
  {
    printf("ERROR: number of blocks must be a positive integer \n");
    return (-1);
  }

  max_iters = static_cast<sunindextype>(atol(argv[++argi]));
  if (max_iters <= 0)
  {
    printf("ERROR: max iterations must be a positive integer \n");
    return (-1);
  }

  print_timing = atoi(argv[++argi]);
  SetTiming(print_timing);

  printf("\nGinkgo linear solver test: method = %s, size = %ld x %ld, blocks = "
         "%ld, condition = %g, max iters. = "
         "%ld\n\n",
         method.c_str(), static_cast<long int>(matrows),
         static_cast<long int>(matcols), static_cast<long int>(num_blocks),
         matcond, static_cast<long int>(max_iters));

  /* Create vectors and matrices */
  std::default_random_engine engine;
  std::uniform_real_distribution<sunrealtype> distribution_real(8, 10);

  x  = OMP_OR_HIP_OR_CUDA(N_VNew_Serial(num_blocks * matcols, sunctx),
                          N_VNew_Hip(num_blocks * matcols, sunctx),
                          N_VNew_Cuda(num_blocks * matcols, sunctx));
  b  = N_VClone(x);
  y  = N_VClone(x);
  s1 = N_VClone(x);
  s2 = N_VClone(x);

  auto common_size    = gko::dim<2>(matrows, matcols);
  auto batch_mat_size = gko::batch_dim<2>(num_blocks, common_size);
  auto batch_vec_size = gko::batch_dim<2>(num_blocks, gko::dim<2>(matrows, 1));

  auto gko_matrixA = GkoMatrixType::create(gko_exec,
                                           batch_mat_size.get_common_size());
  auto gko_matdata =
    gko::matrix_data<sunrealtype>(batch_mat_size.get_common_size());
  fill_matrix_data(gko_matdata);
  gko_matrixA->read(gko_matdata);
  // FIXME: batch::matrix::Csr needs additional param for num_nnz, dynamic_cast
  // check
  auto gko_batch_matrixA =
    gko::share(GkoBatchMatrixType::create(gko_exec, batch_mat_size));
  for (int b = 0; b < num_blocks; ++b)
  {
    gko_batch_matrixA->create_view_for_item(b)->read(gko_matdata);
  }

  auto gko_matrixB = GkoMatrixType::create(gko_exec,
                                           batch_vec_size.get_common_size());
  auto gko_matdata2 =
    gko::matrix_data<sunrealtype>(batch_vec_size.get_common_size(),
                                  distribution_real, engine);
  gko_matrixA->read(gko_matdata2);
  // FIXME: batch::matrix::Csr needs additional param for num_nnz, dynamic_cast
  // check
  auto gko_batch_matrixB =
    gko::share(GkoBatchMatrixType::create(gko_exec, batch_vec_size));
  for (int b = 0; b < num_blocks; ++b)
  {
    gko_batch_matrixB->create_view_for_item(b)->read(gko_matdata2);
  }

  /* Fill x with random data */
  xdata = N_VGetArrayPointer(x);
  for (sunindextype blocki = 0; blocki < num_blocks; blocki++)
  {
    for (sunindextype i = 0; i < matcols; i++)
    {
      xdata[blocki * matcols + i] = distribution_real(engine);
    }
  }
  OMP_OR_HIP_OR_CUDA(, N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x));

  /* Set scaling vectors */
  N_VConst(sunrealtype{1.0}, s1);
  N_VConst(sunrealtype{1.0}, s2);

  /* Wrap ginkgo matrices for SUNDIALS->
     Matrix is overloaded to a SUNMatrix. */
  SUNMatrixType A{gko_batch_matrixA, sunctx};
  SUNMatrixType B{gko_batch_matrixB, sunctx};

  /* Copy A and x into B and y to print in case of solver failure */
  SUNMatCopy(A, B);
  N_VScale(ONE, x, y);

  /* Create right-hand side vector for linear solve */
  fails += SUNMatMatvecSetup(A);
  fails += SUNMatMatvec(A, x, b);
  if (fails)
  {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");

    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);
    N_VDestroy(s1);
    N_VDestroy(s2);

    return (1);
  }

  /*
   * Create linear solver.
   */

  // auto precond_factory = gko::share(
  //   gko::batch::preconditioner::Identity<sunrealtype>::build().on(gko_exec));

  std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;
  if (method == "bicgstab")
  {
    using GkoSolverType = gko::batch::solver::Bicgstab<sunrealtype>;
    using SUNLinearSolverViewType =
      BlockLinearSolver<GkoSolverType, GkoBatchMatrixType>;
    LS =
      std::make_unique<SUNLinearSolverViewType>(gko_exec,
                                                gko::batch::stop::tolerance_type::absolute,
                                                nullptr, num_blocks, sunctx);
  }
  else if (method == "cg")
  {
    // using GkoSolverType = gko::solver::BatchCg<sunrealtype>;
    // using SUNLinearSolverViewType =
    //   BlockLinearSolver<GkoSolverType, GkoBatchMatrixType>;
    // LS =
    //   std::make_unique<SUNLinearSolverViewType>(gko_exec,
    //                                             gko::stop::batch::ToleranceType::absolute,
    //                                             precond_factory, num_blocks,
    //                                             sunctx);
  }
  else if (method == "direct")
  {
    // using SUNLinearSolverViewType = BlockLinearSolver<SUNMatrixType>;
    // LS = std::make_unique<SUNLinearSolverViewType>(gko_exec, num_blocks,
    // sunctx);
  }
  else if (method == "gmres")
  {
    // using GkoSolverType = gko::solver::BatchGmres<sunrealtype>;
    // using SUNLinearSolverViewType =
    //   BlockLinearSolver<GkoSolverType, GkoBatchMatrixType>;
    // LS =
    //   std::make_unique<SUNLinearSolverViewType>(gko_exec,
    //                                             gko::stop::batch::ToleranceType::absolute,
    //                                             precond_factory, num_blocks,
    //                                             sunctx);
  }
  else if (method == "idr")
  {
    // using GkoSolverType = gko::solver::BatchIdr<sunrealtype>;
    // using SUNLinearSolverViewType =
    //   BlockLinearSolver<GkoSolverType, GkoBatchMatrixType>;
    // LS =
    //   std::make_unique<SUNLinearSolverViewType>(gko_exec,
    //                                             gko::stop::batch::ToleranceType::absolute,
    //                                             precond_factory, num_blocks,
    //                                             sunctx);
  }

  /* Run Tests */
  fails += Test_SUNLinSolGetID(LS->Convert(), SUNLINEARSOLVER_GINKGOBLOCK, 0);
  fails += Test_SUNLinSolGetType(LS->Convert(),
                                 SUNLINEARSOLVER_MATRIX_ITERATIVE, 0);
  fails += Test_SUNLinSolInitialize(LS->Convert(), 0);
  fails += Test_SUNLinSolSetScalingVectors(LS->Convert(), s1, s2, 0);
  fails += Test_SUNLinSolSetup(LS->Convert(), A, 0);
  fails += Test_SUNLinSolSolve(LS->Convert(), A, x, b, solve_tolerance, SUNTRUE,
                               0);

  /* Print result */
  if (fails)
  {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
    /* only print if its a small problem */
    if (matcols < 4)
    {
      printf("\nx (original) =\n");
      N_VPrint(y);
      printf("\nx (computed) =\n");
      N_VPrint(x);
      printf("\nb =\n");
      N_VPrint(b);
    }
  }
  else { printf("\nSUCCESS: SUNLinSol module passed all tests \n \n"); }

  /* Free solver, matrix and vectors */
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(b);
  N_VDestroy(s1);
  N_VDestroy(s2);

  return (fails);
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_vector(N_Vector expected, N_Vector actual, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* copy vectors to host */
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(actual),
                     N_VCopyFromDevice_Cuda(actual));
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(expected),
                     N_VCopyFromDevice_Cuda(expected));

  /* get vector data */
  xdata = N_VGetArrayPointer(actual);
  ydata = N_VGetArrayPointer(expected);

  /* check data lengths */
  xldata = N_VGetLength(actual);
  yldata = N_VGetLength(expected);

  if (xldata != yldata)
  {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return (1);
  }

  /* check vector data */
  for (i = 0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO)
  {
    printf("Check_vector failures:\n");
    for (i = 0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        printf("  actual[%ld] = %g != %e (err = %g)\n", (long int)i, xdata[i],
               ydata[i], SUNRabs(xdata[i] - ydata[i]));
  }

  return failure > 0;
}

void sync_device()
{
  OMP_OR_HIP_OR_CUDA(, hipDeviceSynchronize(), cudaDeviceSynchronize());
}
