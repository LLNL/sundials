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
#include <cstdio>
#include <cstdlib>
#include <map>
#include <ginkgo/ginkgo.hpp>
#include <random>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_ginkgo.hpp>
#include <unordered_map>

#include "test_sunlinsol.h"

#if defined(USE_HIP)
#include <nvector/nvector_hip.h>
#define HIP_OR_CUDA(a, b) a
constexpr auto N_VNew = N_VNew_Raja;
#elif defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#define HIP_OR_CUDA(a, b) b
constexpr auto N_VNew = N_VNew_Cuda;
#elif defined(USE_OMP)
#include <nvector/nvector_serial.h>
#define HIP_OR_CUDA(a, b)
constexpr auto N_VNew = N_VNew_Serial;
#else
#include <nvector/nvector_serial.h>
#define HIP_OR_CUDA(a, b)
constexpr auto N_VNew = N_VNew_Serial;
#endif

#if defined(USE_CSR)
using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
#else
using GkoMatrixType = gko::matrix::Dense<sunrealtype>;

#endif

const std::unordered_map<std::string, int> methods{
  {"bicg", 0}, {"bicgstab", 1}, {"cg", 2}, {"cgs", 3},
  {"fcg", 4},  {"gmres", 5},    {"idr", 6}
  // {"multigrid",7} // does not support the combined stopping criteria we use
  // {"cbgmres", 8} // does not support setting stopping criteria
};

template <typename ValueType>
void fill_matrix(gko::matrix::Dense<ValueType>* matrix)
{
  const auto matrows = matrix->get_size()[0];
  const ValueType coefs[] = {-1, 2, -1};
  for (auto i = 0; i < matrows; ++i)
  {
    for (auto dofs : {-1, 0, 1})
    {
      if (0 <= i + dofs && i + dofs < matrows)
      {
        matrix->at(i, i + dofs) = coefs[dofs + 1];
      }
    }
  }
}

template <typename ValueType, typename IndexType>
void fill_matrix(gko::matrix::Csr<ValueType, IndexType>* matrix)
{
  const auto matrows = matrix->get_size()[0];
  auto row_ptrs = matrix->get_row_ptrs();
  auto col_idxs = matrix->get_col_idxs();
  auto values = matrix->get_values();
  int pos = 0;

  const ValueType coefs[] = {-1, 2, -1};
  row_ptrs[0] = pos;
  for (auto i = 0; i < matrows; ++i)
  {
    for (auto ofs : {-1, 0, 1})
    {
      if (0 <= i + ofs && i + ofs < matrows)
      {
        values[pos] = coefs[ofs + 1];
        col_idxs[pos] = i + ofs;
        ++pos;
      }
    }
    row_ptrs[i + 1] = pos;
  }
}


/* -------------------------------------------------------------------------- *
 * SUNLinSol_Ginkgo Testing Routine                                           *
 * -------------------------------------------------------------------------- */

int main(int argc, char* argv[])
{
  int argi{0};
  int fails{0}; /* counter for test failures    */

  sundials::Context sunctx;

#if defined(USE_HIP)
  auto gko_exec{gko::HipExecutor::create(0, gko::OmpExecutor::create(), true)};
#elif defined(USE_CUDA)
  auto gko_exec{gko::CudaExecutor::create(0, gko::OmpExecutor::create(), true)};
#elif defined(USE_OMP)
  auto gko_exec{gko::OmpExecutor::create()};
#else
  auto gko_exec{gko::ReferenceExecutor::create()};
#endif

  /* check input and set matrix dimensions */
  if (argc < 5)
  {
    std::cerr << "ERROR: FOUR (4) Inputs required:\n"
              << "  1) method\n"
              << "  2) number of matrix columns\n"
              << "  3) max iterations\n"
              << "  4) print timing\n";
    return 1;
  }

  std::string method{argv[++argi]};

  std::transform(method.begin(), method.end(), method.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (methods.count(method) == 0)
  {
    std::cerr << "ERROR: method must be one of ";
    for (const auto& m : methods)
    {
      std::cout << m.first << ", ";
    }
    std::cout << std::endl;
    return 1;
  }

  auto matcols{static_cast<const sunindextype>(atoll(argv[++argi]))};
  if (matcols <= 0)
  {
    std::cerr << "ERROR: number of matrix columns must be a positive integer\n";
    return 1;
  }
  auto matrows{matcols};

  auto max_iters{static_cast<const unsigned long>(atoll(argv[++argi]))};
  if (max_iters <= 0)
  {
    std::cerr << "ERROR: max iterations must be a positive integer\n";
    return 1;
  }

  int print_timing{atoi(argv[++argi])};
  SetTiming(print_timing);

  std::cout << "Ginkgo linear solver test:\n"
            << "  method     = " << method.c_str() << "\n"
            << "  size       = " << matrows << " x " << matcols << "\n"
            << "  max iters. = " << max_iters << "\n\n";

  /* Create vectors and matrices */
  std::default_random_engine engine;
  std::uniform_real_distribution<sunrealtype> distribution_real(8, 10);

  N_Vector x{N_VNew(matcols, sunctx)};
  N_Vector b{N_VClone(x)};
  N_Vector y{N_VClone(x)};

  auto matrix_dim{gko::dim<2>(matrows, matcols)};
  auto matrix_nnz{3 * matrows - 2};

#if defined(USE_CSR)
  auto gko_matrix_A = gko::share(GkoMatrixType::create(gko_exec, matrix_dim,
                                                       matrix_nnz));
  auto gko_matrix_B = gko::share(GkoMatrixType::create(gko_exec, matrix_dim,
                                                       matrix_nnz));
#else
  auto gko_matrix_A = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));
  auto gko_matrix_B = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));
#endif

  fill_matrix(gko::lend(gko_matrix_A));

  /* Fill x with random data */
  auto xdata{N_VGetArrayPointer(x)};
  for (sunindextype i = 0; i < matcols; i++)
  {
    xdata[i] = distribution_real(engine);
  }
  HIP_OR_CUDA(N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x));

  /* Wrap ginkgo matrices for SUNDIALS->
     Matrix is overloaded to a SUNMatrix. */
  sundials::ginkgo::Matrix<GkoMatrixType> A{gko_matrix_A, sunctx};
  sundials::ginkgo::Matrix<GkoMatrixType> B{gko_matrix_B, sunctx};

  /* Copy A and x into B and y to print in case of solver failure */
  SUNMatCopy(A, B);
  N_VScale(ONE, x, y);

  /* Create right-hand side vector for linear solve */
  fails += SUNMatMatvecSetup(A);
  fails += SUNMatMatvec(A, x, b);
  if (fails)
  {
    std::cerr << "FAIL: SUNLinSol SUNMatMatvec failure\n";
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);
    return 1;
  }

  /*
   * Create linear solver.
   */

  /* Use default stopping criteria */
  auto crit{sundials::ginkgo::DefaultStop::build()
            .with_max_iters(max_iters)
            .on(gko_exec)};

  /* Use a Jacobi preconditioner */
  auto precon{gko::preconditioner::Jacobi<sunrealtype, sunindextype>::build()
              .on(gko_exec)};


  std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;

  if (method == "bicg")
  {
    using GkoSolverType = gko::solver::Bicg<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, sundials::ginkgo::Matrix<GkoMatrixType>>>
      (std::move(gko_solver_factory), sunctx);
  }
  else if (method == "bicgstab")
  {
    using GkoSolverType = gko::solver::Bicgstab<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, sundials::ginkgo::Matrix<GkoMatrixType>>>
      (std::move(gko_solver_factory), sunctx);
  }
  else if (method == "cg")
  {
    using GkoSolverType = gko::solver::Cg<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, sundials::ginkgo::Matrix<GkoMatrixType>>>
      (std::move(gko_solver_factory), sunctx);
  }
  else if (method == "cgs")
  {
    using GkoSolverType = gko::solver::Cgs<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, sundials::ginkgo::Matrix<GkoMatrixType>>>
      (std::move(gko_solver_factory), sunctx);
  }
  else if (method == "fcg")
  {
    using GkoSolverType = gko::solver::Fcg<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, sundials::ginkgo::Matrix<GkoMatrixType>>>
      (std::move(gko_solver_factory), sunctx);
  }
  else if (method == "gmres")
  {
    using GkoSolverType = gko::solver::Gmres<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, sundials::ginkgo::Matrix<GkoMatrixType>>>
      (std::move(gko_solver_factory), sunctx);
  }
  else if (method == "idr")
  {
    using GkoSolverType = gko::solver::Idr<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, sundials::ginkgo::Matrix<GkoMatrixType>>>
      (std::move(gko_solver_factory), sunctx);
  }

  /* Run Tests */
  fails += Test_SUNLinSolGetID(LS->get(), SUNLINEARSOLVER_GINKGO, 0);
  fails += Test_SUNLinSolGetType(LS->get(), SUNLINEARSOLVER_MATRIX_ITERATIVE, 0);
  fails += Test_SUNLinSolInitialize(LS->get(), 0);
  fails += Test_SUNLinSolSetup(LS->get(), A, 0);
  fails += Test_SUNLinSolSolve(LS->get(), A, x, b, 1000 * std::numeric_limits<sunrealtype>::epsilon(), SUNTRUE, 0);

  /* Print result */
  if (fails)
  {
    std::cerr << "FAIL: SUNLinSol module failed " << fails << " tests\n\n";
    /* only print if its a small problem */
    if (matcols < 4)
    {
      std::cout << "\nx (original) =\n";
      N_VPrint(y);
      std::cout << "\nx (computed) =\n";
      N_VPrint(x);
      std::cout << "\nb =\n";
      N_VPrint(b);
    }
  }
  else
  {
    std::cout << "\nSUCCESS: SUNLinSol module passed all tests\n\n";
  }

  /* Print solve information */
  std::cout << "Number of linear solver iterations: "
            << static_cast<long int>(SUNLinSolNumIters(LS->get()))
            << std::endl;
  std::cout << "Final residual norm: "
            << SUNLinSolResNorm(LS->get())
            << std::endl;

  /* Free solver, matrix and vectors */
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(b);

  return fails;
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/

int check_vector(N_Vector expected, N_Vector actual, realtype tol)
{
  int failure{0};

  /* copy vectors to host */
  HIP_OR_CUDA(N_VCopyFromDevice_Hip(actual), N_VCopyFromDevice_Cuda(actual));
  HIP_OR_CUDA(N_VCopyFromDevice_Hip(expected), N_VCopyFromDevice_Cuda(expected));

  /* get vector data */
  auto xdata{N_VGetArrayPointer(actual)};
  auto ydata{N_VGetArrayPointer(expected)};

  /* check data lengths */
  auto xldata{N_VGetLength(actual)};
  auto yldata{N_VGetLength(expected)};

  if (xldata != yldata)
  {
    std::cerr << ">>> ERROR: check_vector: Different data array lengths\n";
    return 1;
  }

  /* check vector data */
  for (sunindextype i = 0; i < xldata; i++) {
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);
  }

  if (failure > ZERO)
  {
    std::cerr << "Check_vector failures:\n";
    for (sunindextype i = 0; i < xldata; i++)
    {
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
      {
        std::cerr << "  x[" << i << "] = " << xdata[i] << " != " << ydata[i]
                  << " (err = " << abs(xdata[i] - ydata[i]) << ")\n";
      }
    }
  }

  return failure > 0;
}

void sync_device()
{
  HIP_OR_CUDA(hipDeviceSynchronize(), cudaDeviceSynchronize());
}
