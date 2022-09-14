/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details->
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the testing routine to check the Ginkgo SUNLinSol implementation.
 * -------------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------------- *
 * Matrix and solver options                                                  *
 * -------------------------------------------------------------------------- */

// "multigrid" does not support the combined stopping criteria we use
// "cbgmres" does not support setting stopping criteria
const std::unordered_map<std::string, int> methods
  {{"bicg", 0}, {"bicgstab", 1}, {"cg", 2}, {"cgs", 3},
   {"fcg", 4},  {"gmres", 5},    {"idr", 6}};

const std::unordered_map<std::string, int> matrix_types
  {{"csr", 0}, {"dense", 1}};

/* -------------------------------------------------------------------------- *
 * Matrix fill functions                                                      *
 * -------------------------------------------------------------------------- */

void fill_matrix(gko::matrix::Csr<sunrealtype, sunindextype>* matrix)
{
  const auto mat_rows = matrix->get_size()[0];
  const auto mat_cols = matrix->get_size()[1];
  auto row_ptrs = matrix->get_row_ptrs();
  auto col_idxs = matrix->get_col_idxs();
  auto mat_data = matrix->get_values();
  int idx = 0;

  // Matrix entries
  const sunrealtype vals[] = {-1, 2, -1};

  // Fill matrix
  row_ptrs[0] = idx;
  for (auto row = 0; row < mat_rows; ++row)
  {
    for (auto diag_offset : {-1, 0, 1})
    {
      auto col = row + diag_offset;
      if (0 <= col && col < mat_cols)
      {
        mat_data[idx] = vals[diag_offset + 1];
        col_idxs[idx] = col;
        ++idx;
      }
    }
    row_ptrs[row + 1] = idx;
  }
}

void fill_matrix(gko::matrix::Dense<sunrealtype>* matrix)
{
  const auto mat_rows = matrix->get_size()[0];
  const auto mat_cols = matrix->get_size()[1];

  // Matrix entries
  const sunrealtype vals[] = {-1, 2, -1};

  // Fill matrix
  for (auto row = 0; row < mat_rows; ++row)
  {
    for (auto diag_offset : {-1, 0, 1})
    {
      auto col = row + diag_offset;
      if (0 <= col && col < mat_cols)
      {
        matrix->at(row, col) = vals[diag_offset + 1];
      }
    }
  }
}

/* -------------------------------------------------------------------------- *
 * Extra test functions                                                       *
 * -------------------------------------------------------------------------- */

template<class GkoSolverType, class GkoMatrixType>
void Test_Move(std::unique_ptr<typename GkoSolverType::Factory>&& gko_solver_factory,
               sundials::Context& sunctx)
{
  // Move constructor
  sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType> solver{std::move(gko_solver_factory), sunctx};
  sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType> solver2{std::move(solver)};
  assert(solver2.gkofactory());
  assert(solver2.gkoexec());
  assert(SUNLinSolNumIters(solver2) == 0);

  // Move assignment
  sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType> solver3;
  solver3 = std::move(solver2);
  assert(solver3.gkofactory());
  assert(solver3.gkoexec());
  assert(SUNLinSolNumIters(solver3) == 0);

  std::cout << "    PASSED test -- Test_Move\n";
}

/* -------------------------------------------------------------------------- *
 * SUNLinSol_Ginkgo Testing Routine                                           *
 * -------------------------------------------------------------------------- */

int main(int argc, char* argv[])
{
  int argi{0};
  int fails{0}; /* counter for test failures */

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

  /* ------------ *
   * Check inputs *
   * ------------ */

  if (argc < 7)
  {
    std::cerr << "ERROR: SIX (6) inputs required:\n"
              << "  1) method\n"
              << "  2) matrix type\n"
              << "  3) number of matrix columns\n"
              << "  4) condition number\n"
              << "  5) max iterations\n"
              << "  6) print timing\n";
    return 1;
  }

  std::string method{argv[++argi]};
  std::transform(method.begin(), method.end(), method.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (!methods.count(method))
  {
    std::cerr << "ERROR: method must be one of ";
    for (const auto& m : methods)
    {
      std::cout << m.first << ", ";
    }
    std::cout << std::endl;
    return 1;
  }

  std::string matrix_type{argv[++argi]};
  std::transform(matrix_type.begin(), matrix_type.end(), matrix_type.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (!matrix_types.count(matrix_type))
  {
    std::cerr << "ERROR: matrix type must be one of ";
    for (const auto& m : matrix_types)
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

  auto matcond{static_cast<sunrealtype>(atof(argv[++argi]))};
  if (matcond < 0)
  {
    std::cerr << "ERROR: matrix condition number must be positive or 0 (poisson test)\n";
    return 1;
  }

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
            << "  matrix     = " << matrix_type.c_str() << "\n"
            << "  size       = " << matrows << " x " << matcols << "\n"
            << "  cond       = " << matcond << "\n"
            << "  max iters. = " << max_iters << "\n\n";

  /* ------------------------------- *
   * Create solution and RHS vectors *
   * ------------------------------- */

  N_Vector x{N_VNew(matcols, sunctx)};
  N_Vector b{N_VClone(x)};

  /* Fill x with random data */
  std::default_random_engine engine;
  std::uniform_real_distribution<sunrealtype> distribution_real(8, 10);

  auto xdata{N_VGetArrayPointer(x)};
  for (sunindextype i = 0; i < matcols; i++)
  {
    xdata[i] = distribution_real(engine);
  }
  HIP_OR_CUDA(N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x));

  /* -------------------- *
   * Create system matrix *
   * -------------------- */

  /* Wrap ginkgo matrix for SUNDIALS,
     Matrix is overloaded to a SUNMatrix. */

  std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> A;

  auto matrix_dim{gko::dim<2>(matrows, matcols)};

  if (matrix_type == "csr")
  {
    using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
    auto matrix_nnz{3 * matrows - 2};
    auto gko_matrix = gko::share(GkoMatrixType::create(gko_exec, matrix_dim,
                                                       matrix_nnz));

    if (matcond)
    {
      auto gko_matdata{gko::matrix_data<sunrealtype, sunindextype>::cond
                       (matrows, gko::remove_complex<sunrealtype>{matcond},
                        distribution_real, engine)};
      gko_matdata.remove_zeros();
      gko_matrix->read(gko_matdata);
    }
    else
    {
      fill_matrix(gko::lend(gko_matrix));
    }
    A = std::make_unique<sundials::ginkgo::Matrix<GkoMatrixType>>
      (std::move(gko_matrix), sunctx);
  }
  else if (matrix_type == "dense")
  {
    using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
    auto gko_matrix = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));
    if (matcond)
    {
      auto gko_matdata{gko::matrix_data<sunrealtype, sunindextype>::cond
                       (matrows, gko::remove_complex<sunrealtype>{matcond},
                        distribution_real, engine)};
      gko_matdata.remove_zeros();
      gko_matrix->read(gko_matdata);
    }
    else
    {
      fill_matrix(gko::lend(gko_matrix));
    }
    A = std::make_unique<sundials::ginkgo::Matrix<GkoMatrixType>>
      (std::move(gko_matrix), sunctx);
  }

  /* Create right-hand side vector for linear solve */
  fails += SUNMatMatvecSetup(A->get());
  fails += SUNMatMatvec(A->get(), x, b);
  if (fails)
  {
    std::cerr << "FAIL: SUNLinSol SUNMatMatvec failure\n";
    N_VDestroy(x);
    N_VDestroy(b);
    return 1;
  }

  /* -------------------- *
   * Create linear solver *
   * -------------------- */

  /* Use default stopping criteria */
  auto crit{sundials::ginkgo::DefaultStop::build()
            .with_max_iters(max_iters)
            .on(gko_exec)};

  /* Use a Jacobi preconditioner */
  auto precon{gko::preconditioner::Jacobi<sunrealtype, sunindextype>::build()
              .on(gko_exec)};


  /* Wrap ginkgo matrix for SUNDIALS,
     Matrix is overloaded to a SUNLinearSolver. */
  std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;

  if (method == "bicg")
  {
    using GkoSolverType = gko::solver::Bicg<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "bicgstab")
  {
    using GkoSolverType = gko::solver::Bicgstab<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "cg")
  {
    using GkoSolverType = gko::solver::Cg<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "cgs")
  {
    using GkoSolverType = gko::solver::Cgs<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "fcg")
  {
    using GkoSolverType = gko::solver::Fcg<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "gmres")
  {
    using GkoSolverType = gko::solver::Gmres<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "idr")
  {
    using GkoSolverType = gko::solver::Idr<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                            .with_criteria(std::move(crit))
                            .with_preconditioner(std::move(precon))
                            .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(), sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>>
        (std::move(gko_solver_factory), sunctx);
    }
  }

  /* Run Tests */
  fails += Test_SUNLinSolGetID(LS->get(), SUNLINEARSOLVER_GINKGO, 0);
  fails += Test_SUNLinSolGetType(LS->get(), SUNLINEARSOLVER_MATRIX_ITERATIVE, 0);
  fails += Test_SUNLinSolInitialize(LS->get(), 0);
  fails += Test_SUNLinSolSetup(LS->get(), A->get(), 0);
  fails += Test_SUNLinSolSolve(LS->get(), A->get(), x, b,
                               1000 * std::numeric_limits<sunrealtype>::epsilon(),
                               SUNTRUE, 0);

  /* Print result */
  if (fails)
  {
    std::cerr << "FAIL: SUNLinSol module failed " << fails << " tests\n\n";
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
  N_VDestroy(b);

  return fails;
}

/* -------------------------------------------------------------------------- *
 * Implementation-specific 'check' routines                                   *
 * -------------------------------------------------------------------------- */

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
