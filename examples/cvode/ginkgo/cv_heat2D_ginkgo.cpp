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
 * Example problem:
 *
 * The following test simulates a simple anisotropic 2D heat equation,
 *
 *   u_t = kx u_xx + ky u_yy + b,
 *
 * for t in [0, 1] and (x,y) in [0, 1]^2, with initial conditions
 *
 *   u(0,x,y) = sin^2(pi x) sin^2(pi y),
 *
 * stationary boundary conditions
 *
 *   u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
 *
 * and the heat source
 *
 *   b(t,x,y) = -2 pi sin^2(pi x) sin^2(pi y) sin(pi t) cos(pi t)
 *              - kx 2 pi^2 (cos^2(pi x) - sin^2(pi x)) sin^2(pi y) cos^2(pi t)
 *              - ky 2 pi^2 (cos^2(pi y) - sin^2(pi y)) sin^2(pi x) cos^2(pi t).
 *
 * Under this setup, the problem has the analytical solution
 *
 *    u(t,x,y) = sin^2(pi x) sin^2(pi y) cos^2(pi t).
 *
 * The spatial derivatives are computed using second-order centered differences,
 * with the data distributed over nx * ny points on a uniform spatial grid. The
 * problem is advanced in time with BDF methods using an inexact Newton method
 * paired with the CG linear solver from Ginkgo.Several command line options are
 * available to change the problem parameters and CVODE settings. Use the flag
 * --help for more information.
 * ---------------------------------------------------------------------------*/

// Include user data structure and utility functions for this problem
#include "cv_heat2D_ginkgo.hpp"

// Include integrator, vector, matrix, and linear solver headers
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>
#include <sunlinsol/sunlinsol_ginkgo.hpp>

using GkoMatrixType = gko::matrix::Csr<sunrealtype>;
using GkoSolverType = gko::solver::Cg<sunrealtype>;

using SUNGkoMatrixType = sundials::ginkgo::Matrix<GkoMatrixType>;
using SUNGkoSolverType = sundials::ginkgo::LinearSolver<GkoSolverType,
                                                        GkoMatrixType>;

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// ODE right hand side function
int f(sunrealtype t, N_Vector u, N_Vector f, void *user_data);

// Jacobian of RHS function
int J(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // SUNDIALS context
  sundials::Context sunctx;

  // ---------------
  // Setup UserData
  // ---------------

  // Read input options
  UserData udata;
  std::vector<std::string> args(argv + 1, argv + argc);
  if (ReadInputs(args, udata)) return 1;
  PrintUserData(udata);

  // ---------------
  // Create vectors
  // ---------------

  // Create solution vector
  N_Vector u = N_VNew_Serial(udata.nodes, sunctx);
  if (check_ptr(u, "N_VNew_Parallel")) return 1;

  // Set initial condition
  int flag = Solution(ZERO, u, udata);
  if (check_flag(flag, "Solution")) return 1;

  // Create error vector
  N_Vector e = N_VClone(u);
  if (check_ptr(e, "N_VClone")) return 1;

  // ---------------------------------------
  // Create Ginkgo matrix and linear solver
  // ---------------------------------------

#if defined(USE_CUDA)
  auto gko_exec{gko::CudaExecutor::create(0, gko::OmpExecutor::create(), false,
                                          gko::allocation_mode::device)};
#elif defined(USE_HIP)
  auto gko_exec{gko::HipExecutor::create(0, gko::OmpExecutor::create(), false,
                                         gko::allocation_mode::device)};
#elif defined(USE_OMP)
  auto gko_exec{gko::OmpExecutor::create()};
#else
  auto gko_exec{gko::ReferenceExecutor::create()};
#endif

  auto gko_matrix_dim = gko::dim<2>(udata.nodes, udata.nodes);
  auto gko_matrix_nnz{(5 * (udata.nx - 2) + 2) * (udata.ny - 2) + 2 * udata.nx};
  auto gko_matrix = gko::share(GkoMatrixType::create(gko_exec, gko_matrix_dim,
                                                     gko_matrix_nnz));

  SUNGkoMatrixType A{gko_matrix, sunctx};

  // Use default stopping criteria
  auto crit{sundials::ginkgo::DefaultStop::build()
            .with_max_iters(static_cast<gko::uint64>(udata.liniters))
            .on(gko_exec)};

  auto gko_solver_factory = gko::share(GkoSolverType::build()
                                       .with_criteria(std::move(crit))
                                       .on(gko_exec));

  SUNGkoSolverType LS{gko_solver_factory, sunctx};

  // --------------
  // Setup CVODE
  // --------------

  // Create integrator
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_ptr(cvode_mem, "CVodeCreate")) return 1;

  // Initialize integrator
  flag = CVodeInit(cvode_mem, f, ZERO, u);
  if (check_flag(flag, "CVodeInit")) return 1;

  // Specify tolerances
  flag = CVodeSStolerances(cvode_mem, udata.rtol, udata.atol);
  if (check_flag(flag, "CVodeSStolerances")) return 1;

  // Attach user data
  flag = CVodeSetUserData(cvode_mem, (void *) &udata);
  if (check_flag(flag, "CVodeSetUserData")) return 1;

  // Attach linear solver
  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_flag(flag, "CVodeSetLinearSolver")) return 1;

  // Attach Jacobian function
  flag = CVodeSetJacFn(cvode_mem, J);
  if (check_flag(flag, "CVodeSetJacFn")) return 1;

  // Set linear solver tolerance factor
  flag = CVodeSetEpsLin(cvode_mem, udata.epslin);
  if (check_flag(flag, "CVodeSetEpsLin")) return 1;

  // Set max steps between outputs
  flag = CVodeSetMaxNumSteps(cvode_mem, udata.maxsteps);
  if (check_flag(flag, "CVodeSetMaxNumSteps")) return 1;

  // Set stopping time
  flag = CVodeSetStopTime(cvode_mem, udata.tf);
  if (check_flag(flag, "CVodeSetStopTime")) return 1;

  // -----------------------
  // Loop over output times
  // -----------------------

  sunrealtype t     = ZERO;
  sunrealtype dTout = udata.tf / udata.nout;
  sunrealtype tout  = dTout;

  // Inital output
  flag = OpenOutput(udata);
  if (check_flag(flag, "OpenOutput")) return 1;

  flag = WriteOutput(t, u, e, udata);
  if (check_flag(flag, "WriteOutput")) return 1;

  for (int iout = 0; iout < udata.nout; iout++)
  {
    // Evolve in time
    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if (check_flag(flag, "CVode")) break;

    // Output solution and error
    flag = WriteOutput(t, u, e, udata);
    if (check_flag(flag, "WriteOutput")) return 1;

    // Update output time
    tout += dTout;
    tout = (tout > udata.tf) ? udata.tf : tout;
  }

  // Close output
  flag = CloseOutput(udata);
  if (check_flag(flag, "CloseOutput")) return 1;

  // --------------
  // Final outputs
  // --------------

  // Print final integrator stats
  std::cout << "Final integrator statistics:" << std::endl;
  flag = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "CVodePrintAllStats")) return 1;

  // Output final error
  flag = SolutionError(t, u, e, udata);
  if (check_flag(flag, "SolutionError")) return 1;

  sunrealtype maxerr = N_VMaxNorm(e);

  std::cout << std::scientific
            << std::setprecision(std::numeric_limits<sunrealtype>::digits10)
            << "\nMax error = " << maxerr << std::endl;

  // --------------------
  // Clean up and return
  // --------------------

  CVodeFree(&cvode_mem); // Free integrator memory
  N_VDestroy(u);         // Free vectors
  N_VDestroy(e);

  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------

// f routine to compute the ODE RHS function f(t,y).
int f(sunrealtype t, N_Vector u, N_Vector f, void *user_data)
{
  // Access problem data and set shortcuts
  auto udata = static_cast<UserData*>(user_data);
  const auto nx = udata->nx;
  const auto ny = udata->ny;
  const auto dx = udata->dx;
  const auto dy = udata->dy;
  const auto kx = udata->kx;
  const auto ky = udata->ky;

  // Access data arrays
  sunrealtype *uarray = N_VGetArrayPointer(u);
  if (check_ptr(uarray, "N_VGetArrayPointer")) return -1;

  sunrealtype *farray = N_VGetArrayPointer(f);
  if (check_ptr(farray, "N_VGetArrayPointer")) return -1;

  // Constants for computing f(t,y)
  const sunrealtype cx = kx / (dx * dx);
  const sunrealtype cy = ky / (dy * dy);
  const sunrealtype cc = -TWO * (cx + cy);

  const auto bx = kx * TWO * PI * PI;
  const auto by = ky * TWO * PI * PI;

  const auto sin_t_cos_t = sin(PI * t) * cos(PI * t);
  const auto cos_sqr_t   = cos(PI * t) * cos(PI * t);

  // Initialize RHS vector to zero (handles boundary conditions)
  N_VConst(ZERO, f);

#if defined(USE_CUDA) || defined(USE_HIP)


#else
  // Iterate over domain interior and fill the RHS vector
  for (sunindextype j = 1; j < ny - 1; j++)
  {
    for (sunindextype i = 1; i < nx - 1; i++)
    {
      auto x  = i * dx;
      auto y  = j * dy;

      auto sin_sqr_x = sin(PI * x) * sin(PI * x);
      auto sin_sqr_y = sin(PI * y) * sin(PI * y);

      auto cos_sqr_x = cos(PI * x) * cos(PI * x);
      auto cos_sqr_y = cos(PI * y) * cos(PI * y);

      // center, north, south, east, and west indices
      auto idx_c = i + j * nx;
      auto idx_n = i + (j + 1) * nx;
      auto idx_s = i + (j - 1) * nx;
      auto idx_e = (i + 1) + j * nx;
      auto idx_w = (i - 1) + j * nx;

      farray[idx_c] =
        cc * uarray[idx_c]
        + cx * (uarray[idx_w] + uarray[idx_e])
        + cy * (uarray[idx_s] + uarray[idx_n])
        -TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t
        -bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t
        -by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
    }
  }
#endif

  // Return success
  return 0;
}

// J routine to compute the ODE RHS Jacobian function df/dy(t,y).
int J(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  auto udata = static_cast<UserData*>(user_data);
  const auto nx = udata->nx;
  const auto ny = udata->ny;
  const auto dx = udata->dx;
  const auto dy = udata->dy;
  const auto kx = udata->kx;
  const auto ky = udata->ky;

  // Constants for computing J(t,y)
  const sunrealtype cx = kx / (dx * dx);
  const sunrealtype cy = ky / (dy * dy);
  const sunrealtype cc = -TWO * (cx + cy);

  auto J_gko = static_cast<SUNGkoMatrixType*>(J->content)->gkomtx();

  sunindextype* row_ptrs = J_gko->get_row_ptrs();
  sunindextype* col_idxs = J_gko->get_col_idxs();
  sunrealtype*  mat_data = J_gko->get_values();

#if defined(USE_CUDA) || defined(USE_HIP)
#else
  // Explicitly set boundary entries to zero so J(t,y) has the same sparsity
  // pattern as A = I - gamma * J(t,y).

  // Fill southern boundary entries (j = 0)
  for (sunindextype i = 0; i < nx; i++)
  {
    mat_data[i] = ZERO;
    col_idxs[i] = i;
    row_ptrs[i] = i;
  }

  // Fill western boundary entries (i = 0)
  for (sunindextype j = 1; j < ny - 1; j++)
  {
    auto col = j * nx;
    auto idx = (5 * (nx - 2) + 2) * (j - 1) + nx;
    mat_data[idx] = ZERO;
    col_idxs[idx] = col;
    row_ptrs[col] = idx;
  }

  // Fill eastern boundary entries (i = nx - 1)
  for (sunindextype j = 1; j < ny - 1; j++)
  {
    auto col = (nx - 1) + j * nx;
    auto idx = (5 * (nx - 2) + 2) * (j - 1) + nx + 1 + 5 * (nx - 2);
    mat_data[idx] = ZERO;
    col_idxs[idx] = col;
    row_ptrs[col] = idx;
  }

  // Fill northern boundary entries (j = ny - 1)
  for (sunindextype i = 0; i < nx; i++)
  {
    auto col = i + (ny - 1) * nx;
    auto idx = (5 * (nx - 2) + 2) * (ny - 2) + nx + i;
    mat_data[idx] = ZERO;
    col_idxs[idx] = col;
    row_ptrs[col] = idx;
  }
  row_ptrs[nx * ny] = (5 * (nx - 2) + 2) * (ny - 2) + 2 * nx;

  // Fill interior entries
  for (sunindextype j = 1; j < ny - 1; j++)
  {
    for (sunindextype i = 1; i < nx - 1; i++)
    {
      auto row   = i + j * nx;
      auto col_s = row - nx;
      auto col_w = row - 1;
      auto col_c = row;
      auto col_e = row + 1;
      auto col_n = row + nx;

      // Number of non-zero entries from preceding rows
      auto prior_nnz = (5 * (nx - 2) + 2) * (j - 1) + nx;

      // Starting index for this row
      auto idx = prior_nnz + 1 + 5 * (i - 1);

      mat_data[idx]     = cy;
      mat_data[idx + 1] = cx;
      mat_data[idx + 2] = cc;
      mat_data[idx + 3] = cx;
      mat_data[idx + 4] = cy;

      col_idxs[idx]     = col_s;
      col_idxs[idx + 1] = col_w;
      col_idxs[idx + 2] = col_c;
      col_idxs[idx + 3] = col_e;
      col_idxs[idx + 4] = col_n;

      row_ptrs[row] = idx;
    }
  }
#endif

  // auto nnz =  (5 * (nx - 2) + 2) * (ny - 2) + 2 * nx;
  // std::cout << "mat_data" << std::endl;
  // for (auto i = 0; i < nnz; i++)
  //   std::cout << i << " : " << mat_data[i] << std::endl;
  // std::cout << "col_idxs" << std::endl;
  // for (auto i = 0; i < nnz; i++)
  //   std::cout << i << " : " << col_idxs[i] << std::endl;
  // std::cout << "row_ptrs" << std::endl;
  // for (auto i = 0; i <= nx * ny; i++)
  //   std::cout << i << " : " << row_ptrs[i] << std::endl;

  return 0;
}

//---- end of file ----
