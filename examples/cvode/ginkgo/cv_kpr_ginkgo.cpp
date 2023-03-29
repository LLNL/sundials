/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Kvaerno-Prothero-Robinson ODE test problem:
 *
 *   [u]' = [ a  b ] [ (-1 + u^2 - r(t)) / (2u) ] + [ r'(t) / (2u) ]
 *   [v]    [ c  d ] [ (-2 + v^2 - s(t)) / (2v) ]   [ s'(t) / (2v) ]
 *
 * This problem has analytical solution given by
 *
 *   u(t) = sqrt(1 + r(t))
 *   v(t) = sqrt(2 + s(t))
 *
 * where, in this test, we use the functions
 *
 *   r(t) = 0.5 * cos(t)
 *   s(t) = cos(2t)
 * ---------------------------------------------------------------------------*/

#include "cv_kpr_ginkgo.hpp"

// Include integrator, matrix, linear solver, and vector headers
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_ginkgo.hpp>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
using GkoSolverType = gko::solver::Gmres<sunrealtype>;

using SUNGkoMatrixType = sundials::ginkgo::Matrix<GkoMatrixType>;
using SUNGkoSolverType = sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>;

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrators
// -----------------------------------------------------------------------------

// ODE right-hand side function
int f(realtype t, N_Vector y, N_Vector ydot, void* user_data);

// Jacobian of RHS function
int J(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Read input options
  Options opts;
  std::vector<std::string> args(argv + 1, argv + argc);

  int flag = ReadInputs(args, opts);
  if (check_flag(flag, "ReadInputs")) return 1;

  // Create initial condition vector
  N_Vector y = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) return 1;

  realtype utrue, vtrue;
  flag = true_sol(ZERO, &utrue, &vtrue);
  if (check_flag(flag, "true_sol")) return 1;

  realtype* ydata = N_VGetArrayPointer(y);
  ydata[0]        = utrue;
  ydata[1]        = vtrue;

  // Create Ginkgo dense matrix and linear solver
#if defined(USE_OMP)
  auto gko_exec{gko::OmpExecutor::create()};
#else
  auto gko_exec{gko::ReferenceExecutor::create()};
#endif

  auto gko_matrix = gko::share(GkoMatrixType::create(gko_exec, gko::dim<2>(2, 2)));
  SUNGkoMatrixType A{gko_matrix, sunctx};

  // Use default stopping criteria
  auto crit{sundials::ginkgo::DefaultStop::build().with_max_iters(gko::uint64{2}).on(gko_exec)};

  auto gko_solver_factory = gko::share(GkoSolverType::build().with_criteria(std::move(crit)).on(gko_exec));

  SUNGkoSolverType LS{gko_solver_factory, sunctx};

  // Create CVODE memory structure
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_ptr(cvode_mem, "CVodeCreate")) return 1;

  // Attach RHS function and set initial condition
  flag = CVodeInit(cvode_mem, f, ZERO, y);
  if (check_flag(flag, "CVodeInit")) return 1;

  // Set integraton tolerances
  flag = CVodeSStolerances(cvode_mem, opts.rtol, opts.atol);
  if (check_flag(flag, "CVodeSStolerances")) return 1;

  // Attach matrix and linear solver
  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_flag(flag, "CVodeSetLinearSolver")) return 1;

  // Set Jacobian function
  flag = CVodeSetJacFn(cvode_mem, J);
  if (check_flag(flag, "CVodeSetJacFn")) return 1;

  // Attach user data pointer
  realtype udata[4] = {-TWO, HALF, HALF, -ONE};
  flag              = CVodeSetUserData(cvode_mem, udata);
  if (check_flag(flag, "CVodeSetUserData")) return 1;

  // Initial time and fist output time
  realtype tret = ZERO;
  realtype tout = tret + opts.dtout;

  // Output initial contion
  std::cout << std::scientific;
  std::cout << std::setprecision(std::numeric_limits<realtype>::digits10);
  std::cout << "           t              ";
  std::cout << "          u              ";
  std::cout << "          v              ";
  std::cout << "        u err            ";
  std::cout << "        v err" << std::endl;
  for (int i = 0; i < 9; i++) std::cout << "--------------";
  std::cout << std::endl;

  std::cout << std::setw(22) << tret << std::setw(25) << ydata[0] << std::setw(25) << ydata[1] << std::setw(25)
            << std::abs(ydata[0] - utrue) << std::setw(25) << std::abs(ydata[1] - vtrue) << std::endl;

  // Advance in time
  for (int i = 0; i < opts.nout; i++) {
    flag = CVode(cvode_mem, tout, y, &tret, CV_NORMAL);
    if (check_flag(flag, "CVode")) return 1;

    flag = true_sol(tret, &utrue, &vtrue);
    if (check_flag(flag, "true_sol")) return 1;

    std::cout << std::setw(22) << tret << std::setw(25) << ydata[0] << std::setw(25) << ydata[1] << std::setw(25)
              << std::abs(ydata[0] - utrue) << std::setw(25) << std::abs(ydata[1] - vtrue) << std::endl;

    // update output time
    tout += opts.dtout;
  }
  for (int i = 0; i < 9; i++) std::cout << "--------------";
  std::cout << std::endl;

  // Print some final statistics
  flag = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "CVodePrintAllStats")) return 1;

  // Clean up and return with successful completion
  N_VDestroy(y);
  CVodeFree(&cvode_mem);

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the ODE RHS function:
 *   [a  b] * [ (-1 + u^2 - r(t)) / (2*u) ] + [ r'(t) / (2u) ]
 *   [c  d]   [ (-2 + v^2 - s(t)) / (2*v) ]   [ s'(t) / (2v) ]
 * ---------------------------------------------------------------------------*/

int f(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  realtype* udata  = (realtype*)user_data;
  const realtype a = udata[0];
  const realtype b = udata[1];
  const realtype c = udata[2];
  const realtype d = udata[3];

  realtype* ydata  = N_VGetArrayPointer(y);
  const realtype u = ydata[0];
  const realtype v = ydata[1];

  const realtype tmp1 = (-ONE + u * u - r(t)) / (TWO * u);
  const realtype tmp2 = (-TWO + v * v - s(t)) / (TWO * v);

  realtype* fdata = N_VGetArrayPointer(ydot);
  fdata[0]        = a * tmp1 + b * tmp2 + rdot(t) / (TWO * u);
  fdata[1]        = c * tmp1 + d * tmp2 + sdot(t) / (TWO * v);

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the ODE RHS Jacobin:
 *   [a/2 + (a(1+r(t))-rdot(t))/(2u^2)     b/2 + b*(2+s(t))/(2*v^2)         ]
 *   [c/2 + c(1+r(t))/(2u^2)               d/2 + (d(2+s(t))-sdot(t))/(2u^2) ]
 * ---------------------------------------------------------------------------*/

int J(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype* udata  = (realtype*)user_data;
  const realtype a = udata[0];
  const realtype b = udata[1];
  const realtype c = udata[2];
  const realtype d = udata[3];

  realtype* ydata = N_VGetArrayPointer(y);
  auto J_gko      = static_cast<SUNGkoMatrixType*>(J->content)->GkoMtx();
  auto J_data     = J_gko->get_values();

  const realtype u = ydata[0];
  const realtype v = ydata[1];

  // Ginkgo dense matrix is stored in row-major format
  J_data[0] = a / TWO + (a * (ONE + r(t)) - rdot(t)) / (TWO * u * u);
  J_data[1] = b / TWO + b * (TWO + s(t)) / (TWO * v * v);
  J_data[2] = c / TWO + c * (ONE + r(t)) / (TWO * u * u);
  J_data[4] = d / TWO + (d * (TWO + s(t)) - sdot(t)) / (TWO * v * v);

  return 0;
}

/*---- end of file ----*/
