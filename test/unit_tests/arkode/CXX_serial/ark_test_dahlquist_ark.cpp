/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * IMEX multirate Dahlquist problem:
 *
 * y' = lambda_e * y + lambda_i * y
 * ---------------------------------------------------------------------------*/

// Header files
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>

#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

// Constants
#define NEG_ONE SUN_RCONST(-1.0)
#define ZERO    SUN_RCONST(0.0)
#define ONE     SUN_RCONST(1.0)

// Method types
enum class method_type
{
  expl,
  impl,
  imex
};

// Problem parameters
struct ProblemData
{
  sunrealtype lambda_e = NEG_ONE;
  sunrealtype lambda_i = NEG_ONE;
};

// Problem options
struct ProblemOptions
{
  // Initial time
  sunrealtype t0 = ZERO;

  // Number of time steps
  int nsteps = 3;

  // Relative and absolute tolerances
  sunrealtype reltol = SUN_RCONST(1.0e-4);
  sunrealtype abstol = SUN_RCONST(1.0e-6);

  // Step size
  sunrealtype h = SUN_RCONST(0.01);

  // Interpolant type
  // 0 = Hermite
  // 1 = Lagrange
  int interp_type = 0;
};

// User-supplied Functions called by the solver
int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private function to check function return values
int check_flag(void* flagvalue, const std::string funcname, int opt);

// Test drivers
int run_tests(method_type type, ProblemData& prob_data,
              ProblemOptions& prob_opts, sundials::Context& sunctx);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // Problem data and options structures
  ProblemData prob_data;
  ProblemOptions prob_opts;

  // Check for inputs
  if (argc > 1) { prob_opts.interp_type = std::stoi(argv[1]); }

  // Output problem setup
  std::cout << "\nDahlquist ODE test problem:\n"
            << "  lambda expl  = " << prob_data.lambda_e << "\n"
            << "  lambda impl  = " << prob_data.lambda_i << "\n"
            << "  step size    = " << prob_opts.h << "\n"
            << "  num steps    = " << prob_opts.nsteps << "\n"
            << "  relative tol = " << prob_opts.reltol << "\n"
            << "  absolute tol = " << prob_opts.abstol << "\n"
            << "  interp type  = " << prob_opts.interp_type << "\n";

  // Create SUNDIALS context
  sundials::Context sunctx;

  // Test methods
  int numfails = 0;

  // Explicit
  numfails += run_tests(method_type::expl, prob_data, prob_opts, sunctx);

  // Implicit
  numfails += run_tests(method_type::impl, prob_data, prob_opts, sunctx);

  // IMEX
  numfails += run_tests(method_type::imex, prob_data, prob_opts, sunctx);

  if (numfails) { std::cout << "\n\nFailed " << numfails << " tests!\n"; }
  else { std::cout << "\n\nAll tests passed!\n"; }

  // Return test status
  return numfails;
}

// -----------------------------------------------------------------------------
// Test drivers
// -----------------------------------------------------------------------------

int run_tests(method_type type, ProblemData& prob_data,
              ProblemOptions& prob_opts, sundials::Context& sunctx)
{
  // Reusable error-checking flag
  int flag;

  // Test failure counter
  int numfails = 0;

  // Create initial condition vector
  N_Vector y = N_VNew_Serial(1, sunctx);
  if (check_flag((void*)y, "N_VNew_Serial", 0)) return 1;

  N_VConst(SUN_RCONST(1.0), y);

  // Create matrix and linear solver (if necessary)
  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;

  if (type == method_type::impl || type == method_type::imex)
  {
    // Initialize dense matrix data structures and solvers
    A = SUNDenseMatrix(1, 1, sunctx);
    if (check_flag((void*)A, "SUNDenseMatrix", 0)) return 1;

    LS = SUNLinSol_Dense(y, A, sunctx);
    if (check_flag((void*)LS, "SUNLinSol_Dense", 0)) return 1;
  }

  // -----------------
  // Create integrator
  // -----------------

  // Create integrator based on type
  void* arkstep_mem = nullptr;

  if (type == method_type::expl)
  {
    arkstep_mem = ARKStepCreate(fe, nullptr, prob_opts.t0, y, sunctx);
  }
  else if (type == method_type::impl)
  {
    arkstep_mem = ARKStepCreate(nullptr, fi, prob_opts.t0, y, sunctx);
  }
  else if (type == method_type::imex)
  {
    arkstep_mem = ARKStepCreate(fe, fi, prob_opts.t0, y, sunctx);
  }
  else { return 1; }
  if (check_flag((void*)arkstep_mem, "ARKStepCreate", 0)) return 1;

  // Set user data
  flag = ARKStepSetUserData(arkstep_mem, &prob_data);
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

  // Specify tolerances
  flag = ARKStepSStolerances(arkstep_mem, prob_opts.reltol, prob_opts.abstol);
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  // Specify fixed time step size
  flag = ARKStepSetFixedStep(arkstep_mem, prob_opts.h);
  if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;

  // Lagrange interpolant (removes additional RHS evaluation with DIRK methods)
  if (prob_opts.interp_type == 1)
  {
    flag = ARKStepSetInterpolantType(arkstep_mem, ARK_INTERP_LAGRANGE);
    if (check_flag(&flag, "ARKStepSetInterpolantType", 1)) return 1;
  }

  if (type == method_type::impl || type == method_type::imex)
  {
    // Attach linear solver
    flag = ARKStepSetLinearSolver(arkstep_mem, LS, A);
    if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;

    // Set Jacobian function
    flag = ARKStepSetJacFn(arkstep_mem, Ji);
    if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;

    // Specify linearly implicit RHS, with non-time-dependent Jacobian
    flag = ARKStepSetLinear(arkstep_mem, 0);
    if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;
  }

  // ---------------------------
  // Evolve with various methods
  // ---------------------------

  // Methods to test
  int num_methods;

  if (type == method_type::expl)
  {
    std::cout << "\n========================\n"
              << "Test explicit RK methods\n"
              << "========================\n";
    num_methods = 1;
  }
  else if (type == method_type::impl)
  {
    std::cout << "\n========================\n"
              << "Test implicit RK methods\n"
              << "========================\n";
    num_methods = 1;
  }
  else if (type == method_type::imex)
  {
    std::cout << "\n=====================\n"
              << "Test IMEX ARK methods\n"
              << "=====================\n";
    num_methods = 1;
  }
  else { return 1; }

  for (int i = 0; i < num_methods; i++)
  {
    std::cout << "\nTesting method " << i << "\n";

    // -------------
    // Select method
    // -------------

    ARKodeButcherTable Be          = nullptr;
    ARKodeButcherTable Bi          = nullptr;
    int stages                     = 0;
    bool fsal_explicit_first_stage = false;
    bool implicit_first_stage      = false;

    if (type == method_type::expl)
    {
      // Explicit Euler
      stages                    = 1;
      fsal_explicit_first_stage = false;
      implicit_first_stage      = false;

      Be          = ARKodeButcherTable_Alloc(stages, SUNFALSE);
      Be->A[0][0] = ZERO;
      Be->b[0]    = ONE;
      Be->c[0]    = ZERO;
      Be->q       = 1;
      Bi          = nullptr;
    }
    else if (type == method_type::impl)
    {
      // Implicit Euler
      stages                    = 1;
      fsal_explicit_first_stage = false;
      implicit_first_stage      = true;

      Bi          = ARKodeButcherTable_Alloc(stages, SUNFALSE);
      Bi->A[0][0] = ONE;
      Bi->b[0]    = ONE;
      Bi->c[0]    = ONE;
      Bi->q       = 1;
      Be          = nullptr;
    }
    else if (type == method_type::imex)
    {
      // IMEX Euler
      stages                    = 2;
      fsal_explicit_first_stage = true;
      implicit_first_stage      = false;

      Be          = ARKodeButcherTable_Alloc(stages, SUNFALSE);
      Be->A[1][0] = ONE;
      Be->b[0]    = ONE;
      Be->c[1]    = ONE;
      Be->q       = 1;

      Bi          = ARKodeButcherTable_Alloc(stages, SUNFALSE);
      Bi->A[1][1] = ONE;
      Bi->b[1]    = ONE;
      Bi->c[1]    = ONE;
      Bi->q       = 1;
    }
    else { return 1; }

    // Attach Butcher tables
    flag = ARKStepSetTables(arkstep_mem, 1, 0, Bi, Be);
    if (check_flag(&flag, "ARKStepSetTables", 1)) return 1;

    ARKodeButcherTable_Free(Be);
    ARKodeButcherTable_Free(Bi);
    Be = nullptr;
    Bi = nullptr;

    // --------------
    // Evolve in time
    // --------------

    sunrealtype t  = prob_opts.t0;
    sunrealtype tf = prob_opts.nsteps * prob_opts.h;

    for (int i = 0; i < prob_opts.nsteps; i++)
    {
      // Advance in time
      flag = ARKStepEvolve(arkstep_mem, tf, y, &t, ARK_ONE_STEP);
      if (check_flag(&flag, "ARKStepEvolve", 1)) return 1;

      // Update output time
      tf += prob_opts.h;
    }

    // -----------------
    // Output statistics
    // -----------------

    long int nst, nfe, nfi;       // integrator
    long int nni, ncfn;           // nonlinear solver
    long int nsetups, nje, nfeLS; // linear solver

    flag = ARKStepGetNumSteps(arkstep_mem, &nst);
    if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return 1;

    flag = ARKStepGetNumRhsEvals(arkstep_mem, &nfe, &nfi);
    if (check_flag(&flag, "ARKStepGetNumRhsEvals", 1)) return 1;

    if (type == method_type::impl || type == method_type::imex)
    {
      flag = ARKStepGetNumNonlinSolvIters(arkstep_mem, &nni);
      if (check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1)) return 1;

      flag = ARKStepGetNumNonlinSolvConvFails(arkstep_mem, &ncfn);
      if (check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1)) return 1;

      flag = ARKStepGetNumLinSolvSetups(arkstep_mem, &nsetups);
      if (check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1)) return 1;

      flag = ARKStepGetNumJacEvals(arkstep_mem, &nje);
      if (check_flag(&flag, "ARKStepGetNumJacEvals", 1)) return 1;

      flag = ARKStepGetNumLinRhsEvals(arkstep_mem, &nfeLS);
      check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1);
    }

    sunrealtype pow = ZERO;
    if (type == method_type::expl || type == method_type::imex)
    {
      pow += prob_data.lambda_e;
    }
    if (type == method_type::impl || type == method_type::imex)
    {
      pow += prob_data.lambda_i;
    }
    sunrealtype ytrue = exp(pow * t);

    sunrealtype* ydata = N_VGetArrayPointer(y);
    sunrealtype error  = ytrue - ydata[0];

    std::cout << "\nARKStep Statistics:\n"
              << "  Time        = " << t << "\n"
              << "  y(t)        = " << ytrue << "\n"
              << "  y_n         = " << ydata[0] << "\n"
              << "  Error       = " << error << "\n"
              << "  Steps       = " << nst << "\n"
              << "  Fe evals    = " << nfe << "\n"
              << "  Fi evals    = " << nfi << "\n";

    if (type == method_type::impl || type == method_type::imex)
    {
      std::cout << "  NLS iters   = " << nni << "\n"
                << "  NLS fails   = " << ncfn << "\n"
                << "  LS setups   = " << nsetups << "\n"
                << "  LS Fi evals = " << nfeLS << "\n"
                << "  Ji evals    = " << nje << "\n";
    }

    // ----------------
    // Check statistics
    // ----------------

    // expected number of explicit functions evaluations
    if (type == method_type::expl || type == method_type::imex)
    {
      long int expected;

      if (fsal_explicit_first_stage)
      {
        expected = stages + (stages - 1) * (prob_opts.nsteps - 1);
      }
      else { expected = stages * prob_opts.nsteps; }

      if (prob_opts.interp_type == 0 && implicit_first_stage)
      {
        expected += prob_opts.nsteps;
      }

      if (nfe != expected)
      {
        numfails++;
        std::cout << "Fe RHS evals:\n"
                  << "  actual:   " << nfe << "\n"
                  << "  expected: " << expected << "\n";
      }
    }

    // expected number of implicit functions evaluations
    if (type == method_type::impl || type == method_type::imex)
    {
      long int expected;

      if (fsal_explicit_first_stage)
      {
        expected = stages + (stages - 1) * (prob_opts.nsteps - 1) + nni;
      }
      else { expected = stages * prob_opts.nsteps + nni; }

      if (prob_opts.interp_type == 0 && implicit_first_stage)
      {
        expected += prob_opts.nsteps;
      }

      if (nfi != expected)
      {
        numfails++;
        std::cout << "Fi RHS evals:\n"
                  << "  actual:   " << nfi << "\n"
                  << "  expected: " << expected << "\n";
      }
    }

    if (numfails) { std::cout << "Failed " << numfails << " checks\n"; }
    else { std::cout << "All checks passed\n"; }

    // -------------------
    // Setup for next test
    // -------------------

    // Free table(s)

    // Reset state vector to the initial condition
    N_VConst(SUN_RCONST(1.0), y);

    // Re-initialize integrator based on type
    if (type == method_type::expl)
    {
      flag = ARKStepReInit(arkstep_mem, fe, nullptr, prob_opts.t0, y);
    }
    else if (type == method_type::impl)
    {
      flag = ARKStepReInit(arkstep_mem, nullptr, fi, prob_opts.t0, y);
    }
    else if (type == method_type::imex)
    {
      flag = ARKStepReInit(arkstep_mem, fe, fi, prob_opts.t0, y);
    }
    else { return 1; }
    if (check_flag(&flag, "ARKStepReInit", 1)) return 1;
  }

  // Clean up
  ARKStepFree(&arkstep_mem);
  if (type == method_type::impl || type == method_type::imex)
  {
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
  }
  N_VDestroy(y);

  return numfails;
}

// -----------------------------------------------------------------------------
// Functions called by the solver
// -----------------------------------------------------------------------------

// Explicit ODE RHS function fe(t,y)
int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data    = N_VGetArrayPointer(y);
  sunrealtype* yd_data   = N_VGetArrayPointer(ydot);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  yd_data[0] = prob_data->lambda_e * y_data[0];

  return 0;
}

// Implicit ODE RHS function fi(t,y)
int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data    = N_VGetArrayPointer(y);
  sunrealtype* yd_data   = N_VGetArrayPointer(ydot);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  yd_data[0] = prob_data->lambda_i * y_data[0];

  return 0;
}

// Jacobian routine to compute J(t,y) = dfi/dy.
int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* J_data    = SUNDenseMatrix_Data(J);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  J_data[0] = prob_data->lambda_i;

  return 0;
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Check function return value
int check_flag(void* flagvalue, const std::string funcname, int opt)
{
  int* errflag;

  // Check if function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == nullptr)
  {
    std::cerr << "\nMEMORY_ERROR: " << funcname
              << " failed - returned NULL pointer\n\n";
    return 1;
  }
  // Check if flag < 0
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      std::cerr << "\nSUNDIALS_ERROR: " << funcname
                << " failed with flag = " << *errflag << "\n\n";
      return 1;
    }
  }

  return 0;
}
