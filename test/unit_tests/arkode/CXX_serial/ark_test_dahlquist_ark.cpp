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
#define TWO     SUN_RCONST(2.0)

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
int run_tests(ARKodeButcherTable Be, ARKodeButcherTable Bi,
              ProblemData& prob_data, ProblemOptions& prob_opts,
              sundials::Context& sunctx);

int get_method_properties(ARKodeButcherTable Be, ARKodeButcherTable Bi,
                          int& stages, bool& explicit_first_stage,
                          bool& stiffly_accurate, bool& fsal);

int expected_rhs_evals(method_type type, int interp_type, int stages,
                       bool explicit_first_stage, bool stiffly_accurate,
                       bool fsal, void* arkstep_mem, long int& nfe_expected,
                       long int& nfi_expected);

int check_rhs_evals(void* arkstep_mem, long int nfe_expected,
                    long int nfi_expected);

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
            << "  relative tol = " << prob_opts.reltol << "\n"
            << "  absolute tol = " << prob_opts.abstol << "\n"
            << "  interp type  = " << prob_opts.interp_type << "\n";

  // Create SUNDIALS context
  sundials::Context sunctx;

  // Test methods
  int numfails = 0;

  ARKodeButcherTable Be = nullptr;
  ARKodeButcherTable Bi = nullptr;

  // Explicit
  std::cout << "\n========================\n"
            << "Test explicit RK methods\n"
            << "========================\n";

  // Explicit Euler
  Be          = ARKodeButcherTable_Alloc(1, SUNFALSE);
  Be->A[0][0] = ZERO;
  Be->b[0]    = ONE;
  Be->c[0]    = ZERO;
  Be->q       = 1;
  Bi          = nullptr;

  numfails += run_tests(Be, Bi, prob_data, prob_opts, sunctx);

  ARKodeButcherTable_Free(Be);
  ARKodeButcherTable_Free(Bi);
  Be = nullptr;
  Bi = nullptr;

  // Implicit
  std::cout << "\n========================\n"
            << "Test implicit RK methods\n"
            << "========================\n";

  // Implicit Euler
  Bi          = ARKodeButcherTable_Alloc(1, SUNFALSE);
  Bi->A[0][0] = ONE;
  Bi->b[0]    = ONE;
  Bi->c[0]    = ONE;
  Bi->q       = 1;
  Be          = nullptr;

  numfails += run_tests(Be, Bi, prob_data, prob_opts, sunctx);

  ARKodeButcherTable_Free(Be);
  ARKodeButcherTable_Free(Bi);
  Be = nullptr;
  Bi = nullptr;

  // IMEX

  std::cout << "\n=====================\n"
            << "Test IMEX ARK methods\n"
            << "=====================\n";

  // IMEX Euler
  Be          = ARKodeButcherTable_Alloc(2, SUNFALSE);
  Be->A[1][0] = ONE;
  Be->b[0]    = ONE;
  Be->c[1]    = ONE;
  Be->q       = 1;

  Bi          = ARKodeButcherTable_Alloc(2, SUNFALSE);
  Bi->A[1][1] = ONE;
  Bi->b[1]    = ONE;
  Bi->c[1]    = ONE;
  Bi->q       = 1;

  numfails += run_tests(Be, Bi, prob_data, prob_opts, sunctx);

  ARKodeButcherTable_Free(Be);
  ARKodeButcherTable_Free(Bi);
  Be = nullptr;
  Bi = nullptr;

  if (numfails) { std::cout << "\n\nFailed " << numfails << " tests!\n"; }
  else { std::cout << "\n\nAll tests passed!\n"; }

  // Return test status
  return numfails;
}

// -----------------------------------------------------------------------------
// Test drivers
// -----------------------------------------------------------------------------

int run_tests(ARKodeButcherTable Be, ARKodeButcherTable Bi,
              ProblemData& prob_data, ProblemOptions& prob_opts,
              sundials::Context& sunctx)
{
  // Reusable error-checking flag
  int flag;

  // Test failure counter
  int numfails = 0;

  // Determine method type
  method_type type;
  if (Be && !Bi) { type = method_type::expl; }
  else if (!Be && Bi) { type = method_type::impl; }
  else if (Be && Bi) { type = method_type::imex; }
  else
  {
    std::cerr << "ERROR: Both Butcher tables are NULL" << std::endl;
    return 1;
  }

  // Get method properties
  int stages;
  bool explicit_first_stage, stiffly_accurate, fsal;
  flag = get_method_properties(Be, Bi, stages, explicit_first_stage,
                               stiffly_accurate, fsal);
  if (check_flag(&flag, "get_method_properties", 1)) return 1;

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
  else
  {
    arkstep_mem = ARKStepCreate(fe, fi, prob_opts.t0, y, sunctx);
  }
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

  // Attach Butcher tables
  flag = ARKStepSetTables(arkstep_mem, 1, 0, Bi, Be);
  if (check_flag(&flag, "ARKStepSetTables", 1)) return 1;

  // --------------
  // Evolve in time
  // --------------

  sunrealtype t_ret = prob_opts.t0;
  sunrealtype t_out = 3 * prob_opts.h;

  long int nfe_expected, nfi_expected;

  for (int i = 0; i < 3; i++)
  {
    std::cout << "--------------------" << std::endl;

    // Advance in time
    flag = ARKStepEvolve(arkstep_mem, t_out, y, &t_ret, ARK_ONE_STEP);
    if (check_flag(&flag, "ARKStepEvolve", 1)) return 1;

    // Update output time
    t_out += prob_opts.h;

    // Check statistics
    flag = expected_rhs_evals(type, prob_opts.interp_type, stages,
                              explicit_first_stage, stiffly_accurate, fsal,
                              arkstep_mem,nfe_expected, nfi_expected);
    if (check_flag(&flag, "expected_rhs_evals", 1)) return 1;

    numfails += check_rhs_evals(arkstep_mem, nfe_expected, nfi_expected);

    if (numfails)
    {
      std::cout << "Failed " << numfails << " checks\n";
      break;
    }
  }
  std::cout << "--------------------" << std::endl;

  // ----------------
  // Get dense output
  // ----------------

  sunrealtype h_last;
  flag = ARKStepGetLastStep(arkstep_mem, &h_last);
  if (check_flag(&flag, "ARKStepGetLastStep", 1)) return 1;

  flag = ARKStepGetDky(arkstep_mem, t_ret - h_last / TWO, 0, y);
  if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;

  // Stiffly accurate (and FSAL) methods do not require an additional RHS
  // evaluation to get the new RHS value at the end of a step for dense output
  if (prob_opts.interp_type == 0 && !stiffly_accurate)
  {
    if (type == method_type::expl || type == method_type::imex)
    {
      nfe_expected++;
    }
    if (type == method_type::impl || type == method_type::imex)
    {
      nfi_expected++;
    }
  }
  numfails += check_rhs_evals(arkstep_mem, nfe_expected, nfi_expected);

  std::cout << "--------------------" << std::endl;

  // --------------------
  // Additional time step
  // --------------------

  // Advance in time
  flag = ARKStepEvolve(arkstep_mem, t_out, y, &t_ret, ARK_ONE_STEP);
  if (check_flag(&flag, "ARKStepEvolve", 1)) return 1;

  // Update output time
  t_out += prob_opts.h;

  // Check statistics
  flag = expected_rhs_evals(type, prob_opts.interp_type, stages,
                            explicit_first_stage, stiffly_accurate, fsal,
                            arkstep_mem, nfe_expected, nfi_expected);
  if (check_flag(&flag, "expected_rhs_evals", 1)) return 1;

  numfails += check_rhs_evals(arkstep_mem, nfe_expected, nfi_expected);

  std::cout << "--------------------" << std::endl;

  // --------
  // Clean up
  // --------

  ARKStepFree(&arkstep_mem);
  if (type == method_type::impl || type == method_type::imex)
  {
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
  }
  N_VDestroy(y);

  return numfails;
}

int get_method_properties(ARKodeButcherTable Be, ARKodeButcherTable Bi,
                          int& stages, bool& explicit_first_stage,
                          bool& stiffly_accurate, bool& fsal)
{
  stages = 0;
  if (Bi) { stages = Bi->stages; }
  else if (Be) { stages = Be->stages; }
  else
  {
    std::cerr << "ERROR: Both Butcher tables are NULL!" << std::endl;
    return 1;
  }

  // Check for explicit first stage
  explicit_first_stage = true;
  if (Bi)
  {
    if (std::abs(Bi->A[0][0]) > ZERO) { explicit_first_stage = false; }
  }
  if (Be)
  {
    if (std::abs(Be->A[0][0]) > ZERO) { explicit_first_stage = false; }
  }

  // Check for stiffly accurate method
  stiffly_accurate = true;
  if (Bi)
  {
    for (int i = 0; i < stages; i++)
    {
      if (std::abs(Bi->b[i] - Bi->A[stages - 1][i]) > ZERO)
      {
        stiffly_accurate = false;
      }
    }
  }
  if (Be)
  {
    for (int i = 0; i < stages; i++)
    {
      if (std::abs(Be->b[i] - Be->A[stages - 1][i]) > ZERO)
      {
        stiffly_accurate = false;
      }
    }
  }

  // Check for first same as last (FSAL) property
  fsal = explicit_first_stage && stiffly_accurate;

  return 0;
}

int expected_rhs_evals(method_type type, int interp_type, int stages,
                       bool explicit_first_stage, bool stiffly_accurate,
                       bool fsal, void* arkstep_mem,
                       long int& nfe_expected, long int& nfi_expected)
{
  int flag = 0;

  // Get number of steps and nonlinear solver iterations
  long int nst = 0;
  flag         = ARKStepGetNumSteps(arkstep_mem, &nst);
  if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return 1;

  long int nni = 0;
  if (type == method_type::impl || type == method_type::imex)
  {
    flag = ARKStepGetNumNonlinSolvIters(arkstep_mem, &nni);
    if (check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1)) return 1;
  }

  // Expected number of explicit functions evaluations
  nfe_expected = 0;
  if (type == method_type::expl || type == method_type::imex)
  {
    if (fsal)
    {
      // Save one function evaluation after first step
      nfe_expected = stages + (stages - 1) * (nst - 1);
    }
    else
    {
      nfe_expected = stages * nst;
    }

    if (interp_type == 0 && !explicit_first_stage)
    {
      if (stiffly_accurate)
      {
        // One extra evaluation in the first step only
        nfe_expected++;
      }
      else
      {
        // One extra evaluation in each step
        nfe_expected += nst;
      }
    }
  }

  // Expected number of implicit functions evaluations
  nfi_expected = 0;
  if (type == method_type::impl || type == method_type::imex)
  {
    if (fsal)
    {
      // Save one function evaluation after first step
      nfi_expected = stages + (stages - 1) * (nst - 1) + nni;
    }
    else
    {
      nfi_expected = stages * nst + nni;
    }

    if (interp_type == 0 && !explicit_first_stage)
    {
      if (stiffly_accurate)
      {
        // One extra evaluation in the first step only
        nfi_expected++;
      }
      else
      {
        // One extra evaluation in each step
        nfi_expected += nst;
      }
    }
  }

  std::cout << "Steps: " << nst << std::endl;
  std::cout << "NLS iters: " << nni << std::endl;

  return 0;
}

int check_rhs_evals(void* arkstep_mem, long int nfe_expected,
                    long int nfi_expected)
{
  int flag = 0;

  long int nst = 0;
  flag         = ARKStepGetNumSteps(arkstep_mem, &nst);
  if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return 1;

  long int nfe, nfi;
  flag = ARKStepGetNumRhsEvals(arkstep_mem, &nfe, &nfi);
  if (check_flag(&flag, "ARKStepGetNumRhsEvals", 1)) return 1;

  std::cout << "Fe RHS evals:\n"
            << "  actual:   " << nfe << "\n"
            << "  expected: " << nfe_expected << "\n";
  std::cout << "Fi RHS evals:\n"
            << "  actual:   " << nfi << "\n"
            << "  expected: " << nfi_expected << "\n";

  if (nfe != nfe_expected || nfi != nfi_expected)
  {
    std::cerr << ">>> Check failed <<<" << std::endl;
    return 1;
  }

  return 0;
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
