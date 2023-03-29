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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_mristep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

// Macros for problem constants
#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define TWO    SUN_RCONST(2.0)
#define TWENTY SUN_RCONST(20.0)

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute r(t)
static sunrealtype r(sunrealtype t) { return HALF * cos(t); }

// Compute the derivative of r(t)
static sunrealtype rdot(sunrealtype t) { return -HALF * sin(t); }

// Compute s(t)
static sunrealtype s(sunrealtype t) { return cos(TWENTY * t); }

// Compute the derivative of s(t)
static sunrealtype sdot(sunrealtype t) { return -TWENTY * sin(TWENTY * t); }

// Compute the true solution
static int ytrue(sunrealtype t, N_Vector y)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);

  ydata[0] = sqrt(ONE + r(t));
  ydata[1] = sqrt(TWO + s(t));

  return 0;
}

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

/* -----------------------------------------------------------------------------
 * Compute the slow ODE RHS function:
 *   [a  b] * [ (-1 + u^2 - r(t)) / (2*u) ] + [ r'(t) / (2u) ]
 *   [c  d]   [ (-2 + v^2 - s(t)) / (2*v) ]   [ s'(t) / (2v) ]
 * ---------------------------------------------------------------------------*/
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  const sunrealtype tmp1 = (-ONE + u * u - r(t)) / (TWO * u);
  const sunrealtype tmp2 = (-TWO + v * v - s(t)) / (TWO * v);

  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  fdata[0] = a * tmp1 + b * tmp2 + rdot(t) / (TWO * u);
  fdata[1] = c * tmp1 + d * tmp2 + sdot(t) / (TWO * v);

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the fast ODE RHS function
 * ---------------------------------------------------------------------------*/
int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(ZERO, ydot);
  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the slow ODE RHS Jacobin:
 *   [a/2 + (a(1+r(t))-rdot(t))/(2u^2)     b/2 + b*(2+s(t))/(2*v^2)         ]
 *   [c/2 + c(1+r(t))/(2u^2)               d/2 + (d(2+s(t))-sdot(t))/(2u^2) ]
 * ---------------------------------------------------------------------------*/
int J(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  Jdata[0] = a / TWO + (a * (ONE + r(t)) - rdot(t)) / (TWO * u * u);
  Jdata[1] = c / TWO + c * (ONE + r(t)) / (TWO * u * u);
  Jdata[2] = b / TWO + b * (TWO + s(t)) / (TWO * v * v);
  Jdata[3] = d / TWO + (d * (TWO + s(t)) - sdot(t)) / (TWO * v * v);

  return 0;
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Check function return flag
int check_flag(int flag, const std::string funcname)
{
  if (!flag) return 0;
  if (flag < 0) std::cerr << "ERROR: ";
  if (flag > 0) std::cerr << "WARNING: ";
  std::cerr << funcname << " returned " << flag << std::endl;
  return 1;
}

// Check if a function returned a NULL pointer
int check_ptr(void* ptr, const std::string funcname)
{
  if (ptr) return 0;
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}

// -----------------------------------------------------------------------------
// Test main
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Comparison tolerance
  sunrealtype tol = 100 * std::sqrt(SUN_UNIT_ROUNDOFF);
  if (argc > 1)
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    tol = std::stof(argv[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    tol = std::stod(argv[1]);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    tol = std::stold(argv[1]);
#else
#error "SUNDIALS precision macro not defined"
#endif
    if (tol <= ZERO)
    {
      std::cerr << "ERROR: Invalid tolerance, tol = " << tol << std::endl;
      return 1;
    }
  }

  // Integration tolerances
  const sunrealtype atol = 100 * SUN_UNIT_ROUNDOFF;
#if defined(SUNDIALS_SINGLE_PRECISION)
  const sunrealtype rtol = SUN_RCONST(1.0e-3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  const sunrealtype rtol = SUN_RCONST(1.0e-9);
#else
#error "SUNDIALS precision macro not defined"
#endif

  // Create initial condition
  N_Vector y = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) return 1;

  int flag = ytrue(ZERO, y);
  if (check_flag(flag, "ytrue")) return 1;

  // Create the fast integrator memory structure
  void* inner_arkode_mem = ARKStepCreate(f0, NULL, ZERO, y, sunctx);
  if (check_ptr(inner_arkode_mem, "ARKStepCreate")) return 1;

  // Create inner stepper
  MRIStepInnerStepper inner_stepper;
  flag = ARKStepCreateMRIStepInnerStepper(inner_arkode_mem,
                                          &inner_stepper);
  if (check_flag(flag, "ARKStepCreateMRIStepInnerStepper")) return 1;

  // Create MRIStep memory structure
  void* arkode_mem = MRIStepCreate(nullptr, f, ZERO, y, inner_stepper, sunctx);
  if (check_ptr(arkode_mem, "MRIStepCreate")) return 1;

  flag = MRIStepSStolerances(arkode_mem, rtol, atol);
  if (check_flag(flag, "MRIStepSStolerances")) return 1;

  flag = MRIStepSetFixedStep(arkode_mem, SUN_RCONST(1.0e-5));
  if (check_flag(flag, "MRIStepSetFixedStep")) return 1;

  SUNMatrix A = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(A, "SUNDenseMatrix")) return 1;

  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_ptr(LS, "SUNLinSol_Dense")) return 1;

  flag = MRIStepSetLinearSolver(arkode_mem, LS, A);
  if (check_flag(flag, "MRIStepSetLinearSolver")) return 1;

  sunrealtype udata[4] = {-TWO, HALF, HALF, -ONE};

  flag = MRIStepSetUserData(arkode_mem, udata);
  if (check_flag(flag, "MRIStepSetUserData")) return 1;

  // Initial time and fist output time
  sunrealtype tret = ZERO;
  sunrealtype tout = tret + SUN_RCONST(0.1);

  // Advance one step in time
  flag = MRIStepEvolve(arkode_mem, tout, y, &tret, ARK_ONE_STEP);
  if (check_flag(flag, "MRIStep")) return 1;

  // Get the internal finite difference approximation to J
  SUNMatrix Jdq;
  flag = MRIStepGetJac(arkode_mem, &Jdq);
  if (check_flag(flag, "MRIStepGetJac")) return 1;

  // Get the step and time at which the approximation was computed
  long int nst_Jdq;
  flag = MRIStepGetJacNumSteps(arkode_mem, &nst_Jdq);
  if (check_flag(flag, "MRIStepGetJacNumSteps")) return 1;

  sunrealtype t_Jdq;
  flag = MRIStepGetJacTime(arkode_mem, &t_Jdq);
  if (check_flag(flag, "MRIStepGetJacTime")) return 1;

  // Compute the true Jacobian
  SUNMatrix Jtrue = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(Jtrue, "SUNDenseMatrix")) return 1;

  flag = ytrue(t_Jdq, y);
  if (check_flag(flag, "ytrue")) return 1;

  flag = J(t_Jdq, y, nullptr, Jtrue, &udata, nullptr, nullptr, nullptr);
  if (check_flag(flag, "J")) return 1;

  // Compare finite difference and true Jacobian
  sunrealtype* Jdq_data = SUNDenseMatrix_Data(Jdq);
  if (check_ptr(Jdq_data, "SUNDenseMatrix_Data")) return 1;

  sunrealtype* Jtrue_data = SUNDenseMatrix_Data(Jtrue);
  if (check_ptr(Jtrue_data, "SUNDenseMatrix_Data")) return 1;

  // Output Jacobian data
  std::cout << std::scientific;
  std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  std::cout << "Jac nst = " << nst_Jdq << std::endl;
  std::cout << "Jac t   = " << t_Jdq << std::endl;
  std::cout << std::endl;
  std::cout << std::setw(8) << std::right << "Index" << std::setw(25)
            << std::right << "J DQ" << std::setw(25) << std::right << "J true"
            << std::setw(25) << std::right << "absolute difference"
            << std::setw(25) << std::right << "relative difference" << std::endl;
  for (int i = 0; i < 4 * 25 + 8; i++) std::cout << "-";
  std::cout << std::endl;

  int result = 0;
  sunindextype ldata = SUNDenseMatrix_LData(Jtrue);
  for (sunindextype i = 0; i < ldata; i++)
  {
    std::cout << std::setw(8) << std::right << i << std::setw(25) << std::right
              << Jdq_data[i] << std::setw(25) << std::right << Jtrue_data[i]
              << std::setw(25) << std::right
              << std::abs(Jdq_data[i] - Jtrue_data[i])
              << std::setw(25) << std::right
              << std::abs(Jdq_data[i] - Jtrue_data[i])/Jtrue_data[i]
              << std::endl;
    result += SUNRCompareTol(Jdq_data[i], Jtrue_data[i], tol);
  }

  // Clean up and return with successful completion
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNMatDestroy(Jtrue);
  SUNLinSolFree(LS);
  MRIStepInnerStepper_Free(&inner_stepper);
  ARKStepFree(&inner_arkode_mem);
  MRIStepFree(&arkode_mem);

  return result;
}

/*---- end of file ----*/
