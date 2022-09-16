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
 * Kvaerno-Prothero-Robinson ODE test problem, see .cpp file for details
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// Include integrator, matrix, linear solver, and vector headers
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>
#include <sunlinsol/sunlinsol_ginkgo.hpp>

// Macros for problem constants
#define ZERO    RCONST(0.0)
#define HALF    RCONST(0.5)
#define ONE     RCONST(1.0)
#define TWO     RCONST(2.0)
#define TWENTY  RCONST(20.0)

// -----------------------------------------------------------------------------
// Problem options
// -----------------------------------------------------------------------------

struct Options
{
  // Relative and absolute tolerances
  realtype rtol = RCONST(1.0e-6);
  realtype atol = RCONST(1.0e-10);

  // Output options
  realtype dtout = ONE; // output interval
  int      nout  = 10;  // number of outputs
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrators
// -----------------------------------------------------------------------------

// ODE right-hand side function
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

// Jacobian of RHS function
int J(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute r(t)
static realtype r(realtype t)
{
  return HALF * cos(t);
}

// Compute the derivative of r(t)
static realtype rdot(realtype t)
{
  return -HALF * sin(t);
}

// Compute s(t)
static realtype s(realtype t)
{
  return cos(TWENTY * t);
}

// Compute the derivative of s(t)
static realtype sdot(realtype t)
{
  return -TWENTY * sin(TWENTY * t);
}

// Compute the true solution
static int true_sol(realtype t, realtype* u, realtype* v)
{
  *u = sqrt(ONE + r(t));
  *v = sqrt(TWO + s(t));

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
  std::cerr << funcname << " returned " << flag << std::endl;
  return 1;
}

// Check if a function returned a NULL pointer
int check_ptr(void *ptr, const std::string funcname)
{
  if (ptr) return 0;
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}

inline void find_arg(std::vector<std::string> &args, const std::string key, realtype &dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    dest = stof(*(it + 1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    dest = stod(*(it + 1));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    dest = stold(*(it + 1));
#endif
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string> &args, const std::string key, long int &dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoll(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string> &args, const std::string key, int &dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoi(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string> &args, const std::string key, bool &dest,
                     bool store = true)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = store;
    args.erase(it);
  }
}

// Print command line options
void InputHelp()
{
  std::cout << std::endl;
  std::cout << "Command line options:" << std::endl;
  std::cout << "  --help         : print options and exit\n";
  std::cout << "  --rtol         : relative tolerance\n";
  std::cout << "  --atol         : absolute tolerance\n";
  std::cout << "  --dtout        : output interval\n";
  std::cout << "  --nout         : number of outputs\n";
}

int ReadInputs(std::vector<std::string> &args, Options &opts, SUNContext ctx)
{
  if (find(args.begin(), args.end(), "--help") != args.end())
  {
    InputHelp();
    return 1;
  }

  // Problem options
  find_arg(args, "--rtol", opts.rtol);
  find_arg(args, "--atol", opts.atol);
  find_arg(args, "--dtout", opts.dtout);
  find_arg(args, "--nout", opts.nout);

  return 0;
}

/*---- end of file ----*/
