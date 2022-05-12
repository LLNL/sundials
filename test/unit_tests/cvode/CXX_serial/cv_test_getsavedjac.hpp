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

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_dense.h"

// Macros for problem constants
#define ZERO    RCONST(0.0)
#define HALF    RCONST(0.5)
#define ONE     RCONST(1.0)
#define TWO     RCONST(2.0)
#define TWENTY  RCONST(20.0)

using namespace std;

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
static int ytrue(realtype t, N_Vector y)
{
  realtype *ydata = N_VGetArrayPointer(y);

  ydata[0] = sqrt(ONE + r(t));
  ydata[1] = sqrt(TWO + s(t));

  return 0;
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Check function return flag
int check_flag(int flag, const string funcname)
{
  if (!flag) return 0;
  if (flag < 0) cerr << "ERROR: ";
  cerr << funcname << " returned " << flag << endl;
  return 1;
}

// Check if a function returned a NULL pointer
int check_ptr(void *ptr, const string funcname)
{
  if (ptr) return 0;
  cerr << "ERROR: " << funcname << " returned NULL" << endl;
  return 1;
}

/*---- end of file ----*/
