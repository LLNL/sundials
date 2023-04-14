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
 * Prothero-Robinson ODE test problem.
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

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include <cvode/cvode_impl.h>

#include "test_pr.hpp"
#include "test_utilities.hpp"

// ODE Rhs function wrapper
int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  PRData* udata = static_cast<PRData*>(user_data);
  return PR_Rhs(t, y, ydot, udata);
}

// ODE Jacobian function wrapper
int ode_jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  PRData* udata = static_cast<PRData*>(user_data);
  return PR_Jac(t, y, J, udata);
}

// Print CVODE Nordsieck History Array
int PrintNordsieck(int step, sunrealtype* thist, N_Vector tmp, N_Vector* zn,
                   sunrealtype step_size, int order, PRData& udata)
{
  const sunrealtype  tret = thist[step];
  const sunindextype N    = N_VGetLength(tmp);

  sunrealtype* tmpdata = N_VGetArrayPointer(tmp);
  sunrealtype* vdata   = N_VGetArrayPointer(zn[0]);
  sunrealtype  scale   = ONE;

  std::cout << std::setw(4) << "idx"
            << std::setw(25) << "scaled zn"
            << std::setw(25) << "true value"
            << std::setw(25) << "relative error"
            << std::endl;

  PR_true_dot(tret, tmp, 0, udata);
  for (sunindextype i = 0; i < N; i++)
  {
    std::cout << std::setw(4) << 0
              << std::setw(25) << scale * vdata[i]
              << std::setw(25) << tmpdata[i]
              << std::setw(25) << std::abs(scale * vdata[i] - tmpdata[i]) / std::abs(tmpdata[i])
              << std::endl;
  }

  scale = ONE / step_size;
  vdata = N_VGetArrayPointer(zn[1]);
  PR_true_dot(tret, tmp, 1, udata);
  for (sunindextype i = 0; i < N; i++)
  {
    std::cout << std::setw(4) << 1
              << std::setw(25) << scale * vdata[i]
              << std::setw(25) << tmpdata[i]
              << std::setw(25) << std::abs(scale * vdata[i] - tmpdata[i]) / std::abs(tmpdata[i])
              << std::endl;
  }


  if (order > 1)
  {
    scale = TWO / std::pow(step_size, 2);
    vdata = N_VGetArrayPointer(zn[2]);
    PR_true_dot(tret, tmp, 2, udata);
    for (sunindextype i = 0; i < N; i++)
    {
      std::cout << std::setw(4) << 2
                << std::setw(25) << scale * vdata[i]
                << std::setw(25) << tmpdata[i]
                << std::setw(25) << std::abs(scale * vdata[i] - tmpdata[i]) / std::abs(tmpdata[i])
                << std::endl;
    }
  }

  if (order > 2)
  {
    scale = SIX / std::pow(step_size, 3);
    vdata = N_VGetArrayPointer(zn[3]);
    PR_true_dot(tret, tmp, 3, udata);
    for (sunindextype i = 0; i < N; i++)
    {
      std::cout << std::setw(4) << 3
                << std::setw(25) << scale * vdata[i]
                << std::setw(25) << tmpdata[i]
                << std::setw(25) << std::abs(scale * vdata[i] - tmpdata[i]) / std::abs(tmpdata[i])
                << std::endl;
    }
  }

  if (order > 3)
  {
    scale = TWENTYFOUR / std::pow(step_size, 4);
    vdata = N_VGetArrayPointer(zn[4]);
    PR_true_dot(tret, tmp, 4, udata);
    for (sunindextype i = 0; i < N; i++)
    {
      std::cout << std::setw(4) << 4
                << std::setw(25) << scale * vdata[i]
                << std::setw(25) << tmpdata[i]
                << std::setw(25) << std::abs(scale * vdata[i] - tmpdata[i]) / std::abs(tmpdata[i])
                << std::endl;
    }
  }

  if (order > 4)
  {
    scale = ONEHUNDREDTWENTY / std::pow(step_size, 5);
    vdata = N_VGetArrayPointer(zn[5]);
    PR_true_dot(tret, tmp, 5, udata);
    for (sunindextype i = 0; i < N; i++)
    {
      std::cout << std::setw(4) << 5
                << std::setw(25) << scale * vdata[i]
                << std::setw(25) << tmpdata[i]
                << std::setw(25) << std::abs(scale * vdata[i] - tmpdata[i]) / std::abs(tmpdata[i])
                << std::endl;
    }
  }

  std::cout << std::endl;

  return 0;
}


// Test main
int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Create problem data structure
  PRData udata;

  // Create initial condition
  N_Vector y = N_VNew_Serial(1, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) return 1;

  int flag = PR_true(ZERO, y, udata);
  if (check_flag(flag, "PR_true")) return 1;

  // Create CVODE memory structure
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_ptr(cvode_mem, "CVodeCreate")) return 1;

  flag = CVodeInit(cvode_mem, ode_rhs, ZERO, y);
  if (check_flag(flag, "CVodeInit")) return 1;

  // Relative and absolute tolerances
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  const sunrealtype atol = SUN_RCONST(1.0e-10);

  flag = CVodeSStolerances(cvode_mem, rtol, atol);
  if (check_flag(flag, "CVodeSStolerances")) return 1;

  // Use fixed-point nonlinear solver
  SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(y, 2, sunctx);
  if (check_ptr(NLS, "SUNNonlinSol_FixedPoint")) return 1;

  flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if (check_flag(flag, "CVodeSetNonlinearSolver")) return 1;

  flag = CVodeSetMaxNonlinIters(cvode_mem, 10);
  if (check_flag(flag, "CVodeSetMaxNonlinIters")) return 1;

  // flag = CVodeSetNonlinConvCoef(cvode_mem, 0.0001);
  // if (check_flag(flag, "CVodeSetNonlinConvCoef")) return 1;

  // Attach user data
  flag = CVodeSetUserData(cvode_mem, &udata);
  if (check_flag(flag, "CVodeSetUserData")) return 1;

  // Initial time and final times
  sunrealtype tf   = SUN_RCONST(10.0);

  // History of solution times
  sunrealtype tret[31];
  tret[0] = ZERO;

  // Set output formatting
  std::cout << std::scientific;
  std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  std::cout << std::endl;

  // Access private CVODE memory to output Nordsieck array
  CVodeMem cv_mem = (CVodeMem) cvode_mem;

  N_Vector ytmp = N_VClone(y);
  if (check_ptr(ytmp, "N_VClone")) return 1;

  // Advance in time
  // 11 steps - reach 2nd order
  // 14 steps - reach 3rd order
  // 22 steps - reach 4th order
  // 27 steps - reach 5th order
  for (int i = 1; i <= 11; i++)
  {
    flag = CVode(cvode_mem, tf, y, &(tret[i]), CV_ONE_STEP);
    if (check_flag(flag, "CVode")) return 1;

    std::cout << " Step Number: " << std::setw(3) << i
              << " | Time: " << std::setw(21) << tret[i]
              << " | Step Size: " << std::setw(21) << cv_mem->cv_h
              << " | Order: " << cv_mem->cv_q << std::endl;

    // Print Nordsieck array
    PrintNordsieck(i, tret, ytmp, cv_mem->cv_zn, cv_mem->cv_h,
                   cv_mem->cv_q, udata);
    if (check_flag(flag, "PrintNordsieck")) return 1;

    // Test 0/1: Copy the Nordsieck array
    // N_VDestroy(y);
    // N_VDestroy(ytmp);
    // y = N_VNew_Serial(i + 1, sunctx);
    // ytmp = N_VClone(y);
    // N_Vector* znew = N_VCloneVectorArray(6, y);
    // for (int j = 0; j <= 5; j++)
    // {
    //   sunrealtype* zdata    = N_VGetArrayPointer(cv_mem->cv_zn[j]);
    //   sunrealtype* znewdata = N_VGetArrayPointer(znew[j]);
    //   for (int k = 0; k < i + 1; k++)
    //   {
    //     znewdata[k] = zdata[0];
    //   }
    // }
    // flag = CVodeResizeHistory(cvode_mem, &tret[i], znew, 0);
    // if (check_flag(flag, "CVodeResizeHistory")) return 1;
    // N_VDestroyVectorArray(znew, 6);

    // // "Resize" the nonlinear solver
    // SUNNonlinSolFree(NLS);
    // NLS = SUNNonlinSol_FixedPoint(y, 2, sunctx);
    // if (check_ptr(NLS, "SUNNonlinSol_FixedPoint")) return 1;

    // flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
    // if (check_flag(flag, "CVodeSetNonlinearSolver")) return 1;

    // flag = CVodeSetMaxNonlinIters(cvode_mem, 10);
    // if (check_flag(flag, "CVodeSetMaxNonlinIters")) return 1;

    // Test 2: Interpolate past values
    sunrealtype t_hist[4];
    N_Vector y_hist[4];
    int j_limit = std::min(i, 4);
    for (int j = 0; j < j_limit; j++)
    {
      t_hist[0] = tret[i-j];
    }

    flag = CVodeResizeHistory(cvode_mem, t_hist, y_hist, 0);
    if (check_flag(flag, "CVodeResizeHistory")) return 1;

    // "Resize" the nonlinear solver
    SUNNonlinSolFree(NLS);
    NLS = SUNNonlinSol_FixedPoint(y, 2, sunctx);
    if (check_ptr(NLS, "SUNNonlinSol_FixedPoint")) return 1;

    flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
    if (check_flag(flag, "CVodeSetNonlinearSolver")) return 1;

    flag = CVodeSetMaxNonlinIters(cvode_mem, 10);
    if (check_flag(flag, "CVodeSetMaxNonlinIters")) return 1;
  }
  std::cout << std::endl;

  // Print some final statistics
  flag = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "CVodePrintAllStats")) return 1;

  // Clean up and return with successful completion
  N_VDestroy(y);
  N_VDestroy(ytmp);
  SUNNonlinSolFree(NLS);
  CVodeFree(&cvode_mem);

  return 0;
}

/*---- end of file ----*/
