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

int resize_vec(N_Vector v_in, N_Vector v_out, void* user_data)
{
  sunrealtype* v_in_data  = N_VGetArrayPointer(v_in);
  sunrealtype* v_out_data = N_VGetArrayPointer(v_out);

  sunindextype N = N_VGetLocalLength(v_in);

  for (int i = 0; i < N; i++)
  {
    v_out_data[i] = v_in_data[i];
  }
  v_out_data[N] = v_in_data[N - 1];

  return 0;
}

// Print CVODE Nordsieck History Array
int PrintNordsieck(CVodeMem cv_mem)
{
  // const sunindextype N = N_VGetLength(tmp);
  const sunindextype N = 1;
  sunrealtype* vdata   = nullptr;

  std::cout << std::setw(4) << "idx" << std::setw(25) << "zn" << std::endl;
  for (int i = 0; i <= cv_mem->cv_qmax; i++)
  {
    vdata = N_VGetArrayPointer(cv_mem->cv_zn[i]);
    for (sunindextype j = 0; j < N; j++)
    {
      std::cout << std::setw(4) << i << std::setw(25) << vdata[j] << std::endl;
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

  int resize = 0;
  if (argc > 1)
  {
    resize = atoi(argv[1]);
  }

  int steps = 30;
  if (argc > 2)
  {
    steps = atoi(argv[2]);
  }

  if (resize == 0)
  {
    std::cout << "CVODE -- NO RESIZE" << std::endl;
  }
  else
  {
    std::cout << "CVODE -- RESIZE" << std::endl;
  }

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

  // Attach user data
  flag = CVodeSetUserData(cvode_mem, &udata);
  if (check_flag(flag, "CVodeSetUserData")) return 1;

  // Limit max order
  // flag = CVodeSetMaxOrd(cvode_mem, 2);
  // if (check_flag(flag, "CVodeSetMaxOrd")) return 1;

  // Initial time and final times
  sunrealtype tf = SUN_RCONST(10.0);

  // History of solution times
  sunrealtype tret;

  // Set output formatting
  std::cout << std::scientific;
  std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  std::cout << std::endl;

  // Access private CVODE memory to output Nordsieck array
  CVodeMem cv_mem = (CVodeMem) cvode_mem;

  N_Vector ytmp = N_VClone(y);
  if (check_ptr(ytmp, "N_VClone")) return 1;

  N_Vector* yhist = N_VCloneVectorArray(6, y);
  if (check_ptr(yhist, "N_VCloneVectorArray")) return 1;

  sunrealtype thist[6] = {ZERO, ZERO, ZERO, ZERO, ZERO, ZERO};

  int hist_size = 1;

  std::string file_name = "debug_resize_" + std::to_string(resize) + ".txt";
  FILE* debug_file = std::fopen(file_name.c_str(), "w");

  // Initialize saved history
  thist[0] = ZERO;
  N_VScale(ONE, y, yhist[0]);

  // Advance in time
  // 11 steps - reach 2nd order
  // 14 steps - reach 3rd order
  // 22 steps - reach 4th order
  // 27 steps - reach 5th order
  for (int i = 1; i <= steps; i++)
  {
    std::cout << std::flush;
    std::cerr << std::flush;
    std::cout << "\n========== Start Step " << i << " ==========\n\n";

    flag = CVode(cvode_mem, tf, y, &(tret), CV_ONE_STEP);
    if (check_flag(flag, "CVode")) return 1;

    std::cout << " Step Number: " << std::setw(3) << i
              << " | Time: " << std::setw(21) << tret
              << " | Step Size: " << std::setw(21) << cv_mem->cv_h
              << " | Order: " << cv_mem->cv_q << std::endl;

    // Print Nordsieck array (length q_max + 1)
    PrintNordsieck(cv_mem);
    if (check_flag(flag, "PrintNordsieck")) return 1;

    std::cout << "========== End Step " << i << " ==========\n";

    if (resize)
    {
      // std::cout << "\n========== Start Resize " << i << " ==========\n";
      // Test 2: Copy and expand the state

      // Update saved history
      if (i < 6)
      {
        hist_size++;
      }
      else
      {
        hist_size = 6;
      }

      for (int j = 5; j > 0; j--)
      {
        thist[j] = thist[j - 1];
      }
      thist[0] = tret;

      for (int j = 5; j > 0; j--)
      {
        N_VScale(ONE, yhist[j - 1], yhist[j]);
      }
      N_VScale(ONE, y, yhist[0]);

      // Resize all vectors
      for (int j = 0; j < 6; j++)
      {
        sunrealtype* old_data = N_VGetArrayPointer(yhist[j]);

        N_Vector new_vec = N_VNew_Serial(i + 1, sunctx);
        sunrealtype* new_data = N_VGetArrayPointer(new_vec);
        for (int k = 0; k < i + 1; k++)
        {
          new_data[k] = old_data[0];
        }
        N_VDestroy(yhist[j]);
        yhist[j] = new_vec;
      }

      N_VDestroy(y);
      N_VDestroy(ytmp);
      y = N_VNew_Serial(i + 1, sunctx);
      ytmp = N_VClone(y);

      flag = CVodeResizeHistory(cvode_mem, thist, yhist, hist_size, nullptr,
                                resize_vec, resize, debug_file);
      if (check_flag(flag, "CVodeResizeHistory")) return 1;

      // "Resize" the nonlinear solver
      SUNNonlinSolFree(NLS);
      NLS = SUNNonlinSol_FixedPoint(y, 2, sunctx);
      if (check_ptr(NLS, "SUNNonlinSol_FixedPoint")) return 1;

      flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
      if (check_flag(flag, "CVodeSetNonlinearSolver")) return 1;

      flag = CVodeSetMaxNonlinIters(cvode_mem, 10);
      if (check_flag(flag, "CVodeSetMaxNonlinIters")) return 1;

      // // Print Nordsieck array (length q_max + 1)
      // PrintNordsieck(tret, ytmp, cv_mem->cv_zn, cv_mem->cv_hscale, 6, udata);
      // if (check_flag(flag, "PrintNordsieck")) return 1;

      // std::cout << "\n========== End Resize " << i << " ==========\n";
    }
  }
  std::cout << std::endl;

  std::fclose(debug_file);

  // Print some final statistics
  flag = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "CVodePrintAllStats")) return 1;

  // Clean up and return with successful completion
  N_VDestroy(y);
  N_VDestroy(ytmp);
  N_VDestroyVectorArray(yhist, 6);
  SUNNonlinSolFree(NLS);
  CVodeFree(&cvode_mem);

  return 0;
}

/*---- end of file ----*/
