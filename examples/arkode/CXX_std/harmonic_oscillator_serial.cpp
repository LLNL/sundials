/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Solves the 2D harmonic oscillator:
 *
 *   dy1/dt =  y2        y1(0) = 1
 *   dy2/dt = -y1        y2(0) = 0
 *
 * Exact solution: y1(t) = cos(t),  y2(t) = -sin(t).
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>

#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

#include "example_utilities.hpp"

/* --------------------------------------------------------------------------
 * Right-hand-side function
 * -------------------------------------------------------------------------- */

static int rhs(sunrealtype /* t */, N_Vector y, N_Vector ydot,
               void* /* user_data */)
{
  auto yd = N_VGetArrayPointer(y);
  auto fd = N_VGetArrayPointer(ydot);

  /* dy1/dt =  y2
   * dy2/dt = -y1 */
  fd[0] = yd[1];
  fd[1] = -yd[0];

  return 0;
}

/* --------------------------------------------------------------------------
 * Main
 * -------------------------------------------------------------------------- */

int main()
{
  SUNContext sunctx = nullptr;
  if (check_flag(SUNContext_Create(SUN_COMM_NULL, &sunctx), "SUNContext_Create"))
  {
    return 1;
  }

  const sunindextype neq = 2;                 /* number of equations */
  const sunrealtype t0   = SUN_RCONST(0.0);   /* initial time        */
  const sunrealtype tf   = SUN_RCONST(10.0);  /* final time          */
  const sunrealtype dt   = SUN_RCONST(1.0);   /* output interval     */

  N_Vector y = N_VNew_Serial(neq, sunctx);
  if (check_ptr(y, "N_VNew_Serial"))
  {
    SUNContext_Free(&sunctx);
    return 1;
  }

  int ret = 0;

  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0]           = SUN_RCONST(1.0); /* y1(0) = cos(0) = 1 */
  ydata[1]           = SUN_RCONST(0.0); /* y2(0) = sin(0) = 0 */

  void* arkode_mem = ERKStepCreate(rhs, t0, y, sunctx);
  if (check_ptr(arkode_mem, "ERKStepCreate")) { ret = 1; }
  else
  {
    const sunrealtype reltol = SUN_RCONST(1.0e-8);
    const sunrealtype abstol = SUN_RCONST(1.0e-12);
    ret = check_flag(ARKodeSStolerances(arkode_mem, reltol, abstol),
                     "ARKodeSStolerances");

    if (ret == 0)
    {
      std::printf("  %8s  %14s  %14s  %14s  %14s  %10s\n", "t", "y1", "y2",
                  "y1 exact", "y2 exact", "error");
      std::printf("  -------------------------------------------"
                  "---------------------------------------------------\n");

      /* Print initial condition */
      const sunrealtype err0 = std::hypot(ydata[0] - std::cos(t0),
                                          ydata[1] + std::sin(t0));
      std::printf("  %8.4f  %14.10f  %14.10f  %14.10f  %14.10f  %10.2e\n",
                  static_cast<double>(t0), static_cast<double>(ydata[0]),
                  static_cast<double>(ydata[1]),
                  static_cast<double>(std::cos(t0)),
                  static_cast<double>(-std::sin(t0)), static_cast<double>(err0));

      sunrealtype tret = t0;
      sunrealtype tout = t0 + dt;

      while (tout <= tf + SUN_RCONST(0.5) * dt)
      {
        if (check_flag(ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL),
                       "ARKodeEvolve"))
        {
          ret = 1;
          break;
        }

        /* Compare against exact solution */
        const sunrealtype y1_exact = std::cos(tret);
        const sunrealtype y2_exact = -std::sin(tret);
        const sunrealtype error = std::hypot(ydata[0] - y1_exact,
                                             ydata[1] - y2_exact);

        std::printf("  %8.4f  %14.10f  %14.10f  %14.10f  %14.10f  %10.2e\n",
                    static_cast<double>(tret), static_cast<double>(ydata[0]),
                    static_cast<double>(ydata[1]), static_cast<double>(y1_exact),
                    static_cast<double>(y2_exact), static_cast<double>(error));

        tout += dt;
      }

      if (check_flag(ARKodePrintAllStats(arkode_mem, stdout,
                                         SUN_OUTPUTFORMAT_TABLE),
                     "ARKodePrintAllStats"))
      {
        ret = 1;
      }
    }

    ARKodeFree(&arkode_mem);
  }

  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  return ret;
}
