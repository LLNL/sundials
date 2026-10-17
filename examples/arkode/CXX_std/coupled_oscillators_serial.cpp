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
 * Solves a system of N coupled oscillators arranged on a ring:
 *
 *   dq_i/dt =  p_i
 *   dp_i/dt = -q_i + k*(q_{i+1} - 2*q_i + q_{i-1})
 *
 * with periodic boundary conditions (indices mod N) and coupling k.
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>

#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

#include "example_utilities.hpp"

/* --------------------------------------------------------------------------
 * Problem parameters (passed through user_data)
 * -------------------------------------------------------------------------- */

struct OscParams
{
  int num_osc;          /* number of oscillators */
  sunrealtype coupling; /* coupling strength     */
};

/* --------------------------------------------------------------------------
 * Index layout convention
 *
 * The state is a flat array of length 2*N:
 *   y[2*i + 0] = q_i   (position of oscillator i)
 *   y[2*i + 1] = p_i   (momentum of oscillator i)
 * -------------------------------------------------------------------------- */

static inline sunindextype posIndex(int i) { return 2 * i; }

static inline sunindextype momIndex(int i) { return 2 * i + 1; }

/* --------------------------------------------------------------------------
 * RHS: coupled oscillators on a ring
 * -------------------------------------------------------------------------- */

static int rhs(sunrealtype /* t */, N_Vector y, N_Vector ydot, void* user_data)
{
  auto& params         = *static_cast<OscParams*>(user_data);
  int num_osc          = params.num_osc;
  sunrealtype coupling = params.coupling;

  /* Flat pointer access */
  const sunrealtype* ys = N_VGetArrayPointer(y);
  sunrealtype* fs       = N_VGetArrayPointer(ydot);

  for (int i = 0; i < num_osc; ++i)
  {
    int ip = (i + 1) % num_osc;           /* periodic neighbor right */
    int im = (i - 1 + num_osc) % num_osc; /* periodic neighbor left  */

    /* dq_i/dt = p_i */
    fs[posIndex(i)] = ys[momIndex(i)];

    /* dp_i/dt = -q_i + k*(q_{i+1} - 2*q_i + q_{i-1}) */
    fs[momIndex(i)] = -ys[posIndex(i)] +
                      coupling * (ys[posIndex(ip)] -
                                  SUN_RCONST(2.0) * ys[posIndex(i)] +
                                  ys[posIndex(im)]);
  }

  return 0;
}

/* --------------------------------------------------------------------------
 * Compute total energy: E = sum_i [ 0.5*p_i^2 + 0.5*q_i^2
 *                                   + 0.5*k*(q_{i+1} - q_i)^2 ]
 * -------------------------------------------------------------------------- */

static sunrealtype totalEnergy(const sunrealtype* ydata, int num_osc,
                               sunrealtype coupling)
{
  sunrealtype energy = 0;
  for (int i = 0; i < num_osc; ++i)
  {
    int ip = (i + 1) % num_osc;
    energy += SUN_RCONST(0.5) * ydata[momIndex(i)] * ydata[momIndex(i)]; /* kinetic */
    energy += SUN_RCONST(0.5) * ydata[posIndex(i)] * ydata[posIndex(i)]; /* on-site */
    energy += SUN_RCONST(0.5) * coupling *
              (ydata[posIndex(ip)] - ydata[posIndex(i)]) *
              (ydata[posIndex(ip)] - ydata[posIndex(i)]); /* coupling */
  }
  return energy;
}

/* --------------------------------------------------------------------------
 * Main
 * -------------------------------------------------------------------------- */

int main()
{
  /* Parameters */
  OscParams params;
  params.num_osc  = 8;                /* 8 oscillators on a ring */
  params.coupling = SUN_RCONST(0.25); /* moderate coupling       */

  const sunrealtype t0 = SUN_RCONST(0.0);
  const sunrealtype tf = SUN_RCONST(50.0);
  const sunrealtype dt = SUN_RCONST(5.0);

  SUNContext sunctx = nullptr;
  if (check_flag(SUNContext_Create(SUN_COMM_NULL, &sunctx), "SUNContext_Create"))
  {
    return 1;
  }

  const int neq = params.num_osc * 2; /* total number of scalar equations */

  N_Vector y = N_VNew_Serial(neq, sunctx);
  if (check_ptr(y, "N_VNew_Serial"))
  {
    SUNContext_Free(&sunctx);
    return 1;
  }

  int ret = 0;

  sunrealtype* ydata = N_VGetArrayPointer(y);

  /* Excite the first oscillator: q_0 = 1, all others zero. */
  N_VConst(SUN_RCONST(0.0), y);
  ydata[posIndex(0)] = SUN_RCONST(1.0);

  void* ark_mem = ERKStepCreate(rhs, t0, y, sunctx);
  if (check_ptr(ark_mem, "ERKStepCreate")) { ret = 1; }
  else
  {
    ret = check_flag(ARKodeSetUserData(ark_mem, &params), "ARKodeSetUserData");
    if (ret == 0)
    {
      ret = check_flag(ARKodeSStolerances(ark_mem, SUN_RCONST(1.0e-6),
                                          SUN_RCONST(1.0e-14)),
                       "ARKodeSStolerances");
    }

    if (ret == 0)
    {
      const sunrealtype energy0 =
        totalEnergy(ydata, params.num_osc, params.coupling);

      std::printf("  %8s  %14s  %14s  %10s\n", "t", "energy", "E0", "|dE/E0|");
      std::printf("  --------------------------------------------------------\n");
      std::printf("  %8.2f  %14.10f  %14.10f  %10.2e\n", static_cast<double>(t0),
                  static_cast<double>(energy0), static_cast<double>(energy0),
                  0.0);

      sunrealtype tret = t0;
      sunrealtype tout = t0 + dt;

      while (tout <= tf + SUN_RCONST(0.5) * dt)
      {
        if (check_flag(ARKodeEvolve(ark_mem, tout, y, &tret, ARK_NORMAL),
                       "ARKodeEvolve"))
        {
          ret = 1;
          break;
        }

        const sunrealtype energy =
          totalEnergy(ydata, params.num_osc, params.coupling);
        std::printf("  %8.2f  %14.10f  %14.10f  %10.2e\n",
                    static_cast<double>(tret), static_cast<double>(energy),
                    static_cast<double>(energy0),
                    static_cast<double>(std::abs((energy - energy0) / energy0)));

        tout += dt;
      }

      std::printf("\n  Final positions (q_i at t = %.1f):\n",
                  static_cast<double>(tret));
      for (int i = 0; i < params.num_osc; ++i)
      {
        std::printf("    q_%d = %+.8f\n", i,
                    static_cast<double>(ydata[posIndex(i)]));
      }

      if (check_flag(ARKodePrintAllStats(ark_mem, stdout, SUN_OUTPUTFORMAT_TABLE),
                     "ARKodePrintAllStats"))
      {
        ret = 1;
      }
    }

    ARKodeFree(&ark_mem);
  }

  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  return ret;
}
