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
#include <cstddef>
#include <cstdio>
#include <exception>
#include <iostream>
#include <mdspan>
#include <type_traits>
#include <vector>

#include <arkode/arkode_erkstep.h>

#include "example_utilities.hpp"
#include "nvector_stdmdspan.hpp"

/* --------------------------------------------------------------------------
 * Type aliases
 *
 * The state is an N x 2 matrix (N oscillators, 2 components each).
 * layout_right gives row-major storage: state[i, 0] and state[i, 1] are
 * contiguous for each oscillator.
 * -------------------------------------------------------------------------- */

using Span2D    = std::mdspan<sunrealtype, std::dextents<std::size_t, 2>>;
using MdSpanVec = sundials::example::MdSpanVector<Span2D>;

static_assert(std::is_nothrow_move_constructible_v<MdSpanVec>,
              "MdSpanVec must be nothrow move constructible");

/* --------------------------------------------------------------------------
 * Problem parameters (passed through user_data)
 * -------------------------------------------------------------------------- */

struct OscParams
{
  int num_osc;          /* number of oscillators */
  sunrealtype coupling; /* coupling strength     */
};

/* --------------------------------------------------------------------------
 * RHS: coupled oscillators on a ring
 * -------------------------------------------------------------------------- */

static int rhs(sunrealtype /* t */, N_Vector y, N_Vector ydot, void* user_data)
{
  auto& params          = *static_cast<OscParams*>(user_data);
  int num_osc           = params.num_osc;
  sunrealtype coupling  = params.coupling;

  /* Recover the 2D mdspan from the N_Vector. */
  const auto& ys = MdSpanVec::FromNVectorConst(y)->span(); /* N x 2, read-only */
  auto& fs       = MdSpanVec::FromNVector(ydot)->span();   /* N x 2, write     */

  for (int i = 0; i < num_osc; ++i)
  {
    int ip = (i + 1) % num_osc;           /* periodic neighbor right */
    int im = (i - 1 + num_osc) % num_osc; /* periodic neighbor left  */

    /* dq_i/dt = p_i */
    fs[i, 0] = ys[i, 1];

    /* dp_i/dt = -q_i + k*(q_{i+1} - 2*q_i + q_{i-1}) */
    fs[i, 1] = -ys[i, 0] +
               coupling * (ys[ip, 0] - SUN_RCONST(2.0) * ys[i, 0] + ys[im, 0]);
  }

  return 0;
}

/* --------------------------------------------------------------------------
 * Compute total energy: E = sum_i [ 0.5*p_i^2 + 0.5*q_i^2
 *                                   + 0.5*k*(q_{i+1} - q_i)^2 ]
 * -------------------------------------------------------------------------- */

static sunrealtype totalEnergy(const Span2D& state, int num_osc,
                               sunrealtype coupling)
{
  sunrealtype energy = 0;
  for (int i = 0; i < num_osc; ++i)
  {
    int ip = (i + 1) % num_osc;
    energy += SUN_RCONST(0.5) * state[i, 1] * state[i, 1]; /* kinetic */
    energy += SUN_RCONST(0.5) * state[i, 0] * state[i, 0]; /* on-site */
    energy += SUN_RCONST(0.5) * coupling * (state[ip, 0] - state[i, 0]) *
              (state[ip, 0] - state[i, 0]); /* coupling */
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

  int ret = 0;

  /* Scope the vector and integrator so they are destroyed before the context.
   * The storage and vector constructors report failure by throwing. */
  try
  {
    /* Initial condition: N x 2 matrix via mdspan
     *
     * Storage holds the raw data; the mdspan provides 2D (i, component)
     * indexing over it; the MdSpanVector wraps the mdspan as an N_Vector. All
     * three refer to the same memory.
     */

    const int neq = params.num_osc * 2; /* total number of scalar equations */

    std::vector<sunrealtype> storage(neq, SUN_RCONST(0.0));
    Span2D state(storage.data(), params.num_osc, 2);

    /* Excite the first oscillator: q_0 = 1, all others zero. */
    state[0, 0] = SUN_RCONST(1.0);

    /* Wrap the mdspan as a SUNDIALS vector (view - no copy). */
    MdSpanVec y(state, sunctx);

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
          totalEnergy(state, params.num_osc, params.coupling);

        std::printf("  %8s  %14s  %14s  %10s\n", "t", "energy", "E0", "|dE/E0|");
        std::printf("  --------------------------------------------------------\n");
        std::printf("  %8.2f  %14.10f  %14.10f  %10.2e\n",
                    static_cast<double>(t0), static_cast<double>(energy0),
                    static_cast<double>(energy0), 0.0);

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
            totalEnergy(state, params.num_osc, params.coupling);
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
          std::printf("    q_%d = %+.8f\n", i, static_cast<double>(state[i, 0]));
        }

        if (check_flag(ARKodePrintAllStats(ark_mem, stdout,
                                           SUN_OUTPUTFORMAT_TABLE),
                       "ARKodePrintAllStats"))
        {
          ret = 1;
        }
      }

      ARKodeFree(&ark_mem);
    }
  }
  catch (const std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << std::endl;
    ret = 1;
  }

  SUNContext_Free(&sunctx);

  return ret;
}
