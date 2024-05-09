/* -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Program to test the SUNAdjoint capability with ARKODE
 * ---------------------------------------------------------------------------*/

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_core.h>

#include "sunadjoint/sunadjoint_checkpointscheme.h"
#include "sunadjoint/sunadjoint_solver.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"

int lotka_volterra(sunrealtype t, N_Vector uvec, N_Vector udotvec, void* user_data)
{
  sunrealtype* p    = (sunrealtype*)user_data;
  sunrealtype* u    = N_VGetArrayPointer(uvec);
  sunrealtype* udot = N_VGetArrayPointer(udotvec);

  udot[0] = p[0] * u[0] - p[1] * u[0] * u[1];
  udot[1] = -p[2] * u[1] + p[3] * u[0] * u[1];

  return 0;
}

int forward_solution(SUNContext sunctx, void* arkode_mem,
                     SUNAdjointCheckpointScheme checkpoint_scheme,
                     sunrealtype t0, sunrealtype tf, N_Vector u)
{
  sunrealtype params[4] = {1.5, 1.0, 3.0, 1.0};
  ARKStepSetUserData(arkode_mem, (void*)params);

  ARKStepSStolerances(arkode_mem, 1e-4, 1e-10);

  sunrealtype t = t0;
  while (t < tf)
  {
    int flag = ARKStepEvolve(arkode_mem, tf, u, &t, ARK_NORMAL);
    if (flag < 0)
    {
      fprintf(stderr, ">>> ERROR: ARKStepEvolve returned %d\n", flag);
      return -1;
    }
  }

  fprintf(stdout, "Forward Solution:\n");
  N_VPrint(u);

  return 0;
}

int adjoint_solution(SUNContext sunctx, void* arkode_mem,
                     SUNAdjointCheckpointScheme checkpoint_scheme,
                     sunrealtype tf, sunrealtype tout, N_Vector u)
{
  // TODO(CJB): should we use ManyVector to separate IC sensitivities and param sensitivities?
  sunindextype neq        = N_VGetLength(u);
  sunindextype num_cost   = 1;
  sunindextype num_params = 4;
  N_Vector sf             = N_VNew_Serial(neq + num_params, sunctx);

  SUNStepper stepper = NULL;
  ARKStepCreateSUNStepper(arkode_mem, &stepper);
  SUNAdjointSolver adj_solver = NULL;
  SUNAdjointSolver_Create(stepper, num_cost, sf, checkpoint_scheme, sunctx,
                          &adj_solver);
  // SUNAdjointSolver_SetJacFn(adj_solver, );
  // SUNAdjointSolver_SetJacPFn(adj_solver, );

  int flag      = 0;
  sunrealtype t = tf;
  SUNAdjointSolver_Solve(adj_solver, tf, tout, sf, &t, &flag);

  fprintf(stdout, "Adjoint Solution:\n");
  N_VPrint(sf);

  N_VDestroy(sf);
  SUNStepper_Destroy(&stepper);
  SUNAdjointSolver_Destroy(&adj_solver);
}

int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  //
  // Create the initial conditions vector
  //

  sunindextype neq = 2;
  N_Vector u       = N_VNew_Serial(neq, sunctx);
  N_VConst(1.0, u);

  //
  // Create the ARKODE stepper that will be used for both the forward solution and adjoint solution.
  //

  sunrealtype t0   = 0.0;
  sunrealtype tf   = 10.0;
  void* arkode_mem = ARKStepCreate(lotka_volterra, NULL, t0, u, sunctx);

  // Enable checkpointing during the forward solution
  SUNAdjointCheckpointScheme checkpoint_scheme = NULL;
  // SUNAdjointCheckpointScheme_NewEmpty(sunctx, &checkpoint_scheme);
  // ARKStepSetCheckpointScheme(arkode_mem, checkpoint_scheme);

  //
  // Compute the forward solution
  //

  sunrealtype t = t0;
  forward_solution(sunctx, arkode_mem, checkpoint_scheme, t0, tf, u);

  //
  // Now compute the adjoint solution
  //

  adjoint_solution(sunctx, arkode_mem, checkpoint_scheme, tf, t0, u);

  //
  // Cleanup
  //

  N_VDestroy(u);
  ARKStepFree(&arkode_mem);

  return 0;
}