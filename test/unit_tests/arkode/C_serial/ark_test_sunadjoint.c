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
#include <nvector/nvector_manyvector.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_core.h>
#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sunadjoint/sunadjoint_solver.h>
#include <sunmatrix/sunmatrix_dense.h>

static const sunrealtype params[4] = {1.5, 1.0, 3.0, 1.0};

int lotka_volterra(sunrealtype t, N_Vector uvec, N_Vector udotvec, void* user_data)
{
  sunrealtype* p    = (sunrealtype*)user_data;
  sunrealtype* u    = N_VGetArrayPointer(uvec);
  sunrealtype* udot = N_VGetArrayPointer(udotvec);

  udot[0] = p[0] * u[0] - p[1] * u[0] * u[1];
  udot[1] = -p[2] * u[1] + p[3] * u[0] * u[1];

  return 0;
}

int jacobian(sunrealtype t, N_Vector uvec, N_Vector udotvec, SUNMatrix Jac,
             void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* p = (sunrealtype*)user_data;
  sunrealtype* u = N_VGetArrayPointer(uvec);
  sunrealtype* J = SUNDenseMatrix_Data(Jac);

  J[0] = p[0] - p[1] * u[1];
  J[1] = p[3] * u[1];
  J[2] = -p[1] * u[0];
  J[3] = -p[2] * p[3] * u[0];

  return 0;
}

int parameter_jacobian(sunrealtype t, N_Vector uvec, N_Vector udotvec,
                       SUNMatrix Jac, void* user_data, N_Vector tmp1,
                       N_Vector tmp2, N_Vector tmp3)
{
  assert(user_data == params);

  sunrealtype* p = (sunrealtype*)user_data;
  sunrealtype* u = N_VGetArrayPointer(uvec);
  sunrealtype* J = SUNDenseMatrix_Data(Jac);

  J[0] = u[0];
  J[1] = 0.0;
  J[2] = -u[0] * u[1];
  J[3] = 0.0;
  J[4] = 0.0;
  J[5] = -u[1];
  J[6] = 0.0;
  J[7] = u[0] * u[1];

  return 0;
}

int forward_solution(SUNContext sunctx, void* arkode_mem,
                     SUNAdjointCheckpointScheme checkpoint_scheme,
                     sunrealtype t0, sunrealtype tf, N_Vector u)
{
  const sunrealtype dt = 1e-2;
  ARKodeSetUserData(arkode_mem, (void*)params);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetMaxNumSteps(arkode_mem, (tf-t0)/dt+1);

  sunrealtype t = t0;
  while (t < tf)
  {
    int flag = ARKodeEvolve(arkode_mem, tf, u, &t, ARK_NORMAL);
    if (flag < 0)
    {
      fprintf(stderr, ">>> ERROR: ARKodeEvolve returned %d\n", flag);
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
  // TODO(CJB): should we require ManyVector with separate vectors for IC and params to make things cleaner?
  sunindextype neq        = N_VGetLength(u);
  sunindextype num_cost   = 1;
  sunindextype num_params = 4;
  N_Vector sensu0         = N_VClone(u);
  N_Vector sensp          = N_VNew_Serial(num_params, sunctx);
  N_Vector sens[2]        = {sensu0, sensp};
  N_Vector sf             = N_VNew_ManyVector(2, sens, sunctx);

  // TODO(CJB): Load sf with the sensitivity terminal conditions
  N_VConst(0.0, sf);

  SUNAdjointSolver adj_solver;
  ARKStepCreateAdjointSolver(arkode_mem, num_cost, sf, &adj_solver);

  SUNMatrix J  = SUNDenseMatrix(neq, neq, sunctx);
  SUNMatrix Jp = SUNDenseMatrix(num_params, num_params, sunctx);

  SUNAdjointSolver_SetJacFn(adj_solver, jacobian, J);
  SUNAdjointSolver_SetJacPFn(adj_solver, parameter_jacobian, Jp);

  int stop_reason = 0;
  sunrealtype t   = tf;
  SUNAdjointSolver_Solve(adj_solver, tout, sf, &t, &stop_reason);

  fprintf(stdout, "Adjoint Solution:\n");
  N_VPrint(sf);

  N_VDestroy(sf);
  SUNMatDestroy(J);
  SUNMatDestroy(Jp);
  SUNAdjointSolver_Destroy(&adj_solver);

  return 0;
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
  ARKodeSetCheckpointScheme(arkode_mem, checkpoint_scheme);

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
  ARKodeFree(&arkode_mem);

  return 0;
}