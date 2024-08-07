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
 * Program to test the SUNAdjoint capability with ARKODE. The test uses the
 * Lotka-Volterra problem with four parameters as the test case.
 * ---------------------------------------------------------------------------*/

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sunadjoint/sunadjoint_checkpointscheme_basic.h>
#include <sunadjoint/sunadjoint_solver.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "sundials/sundials_context.h"
#include "sundials/sundials_datanode.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"
#include "sunmemory/sunmemory_system.h"

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
  J[2] = -p[1] * u[0];
  J[1] = p[3] * u[1];
  J[3] = p[3] * u[0] - p[2];

  return 0;
}

int jvp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
        N_Vector udotvec, void* user_data, N_Vector tmp)
{
  sunrealtype* p  = (sunrealtype*)user_data;
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  Jv[0] = (p[0] - p[1] * u[1]) * v[0] + (-p[1] * u[0]) * v[1];
  Jv[1] = (p[3] * u[1]) * v[0] + (-p[2] * p[3] * u[0]) * v[1];

  return 0;
}

int parameter_jacobian(sunrealtype t, N_Vector uvec, N_Vector udotvec,
                       SUNMatrix Jac, void* user_data, N_Vector tmp1,
                       N_Vector tmp2, N_Vector tmp3)
{
  if (user_data != params) { return -1; }

  sunrealtype* p = (sunrealtype*)user_data;
  sunrealtype* u = N_VGetArrayPointer(uvec);
  sunrealtype* J = SUNDenseMatrix_Data(Jac);

  // TODO(CJB): this isnt right - matrix is 4x4
  J[0] = u[0];
  J[1] = 0.0;
  J[2] = -u[0] * u[1];
  J[3] = 0.0;
  J[4] = 0.0;
  J[5] = -u[1];
  J[6] = 0.0;
  J[7] = u[0] * u[1];
  J[8] = 0.0;

  return 0;
}

int parameter_jvp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
                  N_Vector udotvec, void* user_data, N_Vector tmp)
{
  if (user_data != params) { return -1; }

  sunrealtype* p  = (sunrealtype*)user_data;
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  // TODO(CJB): this isnt right
  Jv[0] = v[0] * u[0];
  Jv[1] = 0.0;
  Jv[2] = 0.0;
  Jv[3] = 0.0;

  return 0;
}

sunrealtype g(N_Vector u, const sunrealtype* p, sunrealtype t)
{
  /* (sum(u) .^ 2) ./ 2 */
  sunrealtype* uarr = N_VGetArrayPointer(u);
  sunrealtype sum   = SUN_RCONST(0.0);
  for (sunindextype i = 0; i < N_VGetLength(u); i++) { sum += uarr[i]; }
  return (sum * sum) / SUN_RCONST(2.0);
}

void dgdu(N_Vector uvec, N_Vector dgvec, const sunrealtype* p, sunrealtype t)
{
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* dg = N_VGetArrayPointer(dgvec);

  dg[0] = u[0] + u[1];
  dg[1] = u[0] + u[1];
}

int forward_solution(SUNContext sunctx, void* arkode_mem,
                     SUNAdjointCheckpointScheme checkpoint_scheme,
                     sunrealtype t0, sunrealtype tf, sunrealtype dt, N_Vector u)
{
  ARKodeSetUserData(arkode_mem, (void*)params);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetMaxNumSteps(arkode_mem, (tf - t0) / dt + 1);

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
  int retval = 0;

  sunindextype neq        = N_VGetLength(u);
  sunindextype num_cost   = 1;
  sunindextype num_params = 4;
  N_Vector sensu0         = N_VClone(u);
  N_Vector sensp          = N_VNew_Serial(num_params, sunctx);
  N_Vector sens[2]        = {sensu0, sensp};
  N_Vector sf             = N_VNew_ManyVector(2, sens, sunctx);

  // Set the terminal condition for the adjoint system, which
  // should be the the gradient of our cost function at tf.
  dgdu(u, sensu0, params, tf);
  // dgdp()
  N_VConst(0.0, sensp);

  fprintf(stdout, "Adjoint terminal condition:\n");
  N_VPrint(sf);

  SUNAdjointSolver adj_solver;
  retval = ARKStepCreateAdjointSolver(arkode_mem, num_cost, sf, &adj_solver);

  SUNMatrix J  = SUNDenseMatrix(neq, neq, sunctx);
  SUNMatrix Jp = SUNDenseMatrix(num_params, num_params, sunctx);

  retval = SUNAdjointSolver_SetJacFn(adj_solver, jacobian, J,
                                     parameter_jacobian, Jp);

  int stop_reason = 0;
  sunrealtype t   = tf;
  retval = SUNAdjointSolver_Solve(adj_solver, tout, sf, &t, &stop_reason);
  if (stop_reason < 0 || stop_reason > 2)
  {
    fprintf(stderr, "SUNAdjointSolver_Solve stopped with reason %d\n",
            stop_reason);
    return -1;
  }

  fprintf(stdout, "Adjoint Solution:\n");
  N_VPrint(sf);

  N_VDestroy(sf);
  SUNMatDestroy(J);
  SUNMatDestroy(Jp);
  SUNAdjointSolver_Destroy(&adj_solver);

  return 0;
}

int adjoint_solution_jvp(SUNContext sunctx, void* arkode_mem,
                         SUNAdjointCheckpointScheme checkpoint_scheme,
                         sunrealtype tf, sunrealtype tout, N_Vector u)
{
  int retval = 0;

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
  retval = ARKStepCreateAdjointSolver(arkode_mem, num_cost, sf, &adj_solver);

  retval = SUNAdjointSolver_SetJacTimesVecFn(adj_solver, jvp, parameter_jvp);

  int stop_reason = 0;
  sunrealtype t   = tf;
  retval = SUNAdjointSolver_Solve(adj_solver, tout, sf, &t, &stop_reason);

  fprintf(stdout, "Adjoint Solution:\n");
  N_VPrint(sf);

  N_VDestroy(sf);
  SUNAdjointSolver_Destroy(&adj_solver);

  return 0;
}

int adjoint_solution_vjp(SUNContext sunctx, void* arkode_mem,
                         SUNAdjointCheckpointScheme checkpoint_scheme,
                         sunrealtype tf, sunrealtype tout, N_Vector u)
{
  int retval = 0;

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
  retval = ARKStepCreateAdjointSolver(arkode_mem, num_cost, sf, &adj_solver);

  retval = SUNAdjointSolver_SetVecTimesJacFn(adj_solver, jvp, parameter_jvp);

  int stop_reason = 0;
  sunrealtype t   = tf;
  retval = SUNAdjointSolver_Solve(adj_solver, tout, sf, &t, &stop_reason);

  fprintf(stdout, "Adjoint Solution:\n");
  N_VPrint(sf);

  N_VDestroy(sf);
  SUNAdjointSolver_Destroy(&adj_solver);

  return 0;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  // Since this a unit test, we want to abort immediately on any internal error
  SUNContext_PushErrHandler(sunctx, SUNAbortErrHandlerFn, NULL);

  //
  // Create the initial conditions vector
  //

  sunindextype neq = 2;
  N_Vector u       = N_VNew_Serial(neq, sunctx);
  N_VConst(1.0, u);

  //
  // Create the ARKODE stepper that will be used for both the forward solution and adjoint solution.
  //

  const sunrealtype dt = 1e-4;
  sunrealtype t0       = 0.0;
  sunrealtype tf       = 1.0;
  void* arkode_mem     = ARKStepCreate(lotka_volterra, NULL, t0, u, sunctx);

  // Use Forward Euler for debugging
  ARKodeSetOrder(arkode_mem, 1);

  // Enable checkpointing during the forward solution
  SUNAdjointCheckpointScheme checkpoint_scheme = NULL;

  SUNMemoryHelper mem_helper = SUNMemoryHelper_Sys(sunctx);
  SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, mem_helper, 1,
                                          ((tf - t0) / dt + 1) * 6, SUNTRUE,
                                          SUNTRUE, sunctx, &checkpoint_scheme);

  ARKodeSetCheckpointScheme(arkode_mem, checkpoint_scheme);

  //
  // Compute the forward solution
  //

  printf("Initial condition:\n");
  N_VPrint(u);

  sunrealtype t = t0;
  forward_solution(sunctx, arkode_mem, checkpoint_scheme, t0, tf, dt, u);

  //
  // Now compute the adjoint solution
  //

  adjoint_solution(sunctx, arkode_mem, checkpoint_scheme, tf, t0, u);

  // //
  // // Now compute the adjoint solution using Jvp
  // //
  // // TODO(CJB): make sure this reinitializes arkode correctly (probably need SUNAdjointSolver_Reset function)
  // adjoint_solution_jvp(sunctx, arkode_mem, checkpoint_scheme, tf, t0, u);

  // //
  // // Now compute the adjoint solution using vJp
  // //
  // // TODO(CJB): make sure this reinitializes arkode correctly
  // adjoint_solution_vjp(sunctx, arkode_mem, checkpoint_scheme, tf, t0, u);

  //
  // Cleanup
  //

  N_VDestroy(u);
  ARKodeFree(&arkode_mem);

  return 0;
}