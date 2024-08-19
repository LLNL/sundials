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
#include <string.h>
#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sunadjoint/sunadjoint_checkpointscheme_basic.h>
#include <sunadjoint/sunadjoint_solver.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmemory/sunmemory_system.h>
#include "sundials/sundials_types.h"

typedef struct
{
  sunrealtype tf;
  sunrealtype dt;
  int order;
  int check_freq;
  sunbooleantype save_stages;
  sunbooleantype keep_checks;
} ProgramArgs;

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

  Jv[0] = (p[0] - p[1] * u[1]) * v[0] + p[3] * u[1] * v[1];
  Jv[1] = -p[1] * u[0] * v[0] + (-p[2] + p[3] * u[0]) * v[1];

  return 0;
}

int vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
        N_Vector udotvec, void* user_data, N_Vector tmp)
{
  sunrealtype* p  = (sunrealtype*)user_data;
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  Jv[0] = (p[0] - p[1] * u[1]) * v[0] + p[3] * u[1] * v[1];
  Jv[1] = -p[1] * u[0] * v[0] + (-p[2] + p[3] * u[0]) * v[1];

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

int parameter_jvp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
                  N_Vector udotvec, void* user_data, N_Vector tmp)
{
  if (user_data != params) { return -1; }

  sunrealtype* p  = (sunrealtype*)user_data;
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  Jv[0] = u[0] * v[0];
  Jv[1] = -u[0] * u[1] * v[0];
  Jv[2] = -u[1] * v[1];
  Jv[3] = u[0] * u[1] * v[1];

  return 0;
}

int parameter_vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
                  N_Vector udotvec, void* user_data, N_Vector tmp)
{
  if (user_data != params) { return -1; }

  sunrealtype* p  = (sunrealtype*)user_data;
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  Jv[0] = u[0] * v[0];
  Jv[1] = -u[0] * u[1] * v[0];
  Jv[2] = -u[1] * v[1];
  Jv[3] = u[0] * u[1] * v[1];

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

void dgdp(N_Vector uvec, N_Vector dgvec, const sunrealtype* p, sunrealtype t)
{
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* dg = N_VGetArrayPointer(dgvec);

  dg[0] = SUN_RCONST(0.0);
  dg[1] = SUN_RCONST(0.0);
  dg[2] = SUN_RCONST(0.0);
  dg[3] = SUN_RCONST(0.0);
}

int forward_solution(SUNContext sunctx, void* arkode_mem,
                     SUNAdjointCheckpointScheme checkpoint_scheme,
                     sunrealtype t0, sunrealtype tf, sunrealtype dt, N_Vector u)
{
  int retval = 0;

  retval = ARKodeSetUserData(arkode_mem, (void*)params);
  retval = ARKodeSetFixedStep(arkode_mem, dt);

  sunrealtype t = t0;
  while (t < tf)
  {
    retval = ARKodeEvolve(arkode_mem, tf, u, &t, ARK_NORMAL);
    if (retval < 0)
    {
      fprintf(stderr, ">>> ERROR: ARKodeEvolve returned %d\n", retval);
      return -1;
    }
  }

  fprintf(stdout, "Forward Solution:\n");
  N_VPrint(u);

  fprintf(stdout, "ARKODE Stats for Forward Solution:\n");
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  fprintf(stdout, "\n");

  return 0;
}

int adjoint_solution(SUNContext sunctx, void* arkode_mem,
                     SUNAdjointCheckpointScheme checkpoint_scheme,
                     sunrealtype tf, sunrealtype tout, N_Vector u)
{
  int retval = 0;

  sunindextype neq        = N_VGetLength(u);
  sunindextype num_params = 4;
  N_Vector sensu0         = N_VClone(u);
  N_Vector sensp          = N_VNew_Serial(num_params, sunctx);
  N_Vector sens[2]        = {sensu0, sensp};
  N_Vector sf             = N_VNew_ManyVector(2, sens, sunctx);

  // Set the terminal condition for the adjoint system, which
  // should be the the gradient of our cost function at tf.
  dgdu(u, sensu0, params, tf);
  dgdp(u, sensp, params, tf);

  fprintf(stdout, "Adjoint terminal condition:\n");
  N_VPrint(sf);

  SUNAdjointSolver adj_solver;
  retval = ARKStepCreateAdjointSolver(arkode_mem, sf, &adj_solver);

  SUNMatrix J  = SUNDenseMatrix(neq, neq, sunctx);
  SUNMatrix Jp = SUNDenseMatrix(neq, num_params, sunctx);

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

  fprintf(stdout, "\nSUNAdjointSolver Stats:\n");
  SUNAdjointSolver_PrintAllStats(adj_solver, stdout, SUN_OUTPUTFORMAT_TABLE);
  fprintf(stdout, "\n");

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
  sunindextype num_params = 4;
  N_Vector sensu0         = N_VClone(u);
  N_Vector sensp          = N_VNew_Serial(num_params, sunctx);
  N_Vector sens[2]        = {sensu0, sensp};
  N_Vector sf             = N_VNew_ManyVector(2, sens, sunctx);

  // Set the terminal condition for the adjoint system, which
  // should be the the gradient of our cost function at tf.
  dgdu(u, sensu0, params, tf);
  dgdp(u, sensp, params, tf);

  fprintf(stdout, "Adjoint terminal condition:\n");
  N_VPrint(sf);

  SUNAdjointSolver adj_solver;
  retval = ARKStepCreateAdjointSolver(arkode_mem, sf, &adj_solver);

  retval = SUNAdjointSolver_SetJacTimesVecFn(adj_solver, jvp, parameter_jvp);

  int stop_reason = 0;
  sunrealtype t   = tf;
  retval = SUNAdjointSolver_Solve(adj_solver, tout, sf, &t, &stop_reason);

  fprintf(stdout, "Adjoint Solution:\n");
  N_VPrint(sf);

  fprintf(stdout, "\nSUNAdjointSolver Stats:\n");
  SUNAdjointSolver_PrintAllStats(adj_solver, stdout, SUN_OUTPUTFORMAT_TABLE);
  fprintf(stdout, "\n");

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
  sunindextype num_params = 4;
  N_Vector sensu0         = N_VClone(u);
  N_Vector sensp          = N_VNew_Serial(num_params, sunctx);
  N_Vector sens[2]        = {sensu0, sensp};
  N_Vector sf             = N_VNew_ManyVector(2, sens, sunctx);

  // Set the terminal condition for the adjoint system, which
  // should be the the gradient of our cost function at tf.
  dgdu(u, sensu0, params, tf);
  dgdp(u, sensp, params, tf);

  SUNAdjointSolver adj_solver;
  retval = ARKStepCreateAdjointSolver(arkode_mem, sf, &adj_solver);

  retval = SUNAdjointSolver_SetVecTimesJacFn(adj_solver, vjp, parameter_vjp);

  int stop_reason = 0;
  sunrealtype t   = tf;
  retval = SUNAdjointSolver_Solve(adj_solver, tout, sf, &t, &stop_reason);

  fprintf(stdout, "Adjoint Solution:\n");
  N_VPrint(sf);

  fprintf(stdout, "\nSUNAdjointSolver Stats:\n");
  SUNAdjointSolver_PrintAllStats(adj_solver, stdout, SUN_OUTPUTFORMAT_TABLE);
  fprintf(stdout, "\n");

  N_VDestroy(sf);
  SUNAdjointSolver_Destroy(&adj_solver);

  return 0;
}

void print_help(int argc, char* argv[], int exit_code)
{
  if (exit_code)
  {
    fprintf(stderr, "./ark_test_sunadjoint: option not recognized\n");
  }
  else { fprintf(stderr, "./ark_test_sunadjoint "); }
  fprintf(stderr, "options:\n");
  fprintf(stderr, "--tf <real>         the final simulation time\n");
  fprintf(stderr, "--dt <real>         the timestep size\n");
  fprintf(stderr, "--order <int>       the order of the RK method\n");
  fprintf(stderr, "--check-freq <int>  how often to checkpoint (in steps)\n");
  fprintf(stderr, "--no-stages         don't checkpoint stages\n");
  fprintf(stderr,
          "--keep-checks       keep checkpoints around after loading\n");
  fprintf(stderr, "--help              print these options\n");
  exit(exit_code);
}

void parse_args(int argc, char* argv[], ProgramArgs* args)
{
  for (int argi = 1; argi < argc; ++argi)
  {
    const char* arg = argv[argi];
    if (!strcmp(arg, "--tf")) { args->tf = atof(argv[++argi]); }
    else if (!strcmp(arg, "--dt")) { args->dt = atof(argv[++argi]); }
    else if (!strcmp(arg, "--order")) { args->order = atoi(argv[++argi]); }
    else if (!strcmp(arg, "--check-freq"))
    {
      args->check_freq = atoi(argv[++argi]);
    }
    else if (!strcmp(arg, "--no-stages")) { args->save_stages = SUNFALSE; }
    else if (!strcmp(arg, "--keep-checks")) { args->keep_checks = SUNTRUE; }
    else if (!strcmp(arg, "--help")) { print_help(argc, argv, 0); }
    else { print_help(argc, argv, 1); }
  }
}

int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  // Since this a unit test, we want to abort immediately on any internal error
  SUNContext_PushErrHandler(sunctx, SUNAbortErrHandlerFn, NULL);

  ProgramArgs args;
  args.tf          = 10.0;
  args.dt          = 1e-2;
  args.order       = 4;
  args.save_stages = SUNTRUE;
  args.keep_checks = SUNTRUE;
  args.check_freq  = 2;
  parse_args(argc, argv, &args);

  //
  // Create the initial conditions vector
  //

  sunindextype neq = 2;
  N_Vector u       = N_VNew_Serial(neq, sunctx);
  N_VConst(1.0, u);

  //
  // Create the ARKODE stepper that will be used for both the forward solution and adjoint solution.
  //

  const sunrealtype dt = args.dt;
  sunrealtype t0       = 0.0;
  sunrealtype tf       = args.tf;
  const int nsteps     = ((tf - t0) / dt + 1);
  const int order      = args.order;
  void* arkode_mem     = ARKStepCreate(lotka_volterra, NULL, t0, u, sunctx);

  ARKodeSetOrder(arkode_mem, order);
  ARKodeSetMaxNumSteps(arkode_mem, nsteps * 2);

  // Enable checkpointing during the forward solution
  SUNAdjointCheckpointScheme checkpoint_scheme = NULL;

  SUNMemoryHelper mem_helper       = SUNMemoryHelper_Sys(sunctx);
  const int check_interval         = args.check_freq;
  const int ncheck                 = (nsteps * (order + 1));
  const sunbooleantype save_stages = args.save_stages;
  const sunbooleantype keep_check  = args.keep_checks;
  SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, mem_helper,
                                          check_interval, ncheck, save_stages,
                                          keep_check, sunctx, &checkpoint_scheme);

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

  //
  // Now compute the adjoint solution using Jvp
  //
  // TODO(CJB): make sure this reinitializes arkode correctly (probably need SUNAdjointSolver_Reset function)
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