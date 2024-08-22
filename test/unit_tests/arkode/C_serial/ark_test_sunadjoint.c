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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/sundials_core.h>

#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sunadjoint/sunadjoint_checkpointscheme_basic.h>
#include <sunadjoint/sunadjoint_stepper.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmemory/sunmemory_system.h>

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>

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

  printf("Forward Solution:\n");
  N_VPrint(u);

  printf("ARKODE Stats for Forward Solution:\n");
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  printf("\n");

  return 0;
}

int adjoint_solution(SUNContext sunctx, SUNAdjointStepper adj_stepper,
                     SUNAdjointCheckpointScheme checkpoint_scheme,
                     sunrealtype tf, sunrealtype tout, N_Vector sf)
{
  int retval      = 0;
  int stop_reason = 0;
  sunrealtype t   = tf;
  retval = SUNAdjointStepper_Evolve(adj_stepper, tout, sf, &t, &stop_reason);

  printf("Adjoint Solution:\n");
  N_VPrint(sf);

  printf("\nSUNAdjointStepper Stats:\n");
  SUNAdjointStepper_PrintAllStats(adj_stepper, stdout, SUN_OUTPUTFORMAT_TABLE);
  printf("\n");

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
  // Create the ARKODE stepper that will be used for the forward evolution.
  //

  const sunrealtype dt = args.dt;
  sunrealtype t0       = 0.0;
  sunrealtype tf       = args.tf;
  const int nsteps     = ((tf - t0) / dt + 1);
  const int order      = args.order;
  void* arkode_mem     = ARKStepCreate(lotka_volterra, NULL, t0, u, sunctx);

  ARKodeSetOrder(arkode_mem, order);
  ARKodeSetMaxNumSteps(arkode_mem, nsteps * 2);

  // Enable checkpointing during the forward solution.
  // ncheck will be more than nsteps, but for testing purposes we try setting it
  // to nsteps and allow things to be resized automatically.
  const int check_interval                     = args.check_freq;
  const int ncheck                             = nsteps;
  const sunbooleantype save_stages             = args.save_stages;
  const sunbooleantype keep_check              = args.keep_checks;
  SUNAdjointCheckpointScheme checkpoint_scheme = NULL;
  SUNMemoryHelper mem_helper                   = SUNMemoryHelper_Sys(sunctx);
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
  // Create the adjoint stepper
  //

  sunindextype num_params = 4;
  N_Vector sensu0         = N_VClone(u);
  N_Vector sensp          = N_VNew_Serial(num_params, sunctx);
  N_Vector sens[2]        = {sensu0, sensp};
  N_Vector sf             = N_VNew_ManyVector(2, sens, sunctx);

  // Set the terminal condition for the adjoint system, which
  // should be the the gradient of our cost function at tf.
  dgdu(u, sensu0, params, tf);
  dgdp(u, sensp, params, tf);

  printf("Adjoint terminal condition:\n");
  N_VPrint(sf);

  SUNAdjointStepper adj_stepper;
  ARKStepCreateAdjointSolver(arkode_mem, sf, &adj_stepper);

  //
  // Now compute the adjoint solution
  //

  SUNMatrix jac  = SUNDenseMatrix(neq, neq, sunctx);
  SUNMatrix jacp = SUNDenseMatrix(neq, num_params, sunctx);

  SUNAdjointStepper_SetJacFn(adj_stepper, jacobian, jac, parameter_jacobian,
                             jacp);
  adjoint_solution(sunctx, adj_stepper, checkpoint_scheme, tf, t0, sf);

  SUNMatDestroy(jac);
  SUNMatDestroy(jacp);

  //
  // Now compute the adjoint solution using Jvp
  //

  printf("\n-- Redo adjoint problem using JVP --\n\n");
  dgdu(u, sensu0, params, tf);
  dgdp(u, sensp, params, tf);
  SUNAdjointStepper_ReInit(adj_stepper, sf, tf);
  SUNAdjointStepper_SetJacFn(adj_stepper, NULL, NULL, NULL, NULL);
  SUNAdjointStepper_SetJacTimesVecFn(adj_stepper, jvp, parameter_jvp);
  adjoint_solution(sunctx, adj_stepper, checkpoint_scheme, tf, t0, sf);

  //
  // Now compute the adjoint solution using vJp
  //

  printf("\n-- Redo adjoint problem using VJP --\n\n");
  dgdu(u, sensu0, params, tf);
  dgdp(u, sensp, params, tf);
  SUNAdjointStepper_ReInit(adj_stepper, sf, tf);
  SUNAdjointStepper_SetJacTimesVecFn(adj_stepper, NULL, NULL);
  SUNAdjointStepper_SetVecTimesJacFn(adj_stepper, vjp, parameter_vjp);
  adjoint_solution(sunctx, adj_stepper, checkpoint_scheme, tf, t0, sf);

  //
  // Cleanup
  //

  N_VDestroy(u);
  N_VDestroy(sf);
  SUNAdjointCheckpointScheme_Destroy(&checkpoint_scheme);
  SUNAdjointStepper_Destroy(&adj_stepper);
  ARKodeFree(&arkode_mem);

  return 0;
}