/* -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Program to test the SUNAdjoint capability with ARKODE. The test uses the
 * implements the four parameter Lotka-Volterra problem
 *
 *    u' = [dx/dt] = [  p_0*x - p_1*x*y  ]
 *         [dy/dt]   [ -p_2*y + p_3*x*y ].
 *
 * The initial condition is u(t_0) = 1.0 and we use the parameters
 * p  = [1.5, 1.0, 3.0, 1.0]. We compute the sensitivities for the scalar cost
 * function,
 *
 *    g(u(t_f), p) = || 1 - u(t_f, p) ||^2 / 2
 *
 * with respect to the initial condition and the parameters.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/sundials_core.hpp>

#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h>
#include <sundials/sundials_adjointstepper.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmemory/sunmemory_system.h>

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>

#include "problems/lotka_volterra.hpp"

#if defined(SUNDIALS_SINGLE_PRECISION)
#define FWD_TOL SUN_RCONST(1e-2)
#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define FWD_TOL SUN_RCONST(1e-4)
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define FWD_TOL SUN_RCONST(1e-6)
#endif

#define ADJ_TOL SUN_RCONST(1e-2)

using namespace problems::lotka_volterra;

typedef struct
{
  sunrealtype tf;
  sunrealtype dt;
  int order;
  int check_freq;
  sunbooleantype save_stages;
  sunbooleantype keep_checks;
} ProgramArgs;

static sunrealtype params[4] = {SUN_RCONST(1.5), SUN_RCONST(1.0),
                                SUN_RCONST(3.0), SUN_RCONST(1.0)};

static int neg_rhs(sunrealtype t, N_Vector uvec, N_Vector udotvec, void* user_data)
{
  int err = ode_rhs(t, uvec, udotvec, user_data);
  if (err) { return err; }
  else { N_VScale(SUN_RCONST(-1.0), udotvec, udotvec); }
  return 0;
}

static int neg_vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
                   N_Vector udotvec, void* user_data, N_Vector tmp)
{
  int err = ode_vjp(vvec, Jvvec, t, uvec, udotvec, user_data, tmp);
  if (err) { return err; }
  else { N_VScale(SUN_RCONST(-1.0), Jvvec, Jvvec); }
  return 0;
}

static int neg_parameter_vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t,
                             N_Vector uvec, N_Vector udotvec, void* user_data,
                             N_Vector tmp)
{
  int err = parameter_vjp(vvec, Jvvec, t, uvec, udotvec, user_data, tmp);
  if (err) { return err; }
  else { N_VScale(SUN_RCONST(-1.0), Jvvec, Jvvec); }
  return 0;
}

static int check_forward_answer(N_Vector answer)
{
  const sunrealtype u1 = SUN_RCONST(2.77266836);
  const sunrealtype u2 = SUN_RCONST(0.258714765);
  sunrealtype* ans     = N_VGetArrayPointer(answer);

  if (SUNRCompareTol(ans[0], u1, FWD_TOL))
  {
    fprintf(stdout, "\n>>> ans[0] = %g, should be %g\n", ans[0], u1);
    return -1;
  };
  if (SUNRCompareTol(ans[1], u2, FWD_TOL))
  {
    fprintf(stdout, "\n>>> ans[1] = %g, should be %g\n", ans[1], u2);
    return -1;
  };

  return 0;
}

static int check_sensitivities(N_Vector answer)
{
  // The correct answer was generated with the Julia ForwardDiff.jl
  // automatic differentiation package.

  const sunrealtype lambda[2] = {
    SUN_RCONST(3.5202568952661544),
    -SUN_RCONST(2.19271337646507),
  };

  const sunrealtype mu[4] = {SUN_RCONST(4.341147542533404),
                             -SUN_RCONST(2.000933816791803),
                             SUN_RCONST(1.010120676762905),
                             -SUN_RCONST(1.3955943267337996)};

  sunrealtype* ans = N_VGetSubvectorArrayPointer_ManyVector(answer, 0);

  for (sunindextype i = 0; i < 2; ++i)
  {
    if (SUNRCompareTol(ans[i], lambda[i], ADJ_TOL))
    {
      fprintf(stdout, "\n>>> ans[%lld] = %g, should be %g\n", (long long)i,
              ans[i], lambda[i]);
      return -1;
    };
  }

  ans = N_VGetSubvectorArrayPointer_ManyVector(answer, 1);

  for (sunindextype i = 0; i < 4; ++i)
  {
    if (SUNRCompareTol(ans[i], mu[i], ADJ_TOL))
    {
      fprintf(stdout, "\n>>> ans[%lld] = %g, should be %g\n", (long long)i,
              ans[i], mu[i]);
      return -1;
    };
  }

  return 0;
}

static void dgdu(N_Vector uvec, N_Vector dgvec, const sunrealtype* p,
                 sunrealtype t)
{
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* dg = N_VGetArrayPointer(dgvec);

  dg[0] = -SUN_RCONST(1.0) + u[0];
  dg[1] = -SUN_RCONST(1.0) + u[1];
}

static void dgdp(N_Vector uvec, N_Vector dgvec, const sunrealtype* p,
                 sunrealtype t)
{
  sunrealtype* dg = N_VGetArrayPointer(dgvec);

  dg[0] = SUN_RCONST(0.0);
  dg[1] = SUN_RCONST(0.0);
  dg[2] = SUN_RCONST(0.0);
  dg[3] = SUN_RCONST(0.0);
}

static int forward_solution(SUNContext sunctx, void* arkode_mem,
                            SUNAdjointCheckpointScheme checkpoint_scheme,
                            const sunrealtype t0, const sunrealtype tf,
                            const sunrealtype dt, N_Vector u)
{
  int retval = 0;

  retval = ARKodeSetUserData(arkode_mem, (void*)params);
  retval = ARKodeSetFixedStep(arkode_mem, dt);

  sunrealtype t = t0;
  retval        = ARKodeEvolve(arkode_mem, tf, u, &t, ARK_NORMAL);
  if (retval < 0)
  {
    fprintf(stderr, ">>> ERROR: ARKodeEvolve returned %d\n", retval);
    return -1;
  }

  printf("Forward Solution:\n");
  N_VPrint(u);

  printf("ARKODE Stats for Forward Solution:\n");
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  printf("\n");

  return 0;
}

static int adjoint_solution(SUNContext sunctx, SUNAdjointStepper adj_stepper,
                            SUNAdjointCheckpointScheme checkpoint_scheme,
                            const sunrealtype tf, const sunrealtype tout,
                            N_Vector sf)
{
  sunrealtype t = tf;
  SUNAdjointStepper_Evolve(adj_stepper, tout, sf, &t);

  printf("Adjoint Solution:\n");
  N_VPrint(sf);

  printf("\nSUNAdjointStepper Stats:\n");
  SUNAdjointStepper_PrintAllStats(adj_stepper, stdout, SUN_OUTPUTFORMAT_TABLE);
  printf("\n");

  return 0;
}

static void print_help(int argc, char* argv[], int exit_code)
{
  if (exit_code) { fprintf(stderr, "%s: option not recognized\n", argv[0]); }
  else { fprintf(stderr, "%s ", argv[0]); }
  fprintf(stderr, "options:\n");
  fprintf(stderr, "--tf <real>         the final simulation time\n");
  fprintf(stderr, "--dt <real>         the timestep size\n");
  fprintf(stderr, "--order <int>       the order of the RK method\n");
  fprintf(stderr, "--check-freq <int>  how often to checkpoint (in steps)\n");
  fprintf(stderr, "--no-stages         don't checkpoint stages\n");
  fprintf(stderr,
          "--dont-keep         don't keep checkpoints around after loading\n");
  fprintf(stderr, "--help              print these options\n");
  exit(exit_code);
}

static void parse_args(int argc, char* argv[], ProgramArgs* args)
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
    else if (!strcmp(arg, "--dont-keep")) { args->keep_checks = SUNFALSE; }
    else if (!strcmp(arg, "--help")) { print_help(argc, argv, 0); }
    else { print_help(argc, argv, 1); }
  }
}

int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;
  int failcount     = 0;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  // Since this a unit test, we want to abort immediately on any internal error
  SUNContext_PushErrHandler(sunctx, SUNAbortErrHandlerFn, NULL);

  ProgramArgs args;
  args.tf          = SUN_RCONST(1.0);
  args.dt          = SUN_RCONST(1e-4);
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
  N_VConst(SUN_RCONST(1.0), u);

  //
  // Create the ARKODE stepper that will be used for the forward evolution.
  //

  const sunrealtype dt = args.dt;
  sunrealtype t0       = SUN_RCONST(0.0);
  sunrealtype tf       = args.tf;
  const int nsteps     = (int)ceil(((tf - t0) / dt));
  const int order      = args.order;

  void* arkode_mem = ARKStepCreate(ode_rhs, NULL, t0, u, sunctx);
  ARKodeSetOrder(arkode_mem, order);
  // Due to roundoff in the `t` accumulation within the integrator,
  // the integrator may actually use nsteps + 1 time steps to reach tf.
  ARKodeSetMaxNumSteps(arkode_mem, nsteps + 1);

  // Enable checkpointing during the forward solution.
  // ncheck will be more than nsteps, but for testing purposes we try setting it
  // to nsteps and allow things to be resized automatically.
  const int check_interval                     = args.check_freq;
  const int ncheck                             = nsteps;
  const sunbooleantype save_stages             = args.save_stages;
  const sunbooleantype keep_check              = args.keep_checks;
  SUNAdjointCheckpointScheme checkpoint_scheme = NULL;
  SUNMemoryHelper mem_helper                   = SUNMemoryHelper_Sys(sunctx);
  SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                          check_interval, ncheck, save_stages,
                                          keep_check, sunctx, &checkpoint_scheme);
  ARKodeSetAdjointCheckpointScheme(arkode_mem, checkpoint_scheme);

  //
  // Compute the forward solution
  //

  printf("\n-- Do forward problem --\n\n");

  printf("Initial condition:\n");
  N_VPrint(u);

  forward_solution(sunctx, arkode_mem, checkpoint_scheme, t0, tf, dt, u);
  if (check_forward_answer(u))
  {
    ++failcount;
    fprintf(stderr,
            ">>> FAILURE: forward solution does not match correct answer\n");
  }
  else { printf(">>> PASS\n"); }

  //
  // Create the adjoint stepper
  //

  printf("\n-- Do adjoint problem using Jacobian matrix --\n\n");

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
  ARKStepCreateAdjointStepper(arkode_mem, sf, &adj_stepper);

  //
  // Now compute the adjoint solution
  //

  SUNMatrix jac  = SUNDenseMatrix(neq, neq, sunctx);
  SUNMatrix jacp = SUNDenseMatrix(neq, num_params, sunctx);

  SUNAdjointStepper_SetJacFn(adj_stepper, ode_jac, jac, parameter_jacobian, jacp);

  adjoint_solution(sunctx, adj_stepper, checkpoint_scheme, tf, t0, sf);
  if (check_sensitivities(sf))
  {
    ++failcount;
    fprintf(stderr,
            ">>> FAILURE: adjoint solution does not match correct answer\n");
  }
  else { printf("\n>>> PASS\n"); }

  //
  // Now compute the adjoint solution using vjp
  //

  printf("\n-- Redo adjoint problem using VJP --\n\n");
  if (!keep_check)
  {
    N_VConst(SUN_RCONST(1.0), u);
    printf("Initial condition:\n");
    N_VPrint(u);
    ARKStepReInit(arkode_mem, ode_rhs, NULL, t0, u);
    forward_solution(sunctx, arkode_mem, checkpoint_scheme, t0, tf, dt, u);
    if (check_forward_answer(u))
    {
      ++failcount;
      fprintf(stderr,
              ">>> FAILURE: forward solution does not match correct answer\n");
    };
  }
  dgdu(u, sensu0, params, tf);
  dgdp(u, sensp, params, tf);
  SUNAdjointStepper_ReInit(adj_stepper, u, t0, sf, tf);
  SUNAdjointStepper_SetJacFn(adj_stepper, NULL, NULL, NULL, NULL);
  SUNAdjointStepper_SetVecHermitianTransposeJacFn(adj_stepper, ode_vjp,
                                                  parameter_vjp);
  adjoint_solution(sunctx, adj_stepper, checkpoint_scheme, tf, t0, sf);
  if (check_sensitivities(sf))
  {
    ++failcount;
    fprintf(stderr,
            ">>> FAILURE: adjoint solution does not match correct answer\n");
  }
  else { printf(">>> PASS\n"); }

  //
  // Now compute the adjoint solution but using tau = -t so we can test
  // the backwards-in-time code paths.
  //

  printf("\n-- Redo adjoint problem with change of variables tau = -t  --\n\n");

  // Swap the start and end times
  sunrealtype tau0 = tf;
  sunrealtype tauf = t0;

  // Cleanup from the original forward problem
  SUNAdjointCheckpointScheme_Destroy(&checkpoint_scheme);
  SUNAdjointStepper_Destroy(&adj_stepper);
  ARKodeFree(&arkode_mem);

  // Reset the initial condition
  N_VConst(SUN_RCONST(1.0), u);
  printf("Initial condition:\n");
  N_VPrint(u);

  // Recreate the integrator
  arkode_mem = ARKStepCreate(neg_rhs, NULL, tau0, u, sunctx);
  ARKodeSetOrder(arkode_mem, order);
  ARKodeSetMaxNumSteps(arkode_mem, nsteps + 1);
  SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                          check_interval, ncheck, save_stages,
                                          keep_check, sunctx, &checkpoint_scheme);
  ARKodeSetAdjointCheckpointScheme(arkode_mem, checkpoint_scheme);

  forward_solution(sunctx, arkode_mem, checkpoint_scheme, tau0, tauf, -dt, u);
  if (check_forward_answer(u))
  {
    ++failcount;
    fprintf(stderr,
            ">>> FAILURE: forward solution does not match correct answer\n");
  };

  dgdu(u, sensu0, params, tauf);
  dgdp(u, sensp, params, tauf);
  printf("Adjoint terminal condition:\n");
  N_VPrint(sf);

  ARKStepCreateAdjointStepper(arkode_mem, sf, &adj_stepper);
  SUNAdjointStepper_SetVecHermitianTransposeJacFn(adj_stepper, neg_vjp,
                                                  neg_parameter_vjp);

  adjoint_solution(sunctx, adj_stepper, checkpoint_scheme, tauf, tau0, sf);
  if (check_sensitivities(sf))
  {
    ++failcount;
    fprintf(stderr,
            ">>> FAILURE: adjoint solution does not match correct answer\n");
  }
  else { printf(">>> PASS\n"); }

  //
  // Cleanup
  //

  // adjoint related
  SUNMatDestroy(jac);
  SUNMatDestroy(jacp);
  N_VDestroy(sensu0);
  N_VDestroy(sensp);
  N_VDestroy(sf);
  SUNAdjointCheckpointScheme_Destroy(&checkpoint_scheme);
  SUNAdjointStepper_Destroy(&adj_stepper);
  SUNMemoryHelper_Destroy(mem_helper);
  // forward and adjoint related
  N_VDestroy(u);
  ARKodeFree(&arkode_mem);
  SUNContext_Free(&sunctx);

  if (failcount) { fprintf(stderr, ">>> %d TESTS FAILED\n", failcount); }

  return failcount;
}
