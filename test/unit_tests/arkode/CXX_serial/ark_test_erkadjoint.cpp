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
#include <sunadjoint/sunadjoint_checkpointscheme_fixed.h>
#include <sunadjoint/sunadjoint_stepper.h>
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

static int check_forward_backward_answer(N_Vector answer)
{
  const sunrealtype u1 = SUN_RCONST(1.0);
  const sunrealtype u2 = SUN_RCONST(1.0);
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
    SUN_RCONST(-2.19271337646507),
  };

  const sunrealtype mu[4] = {SUN_RCONST(4.341147542533404),
                             SUN_RCONST(-2.000933816791803),
                             SUN_RCONST(1.010120676762905),
                             SUN_RCONST(-1.3955943267337996)};

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

static int check_sensitivities_backward(N_Vector answer)
{
  // The correct answer was generated with the Julia ForwardDiff.jl
  // automatic differentiation package.

  const sunrealtype lambda[2] = {
    SUN_RCONST(1.772850901841113),
    SUN_RCONST(-0.7412891218574361),
  };

  const sunrealtype mu[4] = {SUN_RCONST(0.0), SUN_RCONST(0.0), SUN_RCONST(0.0),
                             SUN_RCONST(0.0)};

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

  dg[0] = SUN_RCONST(-1.0) + u[0];
  dg[1] = SUN_RCONST(-1.0) + u[1];
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
  const int nsteps     = (int)ceil(((tf - t0) / dt + 1));
  const int order      = args.order;

  void* arkode_mem = ARKStepCreate(ode_rhs, NULL, t0, u, sunctx);
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
    fprintf(stderr,
            ">>> FAILURE: forward solution does not match correct answer\n");
    return -1;
  };
  printf(">>> PASS\n");

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
    fprintf(stderr,
            ">>> FAILURE: adjoint solution does not match correct answer\n");
    return -1;
  }
  printf("\n>>> PASS\n");

  //
  // Now compute the adjoint solution using Jvp
  //

  printf("\n-- Redo adjoint problem using JVP --\n\n");
  if (!keep_check)
  {
    N_VConst(SUN_RCONST(1.0), u);
    printf("Initial condition:\n");
    N_VPrint(u);
    ARKStepReInit(arkode_mem, ode_rhs, NULL, t0, u);
    forward_solution(sunctx, arkode_mem, checkpoint_scheme, t0, tf, dt, u);
    if (check_forward_answer(u))
    {
      fprintf(stderr,
              ">>> FAILURE: forward solution does not match correct answer\n");
      return -1;
    }
  }
  dgdu(u, sensu0, params, tf);
  dgdp(u, sensp, params, tf);
  SUNAdjointStepper_ReInit(adj_stepper, u, t0, sf, tf);
  SUNAdjointStepper_SetJacFn(adj_stepper, NULL, NULL, NULL, NULL);
  SUNAdjointStepper_SetJacTimesVecFn(adj_stepper, ode_jvp, parameter_jvp);
  adjoint_solution(sunctx, adj_stepper, checkpoint_scheme, tf, t0, sf);
  if (check_sensitivities(sf))
  {
    fprintf(stderr,
            ">>> FAILURE: adjoint solution does not match correct answer\n");
    return -1;
  };
  printf("\n>>> PASS\n");

  //
  // Now compute the adjoint solution using vJp
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
      fprintf(stderr,
              ">>> FAILURE: forward solution does not match correct answer\n");
      return -1;
    };
  }
  dgdu(u, sensu0, params, tf);
  dgdp(u, sensp, params, tf);
  SUNAdjointStepper_ReInit(adj_stepper, u, t0, sf, tf);
  SUNAdjointStepper_SetJacTimesVecFn(adj_stepper, NULL, NULL);
  SUNAdjointStepper_SetVecTimesJacFn(adj_stepper, ode_vjp, parameter_vjp);
  adjoint_solution(sunctx, adj_stepper, checkpoint_scheme, tf, t0, sf);
  if (check_sensitivities(sf))
  {
    fprintf(stderr,
            ">>> FAILURE: adjoint solution does not match correct answer\n");
    return -1;
  };
  printf(">>> PASS\n");

  //
  // Now compute the adjoint solution but for when forward problem done backwards
  // starting with the forward solution.
  //

  printf("\n-- Redo adjoint problem of forward problem done backwards  --\n\n");

  // Cleanup from the original forward problem and then recreate the integrator
  // for the forward problem done backwards.
  SUNAdjointCheckpointScheme_Destroy(&checkpoint_scheme);
  SUNAdjointStepper_Destroy(&adj_stepper);
  ARKodeFree(&arkode_mem);
  arkode_mem = ARKStepCreate(ode_rhs, NULL, tf, u, sunctx);
  ARKodeSetOrder(arkode_mem, order);
  ARKodeSetMaxNumSteps(arkode_mem, nsteps * 2);
  SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                          check_interval, ncheck, save_stages,
                                          keep_check, sunctx, &checkpoint_scheme);
  ARKodeSetAdjointCheckpointScheme(arkode_mem, checkpoint_scheme);

  printf("Initial condition:\n");
  N_VPrint(u);

  forward_solution(sunctx, arkode_mem, checkpoint_scheme, tf, t0, -dt, u);
  if (check_forward_backward_answer(u))
  {
    fprintf(stderr,
            ">>> FAILURE: forward solution does not match correct answer\n");
    return -1;
  };

  ARKStepCreateAdjointStepper(arkode_mem, sf, &adj_stepper);
  SUNAdjointStepper_SetJacFn(adj_stepper, ode_jac, jac, parameter_jacobian, jacp);
  dgdu(u, sensu0, params, t0);
  dgdp(u, sensp, params, t0);

  adjoint_solution(sunctx, adj_stepper, checkpoint_scheme, t0, tf, sf);
  if (check_sensitivities_backward(sf))
  {
    fprintf(stderr,
            ">>> FAILURE: adjoint solution does not match correct answer\n");
    return -1;
  };
  printf(">>> PASS\n");

  //
  // Cleanup
  //

  SUNMatDestroy(jac);
  SUNMatDestroy(jacp);
  N_VDestroy(u);
  N_VDestroy(sf);
  SUNAdjointCheckpointScheme_Destroy(&checkpoint_scheme);
  SUNAdjointStepper_Destroy(&adj_stepper);
  ARKodeFree(&arkode_mem);

  return 0;
}
