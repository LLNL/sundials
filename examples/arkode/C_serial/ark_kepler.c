/* clang-format off */
/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * We consider the Kepler problem. We choose one body to be the center of our
 * coordinate system and then we use the coordinates q = (q1, q2) to represent
 * the position of the second body relative to the first (center). This yields
 * the ODE:
 *    dq/dt = [ p1 ]
 *            [ p2 ]
 *    dp/dt = [ -q1 / (q1^2 + q2^2)^(3/2) ]
 *          = [ -q2 / (q1^2 + q2^2)^(3/2) ]
 * with the initial conditions
 *    q(0) = [ 1 - e ],  p(0) = [        0          ]
 *           [   0   ]          [ sqrt((1+e)/(1-e)) ]
 * where e = 0.6 is the eccentricity.
 *
 * The Hamiltonian for the system,
 *    H(p,q) = 1/2 * (p1^2 + p2^2) - 1/sqrt(q1^2 + q2^2)
 * is conserved as well as the angular momentum,
 *    L(p,q) = q1*p2 - q2*p1.
 *
 * By default we solve the problem by letting y = [ q, p ]^T then using a 4th
 * order symplectic integrator via the SPRKStep time-stepper of ARKODE with a
 * fixed time-step size.
 *
 * The rootfinding feature of SPRKStep is used to count the number of complete orbits. 
 * This is done by defining the function,
 *    g(q) = q2
 * and providing it to SPRKStep as the function to find the roots for.
 *
 * The program also accepts command line arguments to change the method
 * used and time-stepping strategy. The program has the following CLI arguments:
 * 
 *   --step-mode <fixed, adapt>  should we use a fixed time-step or adaptive time-step (default fixed)
 *   --stepper <SPRK, ERK>       should we use SPRKStep or ARKStep with an ERK method (default SPRK)
 *   --method <string>           which method to use (default ARKODE_SYMPLECTIC_MCLACHLAN_4_4)
 *   --use-compensated-sums      turns on compensated summation in ARKODE where applicable
 *   --dt <Real>                 the fixed-time step size to use if fixed time stepping is turned on (default 0.01)
 *   --tf <Real>                 the final time for the simulation (default 100)
 *   --nout                      number of output times
 *   --count-orbits              use rootfinding to count the number of completed orbits
 *   --check-order               compute the order of the method used and check if it is within range of the expected 
 * 
 * References:
 *    Ernst Hairer, Christain Lubich, Gerhard Wanner
 *    Geometric Numerical Integration: Structure-Preserving
 *    Algorithms for Ordinary Differential Equations
 *    Springer, 2006,
 *    ISSN 0179-3632
 * --------------------------------------------------------------------------*/
/* clang-format on */

#include <arkode/arkode_arkstep.h> /* prototypes for ARKStep fcts., consts */
#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h> /* prototypes for MRIStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector type, fcts., macros  */
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_math.h> /* def. math fcns, 'sunrealtype'           */
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#include "arkode/arkode.h"
#include "sundials/sundials_context.h"

#define NUM_DT 8

typedef struct
{
  sunrealtype ecc;
} * UserData;

typedef struct
{
  N_Vector sol;
  sunrealtype energy_error;
  int method_order;
} ProblemResult;

typedef struct
{
  int step_mode;
  int stepper;
  int num_output_times;
  int use_compsums;
  int count_orbits;
  int check_order;
  sunrealtype dt;
  sunrealtype tf;
  const char* method_name;
} ProgramArgs;

/* Helper functions */
static int ParseArgs(int argc, char* argv[], ProgramArgs* args);
static void PrintArgs(ProgramArgs* args);
static void PrintHelp();
static int ComputeConvergence(int num_dt, sunrealtype* orders,
                              sunrealtype expected_order, sunrealtype a11,
                              sunrealtype a12, sunrealtype a21, sunrealtype a22,
                              sunrealtype b1, sunrealtype b2,
                              sunrealtype* ord_avg, sunrealtype* ord_max,
                              sunrealtype* ord_est);
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* These functions do the interesting work that offer a good example of how to
 * use SPRKStep */
static int SolveProblem(ProgramArgs* args, ProblemResult* result,
                        SUNContext sunctx);
static void InitialConditions(N_Vector y0, sunrealtype ecc);
static sunrealtype Hamiltonian(N_Vector yvec);
static sunrealtype AngularMomentum(N_Vector y);
static int dydt(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int velocity(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int force(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int rootfn(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data);

int main(int argc, char* argv[])
{
  ProgramArgs args;
  ProblemResult result;
  SUNContext sunctx = NULL;
  int retval        = 0;

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /* Parse the command line arguments */
  if (ParseArgs(argc, argv, &args)) { return 1; };

  /* Allocate space for result variables */
  result.sol = N_VNew_Serial(4, sunctx);

  if (!args.check_order)
  {
    /* SolveProblem calls a stepper to evolve the problem to Tf */
    retval = SolveProblem(&args, &result, sunctx);
    if (check_retval(&retval, "SolveProblem", 1)) { return 1; }
  }
  else
  {
    /* Compute the order of accuracy of the method by testing
       it with different step sizes. */
    sunrealtype acc_orders[NUM_DT];
    sunrealtype con_orders[NUM_DT];
    sunrealtype acc_errors[NUM_DT];
    sunrealtype con_errors[NUM_DT];
    int expected_order = ARKodeSPRKStorage_LoadByName(args.method_name)->q;
    N_Vector ref_sol   = N_VClone(result.sol);
    N_Vector error     = N_VClone(result.sol);
    sunrealtype a11 = 0, a12 = 0, a21 = 0, a22 = 0;
    sunrealtype b1 = 0, b2 = 0, b1e = 0, b2e = 0;
    sunrealtype ord_max_acc = 0, ord_max_conv = 0, ord_avg = 0, ord_est = 0;
    sunrealtype refine = SUN_RCONST(.5);
    sunrealtype dt = (expected_order >= 3) ? SUN_RCONST(1e-1) : SUN_RCONST(1e-3);
    sunrealtype dts[NUM_DT] = {dt,
                               dt * refine,
                               dt * refine * refine,
                               dt * pow(refine, 3),
                               dt * pow(refine, 4),
                               dt * pow(refine, 5),
                               dt * pow(refine, 6),
                               dt * pow(refine, 7)};

    /* Create a reference solution using 8th order ERK with a small time step */
    const int old_step_mode     = args.step_mode;
    const int old_stepper       = args.stepper;
    const char* old_method_name = args.method_name;
    args.dt                     = SUN_RCONST(1e-2);
    args.step_mode              = 0;
    args.stepper                = 1;
    args.method_name            = "ARKODE_ARK548L2SAb_ERK_8_4_5";

    /* SolveProblem calls a stepper to evolve the problem to Tf */
    retval = SolveProblem(&args, &result, sunctx);
    if (check_retval(&retval, "SolveProblem", 1)) { return 1; }

    /* Store the reference solution */
    N_VScale(SUN_RCONST(1.0), result.sol, ref_sol);

    /* Restore the program args */
    args.step_mode   = old_step_mode;
    args.stepper     = old_stepper;
    args.method_name = old_method_name;

    /* Compute the error with various step sizes */
    for (int i = 0; i < NUM_DT; i++)
    {
      /* Set the dt to use for this solve */
      args.dt = dts[i];

      /* SolveProblem calls a stepper to evolve the problem to Tf */
      retval = SolveProblem(&args, &result, sunctx);
      if (check_retval(&retval, "SolveProblem", 1)) { return 1; }

      printf("\n");

      /* Compute the error */
      N_VLinearSum(SUN_RCONST(1.0), result.sol, -SUN_RCONST(1.0), ref_sol, error);
      acc_errors[i] = SUNRsqrt(N_VDotProd(error, error)) /
                      ((sunrealtype)N_VGetLength(error));
      con_errors[i] = SUNRabs(result.energy_error);

      a11 += 1;
      a12 += log(dts[i]);
      a21 += log(dts[i]);
      a22 += (log(dts[i]) * log(dts[i]));
      b1 += log(acc_errors[i]);
      b2 += (log(acc_errors[i]) * log(dts[i]));
      b1e += log(con_errors[i]);
      b2e += (log(con_errors[i]) * log(dts[i]));

      if (i >= 1)
      {
        acc_orders[i - 1] = log(acc_errors[i] / acc_errors[i - 1]) /
                            log(dts[i] / dts[i - 1]);
        con_orders[i - 1] = log(con_errors[i] / con_errors[i - 1]) /
                            log(dts[i] / dts[i - 1]);
      }
    }

    /* Compute the order of accuracy */
    retval = ComputeConvergence(NUM_DT, acc_orders, expected_order, a11, a12, a21,
                                a22, b1, b2, &ord_avg, &ord_max_acc, &ord_est);
    printf("Order of accuracy wrt solution:    expected = %d, max = %.4Lf,  "
           "avg "
           "= %.4Lf,  "
           "overall = %.4Lf\n",
           expected_order, (long double)ord_max_acc, (long double)ord_avg,
           (long double)ord_est);

    /* Compute the order of accuracy with respect to conservation */
    retval = ComputeConvergence(NUM_DT, con_orders, expected_order, a11, a12,
                                a21, a22, b1e, b2e, &ord_avg, &ord_max_conv,
                                &ord_est);

    printf("Order of accuracy wrt Hamiltonian: expected = %d, max = %.4Lf,  "
           "avg = %.4Lf,  overall = %.4Lf\n",
           expected_order, (long double)ord_max_conv, (long double)ord_avg,
           (long double)ord_est);

    if (ord_max_acc < (expected_order - RCONST(0.5)))
    {
      printf(">>> FAILURE: computed order of accuracy wrt solution is below "
             "expected (%d)\n",
             expected_order);
      return 1;
    }

    if (ord_max_conv < (expected_order - RCONST(0.5)))
    {
      printf(">>> FAILURE: computed order of accuracy wrt Hamiltonian is below "
             "expected (%d)\n",
             expected_order);
      return 1;
    }
  }

  N_VDestroy(result.sol);
  SUNContext_Free(&sunctx);

  return 0;
}

int SolveProblem(ProgramArgs* args, ProblemResult* result, SUNContext sunctx)
{
  void* arkode_mem       = NULL;
  N_Vector y             = NULL;
  SUNNonlinearSolver NLS = NULL;
  UserData udata         = NULL;
  sunrealtype* ydata     = NULL;
  sunrealtype tout       = NAN;
  sunrealtype tret       = NAN;
  sunrealtype H0         = NAN;
  sunrealtype L0         = NAN;
  sunrealtype num_orbits = 0;
  FILE* conserved_fp     = NULL;
  FILE* solution_fp      = NULL;
  FILE* times_fp         = NULL;
  int rootsfound         = 0;
  int iout               = 0;
  int retval             = 0;

  const int count_orbits     = args->count_orbits;
  const int step_mode        = args->step_mode;
  const int stepper          = args->stepper;
  const int use_compsums     = args->use_compsums;
  const int num_output_times = args->num_output_times;
  const char* method_name    = args->method_name;
  const sunrealtype dt       = args->dt;
  sunrealtype Tf             = args->tf;

  /* Default problem parameters */
  const sunrealtype T0    = SUN_RCONST(0.0);
  const sunrealtype dTout = (Tf - T0) / ((sunrealtype)num_output_times);
  const sunrealtype ecc   = SUN_RCONST(0.6);

  printf("\n   Begin Kepler Problem\n\n");
  PrintArgs(args);

  /* Allocate and fill udata structure */
  udata      = (UserData)malloc(sizeof(*udata));
  udata->ecc = ecc;

  /* Allocate our state vector */
  y     = N_VNew_Serial(4, sunctx);
  ydata = N_VGetArrayPointer(y);

  /* Fill the initial conditions */
  InitialConditions(y, ecc);

  /* Create SPRKStep integrator */
  if (stepper == 0)
  {
    arkode_mem = SPRKStepCreate(force, velocity, T0, y, sunctx);

    /* Optional: enable temporal root-finding */
    if (count_orbits)
    {
      SPRKStepRootInit(arkode_mem, 1, rootfn);
      if (check_retval(&retval, "SPRKStepRootInit", 1)) return 1;
    }

    retval = SPRKStepSetMethodName(arkode_mem, method_name);
    if (check_retval(&retval, "SPRKStepSetMethodName", 1)) return 1;

    retval = SPRKStepSetUseCompensatedSums(arkode_mem, use_compsums);
    if (check_retval(&retval, "SPRKStepSetUseCompensatedSums", 1)) return 1;

    if (step_mode == 0)
    {
      retval = SPRKStepSetFixedStep(arkode_mem, dt);
      if (check_retval(&retval, "SPRKStepSetFixedStep", 1)) return 1;

      retval = SPRKStepSetMaxNumSteps(arkode_mem, ((long int)ceil(Tf / dt)) + 1);
      if (check_retval(&retval, "SPRKStepSetMaxNumSteps", 1)) return 1;
    }
    else
    {
      fprintf(stderr,
              "ERROR: adaptive time-steps are not supported with SPRKStep\n");
      return 1;
    }

    retval = SPRKStepSetUserData(arkode_mem, (void*)udata);
    if (check_retval(&retval, "SPRKStepSetUserData", 1)) return 1;
  }
  else if (stepper == 1)
  {
    arkode_mem = ARKStepCreate(dydt, NULL, T0, y, sunctx);

    retval = ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE", method_name);
    if (check_retval(&retval, "ARKStepSetTableName", 1)) return 1;

    if (count_orbits)
    {
      ARKStepRootInit(arkode_mem, 1, rootfn);
      if (check_retval(&retval, "ARKStepRootInit", 1)) return 1;
    }

    retval = ARKStepSetUserData(arkode_mem, (void*)udata);
    if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;

    retval = ARKStepSetMaxNumSteps(arkode_mem, ((long int)ceil(Tf / dt)) + 1);
    if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return 1;

    if (step_mode == 0) { retval = ARKStepSetFixedStep(arkode_mem, dt); }
    else
    {
      retval = ARKStepSStolerances(arkode_mem, dt, dt);
      if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
    }
  }

  /* Open output files */
  if (stepper == 0)
  {
    const char* fmt1 = "ark_kepler_conserved_%s-dt-%.2e.txt";
    const char* fmt2 = "ark_kepler_solution_%s-dt-%.2e.txt";
    const char* fmt3 = "ark_kepler_times_%s-dt-%.2e.txt";
    char fname[256];
    sprintf(fname, fmt1, method_name, dt);
    conserved_fp = fopen(fname, "w+");
    sprintf(fname, fmt2, method_name, dt);
    solution_fp = fopen(fname, "w+");
    sprintf(fname, fmt3, method_name, dt);
    times_fp = fopen(fname, "w+");
  }
  else
  {
    const char* fmt1 = "ark_kepler_conserved_%s-dt-%.2e.txt";
    const char* fmt2 = "ark_kepler_solution_%s-dt-%.2e.txt";
    const char* fmt3 = "ark_kepler_times_%s-dt-%.2e.txt";
    char fname[256];
    sprintf(fname, fmt1, method_name, dt);
    conserved_fp = fopen(fname, "w+");
    sprintf(fname, fmt2, method_name, dt);
    solution_fp = fopen(fname, "w+");
    sprintf(fname, fmt3, method_name, dt);
    times_fp = fopen(fname, "w+");
  }

  /* Print out starting energy, momentum before integrating */
  tret = T0;
  tout = T0 + dTout;
  H0   = Hamiltonian(y);
  L0   = AngularMomentum(y);
  fprintf(stdout, "t = %.4Lf, H(p,q) = %.16Lf, L(p,q) = %.16Lf\n", tret, H0, L0);
  fprintf(times_fp, "%.16Lf\n", tret);
  fprintf(conserved_fp, "%.16Lf, %.16Lf\n", H0, L0);
  N_VPrintFile(y, solution_fp);

  /* Do integration */
  if (stepper == 0)
  {
    while (iout < num_output_times)
    {
      /* Optional: if the stop time is not set, then its possible that the
         exact requested output time will not be hit (even with a fixed
         time-step due to roundoff error accumulation) and interpolation will be
         used to get the solution at the output time. */
      SPRKStepSetStopTime(arkode_mem, tout);
      retval = SPRKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

      if (retval == ARK_ROOT_RETURN)
      {
        num_orbits += SUN_RCONST(0.5);

        fprintf(stdout, "ROOT RETURN:\t");
        SPRKStepGetRootInfo(arkode_mem, &rootsfound);
        fprintf(stdout, "  g[0] = %3d, y[0] = %3Lg, y[1] = %3Lg, num. orbits is now %.2Lf\n",
                rootsfound, (long double)ydata[0], (long double)ydata[1],
                (long double)num_orbits);
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
                (long double)tret, (long double)(Hamiltonian(y) - H0),
                (long double)(AngularMomentum(y) - L0));
      }
      else if (retval >= 0)
      {
        /* Output current integration status */
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
                (long double)tret, (long double)(Hamiltonian(y) - H0),
                (long double)(AngularMomentum(y) - L0));
        fprintf(times_fp, "%.16Lf\n", (long double)tret);
        fprintf(conserved_fp, "%.16Lf, %.16Lf\n", (long double)Hamiltonian(y),
                (long double)AngularMomentum(y));

        N_VPrintFile(y, solution_fp);

        tout += dTout;
        tout = (tout > Tf) ? Tf : tout;
        iout++;
      }
      else
      {
        fprintf(stderr, "Solver failure, stopping integration\n");
        break;
      }
    }
  }
  else
  {
    while (iout < num_output_times)
    {
      /* Optional: if the stop time is not set, then its possible that the the
         exact requested output time will not be hit (even with a fixed
         time-step due to roundoff error accumulation) and interpolation will be
         used to get the solution at the output time. */
      ARKStepSetStopTime(arkode_mem, tout);
      retval = ARKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

      if (retval == ARK_ROOT_RETURN)
      {
        num_orbits += SUN_RCONST(0.5);

        fprintf(stdout, "ROOT RETURN:\t");
        ARKStepGetRootInfo(arkode_mem, &rootsfound);
        fprintf(stdout, "  g[0] = %3d, y[0] = %3Lg, y[1] = %3Lg, num. orbits is now %.2Lf\n",
                rootsfound, (long double)ydata[0], (long double)ydata[1],
                (long double)num_orbits);
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
                (long double)tret, (long double)(Hamiltonian(y) - H0),
                (long double)(AngularMomentum(y) - L0));
      }
      else if (retval >= 0)
      {
        /* Output current integration status */
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
                (long double)tret, (long double)(Hamiltonian(y) - H0),
                (long double)(AngularMomentum(y) - L0));
        fprintf(times_fp, "%.16Lf\n", (long double)tret);
        fprintf(conserved_fp, "%.16Lf, %.16Lf\n", (long double)Hamiltonian(y),
                (long double)AngularMomentum(y));

        N_VPrintFile(y, solution_fp);

        tout += dTout;
        tout = (tout > Tf) ? Tf : tout;
        iout++;
      }
      else
      {
        fprintf(stderr, "Solver failure, stopping integration\n");
        break;
      }
    }
  }

  /* Copy results */
  N_VScale(SUN_RCONST(1.0), y, result->sol);
  result->energy_error = Hamiltonian(y) - H0;

  free(udata);
  fclose(times_fp);
  fclose(conserved_fp);
  fclose(solution_fp);
  if (NLS) { SUNNonlinSolFree(NLS); }
  N_VDestroy(y);
  if (stepper == 0)
  {
    SPRKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    SPRKStepFree(&arkode_mem);
  }
  else
  {
    ARKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    ARKStepFree(&arkode_mem);
  }

  return 0;
}

void InitialConditions(N_Vector y0vec, sunrealtype ecc)
{
  const sunrealtype zero = SUN_RCONST(0.0);
  const sunrealtype one  = SUN_RCONST(1.0);
  sunrealtype* y0        = N_VGetArrayPointer(y0vec);

  y0[0] = one - ecc;
  y0[1] = zero;
  y0[2] = zero;
  y0[3] = SUNRsqrt((one + ecc) / (one - ecc));
}

sunrealtype Hamiltonian(N_Vector yvec)
{
  sunrealtype H              = 0.0;
  sunrealtype* y             = N_VGetArrayPointer(yvec);
  const sunrealtype sqrt_qTq = SUNRsqrt(y[0] * y[0] + y[1] * y[1]);
  const sunrealtype pTp      = y[2] * y[2] + y[3] * y[3];

  H = SUN_RCONST(0.5) * pTp - SUN_RCONST(1.0) / sqrt_qTq;

  return H;
}

sunrealtype AngularMomentum(N_Vector yvec)
{
  sunrealtype L        = 0.0;
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  L = q1 * p2 - q2 * p1;

  return L;
}

int dydt(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  int retval = 0;

  retval += force(t, yvec, ydotvec, user_data);
  retval += velocity(t, yvec, ydotvec, user_data);

  return retval;
}

int velocity(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  sunrealtype* ydot    = N_VGetArrayPointer(ydotvec);
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  ydot[0] = p1;
  ydot[1] = p2;

  return 0;
}

int force(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y             = N_VGetArrayPointer(yvec);
  sunrealtype* ydot          = N_VGetArrayPointer(ydotvec);
  const sunrealtype q1       = y[0];
  const sunrealtype q2       = y[1];
  const sunrealtype sqrt_qTq = SUNRsqrt(q1 * q1 + q2 * q2);

  ydot[2] = -q1 / SUNRpowerR(sqrt_qTq, SUN_RCONST(3.0));
  ydot[3] = -q2 / SUNRpowerR(sqrt_qTq, SUN_RCONST(3.0));

  return 0;
}

int rootfn(sunrealtype t, N_Vector yvec, sunrealtype* gout, void* user_data)
{
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  const sunrealtype q2 = y[1];

  gout[0] = q2;

  return 0;
}

int ComputeConvergence(int num_dt, sunrealtype* orders,
                       sunrealtype expected_order, sunrealtype a11,
                       sunrealtype a12, sunrealtype a21, sunrealtype a22,
                       sunrealtype b1, sunrealtype b2, sunrealtype* ord_avg,
                       sunrealtype* ord_max, sunrealtype* ord_est)
{
  /* Compute/print overall estimated convergence rate */
  sunrealtype det = 0;
  *ord_avg = 0, *ord_max = 0, *ord_est = 0;
  for (int i = 1; i < num_dt; i++)
  {
    *ord_avg += orders[i - 1];
    *ord_max = SUNMAX(*ord_max, orders[i - 1]);
  }
  *ord_avg = *ord_avg / ((realtype)num_dt - 1);
  det      = a11 * a22 - a12 * a21;
  *ord_est = (a11 * b2 - a21 * b1) / det;
  return 0;
}

int ParseArgs(int argc, char* argv[], ProgramArgs* args)
{
  args->step_mode        = 0;
  args->stepper          = 0;
  args->method_name      = NULL;
  args->count_orbits     = 0;
  args->use_compsums     = 0;
  args->dt               = SUN_RCONST(1e-2);
  args->tf               = SUN_RCONST(100.);
  args->check_order      = 0;
  args->num_output_times = 50;

  for (int argi = 1; argi < argc; argi++)
  {
    if (!strcmp(argv[argi], "--step-mode"))
    {
      argi++;
      if (!strcmp(argv[argi], "fixed")) { args->step_mode = 0; }
      else if (!strcmp(argv[argi], "adapt")) { args->step_mode = 1; }
      else
      {
        fprintf(stderr, "ERROR: --step-mode must be 'fixed' or 'adapt'\n");
        return 1;
      }
    }
    else if (!strcmp(argv[argi], "--stepper"))
    {
      argi++;
      if (!strcmp(argv[argi], "SPRK")) { args->stepper = 0; }
      else if (!strcmp(argv[argi], "ERK")) { args->stepper = 1; }
      else
      {
        fprintf(stderr, "ERROR: --stepper must be 'SPRK' or 'ERK'\n");
        return 1;
      }
    }
    else if (!strcmp(argv[argi], "--method"))
    {
      argi++;
      args->method_name = argv[argi];
    }
    else if (!strcmp(argv[argi], "--dt"))
    {
      argi++;
      args->dt = atof(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--tf"))
    {
      argi++;
      args->tf = atof(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--nout"))
    {
      argi++;
      args->num_output_times = atoi(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--count-orbits")) { args->count_orbits = 1; }
    else if (!strcmp(argv[argi], "--use-compensated-sums"))
    {
      args->use_compsums = 1;
    }
    else if (!strcmp(argv[argi], "--check-order")) { args->check_order = 1; }
    else if (!strcmp(argv[argi], "--help"))
    {
      PrintHelp();
      return 1;
    }
    else
    {
      fprintf(stderr, "ERROR: unrecognized argument %s\n", argv[argi]);
      PrintHelp();
      return 1;
    }
  }

  if (!args->method_name)
  {
    if (args->stepper == 0)
    {
      args->method_name = "ARKODE_SYMPLECTIC_MCLACHLAN_4_4";
    }
    else if (args->stepper == 1)
    {
      args->method_name = "ARKODE_ZONNEVELD_5_3_4";
    }
  }

  return 0;
}

void PrintArgs(ProgramArgs* args)
{
  fprintf(stdout, "Problem Arguments:\n");
  fprintf(stdout, "  stepper:              %d\n", args->stepper);
  fprintf(stdout, "  step mode:            %d\n", args->step_mode);
  fprintf(stdout, "  use compensated sums: %d\n", args->use_compsums);
  fprintf(stdout, "  dt:                   %Lg\n", (long double)args->dt);
  fprintf(stdout, "  Tf:                   %Lg\n", (long double)args->tf);
  fprintf(stdout, "  nout:                 %d\n\n", args->num_output_times);
}

void PrintHelp()
{
  fprintf(stderr, "ark_kepler: an ARKODE example demonstrating the SPRKStep "
                  "time-stepping module solving the Kepler problem\n");
  fprintf(stderr, "  --step-mode <fixed, adapt>  should we use a fixed "
                  "time-step or adaptive time-step (default fixed)\n");
  fprintf(stderr,
          "  --stepper <SPRK, ERK>       should we use SPRKStep or ARKStep "
          "with an ERK method (default SPRK)\n");
  fprintf(stderr, "  --method <string>           which method to use (default "
                  "ARKODE_SYMPLECTIC_MCLACHLAN_4_4)\n");
  fprintf(stderr,
          "  --use-compensated-sums      turns on compensated summation in "
          "ARKODE where applicable\n");
  fprintf(stderr, "  --dt <Real>                 the fixed-time step size to "
                  "use if fixed "
                  "time stepping is turned on (default 0.01)\n");
  fprintf(stderr, "  --tf <Real>                 the final time for the "
                  "simulation (default 100)\n");
  fprintf(stderr, "  --nout <int>                the number of output times "
                  "(default 100)\n");
  fprintf(stderr, "  --count-orbits              use rootfinding to count the "
                  "number of completed orbits\n");
  fprintf(stderr,
          "  --check-order               compute the order of the method used "
          "and check if it is within range of the expected\n");
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval < 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval = NULL;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}
