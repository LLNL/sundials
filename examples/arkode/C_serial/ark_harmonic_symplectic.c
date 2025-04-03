/* clang-format: off */
/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * In this example we consider the simple harmonic oscillator
 *    x''(t) + omega^2*x(t) = 0.
 * We rewrite the second order ODE as the first order ODE model
 *    x'(t) = v(t)
 *    v'(t) = -omega^2*x(t).
 * With the initial conditions x(0) = x0 and v(0) = v0,
 * the analytical solution is
 *    x(t) = A*cos(t*omega + phi),
 *    v(t) = -A*omega*sin(t*omega + phi)
 * where A = sqrt(x0^2 + v0^2/omega) and tan(phi) = v0/(omega*x0).
 * The total energy (potential + kinetic) in this system is
 *    E = (v^2 + omega^2*x^2) / 2
 * E is conserved and is the system Hamiltonian.
 * We simulate the problem on t = [0, 2pi] using the symplectic methods
 * in SPRKStep. Symplectic methods will approximately conserve E.
 *
 * The example has the following command line arguments:
 *   --order <int>               the order of the method to use (default 4)
 *   --dt <Real>                 the fixed-time step size to use (default 0.01)
 *   --nout <int>                the number of output times (default 100)
 *   --use-compensated-sums      turns on compensated summation in ARKODE where
 *                               applicable
 *   --disable-tstop             turns off tstop mode
 * --------------------------------------------------------------------------*/
/* clang-format: on */

#include "ark_harmonic_symplectic.h"

#include <arkode/arkode.h>
#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h> /* prototypes for SPRKStep fcts., consts */
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector type, fcts., macros  */
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_math.h> /* def. math fcns, 'sunrealtype'           */
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

typedef struct
{
  sunrealtype A, phi, omega;
} UserData;

/* RHS functions */
static int xdot(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int vdot(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Other helper functions */
static void Solution(sunrealtype t, N_Vector y, N_Vector solvec, UserData* udata);
static sunrealtype Energy(N_Vector yvec, sunrealtype dt, UserData* udata);

int main(int argc, char* argv[])
{
  ProgramArgs args;
  UserData udata;
  SUNContext sunctx       = NULL;
  N_Vector y              = NULL;
  N_Vector solution       = NULL;
  sunrealtype* ydata      = NULL;
  sunrealtype tout        = NAN;
  sunrealtype tret        = NAN;
  sunrealtype err         = NAN;
  void* arkode_mem        = NULL;
  void* q_mem             = NULL;
  void* p_mem             = NULL;
  int iout                = 0;
  int retval              = 0;
  int order               = 0;
  int use_compsums        = 0;
  int num_output_times    = 0;
  sunrealtype Tf          = SUN_RCONST(0.0);
  sunrealtype dt          = SUN_RCONST(0.0);
  sunrealtype dTout       = SUN_RCONST(0.0);
  const sunrealtype T0    = SUN_RCONST(0.0);
  const sunrealtype A     = SUN_RCONST(10.0);
  const sunrealtype phi   = SUN_RCONST(0.0);
  const sunrealtype omega = SUN_RCONST(1.0);

  /* Parse the command line arguments */
  if (ParseArgs(argc, argv, &args)) { return 1; };

  /* Default integrator options and problem parameters */
  order            = args.order;
  use_compsums     = args.use_compsums;
  num_output_times = args.num_output_times;
  Tf               = args.Tf;
  dt               = args.dt;
  dTout            = (Tf - T0) / ((sunrealtype)num_output_times);

  /* Default problem parameters */

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  printf("\n   Begin simple harmonic oscillator problem\n\n");

  /* Allocate and fill udata structure */
  udata.A     = A;
  udata.phi   = phi;
  udata.omega = omega;

  /* Allocate our state vector [x, v]^T */
  y        = N_VNew_Serial(2, sunctx);
  solution = N_VClone(y);

  /* Fill the initial conditions (x0 then v0) */
  ydata    = N_VGetArrayPointer(y);
  ydata[0] = A * cos(phi);
  ydata[1] = -A * omega * sin(phi);

  if (order > 0) // use SPRKStep
  {
    /* Create SPRKStep integrator */
    arkode_mem = SPRKStepCreate(xdot, vdot, T0, y, sunctx);

    retval = ARKodeSetOrder(arkode_mem, order);
    if (check_retval(&retval, "ARKodeSetOrder", 1)) { return 1; }

    retval = SPRKStepSetUseCompensatedSums(arkode_mem, use_compsums);
    if (check_retval(&retval, "SPRKStepSetUseCompensatedSums", 1)) { return 1; }

    retval = ARKodeSetUserData(arkode_mem, &udata);
    if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }
  }
  else // use SplittingStep
  {
    q_mem = ARKStepCreate(xdot, NULL, T0, y, sunctx);

    retval = ARKodeSetOrder(q_mem, 1);
    if (check_retval(&retval, "ARKodeSetOrder", 1)) { return 1; }

    retval = ARKodeSetFixedStep(q_mem, dt);
    if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

    retval = ARKodeSetUserData(q_mem, &udata);
    if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }

    p_mem = ARKStepCreate(vdot, NULL, T0, y, sunctx);

    retval = ARKodeSetOrder(p_mem, 1);
    if (check_retval(&retval, "ARKodeSetOrder", 1)) { return 1; }

    retval = ARKodeSetFixedStep(p_mem, dt);
    if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

    retval = ARKodeSetUserData(p_mem, &udata);
    if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }

    SUNStepper steppers[2];
    ARKodeCreateSUNStepper(q_mem, &steppers[0]);
    ARKodeCreateSUNStepper(p_mem, &steppers[1]);
    arkode_mem = SplittingStepCreate(steppers, 2, T0, y, sunctx);

    int table_id;
    switch (-order)
    {
    case 1: table_id = SPRKSTEP_DEFAULT_1; break;
    case 2: table_id = SPRKSTEP_DEFAULT_2; break;
    case 3: table_id = SPRKSTEP_DEFAULT_3; break;
    case 4: table_id = SPRKSTEP_DEFAULT_4; break;
    case 5: table_id = SPRKSTEP_DEFAULT_5; break;
    case 6: table_id = SPRKSTEP_DEFAULT_6; break;
    case 8: table_id = SPRKSTEP_DEFAULT_8; break;
    case 10: table_id = SPRKSTEP_DEFAULT_10; break;
    default:
      fprintf(stderr, "Specified order isn't valid SPRKStep order\n");
      break;
    }
    ARKodeSPRKTable table = ARKodeSPRKTable_Load(table_id);

    sunrealtype alpha = 1.0;
    sunrealtype* beta = malloc((2*(table->stages+1))*sizeof(sunrealtype));
    beta[0] = 0.0;
    beta[1] = 0.0;
    for (int i=0; i < table->stages; i++)
    {
      beta[2*(i+1)] = beta[2*i] + table->ahat[i];
      beta[2*(i+1)+1] = beta[2*(i+1)-1] + table->a[i];
    }
    SplittingStepCoefficients coefficients =
      SplittingStepCoefficients_Create(1, table->stages, 2, table->q, &alpha,
        beta);
    retval = SplittingStepSetCoefficients(arkode_mem, coefficients);
    if (check_retval(&retval, "SplittingStepSetCoefficients", 1)) { return 1; }
    free(beta);
    SplittingStepCoefficients_Destroy(&coefficients);
  }

  retval = ARKodeSetFixedStep(arkode_mem, dt);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

  retval = ARKodeSetMaxNumSteps(arkode_mem, ((long int)ceil(Tf / dt)) + 2);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  /* Print out starting energy, momentum before integrating */
  tret = T0;
  tout = T0 + dTout;
  fprintf(stdout, "t = %.6Lf, x(t) = %.6Lf, E = %.6Lf, sol. err = %.6Lf\n",
          (long double)tret, (long double)ydata[0],
          (long double)Energy(y, dt, &udata), (long double)SUN_RCONST(0.0));

  /* Do integration */
  for (iout = 0; iout < num_output_times; iout++)
  {
    if (args.use_tstop) { ARKodeSetStopTime(arkode_mem, tout); }
    retval = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

    /* Compute the analytical solution */
    Solution(tret, y, solution, &udata);

    /* Compute L2 error */
    N_VLinearSum(SUN_RCONST(1.0), y, -SUN_RCONST(1.0), solution, solution);
    err = sqrt(N_VDotProd(solution, solution));

    /* Output current integration status */
    fprintf(stdout, "t = %.6Lf, x(t) = %.6Lf, E = %.6Lf, sol. err = %.16Le\n",
            (long double)tret, (long double)ydata[0],
            (long double)Energy(y, dt, &udata), (long double)err);

    /* Check that solution error is within tolerance */
    if (err > SUNMAX(dt / pow(10, order - 2), 1000 * SUN_UNIT_ROUNDOFF))
    {
      fprintf(stderr, "FAILURE: solution error is too high\n");
      return 1;
    }

    /* Check if the solve was successful, if so, update the time and continue */
    if (retval >= 0)
    {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    else
    {
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
  }

  fprintf(stdout, "\n");
  N_VDestroy(y);
  N_VDestroy(solution);
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  ARKodeFree(&arkode_mem);
  SUNContext_Free(&sunctx);

  return 0;
}

void Solution(sunrealtype t, N_Vector y, N_Vector solvec, UserData* udata)
{
  sunrealtype* sol = N_VGetArrayPointer(solvec);

  /* compute solution */
  sol[0] = udata->A * cos(udata->omega * t + udata->phi);
  sol[1] = -udata->A * udata->omega * sin(udata->omega * t + udata->phi);
}

sunrealtype Energy(N_Vector yvec, sunrealtype dt, UserData* udata)
{
  sunrealtype E            = 0.0;
  sunrealtype* y           = N_VGetArrayPointer(yvec);
  const sunrealtype x      = y[0];
  const sunrealtype v      = y[1];
  const sunrealtype omega2 = udata->omega * udata->omega;

  E = (v * v + omega2 * x * x) / SUN_RCONST(2.0);

  return E;
}

int xdot(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y      = N_VGetArrayPointer(yvec);
  sunrealtype* ydot   = N_VGetArrayPointer(ydotvec);
  const sunrealtype v = y[1];

  ydot[0] = v;

  return 0;
}

int vdot(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  UserData* udata          = (UserData*)user_data;
  sunrealtype* y           = N_VGetArrayPointer(yvec);
  sunrealtype* ydot        = N_VGetArrayPointer(ydotvec);
  const sunrealtype x      = y[0];
  const sunrealtype omega2 = udata->omega * udata->omega;

  ydot[1] = -omega2 * x;

  return 0;
}
