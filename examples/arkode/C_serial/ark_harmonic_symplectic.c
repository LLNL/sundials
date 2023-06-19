/* clang-format: off */
/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
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
 * in SPRKStep. Symplectic methods will approximately conserve U.
 *
 * The problem can be run like so:
 *    ./ark_harmonic_symplectic [order] [dt] [use_compsums]
 *
 * Order sets the order of the method to use, dt is the time step size, and
 * use_compsums turns on (1) or off (0) compensated summation inside SPRKStep.
 * Compensated summation increases accuracy but at increased computational
 * and memory cost.
 * --------------------------------------------------------------------------*/
/* clang-format: on */

#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h> /* prototypes for SPRKStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector type, fcts., macros  */
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_math.h> /* def. math fcns, 'sunrealtype'           */
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "arkode/arkode.h"

typedef struct
{
  sunrealtype A, phi, omega;
} UserData;

typedef struct
{
  int order;
  int num_output_times;
  int use_compsums;
  sunrealtype dt;
} ProgramArgs;

/* RHS functions */
static int Velocity(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Force(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Helper functions */
static void InitialConditions(N_Vector y0);
static void Solution(sunrealtype t, N_Vector y, N_Vector solvec, UserData* udata);
static sunrealtype Energy(N_Vector yvec, sunrealtype dt, UserData* udata);
static int ParseArgs(int argc, char* argv[], ProgramArgs* args);
static void PrintHelp();
static int check_retval(void* returnvalue, const char* funcname, int opt);

int main(int argc, char* argv[])
{
  ProgramArgs args;
  UserData udata;
  SUNContext sunctx      = NULL;
  N_Vector y             = NULL;
  N_Vector solution      = NULL;
  SUNNonlinearSolver NLS = NULL;
  sunrealtype* ydata     = NULL;
  sunrealtype tout       = NAN;
  sunrealtype tret       = NAN;
  sunrealtype err        = NAN;
  void* arkode_mem       = NULL;
  int iout               = 0;
  int retval             = 0;

  /* Parse the command line arguments */
  if (ParseArgs(argc, argv, &args)) { return 1; };

  /* Default integrator options */
  int order                  = args.order;
  int use_compsums           = args.use_compsums;
  const int num_output_times = args.num_output_times;

  /* Default problem parameters */
  const sunrealtype T0    = SUN_RCONST(0.0);
  sunrealtype Tf          = SUN_RCONST(2.0) * M_PI;
  sunrealtype dt          = args.dt;
  const sunrealtype A     = SUN_RCONST(10.0);
  const sunrealtype phi   = SUN_RCONST(0.0);
  const sunrealtype omega = SUN_RCONST(1.0);
  const sunrealtype dTout = (Tf - T0) / ((sunrealtype)num_output_times);

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  printf("\n   Begin simple harmonic oscillator problem\n\n");

  /* Allocate and fill udata structure */
  udata.A     = A;
  udata.phi   = phi;
  udata.omega = omega;

  /* Allocate our state vector */
  y        = N_VNew_Serial(2, sunctx);
  solution = N_VClone(y);

  /* Fill the initial conditions */
  ydata    = N_VGetArrayPointer(y);
  ydata[0] = A*cos(phi);
  ydata[1] = -A*omega*sin(phi);

  /* Create SPRKStep integrator */
  arkode_mem = SPRKStepCreate(Force, Velocity, T0, y, sunctx);

  retval = SPRKStepSetOrder(arkode_mem, order);
  if (check_retval(&retval, "SPRKStepSetOrder", 1)) return 1;

  retval = SPRKStepSetUserData(arkode_mem, &udata);
  if (check_retval(&retval, "SPRKStepSetUserData", 1)) return 1;

  retval = SPRKStepSetUseCompensatedSums(arkode_mem, use_compsums);
  if (check_retval(&retval, "SPRKStepSetUseCompensatedSums", 1)) return 1;

  retval = SPRKStepSetFixedStep(arkode_mem, dt);
  if (check_retval(&retval, "SPRKStepSetFixedStep", 1)) return 1;

  retval = SPRKStepSetMaxNumSteps(arkode_mem, ((long int)ceil(Tf / dt)) + 2);
  if (check_retval(&retval, "SPRKStepSetMaxNumSteps", 1)) return 1;

  /* Print out starting energy, momentum before integrating */
  tret = T0;
  tout = T0 + dTout;
  fprintf(stdout, "t = %.6Lf, x(t) = %.6Lf, E = %.6Lf, sol. err = %.6Lf\n", tret,
          ydata[0], Energy(y, dt, &udata), SUN_RCONST(0.0));

  /* Do integration */
  for (iout = 0; iout < num_output_times; iout++)
  {
    SPRKStepSetStopTime(arkode_mem, tout);
    retval = SPRKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

    /* Compute the anaytical solution */
    Solution(tret, y, solution, &udata);

    /* Compute L2 error */
    N_VLinearSum(SUN_RCONST(1.0), y, -SUN_RCONST(1.0), solution, solution);
    err = sqrt(N_VDotProd(solution, solution));

    /* Output current integration status */
    fprintf(stdout, "t = %.6Lf, x(t) = %.6Lf, E = %.6Lf, sol. err = %.16Lf\n", tret,
            ydata[0], Energy(y, dt, &udata), err);

    /* Check that solution error is within tolerance */
    if (err > SUNMAX(dt / pow(10, order-2), 1000*SUN_UNIT_ROUNDOFF))
    {
      fprintf(stderr, "FAILURE: solution error is too high\n");
      return 1;
    }

    /* Check if the solve was successful, if so, update the time and continue
     */
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

  N_VDestroy(y);
  fprintf(stdout, "\n");
  SPRKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  SPRKStepFree(&arkode_mem);
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

int Velocity(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y    = N_VGetArrayPointer(yvec);
  sunrealtype* ydot = N_VGetArrayPointer(ydotvec);

  ydot[0] = y[1];

  return 0;
}

int Force(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  UserData* udata          = (UserData*)user_data;
  sunrealtype* y           = N_VGetArrayPointer(yvec);
  sunrealtype* ydot        = N_VGetArrayPointer(ydotvec);
  const sunrealtype omega2 = udata->omega * udata->omega;

  ydot[1] = -omega2 * y[0];

  return 0;
}

int ParseArgs(int argc, char* argv[], ProgramArgs* args)
{
  args->order            = 4;
  args->num_output_times = 8;
  args->use_compsums     = 0;
  args->dt               = SUN_RCONST(1e-3);

  for (int argi = 1; argi < argc; argi++)
  {
    if (!strcmp(argv[argi], "--order"))
    {
      argi++;
      args->order = atoi(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--dt"))
    {
      argi++;
      args->dt = atof(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--nout"))
    {
      argi++;
      args->num_output_times = atoi(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--use-compensated-sums"))
    {
      args->use_compsums = 1;
    }
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

  return 0;
}

void PrintHelp()
{
  fprintf(stderr, "ark_harmonic_symplectic: an ARKODE example demonstrating "
                  "the SPRKStep time-stepping module solving a simple harmonic "
                  "oscillator\n");
  fprintf(stderr, "  --order <int>               the order of the method to "
                  "use (default 4)\n");
  fprintf(stderr,
          "  --dt <Real>                 the fixed-time step size to use "
          "(default 0.01)\n");
  fprintf(stderr, "  --nout <int>                the number of output times "
                  "(default 100)\n");
  fprintf(stderr,
          "  --use-compensated-sums      turns on compensated summation in "
          "ARKODE where applicable\n");
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
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return 1;
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nERROR: %s() failed with retval = %d\n\n", funcname,
              *retval);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}
