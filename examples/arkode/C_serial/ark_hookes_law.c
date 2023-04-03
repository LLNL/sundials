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
 * This example considers the motion of a spring with an attached mass using
 * the symplectic integrators available in SPRKStep. Our ODE model is,
 *    x'(t) = v
 *    v'(t) = -k/m x
 * where k is the spring constant and m is the mass. For convenience we choose
 * k = 8 and m = 1 and let omega = sqrt(k/m). We simulate the problem on
 * t = [0, 0.1] with the initial conditions x(0) = 0.2, v(0) = 0.
 * This problem setup has the exact solution
 *    x(t) = 0.2*cos(t*omega),
 *    v(t) = -0.2*omega*sin(t*omega).
 * The potential energy,
 *    U = 1/2*k*x^2
 * is conserved and is the Hamiltonian for the system. The symplectic methods
 * in SPRKStep conserve U provided a sufficiently small time-step size is used.
 *
 * The problem can be run like so:
 *    ./ark_hookes_law [order] [dt] [use_compsums]
 *
 * Order sets the order of the method to use, dt is the time step size, and
 * use_compsums turns on (1) or off (0) compensated summation inside SPRKStep.
 * Compensated summation increases accuracy but at increased cost. 
 * --------------------------------------------------------------------------*/

#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h> /* prototypes for MRIStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector type, fcts., macros  */
#include <stdio.h>
#include <sundials/sundials_math.h> /* def. math fcns, 'sunrealtype'           */
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "arkode/arkode.h"

typedef struct
{
  sunrealtype A, B, omega;
} * UserData;

static int check_retval(void* returnvalue, const char* funcname, int opt);

static void InitialConditions(N_Vector y0);
static sunrealtype Solution(sunrealtype t, N_Vector y, N_Vector solvec,
                            UserData udata);
static sunrealtype Energy(N_Vector yvec, sunrealtype dt, UserData udata);

static int Velocity(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Force(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

int main(int argc, char* argv[])
{
  SUNContext sunctx;
  N_Vector y, solution;
  SUNNonlinearSolver NLS;
  UserData udata;
  sunrealtype tout, tret;
  void* arkode_mem;
  int argi, iout, retval;

  NLS      = NULL;
  y        = NULL;
  solution = NULL;

  /* Default problem parameters */
  const sunrealtype T0    = SUN_RCONST(0.0);
  sunrealtype Tf          = SUN_RCONST(0.1);
  sunrealtype dt          = SUN_RCONST(1e-4);
  const sunrealtype A     = SUN_RCONST(0.2);
  const sunrealtype B     = SUN_RCONST(0.0);
  const sunrealtype omega = SUN_RCONST(64.0);

  /* Default integrator Options */
  int method                 = 0;
  int order                  = 4;
  int use_compsums           = 1;
  const sunrealtype dTout    = SUN_RCONST(0.01);
  const int num_output_times = (int)ceil(Tf / dTout);

  printf("\n   Begin Hooke's Law Problem\n\n");

  /* Parse CLI args */
  argi = 0;
  if (argc > 1) { order = atoi(argv[++argi]); }
  if (argc > 2) { dt = atof(argv[++argi]); }
  if (argc > 3) { use_compsums = atoi(argv[++argi]); }

  /* Allocate and fill udata structure */
  udata        = (UserData)malloc(sizeof(*udata));
  udata->A     = A;
  udata->B     = B;
  udata->omega = omega;

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  /* Allocate our state vector */
  y        = N_VNew_Serial(2, sunctx);
  solution = N_VClone(y);

  /* Fill the initial conditions */
  InitialConditions(y);

  /* Create SPRKStep integrator */
  arkode_mem = SPRKStepCreate(Force, Velocity, T0, y, sunctx);

  retval = SPRKStepSetOrder(arkode_mem, order);
  if (check_retval(&retval, "SPRKStepSetOrder", 1)) return 1;

  retval = SPRKStepSetUserData(arkode_mem, udata);
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
  fprintf(stdout, "t = %.6Lf, energy = %.6Lf\n", tret, Energy(y, dt, udata));

  /* Do integration */
  for (iout = 0; iout < num_output_times; iout++)
  {
    SPRKStepSetStopTime(arkode_mem, tout);
    retval = SPRKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

    /* Output current integration status */
    fprintf(stdout, "t = %.6Lf, sol. err = %.6Lf, energy = %.6Lf\n", tret,
            Solution(tret, y, solution, udata), Energy(y, dt, udata));

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

  free(udata);
  N_VDestroy(y);
  fprintf(stdout, "\n");
  SPRKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  SPRKStepFree(&arkode_mem);
  SUNContext_Free(&sunctx);

  return 0;
}

void InitialConditions(N_Vector y0vec)
{
  sunrealtype* y0 = N_VGetArrayPointer(y0vec);
  y0[0]           = SUN_RCONST(0.2);
  y0[1]           = SUN_RCONST(0.0);
}

sunrealtype Solution(sunrealtype t, N_Vector y, N_Vector solvec, UserData udata)
{
  sunrealtype* sol = N_VGetArrayPointer(solvec);
  /* compute error */
  sol[0] = SUN_RCONST(0.2) * cos(udata->omega * t);
  sol[1] = -SUN_RCONST(0.2) * udata->omega * sin(udata->omega * t);
  N_VLinearSum(SUN_RCONST(1.0), y, -SUN_RCONST(1.0), solvec, solvec);
  sunrealtype err = N_VMaxNorm(solvec);
  /* restore solution vec */
  sol[0] = SUN_RCONST(0.2) * cos(udata->omega * t);
  sol[1] = -SUN_RCONST(0.2) * udata->omega * sin(udata->omega * t);
  return err;
}

sunrealtype Energy(N_Vector yvec, sunrealtype dt, UserData udata)
{
  sunrealtype H            = 0.0;
  sunrealtype* y           = N_VGetArrayPointer(yvec);
  const sunrealtype x      = y[0];
  const sunrealtype v      = y[1];
  const sunrealtype omega2 = udata->omega * udata->omega;

  H = (v * v + omega2 * x * x) / SUN_RCONST(2.0);

  return H;
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
  UserData udata           = (UserData)user_data;
  sunrealtype* y           = N_VGetArrayPointer(yvec);
  sunrealtype* ydot        = N_VGetArrayPointer(ydotvec);
  const sunrealtype omega2 = udata->omega * udata->omega;

  ydot[1] = -omega2 * y[0];

  return 0;
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
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
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
