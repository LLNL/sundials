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
static sunrealtype Hamiltonian(N_Vector yvec, sunrealtype dt, UserData udata);

static int velocity(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int force(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

int main(int argc, char* argv[])
{
  SUNContext sunctx;
  N_Vector y;
  SUNNonlinearSolver NLS;
  UserData udata;
  sunrealtype tout, tret;
  sunrealtype H0;
  void* arkode_mem;
  // FILE *conserved_fp, *solution_fp, *times_fp;
  int argi, iout, retval;

  NLS = NULL;
  y   = NULL;

  /* Default problem parameters */
  const sunrealtype T0    = SUN_RCONST(0.0);
  sunrealtype Tf          = SUN_RCONST(0.1);
  sunrealtype dt          = SUN_RCONST(1e-4);
  const sunrealtype A     = SUN_RCONST(0.2);
  const sunrealtype B     = SUN_RCONST(0.0);
  const sunrealtype omega = SUN_RCONST(64.);

  /* Default integrator Options */
  int step_mode              = 0;
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
  y = N_VNew_Serial(4, sunctx);

  /* Fill the initial conditions */
  InitialConditions(y);

  /* Create SPRKStep integrator */
  arkode_mem = SPRKStepCreate(force, velocity, T0, y, sunctx);

  switch (order)
  {
  case 1:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticEuler());
    fprintf(stdout, "Using symplectic Euler O(h^1)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 2:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticLeapfrog2());
    fprintf(stdout, "Using symplectic Leapfrog O(h^2)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 22:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticPseudoLeapfrog2());
    fprintf(stdout, "Using symplectic Pseudo Leapfrog O(h^2)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 222:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticMcLachlan2());
    fprintf(stdout, "Using symplectic McLachlan O(h^2)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 3:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticRuth3());
    fprintf(stdout, "Using symplectic Ruth O(h^3)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 33:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticMcLachlan3());
    fprintf(stdout, "Using symplectic McLachlan O(h^3)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 4:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticCandyRozmus4());
    fprintf(stdout, "Using symplectic Candy-Rozmus O(h^4)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 44:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticMcLachlan4());
    fprintf(stdout, "Using symplectic McLachlan O(h^4)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 5:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticMcLachlan5());
    fprintf(stdout, "Using symplectic McLachlan O(h^5)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 6:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticYoshida6());
    fprintf(stdout, "Using symplectic Yoshida O(h^6)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 8:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticMcLachlan8());
    fprintf(stdout, "Using symplectic McLachlan O(h^8)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  case 10:
    retval = SPRKStepSetMethod(arkode_mem, ARKodeSymplecticSofroniou10());
    fprintf(stdout, "Using symplectic Sofroniou O(h^10)\n\n");
    if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
    break;
  default: fprintf(stderr, "Not a valid method\n"); return 1;
  }

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
  H0   = Hamiltonian(y, dt, udata);
  fprintf(stdout, "t = %.4Lf, H(p,q) = %.16Lf\n", tret, H0);

  /* Do integration */
  for (iout = 0; iout < num_output_times; iout++)
  {
    SPRKStepSetStopTime(arkode_mem, tout);
    retval = SPRKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

    /* Output current integration status */
    fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf\n", tret, Hamiltonian(y, dt, udata) - H0);

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

void Solution(sunrealtype t, N_Vector solvec, UserData udata)
{
  sunrealtype* sol = N_VGetArrayPointer(solvec);

  sol[0] = udata->A * sin(t * udata->omega) + udata->B * cos(t * udata->omega);

  return;
}

sunrealtype Hamiltonian(N_Vector yvec, sunrealtype dt, UserData udata)
{
  sunrealtype H            = 0.0;
  sunrealtype* y           = N_VGetArrayPointer(yvec);
  const sunrealtype x      = y[0];
  const sunrealtype v      = y[1];
  const sunrealtype omega2 = udata->omega * udata->omega;

  H = (v*v + omega2 * x*x - omega2 * dt * v * x) / SUN_RCONST(2.0);

  return H;
}

int velocity(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  sunrealtype* ydot    = N_VGetArrayPointer(ydotvec);

  ydot[0] = ydot[1];

  return 0;
}

int force(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
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
