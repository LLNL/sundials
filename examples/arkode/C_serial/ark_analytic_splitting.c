/*-----------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem with analytical
 * solution,
 *    dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lambda".  The value of lambda should be negative to
 * result in a well-posed ODE; for values with magnitude larger
 * than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the DIRK method,
 * Newton iteration with the dense SUNLinearSolver, and a
 * user-supplied Jacobian routine.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <nvector/nvector_serial.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

typedef struct
{
  sunrealtype lambda;
} UserData;

/* RHS for f^1(t, y) = -lambda * y */
static int f_linear(const sunrealtype t, const N_Vector y, N_Vector ydot,
                    void* const user_data)
{
  sunrealtype lambda = ((UserData*)user_data)->lambda;
  N_VScale(-lambda, y, ydot);
  return 0;
}

/* RHS for f^2(t, y) = y^2 */
static int f_nonlinear(const sunrealtype t, const N_Vector y, N_Vector ydot,
                       void* const user_data)
{
  N_VProd(y, y, ydot);
  return 0;
}

/* Compute the exact analytic solution */
static N_Vector exact_sol(const N_Vector y0, const sunrealtype tf,
                          const UserData* const user_data)
{
  N_Vector sol             = N_VClone(y0);
  const sunrealtype y0_val = NV_Ith_S(y0, 0);
  const sunrealtype lambda = user_data->lambda;
  NV_Ith_S(sol, 0)         = lambda * y0_val /
                     (y0_val - (y0_val - lambda) * SUNRexp(lambda * tf));
  return sol;
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

int main()
{
  /* Problem parameters */
  const sunrealtype t0 = SUN_RCONST(0.0);   /* initial time */
  const sunrealtype tf = SUN_RCONST(1.0);   /* final time */
  const sunrealtype dt = SUN_RCONST(0.01);  /* operator splitting time step */
  const sunrealtype dt_linear    = dt / 5;  /* linear integrator time step */
  const sunrealtype dt_nonlinear = dt / 10; /* nonlinear integrator time step */

  UserData user_data = {.lambda = SUN_RCONST(2.0)};

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  int flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  /* Initialize vector with initial condition */
  N_Vector y = N_VNew_Serial(1, ctx);
  if (check_flag(y, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(SUN_RCONST(1.0), y);

  N_Vector y_exact = exact_sol(y, tf, &user_data);

  printf("\nAnalytical ODE test problem:\n");
  printf("   lambda = %" GSYM "\n", user_data.lambda);

  /* Create the integrator for the linear partition */
  void* linear_mem = ARKStepCreate(f_linear, NULL, t0, y, ctx);
  if (check_flag(linear_mem, "N_VNew_Serial", 0)) { return 1; }

  flag = ARKodeSetUserData(linear_mem, &user_data);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  flag = ARKodeSetFixedStep(linear_mem, dt_linear);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  /* Create the integrator for the nonlinear partition */
  void* nonlinear_mem = ARKStepCreate(f_nonlinear, NULL, t0, y, ctx);
  if (check_flag(nonlinear_mem, "N_VNew_Serial", 0)) { return 1; }

  flag = ARKodeSetFixedStep(nonlinear_mem, dt_nonlinear);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  /* Create SUNSteppers out of the integrators */
  SUNStepper steppers[2];
  ARKStepCreateSUNStepper(linear_mem, &steppers[0]);
  ARKStepCreateSUNStepper(nonlinear_mem, &steppers[1]);

  /* Create the operator splitting method */
  void* splitting_mem = SplittingStepCreate(steppers, 2, t0, y, ctx);
  if (check_flag(splitting_mem, "SplittingStepCreate", 0)) { return 1; }

  flag = ARKodeSetFixedStep(splitting_mem, dt);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  /* Compute the operator splitting solution */
  sunrealtype tret;
  flag = ARKodeEvolve(splitting_mem, tf, y, &tret, ARK_NORMAL);
  if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }

  /* Print the numerical error */
  N_Vector y_err = N_VClone(y);
  if (check_flag(y_err, "N_VClone", 0)) { return 1; }
  N_VLinearSum(SUN_RCONST(1.0), y, -SUN_RCONST(1.0), y_exact, y_err);
  printf("Error: %" GSYM "\n", N_VMaxNorm(y_err));

  /* Free memory */
  N_VDestroy(y);
  N_VDestroy(y_exact);
  N_VDestroy(y_err);
  ARKodeFree(&linear_mem);
  SUNStepper_Destroy(&steppers[0]);
  ARKodeFree(&nonlinear_mem);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&splitting_mem);
  SUNContext_Free(&ctx);

  return 0;
}
