/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Routine to test that ARKStepReset, ERKStepReset and
 * MRIStepReset function correctly.
 *
 * This runs the same test problem as in
 * examples/arkode/C_serial/ark_analytic.c:
 *    dy/dt = lamda*y + 1/(1+t^2) - lambda*atan(t)
 * for t in various time intervals, with the initial condition
 * y(0)=0, and having analytical solution y(t) = atan(t).
 *
 * To test each time stepper, we:
 * 1. Initialize to analytical solution at 0, evolve over [0,0.1],
 *    and check result
 * 2. Reset to analytical solution at 0.1, evolve over [0.1,0.2],
 *    and check result
 * 3. Reset to analytical solution at 0.3, evolve over [0.3,0.4],
 *    and check result
 * 4. Reset to analytical solution at 0.1, evolve over [0.1,0.2],
 *    and check result
 *-----------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_mristep.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Analytical solution */
static sunrealtype ytrue(sunrealtype t);

/* Private function to check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

/* Private function to check computed solution */
static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol, sunrealtype atol);

/* Main Program */
int main()
{
  /* general problem parameters */
  sunrealtype T0     = SUN_RCONST(0.0);
  sunrealtype dTout  = SUN_RCONST(0.1);
  sunrealtype lambda = SUN_RCONST(-25.0);
  sunrealtype rtol   = SUN_RCONST(1e-3);
  sunrealtype atol   = SUN_RCONST(1e-6);

  /* general problem variables */
  N_Vector y = NULL;
  SUNMatrix A = NULL;
  SUNLinearSolver LS = NULL;
  void *arkode_mem = NULL;
  void *mristep_mem = NULL;
  MRIStepInnerStepper inner_stepper = NULL;
  int retval;
  sunrealtype t;

  SUNContext ctx;

  /* Initial diagnostics output */
  printf("\nARKODE reset functionality tester:\n");

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  /* Initialize vector, matrix, and linaer solver data structures */
  y = N_VNew_Serial(1, ctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return 1;
  A = SUNDenseMatrix(1, 1, ctx);
  if (check_retval((void *)A, "SUNDenseMatrix", 0)) return 1;
  LS = SUNLinSol_Dense(y, A, ctx);
  if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) return 1;


  /******* Part I: ERKStep *******/

  /* Set initial condition, and construct stepper */
  t = T0;
  N_VConst(ytrue(t), y);
  arkode_mem = ERKStepCreate(f, t, y, ctx);
  if (check_retval((void *)arkode_mem, "ERKStepCreate", 0)) return 1;
  retval = ERKStepSetUserData(arkode_mem, (void *) &lambda);
  if (check_retval(&retval, "ERKStepSetUserData", 1)) return 1;
  retval = ERKStepSStolerances(arkode_mem, rtol, atol);
  if (check_retval(&retval, "ERKStepSStolerances", 1)) return 1;
  retval = ERKStepSetMaxNumSteps(arkode_mem, 1000);
  check_retval(&retval, "ERKStepSetMaxNumSteps", 1);

  /* Initially evolve to dTout, and check result */
  retval = ERKStepSetStopTime(arkode_mem, t + dTout);
  check_retval(&retval, "ERKStepSetStopTime", 1);
  retval = ERKStepEvolve(arkode_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "ERKStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Initial ERKStepEvolve had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Initial ERKStepEvolve call successful\n");
  }


  /* Reset state to analytical solution at dTout, evolve to 2*dTout and check result */
  t = T0 + dTout;
  N_VConst(ytrue(t), y);
  retval = ERKStepReset(arkode_mem, t, y);
  check_retval(&retval, "ERKStepReset", 1);
  retval = ERKStepEvolve(arkode_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "ERKStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Second ERKStepEvolve call had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Second ERKStepEvolve call successful\n");
  }

  /* Reset state to analytical solution at 3*dTout, evolve to 4*dTout and check result */
  t = T0 + SUN_RCONST(3.0)*dTout;
  N_VConst(ytrue(t), y);
  retval = ERKStepReset(arkode_mem, t, y);
  check_retval(&retval, "ERKStepReset", 1);
  retval = ERKStepEvolve(arkode_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "ERKStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Third ERKStepEvolve call had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Third ERKStepEvolve call successful\n");
  }

  /* Reset state to analytical solution at dTout, evolve to 2*dTout and check result */
  t = T0 + dTout;
  N_VConst(ytrue(t), y);
  retval = ERKStepReset(arkode_mem, t, y);
  check_retval(&retval, "ERKStepReset", 1);
  retval = ERKStepEvolve(arkode_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "ERKStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Fourth ERKStepEvolve call had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Fourth ERKStepEvolve call successful\n");
  }

  /* Free ERKStep memory structure */
  ERKStepFree(&arkode_mem);
  arkode_mem = NULL;



  /******* Part II: ARKStep *******/

  /* Set initial condition, and construct stepper */
  t = T0;
  N_VConst(ytrue(t), y);
  arkode_mem = ARKStepCreate(NULL, f, t, y, ctx);
  if (check_retval((void *)arkode_mem, "ARKStepCreate", 0)) return 1;
  retval = ARKStepSetUserData(arkode_mem, (void *) &lambda);
  if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;
  retval = ARKStepSStolerances(arkode_mem, rtol, atol);
  if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
  retval = ARKStepSetLinearSolver(arkode_mem, LS, A);
  if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return 1;
  retval = ARKStepSetJacFn(arkode_mem, Jac);
  if (check_retval(&retval, "ARKStepSetJacFn", 1)) return 1;
  retval = ARKStepSetLinear(arkode_mem, 0);
  if (check_retval(&retval, "ARKStepSetLinear", 1)) return 1;
  retval = ARKStepSetMaxNumSteps(arkode_mem, 100);
  check_retval(&retval, "ARKStepSetMaxNumSteps", 1);

  /* Initially evolve to dTout, and check result */
  retval = ARKStepEvolve(arkode_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "ARKStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Initial ARKStepEvolve had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Initial ARKStepEvolve call successful\n");
  }

  /* Reset state to analytical solution at dTout, evolve to 2*dTout and check result */
  t = T0 + dTout;
  N_VConst(ytrue(t), y);
  retval = ARKStepReset(arkode_mem, t, y);
  check_retval(&retval, "ARKStepReset", 1);
  retval = ARKStepEvolve(arkode_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "ARKStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Second ARKStepEvolve call had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Second ARKStepEvolve call successful\n");
  }

  /* Reset state to analytical solution at 3*dTout, evolve to 4*dTout and check result */
  t = T0 + SUN_RCONST(3.0)*dTout;
  N_VConst(ytrue(t), y);
  retval = ARKStepReset(arkode_mem, t, y);
  check_retval(&retval, "ARKStepReset", 1);
  retval = ARKStepEvolve(arkode_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "ARKStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Third ARKStepEvolve call had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Third ARKStepEvolve call successful\n");
  }

  /* Reset state to analytical solution at dTout, evolve to 2*dTout and check result */
  t = T0 + dTout;
  N_VConst(ytrue(t), y);
  retval = ARKStepReset(arkode_mem, t, y);
  check_retval(&retval, "ARKStepReset", 1);
  retval = ARKStepEvolve(arkode_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "ARKStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Fourth ARKStepEvolve call had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Fourth ARKStepEvolve call successful\n");
  }

  /* Free ARKStep memory structure */
  ARKStepFree(&arkode_mem);
  arkode_mem = NULL;



  /******* Part III: MRIStep *******/

  /* Set initial condition, and construct stepper */
  t = T0;
  N_VConst(ytrue(t), y);
  arkode_mem = ARKStepCreate(f0, NULL, t, y, ctx);
  if (check_retval((void *)arkode_mem, "ARKStepCreate", 0)) return 1;
  retval = ARKStepSStolerances(arkode_mem, rtol, atol);
  if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
  retval = ARKStepSetMaxNumSteps(arkode_mem, 100);
  check_retval(&retval, "ARKStepSetMaxNumSteps", 1);
  retval = ARKStepCreateMRIStepInnerStepper(arkode_mem, &inner_stepper);
  if (check_retval(&retval, "ARKStepCreateMRIStepInnerStepper", 1)) return 1;
  mristep_mem = MRIStepCreate(NULL, f, t, y, inner_stepper, ctx);
  if (check_retval((void *)mristep_mem, "MRIStepCreate", 0)) return 1;
  retval = MRIStepSetUserData(mristep_mem, (void *) &lambda);
  if (check_retval(&retval, "MRIStepSetUserData", 1)) return 1;
  retval = MRIStepSetLinearSolver(mristep_mem, LS, A);
  if (check_retval(&retval, "MRIStepSetLinearSolver", 1)) return 1;
  retval = MRIStepSetJacFn(mristep_mem, Jac);
  if (check_retval(&retval, "MRIStepSetJacFn", 1)) return 1;
  retval = MRIStepSetLinear(mristep_mem, 0);
  if (check_retval(&retval, "MRIStepSetLinear", 1)) return 1;
  retval = MRIStepSetMaxNumSteps(mristep_mem, 100);
  check_retval(&retval, "MRIStepSetMaxNumSteps", 1);
  retval = MRIStepSetFixedStep(mristep_mem, dTout*SUN_RCONST(0.105));
  check_retval(&retval, "MRIStepSetFixedStep", 1);

  /* Initially evolve to dTout, and check result */
  retval = MRIStepEvolve(mristep_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "MRIStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Initial MRIStepEvolve had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Initial MRIStepEvolve call successful\n");
  }

  /* Reset state to analytical solution at dTout, evolve to 2*dTout and check result */
  t = T0 + dTout;
  N_VConst(ytrue(t), y);
  retval = MRIStepReset(mristep_mem, t, y);
  check_retval(&retval, "MRIStepReset", 1);
  retval = MRIStepEvolve(mristep_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "MRIStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Second MRIStepEvolve call had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Second MRIStepEvolve call successful\n");
  }

  /* Reset state to analytical solution at 3*dTout, evolve to 4*dTout and check result */
  t = T0 + SUN_RCONST(3.0)*dTout;
  N_VConst(ytrue(t), y);
  retval = MRIStepReset(mristep_mem, t, y);
  check_retval(&retval, "MRIStepReset", 1);
  retval = MRIStepEvolve(mristep_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "MRIStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Third MRIStepEvolve call had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Third MRIStepEvolve call successful\n");
  }

  /* Reset state to analytical solution at dTout, evolve to 2*dTout and check result */
  t = T0 + dTout;
  N_VConst(ytrue(t), y);
  retval = MRIStepReset(mristep_mem, t, y);
  check_retval(&retval, "MRIStepReset", 1);
  retval = MRIStepEvolve(mristep_mem, t + dTout, y, &t, ARK_NORMAL);
  if (check_retval(&retval, "MRIStepEvolve", 1)) return 1;
  if (check_ans(y, t, SUN_RCONST(0.001), SUN_RCONST(0.000001))) {
    printf("  Fourth MRIStepEvolve call had insufficient accuracy\n");
    printf("    t = %" GSYM "\n", t);
    printf("    y = %" GSYM "\n", NV_Ith_S(y,0));
    printf("    ytrue = %" GSYM "\n", ytrue(t));
    printf("    |y-ytrue| = %" GSYM "\n", SUNRabs(ytrue(t)-NV_Ith_S(y,0)));
    return 1;
  } else {
    printf("  Fourth MRIStepEvolve call successful\n");
  }

  /* Free MRIStep and ARKStep memory structures */
  MRIStepFree(&mristep_mem);
  MRIStepInnerStepper_Free(&inner_stepper);
  ARKStepFree(&arkode_mem);
  arkode_mem = NULL;

  /* Clean up and return with success */
  N_VDestroy(y);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNContext_Free(&ctx);
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  sunrealtype *rdata = (sunrealtype *) user_data;   /* cast user_data to sunrealtype */
  sunrealtype lambda = rdata[0];                 /* set shortcut for stiffness parameter */
  sunrealtype u = NV_Ith_S(y,0);                 /* access current solution value */

  /* fill in the RHS function: "NV_Ith_S" accesses the 0th entry of ydot */
  NV_Ith_S(ydot,0) = lambda*u + SUN_RCONST(1.0)/(SUN_RCONST(1.0)+t*t) - lambda*atan(t);

  return 0;                                   /* return with success */
}

/* f0 routine to compute a zero-valued ODE RHS function f(t,y). */
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(SUN_RCONST(0.0), ydot);
  return 0;
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype *rdata = (sunrealtype *) user_data;   /* cast user_data to sunrealtype */
  sunrealtype lambda = rdata[0];                  /* set shortcut for stiffness parameter */
  sunrealtype *Jdata = SUNDenseMatrix_Data(J);

  /* Fill in Jacobian of f: set the first entry of the data array to set the (0,0) entry */
  Jdata[0] = lambda;

  return 0;                                   /* return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a value so check if
             retval >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}

/* analytical solution */
static sunrealtype ytrue(sunrealtype t)
{
  return (atan(t));
}

/* check the computed solution */
static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol, sunrealtype atol)
{
  int      passfail=0;     /* answer pass (0) or fail (1) value    */
  sunrealtype ans, err, ewt;  /* answer data, error, and error weight */

  /* compute solution error */
  ans = ytrue(t);
  ewt = SUN_RCONST(1.0) / (rtol * fabs(ans) + atol);
  err = ewt * fabs(NV_Ith_S(y,0) - ans);

  /* is the solution within the tolerances? */
  passfail = (err < SUN_RCONST(1.0)) ? 0 : 1;

  if (passfail) {
    fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%"GSYM"\n\n", err);
  }

  return(passfail);
}

/*---- end of file ----*/
