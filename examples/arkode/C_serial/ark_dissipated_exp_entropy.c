/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example problem is adapted from:
 *
 * H. Ranocha, M. Sayyari, L. Dalcin, M. Parsani, and D.I. Ketcheson,
 * "Relaxation Runge-Kutta Methods: Fully-Discrete Explicit Entropy-Stable
 * Schemes for the Compressible Euler and Navier-Stokes Equations," SIAM Journal
 * on Scientific Computing, 42(2), 2020, https://doi.org/10.1137/19M1263480.
 * -----------------------------------------------------------------------------
 * This example evolves the equation du/dt = -exp(u) for t in the interval
 * [0, 5] with the initial condition u(0) = 0.5. The equation has the analytic
 * solution u(t) = -log(e^{-0.5} + t) and the dissipated exponential entropy is
 * given by ent(u) = exp(u) with Jacobian ent'(u) = de/du = exp(u).
 *
 * The problem is advanced in time with an explicit or implicit relaxed
 * Runge-Kutta method to ensure dissipation of the entropy.
 * ---------------------------------------------------------------------------*/

/* Header files */
#include <math.h>
#include <stdio.h>

#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

/* Precision-dependent output macros */
#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */

/* ODE RHS function */
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* ODE RHS Jacobian function */
int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Entropy function */
int Ent(N_Vector* y, sunrealtype* e, void* user_data);

/* Entropy Jacobian function */
int JacEnt(N_Vector* y, N_Vector* J, void* user_data);

/* ----------------- *
 * Utility functions *
 * ----------------- */

/* Analytic solution */
int ans(sunrealtype t, N_Vector y);

/* Check return flags and pointers */
int check_flag(int flag, const char* funcname);
int check_ptr(void* ptr, const char* funcname);

/* ------------ *
 * Main Program *
 * ------------ */

int main(int argc, char* argv[])
{
  /* Error-checking flag */
  int flag;

  /* Initial and final times */
  sunrealtype t0 = SUN_RCONST(0.0);
  sunrealtype tf = SUN_RCONST(5.0);

  /* Relative and absolute tolerances */
  sunrealtype reltol = SUN_RCONST(1.0e-6);
  sunrealtype abstol = SUN_RCONST(1.0e-10);

  /* SUNDIALS context, vector, matrix, and linear solver objects */
  SUNContext ctx       = NULL;
  N_Vector y           = NULL;
  N_Vector ytrue       = NULL;
  SUNMatrix A          = NULL;
  SUNLinearSolver LS   = NULL;
  ARKodeButcherTable B = NULL;

  /* Pointer to vector data array */
  sunrealtype* ydata;
  sunrealtype* ytdata;

  /* Initial, current, and change in entropy value */
  sunrealtype ent0, ent, delta_ent;

  /* Solution errors */
  sunrealtype u_err;

  /* ARKODE memory structure */
  void* arkode_mem = NULL;

  /* ARKODE statistics */
  long int nst, nst_a, nfe, nfi, nri, nrf;
  long int nsetups, nje, nfeLS, nni, ncfn, netf;

  /* Number of outputs, output counter, output times, and output file */
  int nout = 10;
  int iout;
  sunrealtype t, tout, dtout;
  FILE* UFID;

  /* Command line options */
  int implicit        = 0;               /* explicit          */
  int relax           = 0;               /* enable relaxation */
  int steps           = 100;             /* number of steps   */
  sunrealtype fixed_h = SUN_RCONST(0.0); /* adaptive steps    */

  /* -------------------- *
   * Output Problem Setup *
   * -------------------- */

  if (argc > 1) implicit = atoi(argv[1]);
  if (argc > 2) relax = atoi(argv[2]);
  if (argc > 3) fixed_h = (sunrealtype)atof(argv[3]);
  if (argc > 4) steps = atoi(argv[4]);

  printf("\nDissipated Exponential Entropy problem:\n");
  if (implicit) printf("   method     = DIRK\n");
  else printf("   method     = ERK\n");
  if (relax) printf("   relaxation = ON\n");
  else printf("   relaxation = OFF\n");
  if (fixed_h > SUN_RCONST(0.0)) printf("   fixed h    = %" GSYM "\n", fixed_h);
  if (implicit || !(fixed_h > SUN_RCONST(0.0)))
  {
    printf("   reltol     = %.1" ESYM "\n", reltol);
    printf("   abstol     = %.1" ESYM "\n", abstol);
  }
  printf("   steps      = %d\n", steps);
  printf("\n");

  /* ------------ *
   * Setup ARKODE *
   * ------------ */

  /* Create the SUNDIALS context object for this simulation */
  flag = SUNContext_Create(NULL, &ctx);
  if (check_flag(flag, "SUNContext_Create")) return 1;

  /* Create serial vector and set the initial condition values */
  y = N_VNew_Serial(2, ctx);
  if (check_ptr(y, "N_VNew_Serial")) return 1;

  ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) return 1;

  ydata[0] = SUN_RCONST(1.0);
  ydata[1] = SUN_RCONST(0.5);

  ytrue = N_VClone(y);
  if (check_ptr(ytrue, "N_VClone")) return 1;

  ytdata = N_VGetArrayPointer(ytrue);
  if (check_ptr(ytdata, "N_VGetArrayPointer")) return 1;

  /* Initialize the ARKStep */
  if (implicit) arkode_mem = ARKStepCreate(NULL, f, t0, y, ctx);
  else arkode_mem = ARKStepCreate(f, NULL, t0, y, ctx);
  if (check_ptr(arkode_mem, "ARKStepCreate")) return 1;

  /* Specify tolerances */
  flag = ARKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(flag, "ARKStepSStolerances")) return 1;

  if (relax)
  {
    /* Enable relaxation methods */
    flag = ARKStepSetRelaxFn(arkode_mem, 1, Ent, JacEnt);
    if (check_flag(flag, "ARKStepSetRelaxFn")) return 1;
  }

  if (fixed_h > 0.0)
  {
    /* Set the step size */
    flag = ARKStepSetFixedStep(arkode_mem, fixed_h);
    if (check_flag(flag, "ARKStepSetFixedStep")) return 1;
  }

  if (implicit)
  {
    /* Create dense matrix and linear solver */
    A = SUNDenseMatrix(2, 2, ctx);
    if (check_ptr(A, "SUNDenseMatrix")) return 1;

    LS = SUNLinSol_Dense(y, A, ctx);
    if (check_ptr(LS, "SUNLinSol_Dense")) return 1;

    /* Attach the matrix and linear solver */
    flag = ARKStepSetLinearSolver(arkode_mem, LS, A);
    if (check_flag(flag, "ARKStepSetLinearSolver")) return 1;

    /* Set Jacobian routine */
    flag = ARKStepSetJacFn(arkode_mem, Jac);
    if (check_flag(flag, "ARKStepSetJacFn")) return 1;
  }
  else if (fixed_h > 0.0)
  {
    /* Use RK4 */
    B = ARKodeButcherTable_Alloc(4, SUNFALSE);

    B->A[1][0] = SUN_RCONST(0.5);
    B->A[2][1] = SUN_RCONST(0.5);
    B->A[3][2] = SUN_RCONST(1.0);

    B->c[1] = SUN_RCONST(0.5);
    B->c[2] = SUN_RCONST(0.5);
    B->c[3] = SUN_RCONST(1.0);

    B->b[0] = SUN_RCONST(1.0) / SUN_RCONST(6.0);
    B->b[1] = SUN_RCONST(1.0) / SUN_RCONST(3.0);
    B->b[2] = SUN_RCONST(1.0) / SUN_RCONST(3.0);
    B->b[3] = SUN_RCONST(1.0) / SUN_RCONST(6.0);

    B->q = 4;
    B->p = 0;

    flag = ARKStepSetTables(arkode_mem, 4, 0, NULL, B);
    if (check_flag(flag, "ARKStepSetTables")) return 1;
  }

  /* Open output stream for results, output comment line */
  UFID = fopen("ark_dissipated_exp_entropy.txt", "w");
  fprintf(UFID, "# vars: t u v entropy u_err v_err entropy_error\n");

  /* --------------- *
   * Advance in Time *
   * --------------- */

  /* Initial time, time between outputs, output time */
  t     = t0;
  dtout = tf / nout;
  tout  = t0 + dtout;

  /* Output the initial condition and entropy */
  flag = Ent(&y, &ent0, NULL);
  if (check_flag(flag, "Ent")) return 1;

  fprintf(UFID,
          "%23.16" ESYM " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM "\n",
          t0, ydata[0], ent0, SUN_RCONST(0.0), SUN_RCONST(0.0));

  printf("        t              u               e                  delta e\n");
  printf("   "
         "---------------------------------------------------------------------"
         "---\n");
  printf(" %14.6" ESYM " %14.6" ESYM " %14.6" ESYM " %14.6" ESYM " %14.6" ESYM "\n",
         t, ydata[0], ent0, SUN_RCONST(0.0), SUN_RCONST(0.0));

  for (iout = 0; iout < steps; iout++)
  {
    /* Evolve in time */
    flag = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_ONE_STEP);
    if (check_flag(flag, "ARKStepEvolve")) break;

    /* Output solution and errors */
    flag = Ent(&y, &ent, NULL);
    if (check_flag(flag, "Ent")) return 1;

    flag = ans(t, ytrue);
    if (check_flag(flag, "ans")) return 1;

    delta_ent = ent - ent0;
    u_err     = ydata[0] - ytdata[0];

    printf(" %14.6" ESYM " %14.6" ESYM " %14.6" ESYM " %14.16" ESYM " %14.6" ESYM "\n",
           t, ydata[0], ent, u_err, delta_ent);

    fprintf(UFID,
            "%23.16" ESYM " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM "\n",
            t, ydata[0], ent, u_err, delta_ent);

    /* Update output time */
    tout += dtout;
    tout = (tout > tf) ? tf : tout;
  }

  printf("   "
         "---------------------------------------------------------------------"
         "---\n");
  fclose(UFID);

  /* ------------ *
   * Output Stats *
   * ------------ */

  /* Get final statistics on how the solve progressed */
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  check_flag(flag, "ARKStepGetNumSteps");

  flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(flag, "ARKStepGetNumStepAttempts");

  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  check_flag(flag, "ARKStepGetNumErrTestFails");

  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(flag, "ARKStepGetNumRhsEvals");

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);

  if (implicit)
  {
    flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
    check_flag(flag, "ARKStepGetNumNonlinSolvIters");

    flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
    check_flag(flag, "ARKStepGetNumNonlinSolvConvFails");

    flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
    check_flag(flag, "ARKStepGetNumLinSolvSetups");

    flag = ARKStepGetNumJacEvals(arkode_mem, &nje);
    check_flag(flag, "ARKStepGetNumJacEvals");

    flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfeLS);
    check_flag(flag, "ARKStepGetNumLinRhsEvals");

    printf("   Total number of Newton iterations = %li\n", nni);
    printf("   Total number of linear solver convergence failures = %li\n", ncfn);
    printf("   Total linear solver setups = %li\n", nsetups);
    printf("   Total number of Jacobian evaluations = %li\n", nje);
    printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  }

  if (relax)
  {
    flag = ARKStepGetNumRelaxSolveIters(arkode_mem, &nri);
    check_flag(flag, "ARKStepGetNumRelaxSolveIters");

    flag = ARKStepGetNumRelaxSolveFails(arkode_mem, &nrf);
    check_flag(flag, "ARKStepGetNumRelaxSolveFails");

    printf("   Total Relaxation iters = %li\n", nri);
    printf("   Total Relaxation fails = %li\n", nrf);
  }
  printf("\n");

  /* -------- *
   * Clean up *
   * -------- */

  /* Free ARKStep integrator and SUNDIALS objects */
  ARKStepFree(&arkode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  N_VDestroy(ytrue);
  SUNContext_Free(&ctx);

  return flag;
}

/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */

/* ODE RHS function f(t,y). */
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);
  fdata[0] = -exp(ydata[0]);
  return 0;
}

/* ODE RHS Jacobian function J(t,y) = df/dy. */
int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);
  Jdata[0] = -exp(ydata[0]);
  return 0;
}

/* Entropy function e(y) */
int Ent(N_Vector* y, sunrealtype* e, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y[0]);
  e[0] = exp(ydata[0]);
  return 0;
}

/* Entropy function Jacobian Je(y) = de/dy */
int JacEnt(N_Vector* y, N_Vector* J, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y[0]);
  sunrealtype* jdata = N_VGetArrayPointer(J[0]);
  jdata[0] = exp(ydata[0]);
  return 0;
}

/* ----------------- *
 * Utility functions *
 * ----------------- */

/* Analytic solution */
int ans(sunrealtype t, N_Vector y)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0] = log(exp(SUN_RCONST(-0.5)) + t);
  return 0;
}

/* Check return flags */
int check_flag(int flag, const char* funcname)
{
  if (flag < 0)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
            funcname, flag);
    return 1;
  }
  return 0;
}

/* Check return pointers */
int check_ptr(void* ptr, const char* funcname)
{
  if (!ptr)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }
  return 0;
}
