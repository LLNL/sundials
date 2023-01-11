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
 * This example evolves system
 *
 *   du/dt = -exp(v)
 *   dv/dt =  exp(u)
 *
 * for t in the interval [0, 5] with the initial condition
 *
 *   u(0) = 1.0
 *   v(0) = 0.5
 *
 * The system has the analytic solution
 *
 *   u = log(e + e^(3/2)) - log(b)
 *   v = log(a * e^(a * t)) - log(b)
 *
 * where log is the natural logarithm, a = sqrt(e) + e, and
 * b = sqrt(e) + e^(a * t).
 *
 * The conserved exponential entropy for the system is given by
 * ent(u,v) = exp(u) + exp(v) with the Jacobian
 * ent'(u,v) = [ de/du de/dv ]^T = [ exp(u) exp(v) ]^T.
 *
 * The problem is advanced in time with an explicit relaxed
 * Runge-Kutta method to ensure conservation of the entropy.
 * ---------------------------------------------------------------------------*/

/* Header files */
#include <math.h>
#include <stdio.h>

#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

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

/* Value of the natural number e */
#define EVAL 2.718281828459045235360287471352662497757247093699959574966

/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */

/* ODE RHS function */
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

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
  ARKodeButcherTable B = NULL;

  /* Pointer to vector data array */
  sunrealtype* ydata;
  sunrealtype* ytdata;

  /* Initial, current, and change in entropy value */
  sunrealtype ent0, ent, ent_err;

  /* Solution errors */
  sunrealtype u_err, v_err;

  /* ARKODE memory structure */
  void* arkode_mem = NULL;

  /* ARKODE statistics */
  long int nst, nst_a, nfe;
  long int nre, nrje, nri, nrf, netf;

  /* Number of outputs, output counter, output times, and output file */
  int nout = 10;
  int iout;
  sunrealtype t, tout, dtout;
  FILE* UFID;

  /* Command line options */
  int relax           = 1;               /* enable relaxation */
  int steps           = 100;             /* number of steps   */
  sunrealtype fixed_h = SUN_RCONST(0.0); /* adaptive steps    */

  /* -------------------- *
   * Output Problem Setup *
   * -------------------- */

  if (argc > 1) relax = atoi(argv[1]);
  if (argc > 2) fixed_h = (sunrealtype)atof(argv[2]);
  if (argc > 3) steps = atoi(argv[3]);

  printf("\nConserved Exponential Entropy problem:\n");
  printf("   method     = ERK\n");
  if (relax) printf("   relaxation = ON\n");
  else printf("   relaxation = OFF\n");
  if (fixed_h > SUN_RCONST(0.0)) printf("   fixed h    = %" GSYM "\n", fixed_h);
  if (!(fixed_h > SUN_RCONST(0.0)))
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

  /* Initial entropy */
  flag = Ent(&y, &ent0, NULL);
  if (check_flag(flag, "Ent")) return 1;

  /* Initialize the ERKStep */
  arkode_mem = ERKStepCreate(f, t0, y, ctx);
  if (check_ptr(arkode_mem, "ERKStepCreate")) return 1;

  /* Specify tolerances */
  flag = ERKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(flag, "ERKStepSStolerances")) return 1;

  if (relax)
  {
    /* Enable relaxation methods */
    flag = ERKStepSetRelaxFn(arkode_mem, 1, Ent, JacEnt);
    if (check_flag(flag, "ERKStepSetRelaxFn")) return 1;
  }

  if (fixed_h > 0.0)
  {
    /* Set the step size */
    flag = ERKStepSetFixedStep(arkode_mem, fixed_h);
    if (check_flag(flag, "ERKStepSetFixedStep")) return 1;

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

    flag = ERKStepSetTable(arkode_mem, B);
    if (check_flag(flag, "ERKStepSetTables")) return 1;
  }

  /* Open output stream for results, output comment line */
  UFID = fopen("ark_conserved_exp_entropy_erk.txt", "w");
  fprintf(UFID, "# vars: t u v entropy u_err v_err entropy_error\n");

  /* --------------- *
   * Advance in Time *
   * --------------- */

  /* Initial time, time between outputs, output time */
  t     = t0;
  dtout = tf / nout;
  tout  = t0 + dtout;

  /* Output the initial condition and entropy */
  fprintf(UFID,
          "%23.16" ESYM " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM
          " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM "\n",
          t0, ydata[0], ydata[1], ent0, SUN_RCONST(0.0), SUN_RCONST(0.0),
          SUN_RCONST(0.0));

  printf("        t              u              v              e            "
         "delta e\n");
  printf("   "
         "---------------------------------------------------------------------"
         "---\n");
  printf(" %14.6" ESYM " %14.6" ESYM " %14.6" ESYM " %14.6" ESYM " %14.6" ESYM
         "\n",
         t, ydata[0], ydata[1], ent0, SUN_RCONST(0.0));

  for (iout = 0; iout < steps; iout++)
  {
    /* Evolve in time */
    flag = ERKStepEvolve(arkode_mem, tout, y, &t, ARK_ONE_STEP);
    if (check_flag(flag, "ERKStepEvolve")) break;

    /* Output solution and errors */
    flag = Ent(&y, &ent, NULL);
    if (check_flag(flag, "Ent")) return 1;

    flag = ans(t, ytrue);
    if (check_flag(flag, "ans")) return 1;

    ent_err = ent - ent0;
    u_err   = ydata[0] - ytdata[0];
    v_err   = ydata[1] - ytdata[1];

    printf(" %14.6" ESYM " %14.6" ESYM " %14.6" ESYM " %14.6" ESYM " %14.6" ESYM
           "\n",
           t, ydata[0], ydata[1], ent, ent_err);

    fprintf(UFID,
            "%23.16" ESYM " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM
            " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM "\n",
            t, ydata[0], ydata[1], ent, u_err, v_err, ent_err);

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
  flag = ERKStepGetNumSteps(arkode_mem, &nst);
  check_flag(flag, "ERKStepGetNumSteps");

  flag = ERKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(flag, "ERKStepGetNumStepAttempts");

  flag = ERKStepGetNumErrTestFails(arkode_mem, &netf);
  check_flag(flag, "ERKStepGetNumErrTestFails");

  flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
  check_flag(flag, "ERKStepGetNumRhsEvals");

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Total RHS evals = %li\n", nfe);

  if (relax)
  {
    flag = ERKStepGetNumRelaxFnEvals(arkode_mem, &nre);
    check_flag(flag, "ERKStepGetNumRelaxFnEvals");

    flag = ERKStepGetNumRelaxJacEvals(arkode_mem, &nrje);
    check_flag(flag, "ERKStepGetNumRelaxJacEvals");

    flag = ERKStepGetNumRelaxSolveIters(arkode_mem, &nri);
    check_flag(flag, "ERKStepGetNumRelaxSolveIters");

    flag = ERKStepGetNumRelaxSolveFails(arkode_mem, &nrf);
    check_flag(flag, "ERKStepGetNumRelaxSolveFails");

    printf("   Total Relaxation Fn evals  = %li\n", nre);
    printf("   Total Relaxation Jac evals = %li\n", nrje);
    printf("   Total Relaxation iters     = %li\n", nri);
    printf("   Total Relaxation fails     = %li\n", nrf);
  }
  printf("\n");

  /* -------- *
   * Clean up *
   * -------- */

  /* Free ERKStep integrator and SUNDIALS objects */
  ERKStepFree(&arkode_mem);
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

  fdata[0] = -exp(ydata[1]);
  fdata[1] = exp(ydata[0]);

  return 0;
}

/* Entropy function e(y) */
int Ent(N_Vector* y, sunrealtype* e, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y[0]);

  e[0] = exp(ydata[0]) + exp(ydata[1]);

  return 0;
}

/* Entropy function Jacobian Je(y) = de/dy */
int JacEnt(N_Vector* y, N_Vector* J, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y[0]);
  sunrealtype* jdata = N_VGetArrayPointer(J[0]);

  jdata[0] = exp(ydata[0]);
  jdata[1] = exp(ydata[1]);

  return 0;
}

/* ----------------- *
 * Utility functions *
 * ----------------- */

/* Analytic solution */
int ans(sunrealtype t, N_Vector y)
{
  sunrealtype a, b;
  sunrealtype* ydata = N_VGetArrayPointer(y);

  a = sqrt(EVAL) + EVAL;
  b = sqrt(EVAL) + exp(a * t);

  ydata[0] = log(EVAL + exp(SUN_RCONST(1.5))) - log(b);
  ydata[1] = log(a * exp(a * t)) - log(b);

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
