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
 * The exponential entropy given by
 *
 *   e = exp(u) + exp(v)
 *
 *   e' = [ de/du; de/dv ] = [ exp(u); exp(v) ]
 *
 * is conserved for the analytical solution
 *
 *   u = log(a * exp^(a * t)) - log(b)
 *   v = log(e + e^(3/2)) - log(b)
 *
 * where log is the natural logarithm, a = sqrt(e) + e, and
 * b = sqrt(e) + e^(a * t).
 *
 * ---------------------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>      /* access to ARKStep               */
#include <nvector/nvector_serial.h>     /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h>  /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense SUNLinearSolver */

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

#define EVAL 2.718281828459045235360287471352662497757247093699959574966

/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */

/* ODE RHS function */
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* ODE RHS Jacobian function */
int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Entropy function */
int Ent(N_Vector y, realtype* e, void* user_data);

/* Entropy Jacobian function */
int JacEnt(N_Vector y, N_Vector J, void* user_data);

/* ----------------- *
 * Utility functions *
 * ----------------- */

/* Analytic solution */
int ans(realtype t, N_Vector y);

/* Check return flags and pointers */
int check_flag(int flag, const char* funcname);
int check_ptr(void* ptr, const char* funcname);

/* ------------ *
 * Main Program *
 * ------------ */

int main(int argc, char *argv[])
{
  /* Error-checking flag */
  int flag;

  /* Initial and final times */
  realtype t0 = RCONST(0.0);
  realtype tf = RCONST(5.0);

  /* Relative and absolute tolerances */
  realtype reltol = RCONST(1.0e-6);
  realtype abstol = RCONST(1.0e-10);

  /* SUNDIALS context, vector, matrix, and linear solver objects */
  SUNContext      ctx = NULL;
  N_Vector        y   = NULL;
  SUNMatrix       A   = NULL;
  SUNLinearSolver LS  = NULL;

  /* Pointer to vector data array */
  realtype* ydata;

  /* Initial, current, and change in entropy value */
  realtype ent0, ent, dent;

  /* ARKODE memory structure */
  void* arkode_mem = NULL;

  /* ARKODE statistics */
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;

  /* Number of outputs, output counter, output times, and output file */
  int      nout = 10;
  int      iout;
  realtype t, tout, dtout;
  FILE*    UFID;

  /* Command line options */
  int      implicit = 0;            /* explicit          */
  int      relax    = 1;            /* enable relaxation */
  realtype fixed_h  = RCONST(0.0);  /* adaptive steps    */

  /* -------------------- *
   * Output Problem Setup *
   * -------------------- */

  if (argc > 1) implicit = atoi(argv[1]);
  if (argc > 2) relax    = atoi(argv[2]);
  if (argc > 3) fixed_h  = (realtype) atof(argv[3]);

  printf("\nConserved Exponential Entropy problem:\n");
  if (implicit)
    printf("   method = ERK\n");
  else
    printf("   method = DIRK\n");
  if (relax)
    printf("   relaxation = ON");
  else
    printf("   relaxation = OFF");
  if (fixed_h > 0.0)
    printf("   fixed h = %.6"ESYM"\n", fixed_h);
  if (implicit || fixed_h > 0.0)
  {
    printf("   reltol = %.1"ESYM"\n",  reltol);
    printf("   abstol = %.1"ESYM"\n\n",abstol);
  }

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

  ydata[0] = RCONST(1.0);
  ydata[1] = RCONST(0.5);

  /* Initial entropy */
  flag = Ent(y, &ent0, NULL);
  if (check_flag(flag, "Ent")) return 1;

  /* Initialize the ARKStep */
  if (implicit)
    arkode_mem = ARKStepCreate(NULL, f, t0, y, ctx);
  else
    arkode_mem = ARKStepCreate(f, NULL, t0, y, ctx);
  if (check_ptr(arkode_mem, "ARKStepCreate")) return 1;

  /* Specify tolerances */
  flag = ARKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(flag, "ARKStepSStolerances")) return 1;

  if (relax)
  {
    /* Enable relaxation methods */
    flag = ARKStepSetRelaxFn(arkode_mem, Ent, JacEnt);
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

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt", "w");
  fprintf(UFID,"# t u v e\n");

  /* --------------- *
   * Advance in Time *
   * --------------- */

  /* Initial time, time between outputs, output time */
  t     = t0;
  dtout = tf / nout;
  tout  = t0 + dtout;

  /* Write initial condition to disk */
  fprintf(UFID," %.16"ESYM" %.16"ESYM" %.16"ESYM"\n", t0, ydata[0], ydata[1]);

  printf("        t              u              v              de\n");
  printf("   ---------------------------------------------------------\n");
  printf(" %14.6"ESYM" %14.6"ESYM" %14.6"ESYM" %14.6"ESYM"\n",
         t, ydata[0], ydata[1], RCONST(0.0));

  for (iout = 0; iout < nout; iout++)
  {
    /* Evolve in time */
    flag = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(flag, "ARKStepEvolve")) break;

    /* Output solution and change in entropy */
    flag = Ent(y, &ent, NULL);
    if (check_flag(flag, "Ent")) return 1;

    dent = fabs(ent - ent0);

    printf(" %14.6"ESYM" %14.6"ESYM" %14.6"ESYM" %14.6"ESYM"\n",
           t, ydata[0], ydata[1], dent);
    fprintf(UFID," %.16"ESYM" %.16"ESYM" %.16"ESYM" %.16"ESYM"\n",
            t, ydata[0], ydata[1], dent);

    /* Update output time */
    tout += dtout;
    tout = (tout > tf) ? tf : tout;
  }

  printf("   ---------------------------------------------------------\n");
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
  printf("\n");

  /* -------- *
   * Clean up *
   * -------- */

  /* Free ARKStep integrator and SUNDIALS objects */
  ARKStepFree(&arkode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  SUNContext_Free(&ctx);

  return flag;
}


/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */


/* ODE RHS function f(t,y). */
int f(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  realtype* ydata = N_VGetArrayPointer(y);
  realtype* fdata = N_VGetArrayPointer(ydot);

  fdata[0] = -exp(ydata[1]);
  fdata[1] =  exp(ydata[0]);

  return 0;
}


/* ODE RHS Jacobian function J(t,y) = df/dy. */
int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype* ydata = N_VGetArrayPointer(y);
  realtype* Jdata = SUNDenseMatrix_Data(J);

  /* column 0 */
  Jdata[0] = 0;
  Jdata[1] = -exp(ydata[1]);

  /* column 1 */
  Jdata[2] = exp(ydata[0]);
  Jdata[3] = 0;

  return 0;
}


/* Entropy function e(y) */
int Ent(N_Vector y, realtype* e, void* user_data)
{
  realtype* ydata = N_VGetArrayPointer(y);

  *e = exp(ydata[0]) + exp(ydata[1]);

  return 0;
}


/* Entropy function Jacobian Je(y) = de/dy */
int JacEnt(N_Vector y, N_Vector J, void* user_data)
{
  realtype* ydata = N_VGetArrayPointer(y);
  realtype* jdata = N_VGetArrayPointer(J);

  jdata[0] = exp(ydata[0]);
  jdata[1] = exp(ydata[1]);

  return 0;
}


/* ----------------- *
 * Utility functions *
 * ----------------- */


/* Analytic solution */
int ans(realtype t, N_Vector y)
{
  realtype a, b;
  realtype* ydata = N_VGetArrayPointer(y);

  a = sqrt(EVAL) + EVAL;
  b = sqrt(EVAL) + pow(EVAL, a * t);

  ydata[0] = log(a * exp(a * t)) - log(b);
  ydata[1] = log(EVAL + exp(a * t)) - log(b);

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
