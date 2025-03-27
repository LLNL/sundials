/*-----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
 * The following test employees a manufactured true solution to
 * find the accuracy of numerical solution. The ODE system with
 * 3 components, X = [u,v,w], satisfies the equations,
 *
 * du/dt = (t - 1)*(t - w)*i
 * dv/dt = 2*u*i - v + 4
 * dw/dt = v + t^2*exp(-t)*i - exp(-t)*i + 1
 *
 * for t in the interval [0.0, 5.0], with initial conditions
 * Y0 = [u0,v0,w0], computed substituting t=0 in the true solution,
 *
 * u = -t*exp(-t) + 2*i
 * v = -t^2*exp(-t)*i
 * w = exp(-t)*i + t
 *
 * This program solves the problem with the IMEX method, using a
 * Newton iteration with the SPTFQMR linear solver.
 *
 * 5 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <arkode/arkode_arkstep.h> /* prototypes for ARKStep fcts., consts */
#include <math.h>
#include <stdio.h>
#include <sundials/sundials_types.h>   /* def. of type 'sunrealtype' */
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to sptfqmr SUNLinearSolver       */

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
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

// Compute the true solution
static int Solution(sunrealtype t, N_Vector u, void* user_data);

// Compute the numerical solution error
static int SolutionError(sunrealtype t, N_Vector u, N_Vector e, void* user_data);

// Preconditioner setup and solve functions (Identity in this testing)
static int PSetup(sunrealtype t, N_Vector y, N_Vector f, sunbooleantype jok,
                  sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data);

static int PSolve(sunrealtype t, N_Vector y, N_Vector f, N_Vector r, N_Vector z,
                  sunrealtype gamma, sunrealtype delta, int lr, void* user_data);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* Main Program */
int main(void)
{
  /* general problem parameters */
  sunrealtype T0     = SUN_RCONST(0.0);       /* initial time */
  sunrealtype Tf     = SUN_RCONST(5.0);       /* final time */
  sunrealtype dTout  = SUN_RCONST(1.0);       /* time between outputs */
  sunindextype NEQ   = 3;                     /* number of dependent vars. */
  int Nt             = (int)ceil(Tf / dTout); /* number of output times */
  sunrealtype reltol = 1.0e-6;                /* tolerances */
  sunrealtype abstol = 1.0e-10;
  sunscalartype a, b, c; // some constants to pass as problem data if needed
  int pretype = 1; // 1 for left and 2 for right preconditioners
  int maxl = 100; // Set the maximum number of linear solver iterations

  /* general problem variables */
  int flag;                  /* reusable error-checking flag */
  N_Vector y         = NULL; /* empty vector for storing solution */
  N_Vector True_Sol  = NULL; // vector for storing true solution
  N_Vector Error     = NULL; // vector for storing the error */
  SUNLinearSolver LS = NULL; /* empty linear solver object */
  void* arkode_mem   = NULL; /* empty ARKode memory structure */
  sunscalartype rdata[3];
  FILE* UFID;
  sunrealtype t, tout;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, nnf, ncfn, netf;

  /* Create the SUNDIALS context object for this testing */
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  y = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }

  True_Sol = N_VClone(y);
  Error = N_VClone(y);

  // set up the problem data (unused in this case)
  rdata[0] = a = 1.0;
  rdata[1] = b = 1.0;
  rdata[2] = c = 1.0;

  // Set initial condition
  flag = Solution(0.0, True_Sol, (void*)rdata);
  if (check_flag(&flag, "Solution", 1)) { return 1; }

  NV_Ith_S(y, 0) = NV_Ith_S(True_Sol, 0); /* Set initial conditions */
  NV_Ith_S(y, 1) = NV_Ith_S(True_Sol, 1);
  NV_Ith_S(y, 2) = NV_Ith_S(True_Sol, 2);

  /* Initial problem output */
  printf("\nAnalytic ODE test problem:\n");
  printf("    initial conditions:  u0 = %10.5"FSYM " + " "%.5" FSYM "i  | v0 = " "%10.5"FSYM " + " "%.5" FSYM "i  | w0 = " "%10.5"FSYM " + " "%.5" FSYM "i  \n",
  creal(NV_Ith_S(True_Sol, 0)), cimag(NV_Ith_S(True_Sol, 0)),
  creal(NV_Ith_S(True_Sol, 1)), cimag(NV_Ith_S(True_Sol, 1)),
  creal(NV_Ith_S(True_Sol, 2)), cimag(NV_Ith_S(True_Sol, 2)));

  printf("    problem parameters:  a = %" GSYM " + %" GSYM "i,  b = %" GSYM " + %" GSYM "i,  c = %" GSYM
         " + %" GSYM "i\n",
  creal(a), cimag(a),
  creal(b), cimag(b),
  creal(c), cimag(c));
  printf("    reltol = %.1" ESYM ",  abstol = %.1" ESYM "\n\n", reltol, abstol);

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y, ctx);
  if (check_flag((void*)arkode_mem, "ARKStepCreate", 0)) { return 1; }

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem, (void*)rdata); /* Pass rdata to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  flag = ARKodeSStolerances(arkode_mem, reltol, abstol); /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  flag = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE); /* Specify stiff interpolant */
  if (check_flag(&flag, "ARKodeSetInterpolantType", 1)) { return 1; }

  flag = ARKodeSetDeduceImplicitRhs(arkode_mem, 1); /* Avoid eval of f after stage */
  if (check_flag(&flag, "ARKodeSetDeduceImplicitRhs", 1)) { return 1; }

  LS = SUNLinSol_SPTFQMR(y, pretype, maxl, ctx);
  if (check_flag((void*)LS, "SUNLinSol_SPTFQMR", 0)) { return 1; }

  /* Linear solver interface */
  flag = ARKodeSetLinearSolver(arkode_mem, LS, NULL);
  if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }

  // Attach preconditioner
  flag = ARKodeSetPreconditioner(arkode_mem, PSetup, PSolve);
  if (check_flag(&flag, "ARKodeSetPreconditioner", 1)) { return 1; }

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt", "w");
  fprintf(UFID, "#     t                  u                         v                          w\n");

  /* output initial condition to disk */
  fprintf(UFID," %10.3"FSYM" | " "%10.5"FSYM " + " "%.5" FSYM "i  |  " "%10.5"FSYM " + " "%.5" FSYM "i  |  " "%10.5"FSYM " + " "%.5" FSYM "i  \n", T0,
  creal(NV_Ith_S(y, 0)), cimag(NV_Ith_S(y, 0)),
  creal(NV_Ith_S(y, 1)), cimag(NV_Ith_S(y, 1)),
  creal(NV_Ith_S(y, 2)), cimag(NV_Ith_S(y, 2)));

  /* Main time-stepping loop: calls ARKodeEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t    = T0;
  tout = T0 + dTout;
  printf("      t                  u                         v                          w\n");
  printf("   --------------------------------------------------------------------------------------\n");
  printf(" %10.3"FSYM" | " "%10.5"FSYM " + " "%.5" FSYM "i  |  " "%10.5"FSYM " + " "%.5" FSYM "i  |  " "%10.5"FSYM " + " "%.5" FSYM "i  \n", t,
  creal(NV_Ith_S(y, 0)), cimag(NV_Ith_S(y, 0)),
  creal(NV_Ith_S(y, 1)), cimag(NV_Ith_S(y, 1)),
  creal(NV_Ith_S(y, 2)), cimag(NV_Ith_S(y, 2)));

  for (iout = 0; iout < Nt; iout++)
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL); /* call integrator */
    if (check_flag(&flag, "ARKodeEvolve", 1)) { break; }
    printf(" %10.3"FSYM" | " "%10.5"FSYM " + " "%.5" FSYM "i  |  " "%10.5"FSYM " + " "%.5" FSYM "i  |  " "%10.5"FSYM " + " "%.5" FSYM "i  \n", t,
    creal(NV_Ith_S(y, 0)), cimag(NV_Ith_S(y, 0)),
    creal(NV_Ith_S(y, 1)), cimag(NV_Ith_S(y, 1)),
    creal(NV_Ith_S(y, 2)), cimag(NV_Ith_S(y, 2)));

    fprintf(UFID, " %10.3"FSYM" | " "%10.5"FSYM " + " "%.5" FSYM "i  |  " "%10.5"FSYM " + " "%.5" FSYM "i  |  " "%10.5"FSYM " + " "%.5" FSYM "i  \n", t,
    creal(NV_Ith_S(y, 0)), cimag(NV_Ith_S(y, 0)),
    creal(NV_Ith_S(y, 1)), cimag(NV_Ith_S(y, 1)),
    creal(NV_Ith_S(y, 2)), cimag(NV_Ith_S(y, 2)));

    if (flag >= 0)
    { /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    else
    { /* unsuccessful solve: break */
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   --------------------------------------------------------------------------------------\n");
  fclose(UFID);

  SolutionError(Tf, y, Error, (void*)rdata);

  printf("   --------------------------------------------------------------------------------------\n");
  /* Print some final statistics */
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumStepAttempts", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1);
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  flag = ARKodeGetNumStepSolveFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKodeGetNumStepSolveFails", 1);
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1);
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &nnf);
  check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1);
  flag = ARKodeGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKodeGetNumJacEvals", 1);
  flag = ARKodeGetNumLinRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKodeGetNumLinRhsEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", nnf);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Total number of failed steps from solver failure = %li\n", ncfn);

  /* Clean up and return with successful completion */
  N_VDestroy(y);            /* Free y vector */
  ARKodeFree(&arkode_mem); /* Free integrator memory */
  SUNLinSolFree(LS);        /* Free linear solver */
  SUNContext_Free(&ctx);    /* Free context */

  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunscalartype* rdata = (sunscalartype*)user_data; /* cast user_data to sunscalartype */
  sunscalartype a = rdata[0]; /* access data entries (unused) */
  sunscalartype b = rdata[1];
  sunscalartype c = rdata[2];
  sunscalartype u  = NV_Ith_S(y, 0); /* access solution values */
  sunscalartype v  = NV_Ith_S(y, 1);
  sunscalartype w  = NV_Ith_S(y, 2);

  /* fill in the RHS function */
  NV_Ith_S(ydot, 0) = (t - 1.0)*(t - w)*I;
  NV_Ith_S(ydot, 1) = 2.0*u*I - v + 4.0;
  NV_Ith_S(ydot, 2) = v + t*t*exp(-t)*I -exp(-t)*I + 1.0;

  return 0; /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

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

// Preconditioner setup routine
static int PSetup(sunrealtype t, N_Vector y, N_Vector f, sunbooleantype jok,
                  sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)
{
  return 0;
}

// Preconditioner solve routine for Pz = r
static int PSolve(sunrealtype t, N_Vector y, N_Vector f, N_Vector r, N_Vector z,
                  sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{
  N_Vector ones;
  ones = N_VClone(r);
  N_VConst(1.0, ones);

  N_VProd(ones, r, z);

  /* return with success */
  return 0;
}


// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the exact solution
static int Solution(sunrealtype t, N_Vector u, void* user_data)
{
  sunscalartype* rdata = (sunscalartype*)user_data; /* cast user_data to sunscalartype */
  sunscalartype a = rdata[0];                    /* access data entries */
  sunscalartype b = rdata[1];
  sunscalartype c = rdata[2];
  // Initialize u to zero (handles boundary conditions)
  N_VConst(0.0, u);

  sunscalartype* uarray = N_VGetArrayPointer(u);
  if (check_flag((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

  uarray[0] = -t*exp(-t) + 2.0*I;
  uarray[1] = -t*t*exp(-t)*I;
  uarray[2] = exp(-t)*I + t;

  return 0;
}

// Compute the solution error
static int SolutionError(sunrealtype t, N_Vector u, N_Vector e, void* user_data)
{
  // Compute true solution
  int flag = Solution(t, e, (void*)user_data);
  if (flag != 0) { return -1; }

  // Compute absolute error
  N_VLinearSum(1.0, u, -1.0, e, e);

  sunrealtype error_norm = SUNSQR(N_VDotProd(e, e));

  printf("     Norm of the error is %10.15"ESYM"\n", error_norm);

  return 0;
}

/*---- end of file ----*/