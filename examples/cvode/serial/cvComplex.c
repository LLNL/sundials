/*-----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
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
 * This program solves the problem with the BDF method, using a
 * Newton iteration with the SPTFQMR linear solver.
 *
 * 5 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <cvode/cvode.h> /* prototypes for CVODE fcts., consts. */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <sundials/sundials_types.h>     /* def. of type sunscalartype */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to sptfqmr SUNLinearSolver */

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
  sunrealtype reltol = SUN_RCONST(1.0e-6);    /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-10);
  int maxl = 10;                              /* max linear solver iterations */

  /* general problem variables */
  int flag;                  /* reusable error-checking flag */
  N_Vector y         = NULL; /* empty vector for storing solution */
  N_Vector True_Sol  = NULL; // vector for storing true solution
  N_Vector Error     = NULL; // vector for storing the error */
  SUNLinearSolver LS = NULL; /* empty linear solver object */
  void* cvode_mem    = NULL; /* empty CVode memory structure */
  sunscalartype rdata[3];
  sunrealtype t, tout;
  int iout;

  /* Create the SUNDIALS context object for this testing */
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  y = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  sunscalartype* yvals = N_VGetArrayPointer(y);
  if (check_flag(yvals, "N_VGetArrayPointer", 0)) { return 1; }

  True_Sol = N_VClone(y);
  Error    = N_VClone(y);

  // set up the problem data (unused in this case)
  rdata[0] = SUN_CCONST(1.0, 0.0);
  rdata[1] = SUN_CCONST(1.0, 0.0);
  rdata[2] = SUN_CCONST(1.0, 0.0);

  // Set initial condition
  flag = Solution(0.0, True_Sol, (void*)rdata);
  if (check_flag(&flag, "Solution", 1)) { return 1; }
  N_VScale(1.0, True_Sol, y); /* Set initial conditions */

  /* Initial problem output */
  printf("\nAnalytic ODE test problem:\n");
  printf("    problem parameters:  a = %" GSYM " + %" GSYM "i,  b = %" GSYM
         " + %" GSYM "i,  c = %" GSYM " + %" GSYM "i\n",
         SUN_REAL(rdata[0]), SUN_IMAG(rdata[0]), SUN_REAL(rdata[1]),
         SUN_IMAG(rdata[2]), SUN_REAL(rdata[2]), SUN_IMAG(rdata[2]));
  printf("    reltol = %.1" ESYM ",  abstol = %.1" ESYM "\n\n", reltol, abstol);

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, ctx);
  if (check_flag((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) { return (1); }

  /* Set routines */
  flag = CVodeSetUserData(cvode_mem,
                          (void*)rdata); /* Pass rdata to user functions */
  if (check_flag(&flag, "CVodeSetUserData", 1)) { return 1; }

  flag = CVodeSStolerances(cvode_mem, reltol, abstol); /* Specify tolerances */
  if (check_flag(&flag, "CVodeSStolerances", 1)) { return 1; }

  LS = SUNLinSol_SPTFQMR(y, 0, maxl, ctx);
  if (check_flag((void*)LS, "SUNLinSol_SPTFQMR", 0)) { return 1; }

  /* Linear solver interface */
  flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if (check_flag(&flag, "CVodeSetLinearSolver", 1)) { return 1; }

  /* Main time-stepping loop: calls CVode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t    = T0;
  tout = T0 + dTout;
  printf(
    "     t               u                      v                      w\n");
  printf("   "
         "---------------------------------------------------------------------"
         "-------\n");
  printf(" %8.3" FSYM " | "
         "%8.5" FSYM " + "
         "%8.5" FSYM "i | "
         "%8.5" FSYM " + "
         "%8.5" FSYM "i | "
         "%8.5" FSYM " + "
         "%8.5" FSYM "i\n",
         t, SUN_REAL(yvals[0]), SUN_IMAG(yvals[0]), SUN_REAL(yvals[1]),
         SUN_IMAG(yvals[1]), SUN_REAL(yvals[2]), SUN_IMAG(yvals[2]));

  for (iout = 0; iout < Nt; iout++)
  {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL); /* call integrator */
    if (check_flag(&flag, "CVode", 1)) { break; }
    printf(" %8.3" FSYM " | "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i | "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i | "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i\n",
           t, SUN_REAL(yvals[0]), SUN_IMAG(yvals[0]), SUN_REAL(yvals[1]),
           SUN_IMAG(yvals[1]), SUN_REAL(yvals[2]), SUN_IMAG(yvals[2]));

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
  printf("   "
         "---------------------------------------------------------------------"
         "-------\n");
  SolutionError(Tf, y, Error, (void*)rdata);
  printf("   "
         "---------------------------------------------------------------------"
         "-------\n");

  /* Print all final statistics */
  printf("\nFinal Solver Statistics:\n");
  flag = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "CVodePrintAllStats", 1)) { return 1; }

  /* Clean up and return with successful completion */
  N_VDestroy(y); /* Free vectors */
  N_VDestroy(True_Sol);
  N_VDestroy(Error);
  CVodeFree(&cvode_mem); /* Free integrator memory */
  SUNLinSolFree(LS);     /* Free linear solver */
  SUNContext_Free(&ctx); /* Free context */

  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunscalartype* rdata =
    (sunscalartype*)user_data; /* cast user_data to sunscalartype */
  sunscalartype* yvals  = N_VGetArrayPointer(y);
  sunscalartype* dyvals = N_VGetArrayPointer(ydot);
  sunscalartype a       = rdata[0]; /* access data entries (unused) */
  sunscalartype b       = rdata[1];
  sunscalartype c       = rdata[2];
  sunscalartype u       = yvals[0]; /* access solution values */
  sunscalartype v       = yvals[1];
  sunscalartype w       = yvals[2];

  /* fill in the RHS function */
  dyvals[0] = (t - 1.0) * (t - w) * SUN_I;
  dyvals[1] = SUN_RCONST(2.0) * u * SUN_I - v + SUN_RCONST(4.0);
  dyvals[2] = v + t * t * SUNRexp(-t) * SUN_I - SUNRexp(-t) * SUN_I +
              SUN_RCONST(1.0);

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

/* Compute the exact solution */
static int Solution(sunrealtype t, N_Vector u, void* user_data)
{
  sunscalartype* rdata =
    (sunscalartype*)user_data;      /* cast user_data to sunscalartype */
  sunscalartype a       = rdata[0]; /* access data entries */
  sunscalartype b       = rdata[1];
  sunscalartype c       = rdata[2];
  sunscalartype* uarray = N_VGetArrayPointer(u);
  if (check_flag((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

  uarray[0] = -t * SUNRexp(-t) + SUN_RCONST(2.0) * SUN_I;
  uarray[1] = -t * t * SUNRexp(-t) * SUN_I;
  uarray[2] = SUNRexp(-t) * SUN_I + t;

  return 0;
}

/* Compute the solution error */
static int SolutionError(sunrealtype t, N_Vector u, N_Vector e, void* user_data)
{
  /* Compute true solution */
  int flag = Solution(t, e, (void*)user_data);
  if (flag != 0) { return -1; }

  /* Compute max-norm of the error */
  N_VLinearSum(1.0, u, -1.0, e, e);
  sunrealtype error_norm = N_VMaxNorm(e);
  printf("    Max-norm of the error is %.5" ESYM "\n", error_norm);
  return 0;
}

/*---- end of file ----*/