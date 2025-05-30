/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is a simple example problem for IDA to show its usage for
 * complex-valued applications, and to assess the accuracy of
 * numerical solution. The ODE system with 3 components,
 * X = [u,v,w], satisfies the equations,
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
 * We trivially convert this to the DAE system:
 *
 * 0 = du/dt - (t - 1)*(t - w)*i
 * 0 = dv/dt - 2*u*i + v - 4
 * 0 = dw/dt - v - t^2*exp(-t)*i + exp(-t)*i - 1
 *
 * The implicit systems at each time step are solved using a
 * Newton iteration with the SPFGMR linear solver.
 *
 * 5 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *
 * -----------------------------------------------------------------*/

#include <ida/ida.h> /* prototypes for IDA fcts., consts. */
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector */
#include <stdio.h>
#include <sundials/sundials_math.h>  /* defs. of SUNRabs, SUNRexp, etc. */
#include <sundials/sundials_types.h> /* defs. of sunscalartype, sunindextype */
#include <sunlinsol/sunlinsol_spbcgs.h> /* access to spbcgs SUNLinearSolver */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* Prototypes of functions called by IDA */

/* DAE residual */
static int res(sunrealtype tres, N_Vector yy, N_Vector yp, N_Vector resval,
               void* user_data);

/* Compute the true solution and derivative */
static int Solution(sunrealtype t, N_Vector u, void* user_data);
static int SolutionDerivative(sunrealtype t, N_Vector up, void* user_data);

/* Compute the numerical solution and derivative error */
static int SolutionError(sunrealtype t, N_Vector u, N_Vector e, void* user_data);
static int DerivativeError(sunrealtype t, N_Vector u_p, N_Vector e_p,
                           void* user_data);

/* Private function to check function return values */
static int check_retval(void* returnvalue, const char* funcname, int opt);

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

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
  int maxl           = 10; /* max linear solver iterations */

  /* general problem variables */
  int retval;                 /* reusable error-checking flag */
  N_Vector y          = NULL; /* empty vector for storing solution */
  N_Vector y_p        = NULL; /* empty vector for storing solution derivative */
  N_Vector True_Sol   = NULL; // vector for storing true solution
  N_Vector True_Sol_p = NULL; // vector for storing true solution derivative
  N_Vector Error      = NULL; // vector for storing the error */
  N_Vector Error_p    = NULL; // vector for storing the derivative error */
  SUNLinearSolver LS  = NULL; /* empty linear solver object */
  void* ida_mem       = NULL; /* empty IDA memory structure */
  sunscalartype rdata[3];
  sunrealtype t, tout;
  int iout;

  /* Create the SUNDIALS context object for this testing */
  SUNContext ctx;
  retval = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /* Create vectors */
  y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return 1; }
  y_p                  = N_VClone(y);
  True_Sol             = N_VClone(y);
  True_Sol_p           = N_VClone(y);
  Error                = N_VClone(y);
  Error_p              = N_VClone(y);
  sunscalartype* yvals = N_VGetArrayPointer(y);
  if (check_retval(yvals, "N_VGetArrayPointer", 0)) { return 1; }
  sunscalartype* ypvals = N_VGetArrayPointer(y_p);
  if (check_retval(yvals, "N_VGetArrayPointer", 0)) { return 1; }

  // set up the problem data (unused in this case)
  rdata[0] = SUN_CCONST(1.0, 0.0);
  rdata[1] = SUN_CCONST(1.0, 0.0);
  rdata[2] = SUN_CCONST(1.0, 0.0);

  // Set initial condition
  retval = Solution(0.0, True_Sol, (void*)rdata);
  if (check_retval(&retval, "Solution", 1)) { return 1; }
  retval = SolutionDerivative(0.0, y_p, (void*)rdata);
  if (check_retval(&retval, "SolutionDerivative", 1)) { return 1; }

  N_VScale(1.0, True_Sol, y); /* Set initial conditions */
  N_VScale(1.0, True_Sol_p, y_p);

  /* Initial problem output */
  printf("\nAnalytic ODE test problem:\n");
  printf("    problem parameters:  a = %" GSYM " + %" GSYM "i,  b = %" GSYM
         " + %" GSYM "i,  c = %" GSYM " + %" GSYM "i\n",
         SUN_REAL(rdata[0]), SUN_IMAG(rdata[0]), SUN_REAL(rdata[1]),
         SUN_IMAG(rdata[2]), SUN_REAL(rdata[2]), SUN_IMAG(rdata[2]));
  printf("    reltol = %.1" ESYM ",  abstol = %.1" ESYM "\n\n", reltol, abstol);

  /* Call IDACreate and IDAInit to initialize IDA memory */
  ida_mem = IDACreate(ctx);
  if (check_retval((void*)ida_mem, "IDACreate", 0)) { return (1); }
  retval = IDAInit(ida_mem, res, T0, y, y_p);
  if (check_retval(&retval, "IDAInit", 1)) { return (1); }

  /* Set routines */
  retval = IDASetUserData(ida_mem, (void*)rdata); /* Pass rdata to user functions */
  if (check_retval(&retval, "IDASetUserData", 1)) { return 1; }

  retval = IDASStolerances(ida_mem, reltol, abstol); /* Specify tolerances */
  if (check_retval(&retval, "IDASStolerances", 1)) { return 1; }

  LS = SUNLinSol_SPBCGS(y, 0, maxl, ctx);
  if (check_retval((void*)LS, "SUNLinSol_SPBCGS", 0)) { return 1; }

  /* Linear solver interface */
  retval = IDASetLinearSolver(ida_mem, LS, NULL);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) { return 1; }

  /* Main time-stepping loop: calls CVode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t    = T0;
  tout = T0 + dTout;
  printf(
    "      t              u                        u'                       v"
    "                        v'                       w                        "
    "w'\n");
  printf("   "
         "---------------------------------------------------------------------"
         "------"
         "---------------------------------------------------------------------"
         "-----------\n");
  printf(" %8.3" FSYM " | "
         "%8.5" FSYM " + "
         "%8.5" FSYM "i  |  "
         "%8.5" FSYM " + "
         "%8.5" FSYM "i  |  "
         "%8.5" FSYM " + "
         "%8.5" FSYM "i  |  "
         "%8.5" FSYM " + "
         "%8.5" FSYM "i  |  "
         "%8.5" FSYM " + "
         "%8.5" FSYM "i  |  "
         "%8.5" FSYM " + "
         "%8.5" FSYM "i  \n",
         t, SUN_REAL(yvals[0]), SUN_IMAG(yvals[0]), SUN_REAL(ypvals[0]),
         SUN_IMAG(ypvals[0]), SUN_REAL(yvals[1]), SUN_IMAG(yvals[1]),
         SUN_REAL(ypvals[1]), SUN_IMAG(ypvals[1]), SUN_REAL(yvals[2]),
         SUN_IMAG(yvals[2]), SUN_REAL(ypvals[2]), SUN_IMAG(ypvals[2]));

  for (iout = 0; iout < Nt; iout++)
  {
    retval = IDASolve(ida_mem, tout, &t, y, y_p, IDA_NORMAL); /* call integrator */
    if (check_retval(&retval, "IDASolve", 1)) { break; }
    printf(" %8.3" FSYM " | "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  \n",
           t, SUN_REAL(yvals[0]), SUN_IMAG(yvals[0]), SUN_REAL(ypvals[0]),
           SUN_IMAG(ypvals[0]), SUN_REAL(yvals[1]), SUN_IMAG(yvals[1]),
           SUN_REAL(ypvals[1]), SUN_IMAG(ypvals[1]), SUN_REAL(yvals[2]),
           SUN_IMAG(yvals[2]), SUN_REAL(ypvals[2]), SUN_IMAG(ypvals[2]));

    if (retval >= 0)
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
         "------"
         "---------------------------------------------------------------------"
         "-----------\n");

  SolutionError(Tf, y, Error, (void*)rdata);
  DerivativeError(Tf, y_p, Error_p, (void*)rdata);

  /* Print all final statistics */
  printf("\nFinal Solver Statistics:\n");
  retval = IDAPrintAllStats(ida_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_retval(&retval, "IDAPrintAllStats", 1)) { return 1; }

  /* Clean up and return with successful completion */
  N_VDestroy(y); /* Free vectors */
  N_VDestroy(y_p);
  N_VDestroy(True_Sol);
  N_VDestroy(True_Sol_p);
  N_VDestroy(Error);
  N_VDestroy(Error_p);
  IDAFree(&ida_mem);     /* Free integrator memory */
  SUNLinSolFree(LS);     /* Free linear solver */
  SUNContext_Free(&ctx); /* Free context */

  return (retval);
}

/*
 *--------------------------------------------------------------------
 * Functions called by IDA
 *--------------------------------------------------------------------
 */

/*
 * Define the system residual function.
 */

int res(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
  sunscalartype* rdata =
    (sunscalartype*)user_data; /* cast user_data to sunscalartype */
  sunscalartype* yvals   = N_VGetArrayPointer(yy);
  sunscalartype* ypvals  = N_VGetArrayPointer(yp);
  sunscalartype* resvals = N_VGetArrayPointer(rr);
  sunscalartype a        = rdata[0]; /* access data entries (unused) */
  sunscalartype b        = rdata[1];
  sunscalartype c        = rdata[2];
  sunscalartype u        = yvals[0]; /* access solution values */
  sunscalartype v        = yvals[1];
  sunscalartype w        = yvals[2];
  sunscalartype u_p      = ypvals[0];
  sunscalartype v_p      = ypvals[1];
  sunscalartype w_p      = ypvals[2];

  resvals[0] = u_p - (t - 1.0) * (t - w) * SUN_I;
  resvals[1] = v_p - SUN_RCONST(2.0) * u * SUN_I + v - SUN_RCONST(4.0);
  resvals[2] = w_p - v - t * t * SUNRexp(-t) * SUN_I + SUNRexp(-t) * SUN_I -
               SUN_RCONST(1.0);

  return (0);
}

/*
 *--------------------------------------------------------------------
 * Private helper functions
 *--------------------------------------------------------------------
 */

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }
  else if (opt == 1)
  {
    /* Check if retval < 0 */
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }
  else if (opt == 2 && returnvalue == NULL)
  {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
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
  if (check_retval((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

  uarray[0] = -t * SUNRexp(-t) + SUN_RCONST(2.0) * SUN_I;
  uarray[1] = -t * t * SUNRexp(-t) * SUN_I;
  uarray[2] = SUNRexp(-t) * SUN_I + t;

  return 0;
}

/* Compute the derivative of the exact solution */
static int SolutionDerivative(sunrealtype t, N_Vector u_p, void* user_data)
{
  sunscalartype* rdata =
    (sunscalartype*)user_data;       /* cast user_data to sunscalartype */
  sunscalartype a        = rdata[0]; /* access data entries */
  sunscalartype b        = rdata[1];
  sunscalartype c        = rdata[2];
  sunscalartype* uparray = N_VGetArrayPointer(u_p);
  if (check_retval((void*)uparray, "N_VGetArrayPointer", 0)) { return -1; }

  uparray[0] = t * SUNRexp(-t) - SUNRexp(-t);
  uparray[1] = t * t * SUNRexp(-t) * SUN_I -
               SUN_RCONST(2.0) * t * SUNRexp(-t) * SUN_I;
  uparray[2] = SUN_RCONST(1.0) - SUNRexp(-t) * SUN_I;

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

/* Compute the solution derivative error */
static int DerivativeError(sunrealtype t, N_Vector u_p, N_Vector e_p,
                           void* user_data)
{
  /* Compute true solution */
  int flag = SolutionDerivative(t, e_p, (void*)user_data);
  if (flag != 0) { return -1; }

  /* Compute max-norm of the error */
  N_VLinearSum(1.0, u_p, -1.0, e_p, e_p);
  sunrealtype error_norm = N_VMaxNorm(e_p);
  printf("    Max-norm of the derivative error is %.5" ESYM "\n", error_norm);
  return 0;
}

/*---- end of file ----*/