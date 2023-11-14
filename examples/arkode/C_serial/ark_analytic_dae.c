/*-----------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
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
 * Example problem:
 *
 * The following is a simple example problem with analytical
 * solution,
 *     dy/dt = (t+1)*exp(-y)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 * This has analytical solution
 *      y(t) = log(0.5*t^2 + t + 1)
 *
 * This program solves the problem with the ERK method.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>            /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_serial.h>           /* serial N_Vector types, fcts., macros */
#include <sunmatrix/sunmatrix_dense.h>        /* dense SUNMatrix for mass and Jacobian matrix */
#include <sunlinsol/sunlinsol_dense.h>        /* dense linear solver */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* Newton solver */
#include <sundials/sundials_types.h>          /* def. of type 'realtype' */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SIN(x) (sin(x))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SIN(x) (sinf(x))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SIN(x) (sinl(x))
#endif

/* User-supplied Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);

static int jacobian(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int mass(sunrealtype t, SUNMatrix mass, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Main Program */
int main()
{
  /* general problem parameters */
  sunrealtype T0 = RCONST(0.0);     /* initial time */
  sunrealtype Tf = RCONST(10.0);    /* final time */
  sunrealtype dTout = RCONST(1.0);  /* time between outputs */
  sunindextype NEQ = 2;             /* number of dependent vars. */
  sunrealtype BETA = 1;             /* DAE parameter */
  sunrealtype reltol = 1.0e-6;      /* tolerances */
  sunrealtype abstol = 1.0e-10;

  /* general problem variables */
  int flag;                         /* reusable error-checking flag */
  N_Vector y = NULL;                /* empty vector for storing solution */
  SUNMatrix jacobian_matrix = NULL; /* empty matrix for storing Jacobian */
  SUNLinearSolver ls = NULL;        /* empty linear solver */
  SUNNonlinearSolver nls = NULL;    /* empty nonlinear solver */
  void *arkode_mem = NULL;          /* empty ARKode memory structure */
  FILE *UFID, *FID;
  realtype t, tout;

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  flag = SUNContext_Create(NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) return 1;

  /* Initial problem output */
  printf("\nAnalytical DAE test problem:\n");
  printf("   reltol = %.1"ESYM"\n",  reltol);
  printf("   abstol = %.1"ESYM"\n\n",abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ, ctx);          /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  NV_Ith_S(y,0) = 1.0;             /* Specify initial condition */
  NV_Ith_S(y,1) = BETA;

  jacobian_matrix = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_flag((void *)jacobian_matrix, "SUNDenseMatrix", 0)) return 1;

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y, ctx);
  if (check_flag((void *)arkode_mem, "ARKStepCreate", 0)) return 1;

  flag = ARKStepSetUserData(arkode_mem, &BETA);
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

  flag = ARKStepSetJacFn(arkode_mem, jacobian);
  if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;

  flag = ARKStepSetMassFn(arkode_mem, mass);
  if (check_flag(&flag, "ARKStepSetMassFn", 1)) return 1;

  /* Specify solvers */
  ls = SUNLinSol_Dense(y, jacobian, ctx);
  if (check_flag((void *)ls, "SUNLinSol_Dense", 0)) return 1;

  nls = SUNNonlinSol_Newton(y, ctx);
  if (check_flag((void *)nls, "SUNNonlinSol_Newton", 0)) return 1;

  flag = ARKStepSetLinear(arkode_mem, 1);
  if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;

  /* Specify tolerances */
  flag = ARKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt","w");
  fprintf(UFID,"# t u\n");

  /* output initial condition to disk */
  fprintf(UFID," %.16"ESYM" %.16"ESYM"\n", T0, NV_Ith_S(y,0));

  /* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  tout = T0+dTout;
  printf("        t           u\n");
  printf("   ---------------------\n");
  while (Tf - t > 1.0e-15) {

    flag = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);       /* call integrator */
    if (check_flag(&flag, "ARKStepEvolve", 1)) break;
    printf("  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"\n",
           t, NV_Ith_S(y,0), NV_Ith_S(y,1));           /* access/print solution */
    fprintf(UFID," %.16"ESYM" %.16"ESYM"\n", t, NV_Ith_S(y,0));
    if (flag >= 0) {                                          /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                  /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   ---------------------\n");
  fclose(UFID);

  /* Print final statistics */
  printf("\nFinal Statistics:\n");
  flag = ARKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  /* Print final statistics to a file in CSV format */
  FID = fopen("ark_analytic_nonlin_stats.csv", "w");
  flag = ARKStepPrintAllStats(arkode_mem, FID, SUN_OUTPUTFORMAT_CSV);
  fclose(FID);

  /* Clean up and return with successful completion */
  N_VDestroy(y);               /* Free y vector */
  ARKStepFree(&arkode_mem);    /* Free integrator memory */
  SUNNonlinSolFree(nls);       /* Free nonlinear solver */
  SUNLinSolFree(ls);           /* Free linear solver */
  SUNMatFree(jacobian_matrix); /* Free Jacobian matrix */
  SUNContext_Free(&ctx);       /* Free context */

  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  sunrealtype BETA = (sunrealtype) *user_data;
  sunrealtype y1 = NV_Ith_S(y, 0);
  sunrealtype y2 = NV_Ith_S(y, 1);
  NV_Ith_S(ydot, 0) = (1.0 + t) * y2 - y1;
  NV_Ith_S(ydot, 1) = BETA * y1 - (1.0 + BETA) * y2 + SIN(t);
  return 0;
}

static int jacobian(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype BETA = (sunrealtype) *user_data;
  SM_ELEMENT_D(Jac, 0, 0) = -1.0;
  SM_ELEMENT_D(Jac, 0, 1) = 1.0 + t;
  SM_ELEMENT_D(Jac, 1, 0) = BETA;
  SM_ELEMENT_D(Jac, 1, 1) = -1.0 - BETA * t;
  return 0;
}

static int mass(sunrealtype t, SUNMatrix mass, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  SM_ELEMENT_D(mass, 0, 0) = 1.0;
  SM_ELEMENT_D(mass, 0, 1) = -t;
  SM_ELEMENT_D(mass, 1, 0) = 0.0;
  SM_ELEMENT_D(mass, 1, 1) = 0.0;
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
static int check_flag(void *flagvalue, const char *funcname, int opt)
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
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


/*---- end of file ----*/
