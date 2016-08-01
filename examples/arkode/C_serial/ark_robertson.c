/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Example problem:
 * 
 * The following test simulates the Robertson problem, 
 * corresponding to the kinetics of an autocatalytic reaction.  
 * This is an ODE system with 3 components, Y = [u,v,w], satisfying
 * the equations,
 *    du/dt = -0.04*u + 1e4*v*w
 *    dv/dt = 0.04*u - 1e4*v*w - 3e7*v^2
 *    dw/dt = 3e7*v^2
 * for t in the interval [0.0, 1e11], with initial conditions 
 * Y0 = [1,0,0]. 
 * 
 * This program solves the problem with one of the solvers, ERK, 
 * DIRK or ARK.  For DIRK and ARK, implicit subsystems are solved 
 * using a Newton iteration with the ARKDENSE dense linear solver, 
 * and a user-supplied Jacobian routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics 
 * are printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode.h>            /* prototypes for ARKode fcts., consts. */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., macros */
#include <arkode/arkode_dense.h>      /* prototype for ARKDense solver */
#include <sundials/sundials_dense.h>  /* defs. of DlsMat and DENSE_ELEM */
#include <sundials/sundials_types.h>  /* def. of type 'realtype' */

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);     /* initial time */
  realtype Tf = RCONST(1.e11);   /* final time */
  realtype dTout = (Tf-T0)/100;  /* time between outputs */
  int Nt = ceil(Tf/dTout);       /* number of output times */
  long int NEQ = 3;              /* number of dependent vars. */

  /* general problem variables */
  int flag;                      /* reusable error-checking flag */
  N_Vector y = NULL;             /* empty vector for storing solution */
  void *arkode_mem = NULL;       /* empty ARKode memory structure */
  FILE *UFID;
  realtype t, tout;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;

  /* set up the initial conditions, tolerances, initial time step size */
  realtype u0 = RCONST(1.0);
  realtype v0 = RCONST(0.0);
  realtype w0 = RCONST(0.0);
  realtype reltol = 1.e-4;
  realtype abstol = 1.e-11;
  realtype h0 = 1.e-4 * reltol;

  /* Initial problem output */
  printf("\nRobertson ODE test problem:\n");
  printf("    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n",u0,v0,w0);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ);         /* Create serial vector for solution */
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
  NV_Ith_S(y,0) = u0;             /* Set initial conditions into y */
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(y,2) = w0;
  arkode_mem = ARKodeCreate();    /* Create the solver memory */
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;

  /* Call ARKodeInit to initialize the integrator memory and specify the
     right-hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  /* Set routines */
  flag = ARKodeSetInitStep(arkode_mem, h0);                /* Set custom initial step */
  if (check_flag(&flag, "ARKodeSetInitStep", 1)) return 1;
  flag = ARKodeSetMaxErrTestFails(arkode_mem, 20);         /* Increase max error test fails */
  if (check_flag(&flag, "ARKodeSetMaxErrTestFails", 1)) return 1;
  flag = ARKodeSetMaxNonlinIters(arkode_mem, 8);           /* Increase max nonlin iters  */
  if (check_flag(&flag, "ARKodeSetMaxNonlinIters", 1)) return 1;
  flag = ARKodeSetNonlinConvCoef(arkode_mem, 1.e-7);       /* set nonlinear convergence coeff. */
  if (check_flag(&flag, "ARKodeSetNonlinConvCoef", 1)) return 1;
  flag = ARKodeSetMaxNumSteps(arkode_mem, 100000);         /* Increase max num steps */
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;
  flag = ARKodeSetPredictorMethod(arkode_mem, 1);         /* Specify maximum-order predictor */
  if (check_flag(&flag, "ARKodeSetPredictorMethod", 1)) return 1;
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);   /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Linear solver specification */
  flag = ARKDense(arkode_mem, NEQ);                        /* Specify dense linear solver */
  if (check_flag(&flag, "ARKDense", 1)) return 1;
  flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);             /* Set the Jacobian routine */
  if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return 1;

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt","w");
  fprintf(UFID,"# t u v w\n");

  /* output initial condition to disk */
  fprintf(UFID," %.16e %.16e %.16e %.16e\n", 
	  T0, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));  

  /* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  tout = T0+dTout;
  printf("        t           u           v           w\n");
  printf("   --------------------------------------------------\n");
  printf("  %10.3e  %12.5e  %12.5e  %12.5e\n",
      t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
  for (iout=0; iout<Nt; iout++) {

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);       /* call integrator */
    if (check_flag(&flag, "ARKode", 1)) break;
    printf("  %10.3e  %12.5e  %12.5e  %12.5e\n",              /* access/print solution */
        t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
    fprintf(UFID," %.16e %.16e %.16e %.16e\n", 
	    t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));  
    if (flag >= 0) {                                          /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                  /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   --------------------------------------------------\n");
  fclose(UFID);

  /* Print some final statistics */
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumStepAttempts", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1);
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1);
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1);
  flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKDlsGetNumJacEvals", 1);
  flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKDlsGetNumRhsEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", 
	 nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy_Serial(y);        /* Free y vector */
  ARKodeFree(&arkode_mem);     /* Free integrator memory */
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype u = NV_Ith_S(y,0);   /* access current solution */
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* Fill in ODE RHS function */
  NV_Ith_S(ydot,0) = -0.04*u + 1.e4*v*w;
  NV_Ith_S(ydot,1) = 0.04*u - 1.e4*v*w - 3.e7*v*v;
  NV_Ith_S(ydot,2) = 3.e7*v*v;

  return 0;                     /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype v = NV_Ith_S(y,1);   /* access current solution */
  realtype w = NV_Ith_S(y,2);
  SetToZero(J);                 /* initialize Jacobian to zero */

  /* Fill in the Jacobian of the ODE RHS function */
  DENSE_ELEM(J,0,0) = -0.04;
  DENSE_ELEM(J,0,1) = 1.e4*w;
  DENSE_ELEM(J,0,2) = 1.e4*v;

  DENSE_ELEM(J,1,0) = 0.04;
  DENSE_ELEM(J,1,1) = -1.e4*w - 6.e7*v;
  DENSE_ELEM(J,1,2) = -1.e4*v;

  DENSE_ELEM(J,2,1) = 6.e7*v;

  return 0;                     /* Return with success */
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
