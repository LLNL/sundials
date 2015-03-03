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
 * The following test simulates a simple 1D heat equation,
 *    u_t = k*u_xx + f
 * for t in [0, 10], x in [0, 1], with initial conditions
 *    u(0,x) =  0
 * Dirichlet boundary conditions, i.e. 
 *    u_t(t,0) = u_t(t,1) = 0,
 * and a point-source heating term,
 *    f = 1 for x=0.5.
 * 
 * The spatial derivatives are computed using second-order 
 * centered differences, with the data distributed over N points 
 * on a uniform spatial grid.
 *
 * This program solves the problem with either an ERK or DIRK
 * method.  For the DIRK method, we use a Newton iteration with 
 * the PCG linear solver, and a user-supplied Jacobian-vector 
 * product routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics 
 * are printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>            /* prototypes for ARKode fcts., consts. */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., macros */
#include <arkode/arkode_pcg.h>        /* prototype for ARKPcg solver */
#include <sundials/sundials_types.h>  /* def. of type 'realtype' */

/* user data structure */
typedef struct {
  long int N;    /* number of intervals   */
  realtype dx;   /* mesh spacing          */
  realtype k;    /* diffusion coefficient */
} *UserData;

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
            N_Vector fy, void *user_data, N_Vector tmp);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);

/* Main Program */
int main() {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);   /* initial time */
  realtype Tf = RCONST(1.0);   /* final time */
  int Nt = 10;                 /* total number of output times */
  realtype rtol = 1.e-6;       /* relative tolerance */
  realtype atol = 1.e-10;      /* absolute tolerance */
  UserData udata = NULL;
  realtype *data;
  long int N = 201;            /* spatial mesh size */
  realtype k = 0.5;            /* heat conductivity */
  long int i;

  /* general problem variables */
  int flag;                 /* reusable error-checking flag */
  N_Vector y = NULL;             /* empty vector for storing solution */
  void *arkode_mem = NULL;        /* empty ARKode memory structure */
  FILE *FID, *UFID;
  realtype t, dTout, tout;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nlcf, nni, ncfn, netf;

  /* allocate and fill udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->N = N;
  udata->k = k;
  udata->dx = RCONST(1.0)/(1.0*N-1.0);     /* mesh spacing */

  /* Initial problem output */
  printf("\n1D Heat PDE test problem:\n");
  printf("  N = %li\n", udata->N);
  printf("  diffusion coefficient:  k = %g\n", udata->k);

  /* Initialize data structures */
  y = N_VNew_Serial(N);            /* Create serial vector for solution */
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
  N_VConst(0.0, y);                /* Set initial conditions */
  arkode_mem = ARKodeCreate();     /* Create the solver memory */
  if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;

  /* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem, (void *) udata);   /* Pass udata to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
  flag = ARKodeSetMaxNumSteps(arkode_mem, 10000);         /* Increase max num steps  */
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;
  flag = ARKodeSStolerances(arkode_mem, rtol, atol);      /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Linear solver specification */
  flag = ARKPcg(arkode_mem, 0, N);                        /* Specify the PCG solver */
  if (check_flag(&flag, "ARKPcg", 1)) return 1;
  flag = ARKSpilsSetJacTimesVecFn(arkode_mem, Jac);       /* Set the Jacobian routine */
  if (check_flag(&flag, "ARKSpilsSetJacTimesVecFn", 1)) return 1;

  /* Specify linearly implicit RHS, with non-time-dependent Jacobian */
  flag = ARKodeSetLinear(arkode_mem, 0);
  if (check_flag(&flag, "ARKodeSetLinear", 1)) return 1;

  /* output mesh to disk */
  FID=fopen("heat_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
  fclose(FID);

  /* Open output stream for results, access data array */
  UFID=fopen("heat1D.txt","w");
  data = N_VGetArrayPointer(y);

  /* output initial condition to disk */
  for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
  fprintf(UFID,"\n");

  /* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  dTout = (Tf-T0)/Nt;
  tout = T0+dTout;
  printf("        t      ||u||_rms\n");
  printf("   -------------------------\n");
  printf("  %10.6f  %10.6f\n", t, sqrt(N_VDotProd(y,y)/N));
  for (iout=0; iout<Nt; iout++) {

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);         /* call integrator */
    if (check_flag(&flag, "ARKode", 1)) break;
    printf("  %10.6f  %10.6f\n", t, sqrt(N_VDotProd(y,y)/N));   /* print solution stats */
    if (flag >= 0) {                                            /* successful solve: update output time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                    /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }

    /* output results to disk */
    for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
    fprintf(UFID,"\n");
  }
  printf("   -------------------------\n");
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
  flag = ARKSpilsGetNumLinIters(arkode_mem, &nli);
  check_flag(&flag, "ARKSpilsGetNumLinIters", 1);
  flag = ARKSpilsGetNumJtimesEvals(arkode_mem, &nJv);
  check_flag(&flag, "ARKSpilsGetNumJtimesEvals", 1);
  flag = ARKSpilsGetNumConvFails(arkode_mem, &nlcf);
  check_flag(&flag, "ARKSpilsGetNumConvFails", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total linear iterations = %li\n", nli);
  printf("   Total number of Jacobian-vector products = %li\n", nJv);
  printf("   Total number of linear solver convergence failures = %li\n", nlcf);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy_Serial(y);        /* Free vectors */
  free(udata);                 /* Free user data */
  ARKodeFree(&arkode_mem);     /* Free integrator memory */
  return 0;
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;    /* access problem data */
  long int N  = udata->N;                   /* set variable shortcuts */
  realtype k  = udata->k;
  realtype dx = udata->dx;
  realtype *Y=NULL, *Ydot=NULL;
  realtype c1, c2;
  long int i, isource;

  Y = N_VGetArrayPointer(y);      /* access data arrays */
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return 1;
  Ydot = N_VGetArrayPointer(ydot);
  if (check_flag((void *) Ydot, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                      /* Initialize ydot to zero */

  /* iterate over domain, computing all equations */
  c1 = k/dx/dx;
  c2 = -RCONST(2.0)*k/dx/dx;
  isource = N/2;
  Ydot[0] = 0.0;                 /* left boundary condition */
  for (i=1; i<N-1; i++)
    Ydot[i] = c1*Y[i-1] + c2*Y[i] + c1*Y[i+1];
  Ydot[N-1] = 0.0;               /* right boundary condition */
  Ydot[isource] += 0.01/dx;      /* source term */

  return 0;                      /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
	       N_Vector fy, void *user_data, N_Vector tmp)
{
  UserData udata = (UserData) user_data;     /* variable shortcuts */
  long int N  = udata->N;
  realtype k  = udata->k;
  realtype dx = udata->dx;
  realtype *V=NULL, *JV=NULL;
  realtype c1, c2;
  long int i;

  V = N_VGetArrayPointer(v);       /* access data arrays */
  if (check_flag((void *) V, "N_VGetArrayPointer", 0)) return 1;
  JV = N_VGetArrayPointer(Jv);
  if (check_flag((void *) JV, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, Jv);                         /* initialize Jv product to zero */

  /* iterate over domain, computing all Jacobian-vector products */
  c1 = k/dx/dx;
  c2 = -RCONST(2.0)*k/dx/dx;
  JV[0] = 0.0;
  for (i=1; i<N-1; i++)
    JV[i] = c1*V[i-1] + c2*V[i] + c1*V[i+1];
  JV[N-1] = 0.0;

  return 0;                                  /* Return with success */
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
static int check_flag(void *flagvalue, char *funcname, int opt)
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
