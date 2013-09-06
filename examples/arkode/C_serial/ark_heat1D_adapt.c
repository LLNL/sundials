/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Example problem:
 
 The following test simulates a simple 1D heat equation,
    u_t = k*u_xx + f
 for t in [0, 10], x in [0, 1], with initial conditions
    u(0,x) =  0
 Dirichlet boundary conditions, i.e. 
    u_t(t,0) = u_t(t,1) = 0,
 and a point-source heating term,
    f = 1 for x=0.5.
 
 The spatial derivatives are computed using a three-point 
 centered stencil (second order for a uniform mesh).  The data
 is initially univormly distributed over N points in the interval
 [0, 1], but as the simulation proceeds the mesh is adapted.

 This program solves the problem with either an ERK or DIRK
 method.  For the DIRK method, we use a Newton iteration with 
 the PCG linear solver, and a user-supplied Jacobian-vector 
 product routine.

 100 outputs are printed at equal intervals, and run statistics 
 are printed at the end.
---------------------------------------------------------------*/

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
  long int N;           /* current number of intervals */
  realtype *x;          /* current mesh */
  realtype k;           /* diffusion coefficient */
  realtype refine_tol;  /* adaptivity tolerance */
} *UserData;

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
            N_Vector fy, void *user_data, N_Vector tmp);

/* Private function to check function return values */
realtype * adapt_mesh(N_Vector y, long int *Nnew, UserData udata);
static int project(long int Nold, realtype *xold, N_Vector yold, 
		   long int Nnew, realtype *xnew, N_Vector ynew);
static int check_flag(void *flagvalue, char *funcname, int opt);

/* Main Program */
int main() {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);   /* initial time */
  realtype Tf = RCONST(1.0);   /* final time */
  int Nt = 50;                 /* total number of output times */
  realtype rtol = 1.e-3;       /* relative tolerance */
  realtype atol = 1.e-10;      /* absolute tolerance */
  realtype hscale = 1.0;       /* time step change factor on resizes */
  UserData udata = NULL;
  realtype *data;
  long int N = 21;             /* initial spatial mesh size */
  realtype refine = 3.e-3;     /* adaptivity refinement tolerance */
  realtype k = 0.5;            /* heat conductivity */
  long int i, nst, nst_cur=0, nli, nli_tot=0;

  /* general problem variables */
  int flag;                    /* reusable error-checking flag */
  N_Vector y  = NULL;          /* empty vector for storing solution */
  N_Vector y2 = NULL;          /* empty vector for storing solution */
  N_Vector yt = NULL;          /* empty vector for swapping */
  void *arkode_mem = NULL;     /* empty ARKode memory structure */

  /* allocate and fill initial udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->N = N;
  udata->k = k;
  udata->refine_tol = refine;
  udata->x = malloc(N * sizeof(realtype));
  for (i=0; i<N; i++)  udata->x[i] = 1.0*i/(N-1);

  /* Initial problem output */
  printf("\n1D adaptive Heat PDE test problem:\n");
  printf("  diffusion coefficient:  k = %g\n", udata->k);
  printf("  initial N = %li\n", udata->N);

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

  /* output mesh to disk */
  FILE *XFID=fopen("heat_mesh.txt","w");

  /* output initial mesh to disk */
  for (i=0; i<udata->N; i++)  fprintf(XFID," %.16e", udata->x[i]);
  fprintf(XFID,"\n");


  /* Open output stream for results, access data array */
  FILE *UFID=fopen("heat1D.txt","w");

  /* output initial condition to disk */
  data = N_VGetArrayPointer(y);
  for (i=0; i<udata->N; i++)  fprintf(UFID," %.16e", data[i]);
  fprintf(UFID,"\n");

  /* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  realtype t = T0;
  realtype dTout = (Tf-T0)/Nt;
  realtype tout = T0;
  printf("        t      ||u||_rms    N    steps\n");
  printf("   ------------------------------------\n");
  printf("  %10.6f  %10.6f    %li     %i\n", 
	 t, sqrt(N_VDotProd(y,y)/udata->N), udata->N, 0);
  int iout;
  realtype *xnew=NULL;
  long int Nnew;
  for (iout=0; iout<Nt; iout++) {

    tout += dTout;                                              /* set next output time*/
    tout = (tout > Tf) ? Tf : tout;

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);         /* call integrator */
    if (check_flag(&flag, "ARKode", 1)) break;
    flag = ARKodeGetNumSteps(arkode_mem, &nst);
    check_flag(&flag, "ARKodeGetNumSteps", 1);
    printf("  %10.6f  %10.6f    %li     %li\n",                   /* print solution stats */
	   t, sqrt(N_VDotProd(y,y)/udata->N), udata->N, nst-nst_cur);
    nst_cur = nst;
    if (flag < 0) {                                             /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }

    /* output results and current mesh to disk */
    data = N_VGetArrayPointer(y);
    for (i=0; i<udata->N; i++)  fprintf(UFID," %.16e", data[i]);
    fprintf(UFID,"\n");
    for (i=0; i<udata->N; i++)  fprintf(XFID," %.16e", udata->x[i]);
    fprintf(XFID,"\n");

    /* accumulate total linear iterations for this step */
    flag = ARKSpilsGetNumLinIters(arkode_mem, &nli);
    check_flag(&flag, "ARKSpilsGetNumLinIters", 1);
    nli_tot += nli;

    /* adapt the spatial mesh */
    xnew = adapt_mesh(y, &Nnew, udata);
    if (check_flag(xnew, "ark_adapt", 0)) break;

    /* create N_Vector of new length */
    y2 = N_VNew_Serial(Nnew);
    if (check_flag((void *) y2, "N_VNew_Serial", 0)) break;
    
    /* project solution onto new mesh */
    flag = project(udata->N, udata->x, y, Nnew, xnew, y2);
    if (check_flag(&flag, "project", 1)) break;

    /* delete old vector, old mesh */
    N_VDestroy_Serial(y);
    free(udata->x);
    
    /* swap x and xnew so that new mesh is stored in udata structure */
    udata->x = xnew;
    xnew = NULL;
    udata->N = Nnew;   /* store size of new mesh */
    
    /* swap y and y2 so that y holds new solution */
    yt = y;
    y  = y2;
    y2 = yt;

    /* call ARKodeResize to notify integrator of change in mesh */
    flag = ARKodeResize(arkode_mem, y, hscale, tout, NULL, NULL);
    if (check_flag(&flag, "ARKodeResize", 1)) break;

    /* destroy and re-allocate linear solver memory */
    flag = ARKPcg(arkode_mem, 0, udata->N);
    if (check_flag(&flag, "ARKPcg", 1)) break;
    flag = ARKSpilsSetJacTimesVecFn(arkode_mem, Jac);
    if (check_flag(&flag, "ARKSpilsSetJacTimesVecFn", 1)) break;

  }
  printf("   ------------------------------------\n");
  fclose(UFID);
  fclose(XFID);

  /* Print some final statistics */
  long int nst_a, nfe, nfi, nsetups, nJv, nlcf, nni, ncfn, netf;
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
  free(udata->x);              /* Free user data */
  free(udata);   
  ARKodeFree(&arkode_mem);     /* Free integrator memory */
  return 0;
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(0.0, ydot);                      /* Initialize ydot to zero */
  UserData udata = (UserData) user_data;    /* access problem data */
  long int N  = udata->N;                   /* set variable shortcuts */
  realtype k  = udata->k;
  realtype *x = udata->x;
  realtype *Y = N_VGetArrayPointer(y);      /* access data arrays */
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return 1;
  realtype *Ydot = N_VGetArrayPointer(ydot);
  if (check_flag((void *) Ydot, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over domain, computing all equations */
  realtype dxL, dxR;
  long int i;
  Ydot[0] = 0.0;                 /* left boundary condition */
  for (i=1; i<N-1; i++) {        /* interior */
    dxL = x[i]-x[i-1];
    dxR = x[i+1]-x[i];
    Ydot[i] = Y[i-1]*k*2.0/(dxL*(dxL+dxR)) 
            - Y[i]*k*2.0/(dxL*dxR)
            + Y[i+1]*k*2.0/(dxR*(dxL+dxR));
  }
  Ydot[N-1] = 0.0;               /* right boundary condition */

  /* source term */
  for (i=0; i<N-1; i++) {
    Ydot[i] += 2.0*exp(-200.0*(x[i]-0.25)*(x[i]-0.25))
                 - exp(-400.0*(x[i]-0.7)*(x[i]-0.7))
                 + exp(-500.0*(x[i]-0.4)*(x[i]-0.4))
             - 2.0*exp(-600.0*(x[i]-0.55)*(x[i]-0.55));
  }

  return 0;                      /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
	       N_Vector fy, void *user_data, N_Vector tmp)
{
  N_VConst(0.0, Jv);                         /* initialize Jv product to zero */
  UserData udata = (UserData) user_data;     /* variable shortcuts */
  long int N  = udata->N;
  realtype k  = udata->k;
  realtype *x = udata->x;
  realtype *V = N_VGetArrayPointer(v);       /* access data arrays */
  if (check_flag((void *) V, "N_VGetArrayPointer", 0)) return 1;
  realtype *JV = N_VGetArrayPointer(Jv);
  if (check_flag((void *) JV, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over domain, computing all Jacobian-vector products */
  realtype dxL, dxR;
  long int i;
  JV[0] = 0.0;
  for (i=1; i<N-1; i++) {
    dxL = x[i]-x[i-1];
    dxR = x[i+1]-x[i];
    JV[i] = V[i-1]*k*2.0/(dxL*(dxL+dxR)) 
          - V[i]*k*2.0/(dxL*dxR)
          + V[i+1]*k*2.0/(dxR*(dxL+dxR));
  }
  JV[N-1] = 0.0;

  return 0;                                  /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Adapts the current mesh, using a simple adaptivity strategy of 
   refining when an approximation of the scaled second-derivative is 
   too large.  We only do this in one sweep, so no attempt is made to 
   ensure the resulting mesh meets these same criteria after adaptivity:
      y [input] -- the current solution vector
      Nnew [output] -- the size of the new mesh
      udata [input] -- the current system information 
   The return for this function is a pointer to the new mesh. */
realtype * adapt_mesh(N_Vector y, long int *Nnew, UserData udata)
{
  int i, j, k;

  /* Access current solution and mesh arrays */
  realtype *Y = N_VGetArrayPointer(y);
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return NULL;
  realtype *xold = udata->x;

  /* create marking array */
  int *marks = calloc(udata->N-1, sizeof(int));

  /* /\* perform marking:  */
  /*     0 -> leave alone */
  /*     1 -> refine */
  /* realtype ymax, ymin; */
  /* for (i=0; i<(udata->N-1); i++) { */

  /*   /\* check for refinement *\/ */
  /*   if (fabs(Y[i+1] - Y[i]) > udata->refine_tol) { */
  /*     marks[i] = 1; */
  /*     continue; */
  /*   } */
  /* } */

  /* perform marking: 
      0 -> leave alone
      1 -> refine */
  realtype ydd;
  for (i=1; i<udata->N-1; i++) {

    /* approximate scaled second-derivative */
    ydd = Y[i-1] - 2.0*Y[i] + Y[i+1];

    /* check for refinement */
    if (fabs(ydd) > udata->refine_tol) {
      marks[i-1] = 1;
      marks[i] = 1;
    }
    
  }

  /* allocate new mesh */
  long int num_refine = 0;
  for (i=0; i<udata->N-1; i++) 
    if (marks[i] == 1)   num_refine++;
  long int N_new = udata->N + num_refine;
  *Nnew = N_new;            /* Store new array length */
  realtype *xnew = malloc((N_new) * sizeof(realtype));
  

  /* fill new mesh */
  xnew[0] = udata->x[0];    /* store endpoints */
  xnew[N_new-1] = udata->x[udata->N-1];
  j=1;
  /* iterate over old intervals */
  for (i=0; i<udata->N-1; i++) {
    /* if mark is 0, reuse old interval */ 
    if (marks[i] == 0) {
      xnew[j++] = xold[i+1];
      continue;
    }
    
    /* if mark is 1, refine old interval */
    if (marks[i] == 1) {
      xnew[j++] = 0.5*(xold[i]+xold[i+1]);
      xnew[j++] = xold[i+1];
      continue;
    }
  }

  /* verify that new mesh is legal */
  for (i=0; i<N_new-1; i++) {
    if (xnew[i+1] <= xnew[i]) {
      fprintf(stderr,"adapt_mesh error: illegal mesh created\n");
      free(xnew);
      return NULL;
    }
  }

  free(marks);              /* Delete marking array */
  return xnew;              /* Return with success */
}


/* Projects one vector onto another:
      Nold [input] -- the size of the old mesh
      xold [input] -- the old mesh
      yold [input] -- the vector defined over the old mesh
      Nnew [input] -- the size of the new mesh
      xnew [input] -- the new mesh
      ynew [output] -- the vector defined over the new mesh
                       (allocated prior to calling project) */
static int project(long int Nold, realtype *xold, N_Vector yold, 
		   long int Nnew, realtype *xnew, N_Vector ynew)
{
  /* Access data arrays */
  realtype *Yold = N_VGetArrayPointer(yold);    /* access data arrays */
  if (check_flag((void *) Yold, "N_VGetArrayPointer", 0)) return 1;
  realtype *Ynew = N_VGetArrayPointer(ynew);
  if (check_flag((void *) Ynew, "N_VGetArrayPointer", 0)) return 1;

  /* loop over new mesh, finding corresponding interval within old mesh, 
     and perform piecewise linear interpolation from yold to ynew */
  int iv=0;
  int i, j;
  for (i=0; i<Nnew; i++) {
    
    /* find old interval, start with previous value since sorted */
    for (j=iv; j<Nold-1; j++) {
      if (xnew[i] >= xold[j] && xnew[i] <= xold[j+1]) {
	iv = j;
	break;
      }
      iv = Nold-1;     /* just in case it wasn't found above */
    }

    /* perform interpolation */ 
    Ynew[i] = Yold[iv]*(xnew[i]-xold[iv+1])/(xold[iv]-xold[iv+1]) 
            + Yold[iv+1]*(xnew[i]-xold[iv])/(xold[iv+1]-xold[iv]);
  }

  return 0;            /* Return with success */
}


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
