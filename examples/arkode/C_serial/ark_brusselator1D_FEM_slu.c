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
 * The following test simulates a brusselator problem from chemical 
 * kinetics.  This is a PDE system with 3 components, Y = [u,v,w], 
 * satisfying the equations,
 *    u_t = du*u_xx + a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx + w*u - v*u^2
 *    w_t = dw*w_xx + (b-w)/ep - w*u
 * for t in [0, 80], x in [0, 1], with initial conditions
 *    u(0,x) =  a  + 0.1*sin(pi*x)
 *    v(0,x) = b/a + 0.1*sin(pi*x)
 *    w(0,x) =  b  + 0.1*sin(pi*x),
 * and with stationary boundary conditions, i.e. 
 *    u_t(t,0) = u_t(t,1) = 0
 *    v_t(t,0) = v_t(t,1) = 0
 *    w_t(t,0) = w_t(t,1) = 0.
 * 
 * Here, we use a piecewise linear Galerkin finite element 
 * discretization in space, where all element-wise integrals are 
 * computed using 3-node Gaussian quadrature (since we will have 
 * quartic polynomials in the reaction terms for the u_t and v_t 
 * equations, including the test function).  The time derivative 
 * terms for this system will include a mass matrix, giving rise 
 * to an ODE system of the form
 *      M y_t = L y + R(y),
 * where M is the block mass matrix for each component, L is 
 * the block Laplace operator for each component, and R(y) is 
 * a 3x3 block comprised of the nonlinear reaction terms for 
 * each component.  Since it it highly inefficient to rewrite 
 * this system as
 *      y_t = M^{-1}(L y + R(y)),
 * we solve this system using ARKode, with a user-supplied mass
 * matrix.  We therefore provide functions to evaluate the ODE RHS 
 *    f(t,y) = L y + R(y),
 * its Jacobian
 *    J(t,y) = L + dR/dy,
 * and the mass matrix, M.
 *
 * This program solves the problem with the DIRK method, using a
 * Newton iteration with the ARKSUPERLUMT sparse linear solver.
 *
 * 100 outputs are printed at equal time intervals, and run 
 * statistics are printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>             /* prototypes for ARKode fcts., consts. */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <arkode/arkode_superlumt.h>   /* prototype for ARKSuperLUMT solver */
#include <sundials/sundials_types.h>   /* def. of type 'realtype' */
#include <sundials/sundials_sparse.h>  /* defs. of SlsMat and SlsMat routines */


/* accessor macros between (x,v) location and 1D NVector array */
/* [variables are grouped according to spatial location] */
#define IDX(x,v) (3*(x)+v)

/* constants */
#define ZERO (RCONST(0.0))
#define ONE  (RCONST(1.0))
#define TWO  (RCONST(2.0))
#define HALF (RCONST(0.5))

/* Gaussian quadrature nodes, weights and formula (3 node, 7th-order accurate) */
#define X1(xl,xr)   (HALF*(xl+xr) - HALF*(xr-xl)*RCONST(0.774596669241483377035853079956))
#define X2(xl,xr)   (HALF*(xl+xr))
#define X3(xl,xr)   (HALF*(xl+xr) + HALF*(xr-xl)*RCONST(0.774596669241483377035853079956))
#define W1          (RCONST(0.55555555555555555555555555555556))
#define W2          (RCONST(0.88888888888888888888888888888889))
#define W3          (RCONST(0.55555555555555555555555555555556))
#define Quad(f1,f2,f3,xl,xr) (HALF*(xr-xl)*(W1*f1 + W2*f2 + W3*f3))

/* evaluation macros for variables, basis functions and basis derivatives */
#define ChiL(xl,xr,x) ((xr-x)/(xr-xl))
#define ChiR(xl,xr,x) ((x-xl)/(xr-xl))
#define ChiL_x(xl,xr) (ONE/(xl-xr))
#define ChiR_x(xl,xr) (ONE/(xr-xl))
#define Eval(ul,ur,xl,xr,x) (ul*ChiL(xl,xr,x) + ur*ChiR(xl,xr,x))
#define Eval_x(ul,ur,xl,xr) (ul*ChiL_x(xl,xr) + ur*ChiR_x(xl,xr))


/* user data structure */
typedef struct {  
  int N;         /* number of intervals     */
  realtype *x;   /* mesh node locations     */
  realtype a;    /* constant forcing on u   */
  realtype b;    /* steady-state value of w */
  realtype du;   /* diffusion coeff for u   */
  realtype dv;   /* diffusion coeff for v   */
  realtype dw;   /* diffusion coeff for w   */
  realtype ep;   /* stiffness parameter     */
  N_Vector tmp;  /* temporary vector        */
  SlsMat R;      /* temporary storage       */
} *UserData;


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f_diff(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f_rx(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int MassMatrix(realtype t, SlsMat M, void *user_data, 
		      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jac(realtype t, N_Vector y, N_Vector fy, 
               SlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private helper functions  */
static int LaplaceMatrix(SlsMat Jac, UserData udata);
static int ReactionJac(N_Vector y, SlsMat Jac, UserData udata);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);

/* Main Program */
int main(int argc, char *argv[]) {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);    /* initial time */
  realtype Tf = RCONST(10.0);   /* final time */
  int Nt = 100;                 /* total number of output times */
  int Nvar = 3;                 /* number of solution fields */
  UserData udata = NULL;
  realtype *data;
  int N = 201;                  /* spatial mesh size */
  realtype a = 0.6;             /* problem parameters */
  realtype b = 2.0;
  realtype du = 0.025;
  realtype dv = 0.025;
  realtype dw = 0.025;
  realtype ep = 1.0e-5;         /* stiffness parameter */
  realtype reltol = 1.0e-6;     /* tolerances */
  realtype abstol = 1.0e-10;
  int i, NEQ, NNZ, num_threads;

  /* general problem variables */
  int flag;                     /* reusable error-checking flag */
  N_Vector y = NULL;
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  void *arkode_mem = NULL;
  FILE *FID, *UFID, *VFID, *WFID;
  realtype h, z, t, dTout, tout, u, v, w, pi;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nni, ncfn;
  long int netf, nms, nMv;

  /* if a command-line argument was supplied, set num_threads */
  num_threads = 1;
  if (argc > 1)
    num_threads = strtol(argv[1], NULL, 0);

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->x = NULL;
  udata->tmp = NULL;
  udata->R = NULL;
  if (check_flag((void *)udata, "malloc", 2)) return 1;

  /* store the inputs in the UserData structure */
  udata->N  = N;
  udata->a  = a;
  udata->b  = b;
  udata->du = du;
  udata->dv = dv;
  udata->dw = dw;
  udata->ep = ep;

  /* set total allocated vector length (N-1 intervals, Dirichlet end points) */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D FEM Brusselator PDE test problem:\n");
  printf("    N = %i,  NEQ = %i\n", udata->N, NEQ);
  printf("    num_threads = %i\n", num_threads);
  printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",
	 udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n", 
	 udata->du, udata->dv, udata->dw);
  printf("    reltol = %.1e,  abstol = %.1e\n\n", reltol, abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ);           /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  data = N_VGetArrayPointer(y);     /* Access data array for new NVector y */
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  umask = N_VNew_Serial(NEQ);       /* Create serial vector masks */
  if (check_flag((void *)umask, "N_VNew_Serial", 0)) return 1;
  vmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)vmask, "N_VNew_Serial", 0)) return 1;
  wmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)wmask, "N_VNew_Serial", 0)) return 1;
  udata->tmp = N_VNew_Serial(NEQ);  /* temporary N_Vector inside udata */
  if (check_flag((void *) udata->tmp, "N_VNew_Serial", 0)) return 1;

  /* allocate and set up spatial mesh; this [arbitrarily] clusters 
     more intervals near the end points of the interval */
  udata->x = (realtype *) malloc(N*sizeof(realtype));
  if (check_flag((void *)udata->x, "malloc", 2)) return 1;
  h = 10.0/(N-1);
  for (i=0; i<N; i++) {
    z = -5.0 + h*i;
    udata->x[i] = 0.5/atan(5.0)*atan(z) + 0.5;
  }

  /* Set initial conditions into y */
  pi = RCONST(4.0)*atan(RCONST(1.0));
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*udata->x[i]);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*udata->x[i]);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*udata->x[i]);  /* w */
  }

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  data = N_VGetArrayPointer(umask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,0)] = ONE;

  N_VConst(0.0, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,1)] = ONE;

  N_VConst(0.0, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,2)] = ONE;

  
  /* Create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;

  /* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem, (void *) udata);     /* Pass udata to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);    /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;
  flag = ARKodeResStolerance(arkode_mem, abstol);           /* Specify residual tolerance */
  if (check_flag(&flag, "ARKodeResStolerance", 1)) return 1;

  /* Linear solver specification */
  NNZ = 15*NEQ;
  flag = ARKSuperLUMT(arkode_mem, num_threads, NEQ, NNZ);   /* SuperLUMT sparse solver */
  if (check_flag(&flag, "ARKSuperLUMT", 1)) return 1;
  flag = ARKSlsSetSparseJacFn(arkode_mem, Jac);             /* Set the Jacobian routine */
  if (check_flag(&flag, "ARKSlsSetSparseJacFn", 1)) return 1;
  flag = ARKMassSuperLUMT(arkode_mem, num_threads,          /* Mass matrix linear solver */
			  NEQ, NNZ, MassMatrix);
  if (check_flag(&flag, "ARKMassSuperLUMT", 1)) return 1;

  /* output mesh to disk */
  FID=fopen("bruss_FEM_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->x[i]);
  fclose(FID);

  /* Open output stream for results, access data arrays */
  UFID = fopen("bruss_FEM_u.txt","w");
  VFID = fopen("bruss_FEM_v.txt","w");
  WFID = fopen("bruss_FEM_w.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;

  /* output initial condition to disk */
  for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
  for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
  for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
  fprintf(UFID,"\n");
  fprintf(VFID,"\n");
  fprintf(WFID,"\n");

  /* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t  = T0;
  dTout = Tf/Nt;
  tout = T0+dTout;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms\n");
  printf("   ----------------------------------------------\n");
  for (iout=0; iout<Nt; iout++) {

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);    /* call integrator */
    if (check_flag(&flag, "ARKode", 1)) break;
    u = N_VWL2Norm(y,umask);                               /* access/print solution statistics */
    u = sqrt(u*u/N);
    v = N_VWL2Norm(y,vmask);
    v = sqrt(v*v/N);
    w = N_VWL2Norm(y,wmask);
    w = sqrt(w*w/N);
    printf("  %10.6f  %10.6f  %10.6f  %10.6f\n", t, u, v, w);
    if (flag >= 0) {                                       /* successful solve: update output time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                               /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }

    /* output results to disk */
    for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
    for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
    for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
    fprintf(UFID,"\n");
    fprintf(VFID,"\n");
    fprintf(WFID,"\n");
  }
  printf("   ----------------------------------------------\n");
  fclose(UFID);
  fclose(VFID);
  fclose(WFID);

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
  flag = ARKodeGetNumMassSolves(arkode_mem, &nms);
  check_flag(&flag, "ARKodeGetNumMassSolves", 1);
  flag = ARKodeGetNumMassMultiplies(arkode_mem, &nMv);
  check_flag(&flag, "ARKodeGetNumMassMultiplies", 1);
  flag = ARKSlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKSlsGetNumJacEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total mass matrix solves = %li\n", nms);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy_Serial(y);            /* Free vectors */
  N_VDestroy_Serial(umask);
  N_VDestroy_Serial(vmask);
  N_VDestroy_Serial(wmask);
  DestroySparseMat(udata->R);      /* Free user data */
  N_VDestroy_Serial(udata->tmp);
  free(udata->x);
  free(udata);
  ARKodeFree(&arkode_mem);         /* Free integrator memory */
  return 0;
}


/*------------------------------
  Functions called by the solver
 *------------------------------*/


/* Routine to compute the ODE RHS function f(t,y), where system is of the form
        M y_t = f(t,y) := Ly + R(y) 
   This routine only computes the f(t,y), leaving (M y_t) alone. */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* local data */
  int ier;

  /* clear out RHS (to be careful) */
  N_VConst(0.0, ydot);

  /* add reaction terms to RHS */
  ier = f_rx(t, y, ydot, user_data);
  if (ier != 0)  return ier;
  
  /* add diffusion terms to RHS */
  ier = f_diff(t, y, ydot, user_data);
  if (ier != 0)  return ier;
  
  return 0;
}


/* Routine to compute the diffusion portion of the ODE RHS function f(t,y). */
static int f_diff(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  int N = udata->N;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;

  /* local variables */
  long int i;
  realtype ul, ur, vl, vr, wl, wr;
  realtype xl, xr, f1;
  booleantype left, right;
  realtype *Ydata, *RHSdata;
  
  /* access data arrays */
  Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  RHSdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *)RHSdata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over intervals, filling in residual function */
  for (i=0; i<N-1; i++) {

    /* set booleans to determine whether equations exist on the left/right */
    left  = (i==0)     ? FALSE : TRUE;
    right = (i==(N-2)) ? FALSE : TRUE;

    /* set nodal value shortcuts (interval index aligns with left node) */
    ul = Ydata[IDX(i,0)];
    vl = Ydata[IDX(i,1)];
    wl = Ydata[IDX(i,2)];
    ur = Ydata[IDX(i+1,0)];
    vr = Ydata[IDX(i+1,1)];
    wr = Ydata[IDX(i+1,2)];

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* evaluate L*y on this subinterval
       NOTE: all f values are the same since constant on interval */
    /*    left test function */
    if (left) {
      /*  u */
      f1 = -du * Eval_x(ul,ur,xl,xr) * ChiL_x(xl,xr);
      RHSdata[IDX(i,0)] += Quad(f1,f1,f1,xl,xr);

      /*  v */
      f1 = -dv * Eval_x(vl,vr,xl,xr) * ChiL_x(xl,xr);
      RHSdata[IDX(i,1)] += Quad(f1,f1,f1,xl,xr);
      
      /*  w */
      f1 = -dw * Eval_x(wl,wr,xl,xr) * ChiL_x(xl,xr);
      RHSdata[IDX(i,2)] += Quad(f1,f1,f1,xl,xr);
    }
    /*    right test function */
    if (right) {
      /*  u */
      f1 = -du * Eval_x(ul,ur,xl,xr) * ChiR_x(xl,xr);
      RHSdata[IDX(i+1,0)] += Quad(f1,f1,f1,xl,xr);

      /*  v */
      f1 = -dv * Eval_x(vl,vr,xl,xr) * ChiR_x(xl,xr);
      RHSdata[IDX(i+1,1)] += Quad(f1,f1,f1,xl,xr);

      /*  w */
      f1 = -dw * Eval_x(wl,wr,xl,xr) * ChiR_x(xl,xr);
      RHSdata[IDX(i+1,2)] += Quad(f1,f1,f1,xl,xr);
    }
  }

  return 0;
}



/* Routine to compute the reaction portion of the ODE RHS function f(t,y). */
static int f_rx(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  int N = udata->N;
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;

  /* local variables */
  long int i;
  realtype ul, ur, vl, vr, wl, wr;
  realtype u, v, w, xl, xr, f1, f2, f3;
  booleantype left, right;
  realtype *Ydata, *RHSdata;

  /* access data arrays */
  Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  RHSdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *)RHSdata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over intervals, filling in residual function */
  for (i=0; i<N-1; i++) {

    /* set booleans to determine whether equations exist on the left/right */
    left  = (i==0)     ? FALSE : TRUE;
    right = (i==(N-2)) ? FALSE : TRUE;

    /* set nodal value shortcuts (interval index aligns with left node) */
    ul = Ydata[IDX(i,0)];
    vl = Ydata[IDX(i,1)];
    wl = Ydata[IDX(i,2)];
    ur = Ydata[IDX(i+1,0)];
    vr = Ydata[IDX(i+1,1)];
    wr = Ydata[IDX(i+1,2)];

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* evaluate R(y) on this subinterval */
    /*    left test function */
    if (left) {
      /*  u */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (a - (w+ONE)*u + v*u*u) * ChiL(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (a - (w+ONE)*u + v*u*u) * ChiL(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (a - (w+ONE)*u + v*u*u) * ChiL(xl,xr,X3(xl,xr));
      RHSdata[IDX(i,0)] += Quad(f1,f2,f3,xl,xr);
    
      /*  v */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (w*u - v*u*u) * ChiL(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (w*u - v*u*u) * ChiL(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (w*u - v*u*u) * ChiL(xl,xr,X3(xl,xr));
      RHSdata[IDX(i,1)] += Quad(f1,f2,f3,xl,xr);
    
      /*  w */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = ((b-w)/ep - w*u) * ChiL(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = ((b-w)/ep - w*u) * ChiL(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = ((b-w)/ep - w*u) * ChiL(xl,xr,X3(xl,xr));
      RHSdata[IDX(i,2)] += Quad(f1,f2,f3,xl,xr);
    }
    /*    right test function */
    if (right) {
      /*  u */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (a - (w+ONE)*u + v*u*u) * ChiR(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (a - (w+ONE)*u + v*u*u) * ChiR(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (a - (w+ONE)*u + v*u*u) * ChiR(xl,xr,X3(xl,xr));
      RHSdata[IDX(i+1,0)] += Quad(f1,f2,f3,xl,xr);
    
      /*  v */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (w*u - v*u*u) * ChiR(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (w*u - v*u*u) * ChiR(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (w*u - v*u*u) * ChiR(xl,xr,X3(xl,xr));
      RHSdata[IDX(i+1,1)] += Quad(f1,f2,f3,xl,xr);
    
      /*  w */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = ((b-w)/ep - w*u) * ChiR(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = ((b-w)/ep - w*u) * ChiR(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = ((b-w)/ep - w*u) * ChiR(xl,xr,X3(xl,xr));
      RHSdata[IDX(i+1,2)] += Quad(f1,f2,f3,xl,xr);
    }
  }

  return 0;
}



/* Interface routine to compute the Jacobian of the full RHS function, f(y) */
static int Jac(realtype t, N_Vector y, N_Vector fy, 
               SlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  /* temporary variables */
  int ier;
  UserData udata = (UserData) user_data;
  int N = udata->N;

  /* ensure that Jac is the correct size */
  if ((J->M != N*3) || (J->N != N*3)) {
    printf("Jacobian calculation error: matrix is the wrong size!\n");
    return 1;
  }
  
  /* Fill in the Laplace matrix */
  ier = LaplaceMatrix(J, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling Laplace matrix = %i\n",ier);
    return 1;
  }

  /* Create empty reaction Jacobian matrix (if not done already) */
  if (udata->R == NULL) {
    udata->R = NewSparseMat(J->M, J->N, J->NNZ);
    if (udata->R == NULL) {
      printf("Jac: error in allocating R matrix!\n");
      return 1;
    }
  }
      
  /* Add in the Jacobian of the reaction terms matrix */
  ier = ReactionJac(y, udata->R, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling reaction Jacobian = %i\n",ier);
    return 1;
  }

  /* Add R to J */
  ier = SlsAddMat(J,udata->R);
  if (ier != 0) {
    printf("Jac: error in adding sparse matrices = %i!\n",ier);
    return 1;
  }

  return 0;
}



/* Routine to compute the mass matrix multiplying y_t. */
static int MassMatrix(realtype t, SlsMat M, void *user_data, 
		      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  /* user data structure */
  UserData udata = (UserData) user_data;

  /* set shortcuts */
  int N = udata->N;
  int i, nz=0;

  /* local data */
  realtype xl, xr, f1, f2, f3, dtmp;

  /* clear out mass matrix */
  SlsSetToZero(M);

  /* iterate over columns, filling in matrix entries */
  for (i=0; i<N; i++) {

    /* dependence on u at this node */
    M->colptrs[IDX(i,0)] = nz;

    /*    left u trial function */
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i-1,0);
    }
    /*    this u trial function */
    dtmp = ZERO;
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiR(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiR(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiR(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    M->data[nz] = dtmp;
    M->rowvals[nz++] = IDX(i,0);
    /*    right u trial function */
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i+1,0);
    }


    /* dependence on v at this node */
    M->colptrs[IDX(i,1)] = nz;

    /*    left v trial function */
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i-1,1);
    }
    /*    this v trial function */
    dtmp = ZERO;
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiR(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiR(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiR(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    M->data[nz] = dtmp;
    M->rowvals[nz++] = IDX(i,1);
    /*    right v trial function */
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i+1,1);
    }


    /* dependence on w at this node */
    M->colptrs[IDX(i,2)] = nz;

    /*    left w trial function */
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i-1,2);
    }
    /*    this w trial function */
    dtmp = ZERO;
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiR(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiR(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiR(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    M->data[nz] = dtmp;
    M->rowvals[nz++] = IDX(i,2);
    /*    right w trial function */
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i+1,2);
    }

  }

  /* signal end of data */
  M->colptrs[IDX(N-1,2)+1] = nz;

  return 0;
}





/*-------------------------------
 * Private helper functions
 *-------------------------------*/



/* Routine to compute the Laplace matrix */
static int LaplaceMatrix(SlsMat L, UserData udata)
{

  /* set shortcuts, local variables */
  int N = udata->N;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  int i, nz=0;
  realtype xl, xr;
  
  /* clear out matrix */
  SlsSetToZero(L);

  /* iterate over columns, filling in Laplace matrix entries */
  for (i=0; i<N; i++) {

    /* dependence on u at this node */
    L->colptrs[IDX(i,0)] = nz;

    if (i>1) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] = (-du) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i-1,0);
    }
    if (i<N-1 && i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] = (-du) * Quad(ONE,ONE,ONE,xl,xr) * ChiR_x(xl,xr) * ChiR_x(xl,xr);
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] += (-du) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiL_x(xl,xr);
      L->rowvals[nz++] = IDX(i,0);
    }
    if (i<N-2) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-du) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i+1,0);
    }

    /* dependence on v at this node */
    L->colptrs[IDX(i,1)] = nz;

    if (i>1) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] = (-dv) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i-1,1);
    }
    if (i>0 && i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-dv) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiL_x(xl,xr);
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] += (-dv) * Quad(ONE,ONE,ONE,xl,xr) * ChiR_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i,1);
    }
    if (i<N-2) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-dv) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i+1,1);
    }

    /* dependence on w at this node */
    L->colptrs[IDX(i,2)] = nz;

    if (i>1) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] = (-dw) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i-1,2);
    }
    if (i>0 && i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-dw) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiL_x(xl,xr);
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] += (-dw) * Quad(ONE,ONE,ONE,xl,xr) * ChiR_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i,2);
    }
    if (i<N-2) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-dw) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i+1,2);
    }

  }

  /* signal end of data */
  L->colptrs[IDX(N-1,2)+1] = nz;

  return 0;
}



/* Routine to compute the Jacobian matrix from R(y) */
static int ReactionJac(N_Vector y, SlsMat Jac, UserData udata)
{
  /* set shortcuts, local variables */
  int N = udata->N;
  int i, nz=0;
  realtype ep = udata->ep;
  realtype ul, uc, ur, vl, vc, vr, wl, wc, wr;
  realtype u1l, u2l, u3l, v1l, v2l, v3l, w1l, w2l, w3l;
  realtype u1r, u2r, u3r, v1r, v2r, v3r, w1r, w2r, w3r;
  realtype xl, xc, xr, df1, df2, df3;
  realtype dQdf1l, dQdf2l, dQdf3l, ChiL1l, ChiL2l, ChiL3l, ChiR1l, ChiR2l, ChiR3l;
  realtype dQdf1r, dQdf2r, dQdf3r, ChiL1r, ChiL2r, ChiL3r, ChiR1r, ChiR2r, ChiR3r;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;

  /* initialize all local variables to zero (to avoid uninitialized variable warnings) */
  ul = uc = ur = vl = vc = vr = wl = wc = wr = 0.0;
  u1l = u2l = u3l = v1l = v2l = v3l = w1l = w2l = w3l = 0.0;
  u1r = u2r = u3r = v1r = v2r = v3r = w1r = w2r = w3r = 0.0;
  xl = xc = xr = df1 = df2 = df3 = 0.0;
  dQdf1l = dQdf2l = dQdf3l = ChiL1l = ChiL2l = ChiL3l = ChiR1l = ChiR2l = ChiR3l = 0.0;
  dQdf1r = dQdf2r = dQdf3r = ChiL1r = ChiL2r = ChiL3r = ChiR1r = ChiR2r = ChiR3r = 0.0;

  /* clear out matrix */
  SlsSetToZero(Jac);

  /* iterate over columns, filling in reaction Jacobian */
  for (i=0; i<N; i++) {

    /* set mesh shortcuts */
    if (i>0)
      xl = udata->x[i-1];
    xc = udata->x[i];
    if (i<N-1)
      xr = udata->x[i+1];

    /* set nodal value shortcuts */
    if (i>0) {
      ul = Ydata[IDX(i-1,0)];
      vl = Ydata[IDX(i-1,1)];
      wl = Ydata[IDX(i-1,2)];
    }
    uc = Ydata[IDX(i,0)];
    vc = Ydata[IDX(i,1)];
    wc = Ydata[IDX(i,2)];
    if (i<N-1) {
      ur = Ydata[IDX(i+1,0)];
      vr = Ydata[IDX(i+1,1)];
      wr = Ydata[IDX(i+1,2)];
    }
    if (i>0) {
      u1l = Eval(ul,uc,xl,xc,X1(xl,xc));
      v1l = Eval(vl,vc,xl,xc,X1(xl,xc));
      w1l = Eval(wl,wc,xl,xc,X1(xl,xc));
      u2l = Eval(ul,uc,xl,xc,X2(xl,xc));
      v2l = Eval(vl,vc,xl,xc,X2(xl,xc));
      w2l = Eval(wl,wc,xl,xc,X2(xl,xc));
      u3l = Eval(ul,uc,xl,xc,X3(xl,xc));
      v3l = Eval(vl,vc,xl,xc,X3(xl,xc));
      w3l = Eval(wl,wc,xl,xc,X3(xl,xc));
    }
    if (i<N-1) {
      u1r = Eval(uc,ur,xc,xr,X1(xc,xr));
      v1r = Eval(vc,vr,xc,xr,X1(xc,xr));
      w1r = Eval(wc,wr,xc,xr,X1(xc,xr));
      u2r = Eval(uc,ur,xc,xr,X2(xc,xr));
      v2r = Eval(vc,vr,xc,xr,X2(xc,xr));
      w2r = Eval(wc,wr,xc,xr,X2(xc,xr));
      u3r = Eval(uc,ur,xc,xr,X3(xc,xr));
      v3r = Eval(vc,vr,xc,xr,X3(xc,xr));
      w3r = Eval(wc,wr,xc,xr,X3(xc,xr));
    }

    /* set partial derivative shortcuts */
    if (i>0) {
      dQdf1l = Quad(ONE, ZERO, ZERO, xl, xc);
      dQdf2l = Quad(ZERO, ONE, ZERO, xl, xc);
      dQdf3l = Quad(ZERO, ZERO, ONE, xl, xc);
      ChiL1l = ChiL(xl,xc,X1(xl,xc));
      ChiL2l = ChiL(xl,xc,X2(xl,xc));
      ChiL3l = ChiL(xl,xc,X3(xl,xc));
      ChiR1l = ChiR(xl,xc,X1(xl,xc));
      ChiR2l = ChiR(xl,xc,X2(xl,xc));
      ChiR3l = ChiR(xl,xc,X3(xl,xc));
    }
    if (i<N-1) {
      dQdf1r = Quad(ONE, ZERO, ZERO, xc, xr);
      dQdf2r = Quad(ZERO, ONE, ZERO, xc, xr);
      dQdf3r = Quad(ZERO, ZERO, ONE, xc, xr);
      ChiL1r = ChiL(xc,xr,X1(xc,xr));
      ChiL2r = ChiL(xc,xr,X2(xc,xr));
      ChiL3r = ChiL(xc,xr,X3(xc,xr));
      ChiR1r = ChiR(xc,xr,X1(xc,xr));
      ChiR2r = ChiR(xc,xr,X2(xc,xr));
      ChiR3r = ChiR(xc,xr,X3(xc,xr));
    }


    /*** evaluate dR/dy at this node ***/


    /* dependence on u at this node */
    Jac->colptrs[IDX(i,0)] = nz;

    if (i>1) {
      /*  dR_ul/duc */
      df1 = (-(w1l+ONE) + TWO*v1l*u1l) * ChiL1l * ChiR1l;
      df2 = (-(w2l+ONE) + TWO*v2l*u2l) * ChiL2l * ChiR2l;
      df3 = (-(w3l+ONE) + TWO*v3l*u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,0);

      /*  dR_vl/duc */
      df1 = (w1l - TWO*v1l*u1l) * ChiL1l * ChiR1l;
      df2 = (w2l - TWO*v2l*u2l) * ChiL2l * ChiR2l;
      df3 = (w3l - TWO*v3l*u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,1);

      /*  dR_wl/duc */
      df1 = (-w1l) * ChiL1l * ChiR1l;
      df2 = (-w2l) * ChiL2l * ChiR2l;
      df3 = (-w3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,2);
    }
    if (i>0 && i<N-1) {
      /*  dR_uc/duc */
      df1 = (-(w1r+ONE) + TWO*v1r*u1r) * ChiL1r * ChiL1r;
      df2 = (-(w2r+ONE) + TWO*v2r*u2r) * ChiL2r * ChiL2r;
      df3 = (-(w3r+ONE) + TWO*v3r*u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;

      df1 = (-(w1l+ONE) + TWO*v1l*u1l) * ChiR1l * ChiR1l;
      df2 = (-(w2l+ONE) + TWO*v2l*u2l) * ChiR2l * ChiR2l;
      df3 = (-(w3l+ONE) + TWO*v3l*u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] += dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i,0);

      /*  dR_vc/duc */
      df1 = (w1l - TWO*v1l*u1l) * ChiR1l * ChiR1l;
      df2 = (w2l - TWO*v2l*u2l) * ChiR2l * ChiR2l;
      df3 = (w3l - TWO*v3l*u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (w1r - TWO*v1r*u1r) * ChiL1r * ChiL1r;
      df2 = (w2r - TWO*v2r*u2r) * ChiL2r * ChiL2r;
      df3 = (w3r - TWO*v3r*u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,1);

      /*  dR_wc/duc */
      df1 = (-w1r) * ChiL1r * ChiL1r;
      df2 = (-w2r) * ChiL2r * ChiL2r;
      df3 = (-w3r) * ChiL3r * ChiL3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;

      df1 = (-w1l) * ChiR1l * ChiR1l;
      df2 = (-w2l) * ChiR2l * ChiR2l;
      df3 = (-w3l) * ChiR3l * ChiR3l;
      Jac->data[nz] += dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i,2);
    }
    if (i<N-2) {
      /*  dR_ur/duc */
      df1 = (-(w1r+ONE) + TWO*v1r*u1r) * ChiL1r * ChiR1r;
      df2 = (-(w2r+ONE) + TWO*v2r*u2r) * ChiL2r * ChiR2r;
      df3 = (-(w3r+ONE) + TWO*v3r*u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,0);

      /*  dR_vr/duc */
      df1 = (w1r - TWO*v1r*u1r) * ChiL1r * ChiR1r;
      df2 = (w2r - TWO*v2r*u2r) * ChiL2r * ChiR2r;
      df3 = (w3r - TWO*v3r*u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,1);

      /*  dR_wr/duc */
      df1 = (-w1r) * ChiL1r * ChiR1r;
      df2 = (-w2r) * ChiL2r * ChiR2r;
      df3 = (-w3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,2);
    }


    /* dependence on v at this node */
    Jac->colptrs[IDX(i,1)] = nz;

    if (i>1) {
      /*  dR_ul/dvc */
      df1 = (u1l*u1l) * ChiL1l * ChiR1l;
      df2 = (u2l*u2l) * ChiL2l * ChiR2l;
      df3 = (u3l*u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,0);

      /*  dR_vl/dvc */
      df1 = (-u1l*u1l) * ChiL1l * ChiR1l;
      df2 = (-u2l*u2l) * ChiL2l * ChiR2l;
      df3 = (-u3l*u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,1);
    }
    if (i>0 && i<N-1) {
      /*  dR_uc/dvc */
      df1 = (u1l*u1l) * ChiR1l * ChiR1l;
      df2 = (u2l*u2l) * ChiR2l * ChiR2l;
      df3 = (u3l*u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (u1r*u1r) * ChiL1r * ChiL1r;
      df2 = (u2r*u2r) * ChiL2r * ChiL2r;
      df3 = (u3r*u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,0);

      /*  dR_vc/dvc */
      df1 = (-u1l*u1l) * ChiR1l * ChiR1l;
      df2 = (-u2l*u2l) * ChiR2l * ChiR2l;
      df3 = (-u3l*u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (-u1r*u1r) * ChiL1r * ChiL1r;
      df2 = (-u2r*u2r) * ChiL2r * ChiL2r;
      df3 = (-u3r*u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,1);
    }
    if (i<N-2) {
      /*  dR_ur/dvc */
      df1 = (u1r*u1r) * ChiL1r * ChiR1r;
      df2 = (u2r*u2r) * ChiL2r * ChiR2r;
      df3 = (u3r*u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,0);

      /*  dR_vr/dvc */
      df1 = (-u1r*u1r) * ChiL1r * ChiR1r;
      df2 = (-u2r*u2r) * ChiL2r * ChiR2r;
      df3 = (-u3r*u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,1);
    }


    /* dependence on w at this node */
    Jac->colptrs[IDX(i,2)] = nz;

    if (i>1) {
      /*  dR_ul/dwc */
      df1 = (-u1l) * ChiL1l * ChiR1l;
      df2 = (-u2l) * ChiL2l * ChiR2l;
      df3 = (-u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,0);

      /*  dR_vl/dwc */
      df1 = (u1l) * ChiL1l * ChiR1l;
      df2 = (u2l) * ChiL2l * ChiR2l;
      df3 = (u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,1);

      /*  dR_wl/dwc */
      df1 = (-ONE/ep - u1l) * ChiL1l * ChiR1l;
      df2 = (-ONE/ep - u2l) * ChiL2l * ChiR2l;
      df3 = (-ONE/ep - u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,2);
    }
    if (i>0 && i<N-1) {
      /*  dR_uc/dwc */
      df1 = (-u1l) * ChiR1l * ChiR1l;
      df2 = (-u2l) * ChiR2l * ChiR2l;
      df3 = (-u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (-u1r) * ChiL1r * ChiL1r;
      df2 = (-u2r) * ChiL2r * ChiL2r;
      df3 = (-u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,0);

      /*  dR_vc/dwc */
      df1 = (u1l) * ChiR1l * ChiR1l;
      df2 = (u2l) * ChiR2l * ChiR2l;
      df3 = (u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (u1r) * ChiL1r * ChiL1r;
      df2 = (u2r) * ChiL2r * ChiL2r;
      df3 = (u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,1);

      /*  dR_wc/dwc */
      df1 = (-ONE/ep - u1l) * ChiR1l * ChiR1l;
      df2 = (-ONE/ep - u2l) * ChiR2l * ChiR2l;
      df3 = (-ONE/ep - u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (-ONE/ep - u1r) * ChiL1r * ChiL1r;
      df2 = (-ONE/ep - u2r) * ChiL2r * ChiL2r;
      df3 = (-ONE/ep - u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,2);
    }
    if (i<N-2) {
      /*  dR_ur/dwc */
      df1 = (-u1r) * ChiL1r * ChiR1r;
      df2 = (-u2r) * ChiL2r * ChiR2r;
      df3 = (-u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,0);

      /*  dR_vr/dwc */
      df1 = (u1r) * ChiL1r * ChiR1r;
      df2 = (u2r) * ChiL2r * ChiR2r;
      df3 = (u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,1);

      /*  dR_wr/dwc */
      df1 = (-ONE/ep - u1r) * ChiL1r * ChiR1r;
      df2 = (-ONE/ep - u2r) * ChiL2r * ChiR2r;
      df3 = (-ONE/ep - u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,2);
    }

  }

  /* signal end of data */
  Jac->colptrs[IDX(N-1,2)+1] = nz;

  return 0;
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
