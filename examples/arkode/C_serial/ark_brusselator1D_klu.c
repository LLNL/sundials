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
 * kinetics.  This is n PDE system with 3 components, Y = [u,v,w], 
 * satisfying the equations,
 *    u_t = du*u_xx + a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx + w*u - v*u^2
 *    w_t = dw*w_xx + (b-w)/ep - w*u
 * for t in [0, 80], x in [0, 1], with initial conditions
 *    u(0,x) =  a  + 0.1*sin(pi*x)
 *    v(0,x) = b/a + 0.1*sin(pi*x)
 *    w(0,x) =  b  + 0.1*sin(pi*x),
 * and with stationary boundary conditions, i.e. 
 *    u_t(t,0) = u_t(t,1) = 0,
 *    v_t(t,0) = v_t(t,1) = 0,
 *    w_t(t,0) = w_t(t,1) = 0.
 * Note: these can also be implemented as Dirichlet boundary 
 * conditions with values identical to the initial conditions.
 * 
 * The spatial derivatives are computed using second-order 
 * centered differences, with the data distributed over N points 
 * on a uniform spatial grid.
 *
 * The number of spatial points N, the parameters a, b, du, dv, 
 * dw and ep, as well as the desired relative and absolute solver 
 * tolerances, are provided in the input file 
 * input_brusselator1D.txt.
 * 
 * This program solves the problem with the DIRK method, using a
 * Newton iteration.  The inner linear systems are solved using 
 * the ARKKLU linear solver.
 *
 * 100 outputs are printed at equal intervals, and run statistics 
 * are printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>            /* prototypes for ARKode fcts., consts. */
#include <arkode/arkode_klu.h>        /* prototype for ARKKLU solver */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>  /* def. of type 'realtype' */

/* accessor macros between (x,v) location and 1D NVector array */
#define IDX(x,v) (3*(x)+v)

/* constants */
#define ONE (RCONST(1.0))
#define TWO (RCONST(2.0))

/* user data structure */
typedef struct {  
  int N;         /* number of intervals     */
  realtype dx;   /* mesh spacing            */
  realtype a;    /* constant forcing on u   */
  realtype b;    /* steady-state value of w */
  realtype du;   /* diffusion coeff for u   */
  realtype dv;   /* diffusion coeff for v   */
  realtype dw;   /* diffusion coeff for w   */
  realtype ep;   /* stiffness parameter     */
  SlsMat R;      /* temporary storage       */
} *UserData;


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy, 
	       SlsMat J, void *user_data, 
	       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int LaplaceMatrix(SlsMat Jac, UserData udata);
static int ReactionJac(N_Vector y, SlsMat Jac, UserData udata);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);


/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);    /* initial time */
  realtype Tf = RCONST(10.0);   /* final time */
  int Nt = 10;                  /* total number of output times */
  int Nvar = 3;
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
  int i;
  long int NEQ, NNZ;

  /* general problem variables */
  int flag;                     /* reusable error-checking flag */
  N_Vector y = NULL;
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  void *arkode_mem = NULL;
  realtype pi;
  FILE *FID, *UFID, *VFID, *WFID;
  realtype t = T0;
  realtype dTout = (Tf-T0)/Nt;
  realtype tout = T0+dTout;
  realtype u, v, w;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nni, ncfn, netf;

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  if (check_flag((void *) udata, "malloc", 2)) return 1;

  /* store the inputs in the UserData structure */
  udata->N  = N;
  udata->a  = a;
  udata->b  = b;
  udata->du = du;
  udata->dv = dv;
  udata->dw = dw;
  udata->ep = ep;
  udata->R  = NULL;

  /* set total allocated vector length */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D Brusselator PDE test problem (KLU solver):\n");
  printf("    N = %i,  NEQ = %li\n", udata->N, NEQ);
  printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",
	 udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n", 
	 udata->du, udata->dv, udata->dw);
  printf("    reltol = %.1e,  abstol = %.1e\n\n", reltol, abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ);           /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  udata->dx = RCONST(1.0)/(N-1);    /* set spatial mesh spacing */
  data = N_VGetArrayPointer(y);     /* Access data array for new NVector y */
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  umask = N_VNew_Serial(NEQ);       /* Create serial vector masks */
  if (check_flag((void *)umask, "N_VNew_Serial", 0)) return 1;
  vmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)vmask, "N_VNew_Serial", 0)) return 1;
  wmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)wmask, "N_VNew_Serial", 0)) return 1;

  /* Set initial conditions into y */
  pi = RCONST(4.0)*atan(ONE);
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*i*udata->dx);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*i*udata->dx);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*i*udata->dx);  /* w */
  }

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  data = N_VGetArrayPointer(umask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,0)] = ONE;

  N_VConst(0.0, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,1)] = ONE;

  N_VConst(0.0, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,2)] = ONE;


  /* Create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;
  
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

  /* Specify the KLU sparse linear solver and Jacobian function */
  NNZ = 5*NEQ;
  flag = ARKKLU(arkode_mem, NEQ, NNZ);
  if (check_flag(&flag, "ARKKLU", 1)) return 1;
  flag = ARKSlsSetSparseJacFn(arkode_mem, Jac);
  if (check_flag(&flag, "ARKSlsSetSparseJacFn", 1)) return 1;
 
   /* output spatial mesh to disk */
  FID = fopen("bruss_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
  fclose(FID);

  /* Open output stream for results, access data arrays */
  UFID=fopen("bruss_u.txt","w");
  VFID=fopen("bruss_v.txt","w");
  WFID=fopen("bruss_w.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;

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
    u = N_VWL2Norm(y,umask);
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
  flag = ARKSlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKSlsGetNumJacEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of nonlinear iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy_Serial(y);         /* Free vectors */
  N_VDestroy_Serial(umask);
  N_VDestroy_Serial(vmask);
  N_VDestroy_Serial(wmask);
  DestroySparseMat(udata->R);   /* Free user data */
  free(udata);
  ARKodeFree(&arkode_mem);
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;      /* access problem data */
  int N       = udata->N;                     /* set variable shortcuts */
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  realtype dx = udata->dx;
  realtype *Ydata=NULL, *dYdata=NULL;
  realtype uconst, vconst, wconst, u, ul, ur, v, vl, vr, w, wl, wr;
  int i;

  Ydata = N_VGetArrayPointer(y);     /* access data arrays */
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
  dYdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                        /* initialize ydot to zero */

  /* iterate over domain, computing all equations */
  uconst = du/dx/dx;
  vconst = dv/dx/dx;
  wconst = dw/dx/dx;
  for (i=1; i<N-1; i++) {

    /* set shortcuts */
    u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* u_t = du*u_xx + a - (w+1)*u + v*u^2 */
    dYdata[IDX(i,0)] = (ul - TWO*u + ur)*uconst + a - (w+ONE)*u + v*u*u;

    /* v_t = dv*v_xx + w*u - v*u^2 */
    dYdata[IDX(i,1)] = (vl - TWO*v + vr)*vconst + w*u - v*u*u;

    /* w_t = dw*w_xx + (b-w)/ep - w*u */
    dYdata[IDX(i,2)] = (wl - TWO*w + wr)*wconst + (b-w)/ep - w*u;

  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  return 0;
}


/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(realtype t, N_Vector y, N_Vector fy, 
	       SlsMat J, void *user_data, 
	       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  /* problem data */
  UserData udata = (UserData) user_data;
  int N = udata->N;

  /* ensure that Jac is the correct size */
  if ((J->M != N*3) || (J->N != N*3)) {
    printf("Jacobian calculation error: matrix is the wrong size!\n");
    return 1;
  }
  
  /* Fill in the Laplace matrix */
  if (LaplaceMatrix(J, udata)) {
    printf("Jacobian calculation error in calling LaplaceMatrix!\n");
    return 1;
  }

  /* Add in the Jacobian of the reaction terms matrix */
  if (udata->R == NULL) {
    udata->R = NewSparseMat(J->M, J->N, J->NNZ);
    if (udata->R == NULL) {
      printf("Jacobian calculation error in allocating R matrix!\n");
      return 1;
    }
  }
      
  /* Add in the Jacobian of the reaction terms matrix */
  if (ReactionJac(y, udata->R, udata)) {
    printf("Jacobian calculation error in calling ReactionJac!\n");
    return 1;
  }

  /* Add R to J */
  if (SlsAddMat(J,udata->R) != 0) {
    printf("Jacobian calculation error in adding sparse matrices!\n");
    return 1;
  }

  return 0;
}




/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Routine to compute the stiffness matrix from (L*y) */
static int LaplaceMatrix(SlsMat Lap, UserData udata)
{
  int N = udata->N;  /* set shortcuts */
  int i, nz=0;
  realtype uconst, uconst2, vconst, vconst2, wconst, wconst2;

  /* clear out matrix */
  SlsSetToZero(Lap);

  /* set first column to zero */
  Lap->colptrs[IDX(0,0)] = nz;
  Lap->colptrs[IDX(0,1)] = nz;
  Lap->colptrs[IDX(0,2)] = nz;
  
  /* iterate over nodes, filling in Laplacian entries depending on these */
  uconst  = (udata->du)/(udata->dx)/(udata->dx);
  uconst2 = -TWO*uconst;
  vconst  = (udata->dv)/(udata->dx)/(udata->dx);
  vconst2 = -TWO*vconst;
  wconst  = (udata->dw)/(udata->dx)/(udata->dx);
  wconst2 = -TWO*wconst;
  for (i=1; i<N-1; i++) {

    /* dependence on u at this node */
    Lap->colptrs[IDX(i,0)] = nz;
    if (i>1) {                /* node to left */
      Lap->data[nz] = uconst;
      Lap->rowvals[nz++] = IDX(i-1,0);
    }

    Lap->data[nz] = uconst2;  /* self */
    Lap->rowvals[nz++] = IDX(i,0);

    if (i<N-2) {              /* node to right */
      Lap->data[nz] = uconst;
      Lap->rowvals[nz++] = IDX(i+1,0);
    }

    /* dependence on v at this node */
    Lap->colptrs[IDX(i,1)] = nz;
    if (i>1) {                /* node to left */
      Lap->data[nz] = vconst;
      Lap->rowvals[nz++] = IDX(i-1,1);
    }

    Lap->data[nz] = vconst2;  /* self */
    Lap->rowvals[nz++] = IDX(i,1);

    if (i<N-2) {              /* node to right */
      Lap->data[nz] = vconst;
      Lap->rowvals[nz++] = IDX(i+1,1);
    }

    /* dependence on w at this node */
    Lap->colptrs[IDX(i,2)] = nz;
    if (i>1) {                /* node to left */
      Lap->data[nz] = wconst;
      Lap->rowvals[nz++] = IDX(i-1,2);
    }

    Lap->data[nz] = wconst2;  /* self */
    Lap->rowvals[nz++] = IDX(i,2);

    if (i<N-2) {              /* node to right */
      Lap->data[nz] = wconst;
      Lap->rowvals[nz++] = IDX(i+1,2);
    }

  }

  /* set last column to zero */
  Lap->colptrs[IDX(N-1,0)] = nz;
  Lap->colptrs[IDX(N-1,1)] = nz;
  Lap->colptrs[IDX(N-1,2)] = nz;
  
  /* end of data */
  Lap->colptrs[IDX(N-1,2)+1] = nz;

  return 0;
}



/* Routine to compute the Jacobian matrix from R(y) */
static int ReactionJac(N_Vector y, SlsMat Jac, UserData udata)
{
  int N = udata->N;                            /* set shortcuts */
  int i, nz=0;
  realtype u, v, w;
  realtype ep = udata->ep;
  realtype *Ydata = N_VGetArrayPointer(y);     /* access solution array */
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;

  /* clear out matrix */
  SlsSetToZero(Jac);

  /* set first matrix column to zero */
  Jac->colptrs[IDX(0,0)] = 0;
  Jac->colptrs[IDX(0,1)] = 0;
  Jac->colptrs[IDX(0,2)] = 0;
  
  /* iterate over interior nodes, filling in Jacobian entries */
  for (i=1; i<N-1; i++) {

    /* set nodal value shortcuts */
    u = Ydata[IDX(i,0)];
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* dependence on u at this node */
    Jac->colptrs[IDX(i,0)] = nz;

    Jac->rowvals[nz] = IDX(i,0);        /* fu wrt u */
    Jac->data[nz++] = TWO*u*v - w - ONE;

    Jac->rowvals[nz] = IDX(i,1);        /* fv wrt u */
    Jac->data[nz++] = w - TWO*u*v;

    Jac->rowvals[nz] = IDX(i,2);        /* fw wrt u */
    Jac->data[nz++] = -w;

    /* dependence on v at this node */
    Jac->colptrs[IDX(i,1)] = nz;

    Jac->rowvals[nz] = IDX(i,0);        /* fu wrt v */
    Jac->data[nz++] = u*u;

    Jac->rowvals[nz] = IDX(i,1);        /* fv wrt v */
    Jac->data[nz++] = -u*u;

    /* dependence on w at this node */
    Jac->colptrs[IDX(i,2)] = nz;

    Jac->rowvals[nz] = IDX(i,0);        /* fu wrt w */
    Jac->data[nz++] = -u;

    Jac->rowvals[nz] = IDX(i,1);        /* fv wrt w */
    Jac->data[nz++] = u;

    Jac->rowvals[nz] = IDX(i,2);        /* fw wrt w */
    Jac->data[nz++] = -ONE/ep - u;

  }

  /* set last matrix column to zero */
  Jac->colptrs[IDX(N-1,0)] = nz;
  Jac->colptrs[IDX(N-1,1)] = nz;
  Jac->colptrs[IDX(N-1,2)] = nz;

  /* end of data */
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
