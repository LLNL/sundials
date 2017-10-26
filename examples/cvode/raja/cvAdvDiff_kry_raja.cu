/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This example is based on cvAdvDiff_bnd 
 *                   example by Scott D. Cohen, Alan C. 
 *                   Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem with a banded Jacobian,
 * with the program for its solution by CVODE.
 * The problem is the semi-discrete form of the advection-diffusion
 * equation in 2-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
 * interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
 * are posed, and the initial condition is
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX*MY.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVBAND band linear solver, and a user-supplied
 * Jacobian routine.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .1, .2, ..., 1.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvode/cvode.h>
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver        */
#include <cvode/cvode_spils.h>         /* access to CVSpils interface            */
#include <nvector/nvector_raja.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <RAJA/RAJA.hpp>


/* Real Constants */

#define ATOL  RCONST(1.0e-5) /* scalar absolute tolerance */
#define T0    RCONST(0.0)    /* initial time              */
#define T1    RCONST(0.1)    /* first output time         */
#define DTOUT RCONST(0.1)    /* output time increment     */
#define NOUT  10             /* number of output times    */

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)
#define FIVE RCONST(5.0)


/* Type : _UserData (contains model and discretization parameters) */
struct _UserData {
  sunindextype MX, MY, NEQ;
  realtype dx, dy, XMAX, YMAX;
  realtype hdcoef, hacoef, vdcoef;
};

typedef _UserData *UserData;

/* Problem setup and initialization functions */
static UserData SetUserData(int argc, char** argv);
static void SetIC(N_Vector u, UserData data);

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);
static int jtv(N_Vector v, N_Vector Jv, realtype t,
               N_Vector u, N_Vector fu,
               void *user_data, N_Vector tmp);

/* Private Helper Functions */
static void PrintHeader(realtype reltol, realtype abstol, realtype umax, UserData data);
static void PrintOutput(realtype t, realtype umax, long int nst);
static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(int argc, char** argv)
{
  realtype reltol, abstol, t, tout, umax;
  N_Vector u;
  UserData data;
  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, flag;
  long int nst;

  u = NULL;
  data = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Set model parameters */
  data = SetUserData(argc, argv);
  if(check_flag((void *)data, "malloc", 2)) return(1);

  /* Set the tolerances */
  reltol = ZERO;
  abstol = ATOL;

  /* Create a RAJA vector with initial values */
  u = N_VNew_Raja(data->NEQ);  /* Allocate u vector */
  if(check_flag((void*)u, "N_VNew_Raja", 0)) return(1);

  SetIC(u, data);  /* Initialize u vector */

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the initial time T0, and
   * the initial dependent variable vector u. */
  flag = CVodeInit(cvode_mem, f, T0, u);
  if(check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerance */
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  /* Set the pointer to user-defined data */
  flag = CVodeSetUserData(cvode_mem, data);
  if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Create SPGMR solver structure without preconditioning
   * and the maximum Krylov dimension maxl */
  LS = SUNSPGMR(u, PREC_NONE, 0);
  if(check_flag(&flag, "SUNSPGMR", 1)) return(1);

  /* Set CVSpils linear solver to LS */
  flag = CVSpilsSetLinearSolver(cvode_mem, LS);
  if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);

  /* Set the Jacobian-times-vector function */
  flag = CVSpilsSetJacTimes(cvode_mem, NULL, jtv);
  if(check_flag(&flag, "CVSpilsSetJacTimesVecFn", 1)) return(1);

  /* In loop over output points: call CVode, print results, test for errors */

  umax = N_VMaxNorm(u);
  PrintHeader(reltol, abstol, umax, data);
  for(iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_flag(&flag, "CVode", 1)) break;
    umax = N_VMaxNorm(u);
    flag = CVodeGetNumSteps(cvode_mem, &nst);
    check_flag(&flag, "CVodeGetNumSteps", 1);
    PrintOutput(t, umax, nst);
  }

  PrintFinalStats(cvode_mem);  /* Print some final statistics   */

  N_VDestroy(u);          /* Free the u vector */
  CVodeFree(&cvode_mem);  /* Free the integrator memory */
  free(data);             /* Free the user data */

  return(0);
}

/*
 *-------------------------------------------
 * Problem setup and initialization functions
 *-------------------------------------------
 */

/* Set model and discretization parameters */

UserData SetUserData(int argc, char *argv[])
{
  const sunindextype MX = 10;
  const sunindextype MY = 5;
  const realtype XMAX = RCONST(2.0);    /* domain boundaries         */
  const realtype YMAX = RCONST(1.0);

  /* Allocate user data structure */
  UserData ud = (UserData) malloc(sizeof *ud);
  if(check_flag((void*) ud, "AllocUserData", 2)) return(NULL);

  ud->MX  = MX;
  ud->MY  = MY;
  ud->NEQ = MX*MY;
  ud->XMAX = XMAX;
  ud->YMAX = YMAX;
  ud->dx = XMAX/(MX+1);  /* Set grid coefficients in data */
  ud->dy = YMAX/(MY+1);
  ud->hdcoef = ONE/(ud->dx*ud->dx);
  ud->hacoef = HALF/(TWO*ud->dx);
  ud->vdcoef = ONE/(ud->dy*ud->dy);

  return ud;
}


/* Set initial conditions in u vector */

static void SetIC(N_Vector u, UserData data)
{
  /* Extract needed constants from data */

  const realtype dx = data->dx;
  const realtype dy = data->dy;
  const realtype xmax = data->XMAX;
  const realtype ymax = data->YMAX;
  const sunindextype MY = data->MY;
  const sunindextype NEQ = data->NEQ;

  /* Extract pointer to solution vector data on the host */
  realtype *udata = N_VGetHostArrayPointer_Raja(u);

  sunindextype i, j, tid;
  realtype x, y;


  /* Load initial profile into u vector */

  for (tid=0; tid < NEQ; tid++) {
    i = tid / MY;
    j = tid % MY;

    x = (i+1)*dx;
    y = (j+1)*dy;

    udata[tid] = x*(xmax - x)*y*(ymax - y)*SUNRexp(FIVE*x*y);
  }
  N_VCopyToDevice_Raja(u);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* f routine. Compute f(t,u). */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  UserData data = (UserData) user_data;

  /* Extract needed constants from data */
  const sunindextype MX  = data->MX;
  const sunindextype MY  = data->MY;
  const sunindextype NEQ = data->NEQ;
  const realtype hordc   = data->hdcoef;
  const realtype horac   = data->hacoef;
  const realtype verdc   = data->vdcoef;

  /* Extract pointers to vector data */
  const realtype *udata = N_VGetDeviceArrayPointer_Raja(u);
  realtype *dudata      = N_VGetDeviceArrayPointer_Raja(udot);

  const sunindextype zero = 0;

  RAJA::forall<RAJA::cuda_exec<256> >(zero, NEQ, [=] __device__(sunindextype index) {
    sunindextype i = index/MY;
    sunindextype j = index%MY;

    realtype uab = udata[index];

    realtype udn = (j == 0)    ? ZERO : udata[index - 1];
    realtype uup = (j == MY-1) ? ZERO : udata[index + 1];
    realtype ult = (i == 0)    ? ZERO : udata[index - MY];
    realtype urt = (i == MX-1) ? ZERO : udata[index + MY];

    /* Set diffusion and advection terms and load into udot */

    realtype hdiff = hordc*(ult -TWO*uab + urt);
    realtype hadv  = horac*(urt - ult);
    realtype vdiff = verdc*(udn -TWO*uab + uup);

    dudata[index] = hdiff + hadv + vdiff;
  });

  return(0);
}


/* Jacobian-times-vector routine. */

static int jtv(N_Vector v, N_Vector Jv, realtype t,
               N_Vector u, N_Vector fu,
               void *user_data, N_Vector tmp)
{
  UserData data = (UserData) user_data;

  /* Extract needed constants from data */
  const sunindextype MX  = data->MX;
  const sunindextype MY  = data->MY;
  const sunindextype NEQ = data->NEQ;
  const realtype hordc   = data->hdcoef;
  const realtype horac   = data->hacoef;
  const realtype verdc   = data->vdcoef;

  /* Extract pointers to vector data */
  const realtype *vdata = N_VGetDeviceArrayPointer_Raja(v);
  realtype *Jvdata      = N_VGetDeviceArrayPointer_Raja(Jv);

  const sunindextype zero = 0;

  N_VConst(ZERO, Jv);

  RAJA::forall<RAJA::cuda_exec<256> >(zero, NEQ, [=] __device__(sunindextype index) {
    sunindextype i = index/MY;
    sunindextype j = index%MY;

    Jvdata[index] = -TWO*(verdc+hordc) * vdata[index];
    if (i !=    0) Jvdata[index] += (hordc - horac) * vdata[index-MY];
    if (i != MX-1) Jvdata[index] += (hordc + horac) * vdata[index+MY];
    if (j !=    0) Jvdata[index] += verdc * vdata[index-1];
    if (j != MY-1) Jvdata[index] += verdc * vdata[index+1];
  });
  
  return(0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

/* Print first lines of output (problem description) */

static void PrintHeader(realtype reltol, realtype abstol, realtype umax, UserData data)
{
  printf("\n2-D Advection-Diffusion Equation\n");
  printf("Mesh dimensions = %d X %d\n", data->MX, data->MY);
  printf("Total system size = %d\n", data->NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters: reltol = %Lg   abstol = %Lg\n\n",
         reltol, abstol);
  printf("At t = %Lg      max.norm(u) =%14.6Le \n", T0, umax);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters: reltol = %g   abstol = %g\n\n",
         reltol, abstol);
  printf("At t = %g      max.norm(u) =%14.6e \n", T0, umax);
#else
  printf("Tolerance parameters: reltol = %g   abstol = %g\n\n", reltol, abstol);
  printf("At t = %g      max.norm(u) =%14.6e \n", T0, umax);
#endif

  return;
}

/* Print current value */

static void PrintOutput(realtype t, realtype umax, long int nst)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %4.2Lf   max.norm(u) =%14.6Le   nst = %4ld\n", t, umax, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %4.2f   max.norm(u) =%14.6e   nst = %4ld\n", t, umax, nst);
#else
  printf("At t = %4.2f   max.norm(u) =%14.6e   nst = %4ld\n", t, umax, nst);
#endif

  return;
}

/* Get and print some final statistics */

static void PrintFinalStats(void *cvode_mem)
{
  long lenrw, leniw ;
  long lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int flag;

  flag = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_flag(&flag, "CVodeGetWorkSpace", 1);
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVSpilsGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
  check_flag(&flag, "CVSpilsGetWorkSpace", 1);
  flag = CVSpilsGetNumLinIters(cvode_mem, &nli);
  check_flag(&flag, "CVSpilsGetNumLinIters", 1);
  flag = CVSpilsGetNumPrecEvals(cvode_mem, &npe);
  check_flag(&flag, "CVSpilsGetNumPrecEvals", 1);
  flag = CVSpilsGetNumPrecSolves(cvode_mem, &nps);
  check_flag(&flag, "CVSpilsGetNumPrecSolves", 1);
  flag = CVSpilsGetNumConvFails(cvode_mem, &ncfl);
  check_flag(&flag, "CVSpilsGetNumConvFails", 1);
  flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVSpilsGetNumRhsEvals", 1);

  printf("\nFinal Statistics.. \n\n");
  printf("lenrw   = %5ld     leniw   = %5ld\n", lenrw, leniw);
  printf("lenrwLS = %5ld     leniwLS = %5ld\n", lenrwLS, leniwLS);
  printf("nst     = %5ld\n"                  , nst);
  printf("nfe     = %5ld     nfeLS   = %5ld\n"  , nfe, nfeLS);
  printf("nni     = %5ld     nli     = %5ld\n"  , nni, nli);
  printf("nsetups = %5ld     netf    = %5ld\n"  , nsetups, netf);
  printf("npe     = %5ld     nps     = %5ld\n"  , npe, nps);
  printf("ncfn    = %5ld     ncfl    = %5ld\n\n", ncfn, ncfl);

  return;
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */

  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */

  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */

  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}
