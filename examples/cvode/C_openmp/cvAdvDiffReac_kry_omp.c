/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 *
 * -----------------------------------------------------------------
 */

/* Minimum POSIX version needed for struct timespec and clock_monotonic */
#if !defined(_POSIX_C_SOURCE) || (_POSIX_C_SOURCE < 199309L)
#define _POSIX_C_SOURCE 199309L
#endif

/* POSIX timers */
#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
#include <time.h>
#include <stddef.h>
#include <unistd.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#include <cvode/cvode.h>
#include <cvode/cvode_spils.h>         /* access to CVSpils interface            */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver        */
#include <nvector/nvector_openmp.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

typedef struct
{
  sunindextype Nx;
  sunindextype Ny;
  sunindextype NEQ;

  int num_threads;

  realtype hx;
  realtype hy;

  realtype hordc;
  realtype verdc;
  realtype horac;
  realtype verac;
  realtype reacc;

} *UserData;

/* User defined functions */

static N_Vector SetIC(UserData data, SUNContext sunctx);
static UserData SetUserData(int argc, char *argv[]);
static void Phiu(N_Vector u, N_Vector result, sunindextype NEQ, sunindextype Nx, sunindextype Ny,
                 realtype hordc, realtype verdc, realtype horac, realtype verac);
static int RHS(realtype t, N_Vector u, N_Vector udot, void *userData);
static int Jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *userData, N_Vector tmp);


/* Private Helper Functions */

static void PrintOutput(void *cvode_mem, N_Vector u, realtype t);
static void PrintFinalStats(void *cvode_mem);
static int check_retval(void *retvalvalue, const char *funcname, int opt);


/* private functions */
static double get_time();

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(int argc, char *argv[])
{
  SUNContext sunctx;
  realtype abstol, reltol, t, tout;
  const realtype t_in = 0.0;
  const realtype t_fi = 0.1;
  N_Vector u;
  UserData data;
  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, retval;

  u = NULL;
  data = NULL;
  cvode_mem = NULL;

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if(check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* Allocate memory, set problem data and initial values */
  data = SetUserData(argc, argv);
  u = SetIC(data, sunctx);

  reltol = RCONST(1.0e-5);         /* scalar relative tolerance */
  abstol = reltol * RCONST(100.0); /* scalar absolute tolerance */

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Set the pointer to user-defined data */
  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, RHS, t_in, u);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* Create SPGMR solver structure without preconditioning
   * and the maximum Krylov dimension maxl */
  LS = SUNSPGMR(u, SUN_PREC_NONE, 0);
  if(check_retval(&retval, "SUNSPGMR", 1)) return(1);

  /* Set CVSpils linear solver to LS */
  retval = CVodeSetLinearSolver(cvode_mem, LS);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* set the JAcobian-times-vector function */
  retval = CVodeSetJacTimes(cvode_mem, NULL, Jtv);
  if(check_retval(&retval, "CVodeSetJacTimes", 1)) return(1);


  printf("Solving diffusion-advection-reaction problem with %d unknowns...\n", data->NEQ);

  double start_time, stop_time;
  start_time = get_time();
  retval = CVode(cvode_mem, t_fi, u, &t, CV_NORMAL);
  stop_time = get_time();
  PrintOutput(cvode_mem, u, t);
  if(check_retval(&retval, "CVode", 1))
    return (-1);

  printf("Computation successful!\n");
  /* printf("Execution time on %d threads = %g\n", data->num_threads, stop_time - start_time); */
  printf("L2 norm = %14.6e\n", SUNRsqrt(N_VDotProd(u,u)));

  PrintFinalStats(cvode_mem);

  /* Free memory */
  N_VDestroy(u);
  free(data);
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);

  return(0);
}


/*
 *-------------------------------
 * User defined functions
 *-------------------------------
 */

N_Vector SetIC(UserData data, SUNContext sunctx)
{
  const sunindextype Nx = data->Nx;
  const realtype hx = data->hx;
  const realtype hy = data->hy;

  N_Vector y = N_VNew_OpenMP(data->NEQ, 4, sunctx);
  realtype *ydat = N_VGetArrayPointer(y);
  sunindextype i, j, index;

  for (index = 0; index < data->NEQ; ++index)
  {
    j = index/Nx;
    i = index%Nx;

    realtype y = j * hy;
    realtype x = i * hx;
    realtype tmp = (1 - x) * x * (1 - y) * y;
    ydat[index] = (256.0 * tmp * tmp) + 0.3;
  }
  return y;
}

UserData SetUserData(int argc, char *argv[])
{
  sunindextype dimX = 70; /* Default grid size */
  sunindextype dimY = 80;
  const realtype diffusionConst =  0.01;
  const realtype advectionConst = -10.0;
  const realtype reactionConst  = 100.0;

  /* Allocate user data structure */
  UserData ud = (UserData) malloc(sizeof *ud);
  if(check_retval((void*) ud, "AllocUserData", 2)) return(NULL);

  /* Set the number of threads to use */
  ud->num_threads = 1;     /* default value */
#ifdef _OPENMP
  ud->num_threads = omp_get_max_threads();  /* Overwrite with OMP_NUM_THREADS environment variable */
#endif
  if (argc == 2) {  /* Overwrite with command line values, if supplied */
    ud->num_threads = strtol(argv[1], (char**) NULL, 10);
  }

  /* Set grid size */
  if (argc == 3) {   /* Overwrite grid size with command line values, if supplied */
    dimX = strtol(argv[1], (char**) NULL, 10);
    dimY = strtol(argv[2], (char**) NULL, 10);
  }
  ud->Nx = dimX + 1;
  ud->Ny = dimY + 1;
  ud->NEQ = ud->Nx * ud->Ny;

  /* Compute cell sizes */
  ud->hx = 1.0/((realtype) dimX);
  ud->hy = 1.0/((realtype) dimY);

  /* Compute diffusion coefficients */
  ud->hordc = diffusionConst/(ud->hx * ud->hx);
  ud->verdc = diffusionConst/(ud->hy * ud->hy);

  /* Compute advection coefficient */
  ud->horac = advectionConst/(2.0 * ud->hx);
  ud->verac = advectionConst/(2.0 * ud->hy);

  /* Set reaction coefficient */
  ud->reacc = reactionConst;

  return ud;
}


void Phiu(N_Vector u, N_Vector result, sunindextype NEQ, sunindextype Nx, sunindextype Ny,
          realtype hordc, realtype verdc, realtype horac, realtype verac)
{
  const realtype *uData = N_VGetArrayPointer(u);
  realtype *resultData = N_VGetArrayPointer(result);

  sunindextype i, j, index;

  realtype uij;
  realtype ult;
  realtype urt;
  realtype uup;
  realtype udn;

  realtype hdiff;
  realtype vdiff;
  realtype hadv;
  realtype vadv;

#pragma omp parallel for default(shared) private(index, i, j, uij, ult, urt, uup, udn, hdiff, vdiff, hadv, vadv)
  for (index = 0; index < NEQ; ++index)
  {
    i = index%Nx;
    j = index/Nx;

    uij = uData[index];

    ult = (i == 0)    ? uData[index + 1]  : uData[index - 1];
    urt = (i == Nx-1) ? uData[index - 1]  : uData[index + 1];
    udn = (j == 0)    ? uData[index + Nx] : uData[index - Nx];
    uup = (j == Ny-1) ? uData[index - Nx] : uData[index + Nx];

    hdiff =  hordc*(ult -2.0*uij + urt);
    vdiff =  verdc*(udn -2.0*uij + uup);
    hadv  = -horac*(urt - ult);
    vadv  = -verac*(uup - udn);

    resultData[index] = hdiff + vdiff + hadv + vadv;
  }

}

int RHS(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  UserData data = (UserData) user_data;
  const realtype *udata = N_VGetArrayPointer(u);
  realtype *udotdata = N_VGetArrayPointer(udot);

  const realtype a = -1.0 / 2.0;
  sunindextype i;

  Phiu(u, udot, data->NEQ, data->Nx, data->Ny, data->hordc, data->verdc, data->horac, data->verac);

#pragma omp parallel for default(shared) private(i)
  for(i=0; i<data->NEQ; ++i)
  {
    udotdata[i] += (data->reacc*(udata[i] + a)*(1.0 - udata[i])*udata[i]);
  }

  return 0;
}

int Jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *user_data, N_Vector tmp)
{
  UserData data = (UserData) user_data;
  realtype *udata = N_VGetArrayPointer(u);
  realtype *vdata = N_VGetArrayPointer(v);
  realtype *Jvdata = N_VGetArrayPointer(Jv);

  const realtype a = -1.0 / 2.0;
  sunindextype i;

  Phiu(v, Jv, data->NEQ, data->Nx, data->Ny, data->hordc, data->verdc, data->horac, data->verac);

#pragma omp parallel for default(shared) private(i)
  for(i=0; i<data->NEQ; ++i)
  {
    Jvdata[i] += data->reacc*(3.0*udata[i] + a - 3.0*udata[i]*udata[i])*vdata[i];
  }

  return 0;
}



/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */


/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(void *cvode_mem, N_Vector u, realtype t)
{
  long int nst;
  int qu, retval;
  realtype hu, *udata;

  udata = N_VGetArrayPointer(u);

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetLastOrder(cvode_mem, &qu);
  check_retval(&retval, "CVodeGetLastOrder", 1);
  retval = CVodeGetLastStep(cvode_mem, &hu);
  check_retval(&retval, "CVodeGetLastStep", 1);

}

/* Get and print final statistics */

static void PrintFinalStats(void *cvode_mem)
{
  sunindextype lenrw, leniw ;
  sunindextype lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int retval;

  retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_retval(&retval, "CVodeGetWorkSpace", 1);
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  retval = CVodeGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
  check_retval(&retval, "CVodeGetWorkSpace", 1);
  retval = CVodeGetNumLinIters(cvode_mem, &nli);
  check_retval(&retval, "CVodeGetNumLinIters", 1);
  retval = CVodeGetNumPrecEvals(cvode_mem, &npe);
  check_retval(&retval, "CVodeGetNumPrecEvals", 1);
  retval = CVodeGetNumPrecSolves(cvode_mem, &nps);
  check_retval(&retval, "CVodeGetNumPrecSolves", 1);
  retval = CVodeGetNumConvFails(cvode_mem, &ncfl);
  check_retval(&retval, "CVodeGetNumConvFails", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);

  printf("\nFinal Statistics.. \n\n");
  printf("lenrw   = %5ld     leniw   = %5ld\n", lenrw, leniw);
  printf("lenrwLS = %5ld     leniwLS = %5ld\n", lenrwLS, leniwLS);
  printf("nst     = %5ld\n"                  , nst);
  printf("nfe     = %5ld     nfeLS   = %5ld\n"  , nfe, nfeLS);
  printf("nni     = %5ld     nli     = %5ld\n"  , nni, nli);
  printf("nsetups = %5ld     netf    = %5ld\n"  , nsetups, netf);
  printf("npe     = %5ld     nps     = %5ld\n"  , npe, nps);
  printf("ncfn    = %5ld     ncfl    = %5ld\n\n", ncfn, ncfl);
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a retval so check if
              retval >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_retval(void *retvalvalue, const char *funcname, int opt)
{
  int *errretval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && retvalvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    errretval = (int *) retvalvalue;
    if (*errretval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *errretval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && retvalvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}

/* ======================================================================
 * Timing functions
 * ====================================================================*/

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
time_t base_time_tv_sec = 0; /* Base time; makes time values returned
                                by get_time easier to read when
                                printed since they will be zero
                                based.
                              */
#else
#warning "No posix timers!\n"
#endif

void SetTiming(int onoff)
{
   /* print_time = onoff; */

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime( CLOCK_MONOTONIC, &spec );
  base_time_tv_sec = spec.tv_sec;
#endif
}

/* ----------------------------------------------------------------------
 * Timer
 * --------------------------------------------------------------------*/
static double get_time()
{
#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime( CLOCK_MONOTONIC, &spec );
  double time = (double)(spec.tv_sec - base_time_tv_sec) + ((double)(spec.tv_nsec) / 1E9);
#else
  double time = 0;
#endif
  return time;
}
