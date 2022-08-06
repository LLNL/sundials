/*
 * -----------------------------------------------------------------
 * Programmer(s): Jean M. Sexton @ SMU
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Based on work by Scott D. Cohen, Alan C. Hindmarsh, George Byrne,
 *                  and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the program for
 * its solution by CVODE. The problem is the semi-discrete
 * form of the advection-diffusion equation in 1-D:
 *   du/dt = D d^2 u / dx^2 + c du/dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is the following:
 *   u(x,t=0) = x(2-x)exp(2x) .
 * The PDE is discretized on a uniform grid of size MX+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX.
 * This program solves the problem with the option for nonstiff
 * systems: ADAMS method and fixed-point iteration.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .5, 1.0, ..., 5.
 * Run statistics (optional outputs) are printed at the end.
 *
 * This example uses Hypre vector with "IJ" interface and MPI
 * parallelization. User is expected to be familiar with the Hypre
 * library.
 *
 * Execute with Number of Processors = N,  with 1 <= N <= MX.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cvode/cvode.h>                          /* prototypes for CVODE fcts.                   */
#include <sundials/sundials_types.h>              /* definition of realtype                       */
#include <sundials/sundials_math.h>               /* definition of EXP                            */
#include <nvector/nvector_parhyp.h>               /* nvector implementation                       */
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to the fixed point SUNNonlinearSolver */

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>

#include <mpi.h>                      /* MPI constants and types */

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
#include <cuda_runtime.h>
#endif

/* Problem Constants */

#define ZERO  RCONST(0.0)

#define RTOL  RCONST(1.0e-4) /* scalar absolute tolerance */
#define ATOL  RCONST(1.0e-9) /* scalar absolute tolerance */
#define T0    ZERO           /* initial time              */
#define T1    RCONST(0.5)    /* first output time         */
#define DTOUT RCONST(0.5)    /* output time increment     */
#define NOUT  10             /* number of output times    */

/* Type : UserData
   contains grid constants, parhyp machine parameters, work array. */

typedef struct {
  realtype dx, hdcoef, hacoef;
  int npes, my_pe;
  MPI_Comm comm;
  realtype* z;
} *UserData;

/* Private Helper Functions */

static void SetIC(HYPRE_IJVector Uij, realtype XMAX, realtype dx,
                  sunindextype my_length, sunindextype my_base);

static void PrintIntro(int npes, sunindextype MX);

static void PrintData(realtype t, realtype umax, long int nst);

/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

// CUDA helper functions and kernles
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
static size_t CUDAConfigBlockSize()
{
  size_t blocksize = 256;
  return blocksize;
}

static size_t CUDAConfigGridSize(sunindextype N, size_t &blocksize)
{
  return ((N + blocksize - 1) / blocksize);
}

#define GRID_STRIDE_XLOOP(type, iter, max)                      \
  for (type iter = blockDim.x * blockIdx.x + threadIdx.x;       \
       iter < max;                                              \
       iter += blockDim.x * gridDim.x)

__global__ void copyWithinDeviceKernel(realtype *z, realtype *udata, int n)
{
  GRID_STRIDE_XLOOP(int, i, n)
  {
    z[i+1] = udata[i];
  }
}

__global__ void loadingInitialProfileKernel(sunindextype my_base, HYPRE_Int *iglobal, realtype *udata,
                                            sunindextype my_length, realtype XMAX, realtype dx)
{
  realtype x;
  GRID_STRIDE_XLOOP(int, i, my_length)
  {
    iglobal[i] = my_base + i;
    x=(iglobal[i] + 1)*dx;
    udata[i] = x*(XMAX - x)*SUNRexp(RCONST(4.0) * x / XMAX );
  }
}

__global__ void loadDiffAdvTermsKernel(realtype hordc, realtype horac, realtype *udotdata, realtype *z, sunindextype my_length)
{
  realtype hdiff, hadv;

  GRID_STRIDE_XLOOP(int, i, my_length)
  {
    hdiff = hordc*(z[i] - RCONST(2.0)*z[i+1] + z[i+2]);
    hadv = horac*(z[i+2] - z[i]);
    udotdata[i] = hdiff + hadv;
  }
}
#endif

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  realtype dx, reltol, abstol, t, tout, umax, *z;
  N_Vector u;
  UserData data;
  void *cvode_mem;
  int iout, retval, my_pe, npes;
  long int nst;
  HYPRE_Int local_N, nperpe, nrem, my_base;
  HYPRE_ParVector Upar; /* Declare HYPRE parallel vector */
  HYPRE_IJVector  Uij;  /* Declare "IJ" interface to HYPRE vector */
  SUNNonlinearSolver NLS;
  SUNContext sunctx;

  sunindextype MX  = 10;
  realtype XMAX    = 2;
  realtype dcoef   = RCONST(1.0);
  realtype acoef   = RCONST(0.5);
  int print_status = 1;

  MPI_Comm comm;

  u = NULL;
  data = NULL;
  cvode_mem = NULL;

  /* Get processor number, total number of pe's, and my_pe. */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  /* domain boundary */
  /* mesh dimension (number of equations) */
  /* print integration status */
  if (argc > 2) MX = atoi(argv[1]);
  if (argc > 1) XMAX = atof(argv[2]);
  if (argc > 3) dcoef = atof(argv[3]);
  if (argc > 4) acoef = atof(argv[4]);
  if (argc > 5) print_status = atoi(argv[5]);

  /* Create SUNDIALS context */
  retval = SUNContext_Create(&comm, &sunctx);
  if (check_retval(&retval, "SUNContex_Create", 1, my_pe)) MPI_Abort(comm, 1);

  /* Set partitioning. */
  nperpe = MX / npes;
  nrem = MX - npes*nperpe;
  local_N = (my_pe < nrem) ? nperpe+1 : nperpe;
  my_base = (my_pe < nrem) ? my_pe*local_N : my_pe*nperpe + nrem;

  /* Allocate hypre vector */
  HYPRE_IJVectorCreate(comm, my_base, my_base + local_N - 1, &Uij);
  HYPRE_IJVectorSetObjectType(Uij, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(Uij);

  /* Allocate user defined data */
  data = (UserData) malloc(sizeof *data);  /* Allocate data memory */
  if(check_retval((void *)data, "malloc", 2, my_pe)) MPI_Abort(comm, 1);

  data->comm = comm;
  data->npes = npes;
  data->my_pe = my_pe;

#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  z = (realtype *) malloc((local_N + 2) * sizeof(realtype));
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
  cudaMalloc(&z, (local_N + 2) * sizeof(realtype));
#endif

  reltol = RTOL;  /* Set the tolerances */
  abstol = ATOL;

  dx = data->dx = XMAX / ((realtype)(MX+1));  /* Set grid coefficients in data */
  data->hdcoef = dcoef / (dx * dx);
  data->hacoef = acoef / (RCONST(2.0) * dx);
  data->z = z;

  /* Initialize solution vector. */
  SetIC(Uij, XMAX, dx, local_N, my_base);
  HYPRE_IJVectorAssemble(Uij);
  HYPRE_IJVectorGetObject(Uij, (void**) &Upar);

  u = N_VMake_ParHyp(Upar, sunctx);  /* Create wrapper u around hypre vector */
  if(check_retval((void *)u, "N_VNew", 0, my_pe)) MPI_Abort(comm, 1);

  /* Call CVodeCreate to create the solver memory and specify the
   * Adams-Moulton LMM */
  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0, my_pe)) MPI_Abort(comm, 1);

  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1, my_pe)) MPI_Abort(comm, 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, f, T0, u);
  CVodeSetMaxNumSteps(cvode_mem, 1000000); 
  if(check_retval(&retval, "CVodeInit", 1, my_pe)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1, my_pe)) return(1);

  retval = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1, my_pe)) return(1);


  /* create fixed point nonlinear solver object */
  NLS = SUNNonlinSol_FixedPoint(u, 0, sunctx);
  if(check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0, my_pe)) return(1);

  /* attach nonlinear solver object to CVode */
  retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if(check_retval(&retval, "CVodeSetNonlinearSolver", 1, my_pe)) return(1);

  if (my_pe == 0) PrintIntro(npes, MX);

  if (print_status)
  {
    umax = N_VMaxNorm(u);
    if (my_pe == 0)
    {
      t = T0;
      PrintData(t, umax, 0);
    }
  }

  /* In loop over output points, call CVode, print results, test for error */
  double tic, toc;
  double totaltime = 0.0;

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT)
  {
    tic = MPI_Wtime(); /* start timer */
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    toc = MPI_Wtime(); /* stop timer */
    totaltime += (toc - tic); /* update total time */

    if(check_retval(&retval, "CVode", 1, my_pe)) break;

    if (print_status)
    {
      umax = N_VMaxNorm(u);

      retval = CVodeGetNumSteps(cvode_mem, &nst);
      check_retval(&retval, "CVodeGetNumSteps", 1, my_pe);

      if (my_pe == 0) PrintData(t, umax, nst);
    }
  }

  /* Get the max time across processes */
  double maxtime = 0.0;
  MPI_Reduce(&(totaltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

  if (my_pe == 0)
  {
    /* Print some final statistics */
    CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    printf("total time = %g\n", totaltime);
  }

  N_VDestroy(u);              /* Free hypre vector wrapper */
  HYPRE_IJVectorDestroy(Uij); /* Free the underlying hypre vector */
  CVodeFree(&cvode_mem);      /* Free the integrator memory */
  SUNNonlinSolFree(NLS);      /* Free the nonlinear solver */
  free(data);                 /* Free user data */
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
  cudaFree(z);                /* Free buffer data */
#endif
  SUNContext_Free(&sunctx);   /* Free context */

  MPI_Finalize();

  return(0);
}

/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */

static void SetIC(HYPRE_IJVector Uij, realtype XMAX, realtype dx,
                  sunindextype my_length, sunindextype my_base)
{
  int i;
  HYPRE_Int *iglobal;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u. */
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  udata   = (realtype*) malloc(my_length * sizeof(realtype));
  iglobal = (HYPRE_Int*) malloc(my_length * sizeof(HYPRE_Int));
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
  cudaMalloc(&udata, my_length * sizeof(realtype));
  cudaMalloc(&iglobal, my_length * sizeof(HYPRE_Int));
#endif

  /* Load initial profile into u vector */

#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  for (i = 0; i < my_length; i++) {
    iglobal[i] = my_base + i;
    x = (iglobal[i] + 1)*dx;
    udata[i] = x*(XMAX - x)*SUNRexp(RCONST(2.0)*x);
  }
  HYPRE_IJVectorSetValues(Uij, my_length, iglobal, udata);
  free(iglobal);
  free(udata);

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
  size_t blocksize =  CUDAConfigBlockSize();
  size_t gridsize = CUDAConfigGridSize(my_length, blocksize);
  loadingInitialProfileKernel<<<gridsize, blocksize, 0, 0>>>(my_base, iglobal,
                                                             udata, my_length,
                                                             XMAX, dx);
  HYPRE_IJVectorSetValues(Uij, my_length, iglobal, udata);
  cudaFree(iglobal);
  cudaFree(udata);
#endif
}

/* Print problem introduction */

static void PrintIntro(int npes, sunindextype MX)
{
  printf("\n 1-D advection-diffusion equation, mesh size =%3d \n", MX);
  printf("\n Number of PEs = %3d \n\n", npes);

  return;
}

/* Print data */

static void PrintData(realtype t, realtype umax, long int nst)
{

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %4.2Lf  max.norm(u) =%14.6Le  nst =%4ld \n", t, umax, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %4.2f  max.norm(u) =%14.6e  nst =%4ld \n", t, umax, nst);
#else
  printf("At t = %4.2f  max.norm(u) =%14.6e  nst =%4ld \n", t, umax, nst);
#endif

  return;
}

/***************** Function Called by the Solver ***********************/

/* f routine. Compute f(t,u). */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype ui, ult, urt, hordc, horac, hdiff, hadv;
  realtype *udata, *udotdata, *z;
  int i;
  int npes, my_pe, my_length, my_pe_m1, my_pe_p1, last_pe;
  UserData data;
  MPI_Status status;
  MPI_Comm comm;
  HYPRE_ParVector uhyp;
  HYPRE_ParVector udothyp;

  /* Extract hypre vectors */
  uhyp  = N_VGetVector_ParHyp(u);
  udothyp  = N_VGetVector_ParHyp(udot);

  /* Access hypre vectors local data */
  udata = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));
  udotdata = hypre_VectorData(hypre_ParVectorLocalVector(udothyp));

  /* Extract needed problem constants from data */
  data = (UserData) user_data;
  hordc = data->hdcoef;
  horac = data->hacoef;

  /* Extract parameters for parhyp computation. */
  comm = data->comm;
  npes = data->npes;                           /* Number of processes    */
  my_pe = data->my_pe;                         /* Current process number */
  my_length =  hypre_ParVectorLastIndex(uhyp)  /* Local length of uhyp   */
    - hypre_ParVectorFirstIndex(uhyp) + 1;
  z = data->z;

  /* Compute related parameters. */
  my_pe_m1 = my_pe - 1;
  my_pe_p1 = my_pe + 1;
  last_pe = npes - 1;

  /* Store local segment of u in the working array z. */
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)

  z[0] = ZERO;
  for (i = 1; i <= my_length; i++)
    z[i] = udata[i - 1];
  z[my_length + 1] = ZERO;

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)

  cudaMemset(z, ZERO, (my_length + 2) * sizeof(realtype));
  size_t blocksize =  CUDAConfigBlockSize();
  size_t gridsize = CUDAConfigGridSize(my_length, blocksize);
  copyWithinDeviceKernel<<<gridsize, blocksize, 0, 0>>>(z, udata, my_length);

#endif

  /* Pass needed data to processes before and after current process. */
  if (my_pe != 0)
    MPI_Send(z+1, 1, MPI_SUNREALTYPE, my_pe_m1, 0, comm);
  if (my_pe != last_pe)
    MPI_Send(z+my_length, 1, MPI_SUNREALTYPE, my_pe_p1, 0, comm);

  /* Receive needed data from processes before and after current process. */
  if (my_pe != 0)
    MPI_Recv(z, 1, MPI_SUNREALTYPE, my_pe_m1, 0, comm, &status);
  if (my_pe != last_pe)
    MPI_Recv(z+(my_length+1), 1, MPI_SUNREALTYPE, my_pe_p1, 0, comm, &status);

#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)

  /* Loop over all grid points in current process. */
  for (i=1; i<=my_length; i++) {

    /* Extract u at x_i and two neighboring points */
    ui = z[i];
    ult = z[i-1];
    urt = z[i+1];

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - RCONST(2.0)*ui + urt);
    hadv = horac*(urt - ult);
    udotdata[i-1] = hdiff + hadv;
  }

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA)

  blocksize = CUDAConfigBlockSize();
  gridsize = CUDAConfigGridSize(my_length, blocksize);
  loadDiffAdvTermsKernel<<<gridsize, blocksize, 0, 0>>>(hordc, horac, udotdata,
                                                        z, my_length);

#endif

  return(0);
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval < 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed with retval = %d\n\n",
              id, funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1); }

  return(0);
}
