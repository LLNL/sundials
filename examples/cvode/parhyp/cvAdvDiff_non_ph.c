/*
 * -----------------------------------------------------------------
 * Programmer(s): Jean M. Sexton @ SMU
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Based on work by Scott D. Cohen, Alan C. Hindmarsh, George Byrne,
 *                  and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
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
 *   du/dt = d^2 u / dx^2 + .5 du/dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is the following:
 *   u(x,t=0) = x(2-x)exp(2x) .
 * The PDE is discretized on a uniform grid of size M+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = M.
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
 * Execute with Number of Processors = N,  with 1 <= N <= M.
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

/* Problem Constants */

#define ZERO  RCONST(0.0)
#define HALF  RCONST(0.5)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define XMAX  RCONST(2.0)    /* domain boundary           */
// #define global_M    10             /* mesh dimension            */
// #define NEQ   global_M             /* number of equations       */
#define ATOL  RCONST(1.0e-5) /* scalar absolute tolerance */
#define T0    ZERO           /* initial time              */
#define T1    RCONST(0.5)    /* first output time         */
#define DTOUT RCONST(0.5)    /* output time increment     */
#define NOUT  10             /* number of output times    */

#define MPI_TEST(cond,d,msg,code)  \
  if(cond) {                       \
    if (d->myproc==0) printf(msg); \
    MPI_Abort(d->comm,code);      \
  }

/* Type : UserData
   contains grid constants, parhyp machine parameters, work array. */

typedef struct {
  MPI_Comm  comm;               // MPI communicator
  int       nprocs, myproc;     // # of processes and my process id
  HYPRE_Int global_M, local_M;  // global/local problem size
  realtype  dx, hdcoef, hacoef; // gridpoint spacing and problem coefficients
  realtype  bufs[4];            // Send/Recv buffers: [Lrecvbuf,Lsendbuf,Rsendbuf,Rrecvbuf]
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  realtype* bufs_dev;           // Buffer space (4 realtypes) to be allocated on device
#endif
} *UserData;

/* Private Helper Functions */

static void SetIC(HYPRE_IJVector Uij, realtype dx, sunindextype local_M, sunindextype mybase);

static void PrintIntro(int M, int nprocs);

static void PrintData(realtype t, realtype umax, long int nst);

static void PrintFinalStats(void *cvode_mem);

/* GPU Kernels */

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
__global__ void fKernel(realtype *u, realtype *udot,
                        sunindextype local_M, realtype hdcoef, realtype hacoef);
#endif

/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  realtype  reltol, abstol, t, tout, umax;
  void*     cvode_mem;
  int       iout, retval;
  long int  nst;
  
  UserData           d;
  N_Vector           u;
  HYPRE_Int          M, nperproc, nrem, mybase;
  HYPRE_ParVector    Upar; /* Declare HYPRE parallel vector */
  HYPRE_IJVector     Uij;  /* Declare "IJ" interface to HYPRE vector */
  SUNNonlinearSolver NLS;
  SUNContext         sunctx;

/* Optional: create stream and choose policies to use instead of default */
// #if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
//   NV_ADD_LANG_PREFIX_PH(Stream_t) stream;
//   NV_ADD_LANG_PREFIX_PH(Error_t)  gpuerr;
// #endif

  d         = NULL;
  u         = NULL;
  cvode_mem = NULL;

  /* MPI Init */
  MPI_Init(&argc, &argv);

  /* Allocate UserData */
  d = (UserData) malloc(sizeof *d);
  if(check_retval((void *)d, "malloc", 2, d->myproc)) MPI_Abort(d->comm, 1);

  /* Get MPI info */
  d->comm = MPI_COMM_WORLD;
  MPI_Comm_size(d->comm, &d->nprocs);
  MPI_Comm_rank(d->comm, &d->myproc);

  /* Parse inputs (required: M) */
  MPI_TEST(argc<2, d, "ERROR: ONE (1) input required: M (# of interior gridpoints)\n", 1)
  M = (HYPRE_Int) atoi(argv[1]);
  MPI_TEST(d->nprocs>M, d, "ERROR: There are more MPI processes than the number of interior gridpoints (M)\n", 1)

  /* Set partitioning */
  d->global_M = M;
  nperproc    = M / d->nprocs;
  nrem        = M - d->nprocs*nperproc;
  d->local_M  = (d->myproc < nrem) ? nperproc+1 : nperproc;
  mybase      = (d->myproc < nrem) ? d->myproc*d->local_M : d->myproc*nperproc + nrem;

  /* Create SUNDIALS context */
  retval = SUNContext_Create(&d->comm, &sunctx);
  if (check_retval(&retval, "SUNContex_Create", 1, d->myproc)) MPI_Abort(d->comm, 1);

  /* Populate remaining constants */
  d->dx     = XMAX/((realtype)(M+1));
  d->hdcoef = ONE/(d->dx*d->dx);
  d->hacoef = HALF/(TWO*d->dx);

  /* Initialize Send/Recv buffers */
  // Note: We have zero Dirichlet b.c.; Bdry procs will always have zero in their outer Recv buffer
  memset(d->bufs, 0, 4*sizeof(realtype)); // 
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Malloc)(&d->bufs_dev,    4*sizeof(realtype)); // e.g. cudaMalloc
  NV_ADD_LANG_PREFIX_PH(Memset)( d->bufs_dev, 0, 4*sizeof(realtype)); // e.g. cudaMemset
#endif

  /* Allocate hypre vector */
  HYPRE_IJVectorCreate(d->comm, mybase, mybase + d->local_M - 1, &Uij);
  HYPRE_IJVectorSetObjectType(Uij, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(Uij);

  /* Set tolerances */
  reltol = ZERO;
  abstol = ATOL;

// #if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
//   data_host = d;
//   d = NULL;
//   /* Allocate device UserData */
//   NV_ADD_LANG_PREFIX_PH(Malloc)(&d, sizeof(UserData));
  
//   /* Copy host data to device */
//   NV_ADD_LANG_PREFIX_PH(Memcpy)(d, data_host, sizeof(UserData), NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
// #endif

  /* Initialize solution vector. */
  SetIC(Uij, d->dx, d->local_M, mybase);
  HYPRE_IJVectorAssemble(Uij);
  HYPRE_IJVectorGetObject(Uij, (void**) &Upar);

  u = N_VMake_ParHyp(Upar, sunctx);  /* Create wrapper u around hypre vector */
  if(check_retval((void *)u, "N_VNew", 0, d->myproc)) MPI_Abort(d->comm, 1);

  /* Call CVodeCreate to create the solver memory and specify the
   * Adams-Moulton LMM */
  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0, d->myproc)) MPI_Abort(d->comm, 1);

  retval = CVodeSetUserData(cvode_mem, d);
  if(check_retval(&retval, "CVodeSetUserData", 1, d->myproc)) MPI_Abort(d->comm, 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1, d->myproc)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1, d->myproc)) return(1);

  /* Create fixed point nonlinear solver object */
  NLS = SUNNonlinSol_FixedPoint(u, 0, sunctx);
  if(check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0, d->myproc)) return(1);

  /* Attach nonlinear solver object to CVode */
  retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if(check_retval(&retval, "CVodeSetNonlinearSolver", 1, d->myproc)) return(1);

  if (d->myproc == 0) PrintIntro(M, d->nprocs);

  /* Print initial statistics */
  umax = N_VMaxNorm(u);
  if (d->myproc == 0) {
    t = T0;
    PrintData(t, umax, 0);
  }

  /* In loop over output points, call CVode, print results, test for error */
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    printf("checkpoint 1\n");
    /* Advance to next output time */
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_retval(&retval, "CVode", 1, d->myproc)) break;

    printf("checkpoint 2\n");
    /* Get number of steps taken */
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    if(check_retval(&retval, "CVodeGetNumSteps", 1, d->myproc)) break;

    printf("checkpoint 3\n");
    /* Print time, current max norm, and # steps taken since last output */
    umax = N_VMaxNorm(u);
    if (d->myproc == 0) PrintData(t, umax, nst);
  }

  /* Print some final statistics */
  if (d->myproc == 0) PrintFinalStats(cvode_mem);  

  /* Free memory */
  N_VDestroy(u);              /* Free hypre vector wrapper */
  HYPRE_IJVectorDestroy(Uij); /* Free the underlying hypre vector */
  CVodeFree(&cvode_mem);      /* Free the integrator memory */
  SUNNonlinSolFree(NLS);      /* Free the nonlinear solver */
  free(d->bufs);           /* Free Send/Recv buffers */
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Free)(d->bufs_dev);
#endif
  free(d);                 /* Free user data */
  SUNContext_Free(&sunctx);   /* Free context */

  /* MPI Finalize */
  MPI_Finalize();

  return(0);
}

/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */

static void SetIC(HYPRE_IJVector Uij, realtype dx, sunindextype local_M,
                  sunindextype mybase)
{
  int i;
  HYPRE_Int *iglobal;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u. */
  iglobal = (HYPRE_Int*) malloc(local_M*sizeof(HYPRE_Int));
  udata   = (realtype*)  malloc(local_M*sizeof(realtype));

  /* Load initial profile into u vector */
  for (i = 0; i < local_M; i++) {
    iglobal[i] = mybase + i;
    x = (iglobal[i] + 1)*dx;
    udata[i] = x*(XMAX - x)*SUNRexp(TWO*x);
  }

#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  /* For serial backend, HYPRE_IJVectorSetValues expects host pointers */
  HYPRE_IJVectorSetValues(Uij, local_M, iglobal, udata);

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  /* With GPU backend, HYPRE_IJVectorSetValues expects device pointers */
  HYPRE_Int* iglobal_dev;
  realtype*  udata_dev;

  /* NOTE: The macro NV_ADD_LANG_PREFIX_PH(...) adds the prefix "cuda" or "hip" depending on the backend.
   *       The user may use this macro for portability, or explicitly write cuda/hip calls.
   *       Example: "NV_ADD_LANG_PREFIX_PH(Free)" expands to "cudaFree" when using the CUDA hypre backend.
   */

  /* Allocate device arrays */
  NV_ADD_LANG_PREFIX_PH(Malloc)(&iglobal_dev, local_M*sizeof(HYPRE_Int));
  NV_ADD_LANG_PREFIX_PH(Malloc)(&udata_dev, local_M*sizeof(realtype));
  
  /* Copy host data to device */
  NV_ADD_LANG_PREFIX_PH(Memcpy)(iglobal_dev, iglobal, local_M*sizeof(HYPRE_Int), NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
  NV_ADD_LANG_PREFIX_PH(Memcpy)(udata_dev, udata, local_M*sizeof(realtype), NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));

  /* Set hypre vector values from device arrays, then free the device arrays */
  HYPRE_IJVectorSetValues(Uij, local_M, iglobal_dev, udata_dev);
  NV_ADD_LANG_PREFIX_PH(Free)(iglobal_dev);
  NV_ADD_LANG_PREFIX_PH(Free)(udata_dev);
#endif

  free(iglobal);
  free(udata);
}

/* Print problem introduction */

static void PrintIntro(int M, int nprocs)
{
  printf("\n 1-D advection-diffusion equation, mesh size =%3d \n", M);
  printf("\n Number of PEs = %3d \n\n", nprocs);

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

/* Print some final statistics located in the iopt array */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nni, ncfn, netf;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1, 0);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1, 0);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1, 0);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1, 0);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1, 0);

  printf("\nFinal Statistics: \n\n");
  printf("nst = %-6ld  nfe  = %-6ld  ", nst, nfe);
  printf("nni = %-6ld  ncfn = %-6ld  netf = %ld\n \n", nni, ncfn, netf);
}

/****************************** GPU Kernels ******************************/

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
__global__ void fKernel(realtype* udata, realtype* udotdata, realtype* ubufs,
                        sunindextype local_M, realtype hdcoef, realtype hacoef)
{
  realtype uleft, uright, ucenter, hdiff, hadv;
  sunindextype tid;

  /* Get tid */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  /* Get u values (take from Recv buffers as needed) */
  if (tid < local_M) {
    uleft   = (tid==0        ) ? ubufs[0] : udata[tid-1];
    uright  = (tid==local_M-1) ? ubufs[3] : udata[tid+1];
    ucenter = udata[tid];

    /* Set diffusion and advection terms and load into udot */
    hdiff         = hdcoef*(uleft  - TWO*ucenter + uright);
    hadv          = hacoef*(uright - uleft);
    udotdata[tid] = hdiff + hadv;
  }
}
#endif

/********************* Function Called by the Solver *********************/

/* f routine. Compute f(t,u). */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  printf("checkpoint 1a\n");
  // int      nprocs, myproc, local_M, global_M;
  int      procfirst, proclast, procleft, procright;
  realtype *udata, *udotdata; // Both are length global_M (# of gridpoints excluding bdry)

  UserData        d;
  MPI_Status      status;
  MPI_Request     request;
  // MPI_Status      status_leftward_flow, status_rightward_flow;
  // MPI_Comm        comm;
  HYPRE_ParVector uhyp;
  HYPRE_ParVector udothyp;

  printf("checkpoint 1b\n");
  /* Get reference to UserData */
  d        = (UserData) user_data;

  /* Extract hypre vectors */
  uhyp     = N_VGetVector_ParHyp(u);
  udothyp  = N_VGetVector_ParHyp(udot);

  /* Access hypre vectors local data */
  udata    = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));
  udotdata = hypre_VectorData(hypre_ParVectorLocalVector(udothyp));

//   comm     = d->comm;
//   nprocs   = d->nprocs;
//   myproc   = d->myproc;
//   global_M = d->global_M;
//   local_M  = d->local_M;
//   hdcoef   = d->hdcoef;
//   hacoef   = d->hacoef;
//   bufs     = d->bufs
// #if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
//   bufs_dev = d->bufs_dev;
// #endif
//     MPI_Comm  comm;               // MPI communicator
//   int       nprocs, myproc;     // # of processes and my process id
//   HYPRE_Int global_M, local_M;  // global/local problem size
//   realtype  dx, hdcoef, hacoef; // gridpoint spacing and problem coefficients
//   realtype  bufs[4];            // Send/Recv buffers: [Lrecvbuf,Lsendbuf,Rsendbuf,Rrecvbuf]
// #if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
//   realtype* bufs_dev;           // Buffer space (4 realtypes) to be allocated on device
// #endif

  // local_M =  hypre_ParVectorLastIndex(uhyp)  /* Local length of uhyp   */
  //            - hypre_ParVectorFirstIndex(uhyp) + 1;

  printf("checkpoint 1c\n");
  /* Relevant process ids (utilize behavior of MPI_PROC_NULL in Send/Recv's) */
  // Note: Send/Recv's with MPI_PROC_NULL immediately return successful without modifying buffers.
  procfirst = 0;
  proclast  = d->nprocs-1;
  procleft  = (d->myproc==procfirst) ? MPI_PROC_NULL : d->myproc-1;
  procright = (d->myproc==proclast ) ? MPI_PROC_NULL : d->myproc+1;

  // /* Store local segment of u in the working array bufs. */
  // for (i = 1; i <= d->local_M; i++)
  //   d->bufs[i] = udata[i - 1];

  /* Populate Send buffers */
  // Note: bufs order is [Lrecvbuf,Lsendbuf,Rsendbuf,Rrecvbuf]
  // Note: MPI transfers are host-to-host (if MPI is CUDA/HIP-aware, could do device-to-device)
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  d->bufs[1] = udata[0];          //  leftmost local gridpoint is sent left
  d->bufs[2] = udata[d->local_M]; // rightmost local gridpoint is sent right
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Memcpy)(&d->bufs[1],&udata[0],         sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyDeviceToHost));
  NV_ADD_LANG_PREFIX_PH(Memcpy)(&d->bufs[2],&udata[d->local_M],sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyDeviceToHost));
#endif

  printf("checkpoint 1d\n");
  /* Nonblocking leftward flow of data */
  MPI_Irecv(&d->bufs[3],1,MPI_SUNREALTYPE,procright,0,d->comm,&request); // Receive from procright
  MPI_Send( &d->bufs[1],1,MPI_SUNREALTYPE,procleft ,0,d->comm);          // Send to procleft
  MPI_Wait( &request,&status);

  /* Nonblocking rightward flow of data */
  MPI_Irecv(&d->bufs[0],1,MPI_SUNREALTYPE,procleft ,0,d->comm,&request); // Receive from procleft
  MPI_Send( &d->bufs[2],1,MPI_SUNREALTYPE,procright,0,d->comm);          // Send to procright
  MPI_Wait( &request,&status);
  // MPI_Isendrecv(&d->bufs[2],1,MPI_SUNREALTYPE,procright,0, // Send to procright
  //               &d->bufs[0],1,MPI_SUNREALTYPE,procleft,0,  // Receive from procleft
  //               d->comm,&request);
  // MPI_Wait(&request,&status);

  printf("checkpoint 1e\n");
  /* Move received data from host to GPU */
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Memcpy)(&d->bufs_dev[0],&d->bufs[0],sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
  NV_ADD_LANG_PREFIX_PH(Memcpy)(&d->bufs_dev[3],&d->bufs[3],sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
#endif

  // /* Pass needed data to processes before and after current process. */
  // if (d->myproc != procfirst)
  //   MPI_Send(&d->bufs[1], 1, MPI_SUNREALTYPE, procleft, 0, d->comm);           // Send to left
  // if (d->myproc != proclast)
  //   MPI_Send(&d->bufs[local_M], 1, MPI_SUNREALTYPE, procright, 0, d->comm);  // Send to right

  // /* Receive needed data from processes before and after current process. */
  // if (d->myproc != procfirst)
  //   MPI_Recv(&d->bufs[0], 1, MPI_SUNREALTYPE, procleft, 0, d->comm, &status);  // Recv from left
  // else
  //   d->bufs[0] = ZERO;
  // if (d->myproc != proclast)
  //   MPI_Recv(&d->bufs[d->local_M+1], 1, MPI_SUNREALTYPE, procright, 0, d->comm, // Recv from right
  //            &status);
  // else
  //   d->bufs[d->local_M + 1] = ZERO;

  printf("checkpoint 1f\n");
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  int      i;
  realtype uleft, uright, ucenter, hdiff, hadv;
  /* Loop over all grid points in current process. */
  for (i=1; i<=d->local_M; i++) {
    /* Extract u at x_i and two neighboring points */
    uleft   = (i==0           ) ? ubufs[0] : udata[i-1];
    uright  = (i==d->local_M-1) ? ubufs[3] : udata[i+1];
    ucenter = udata[i];

    /* Set diffusion and advection terms and load into udot */
    hdiff       = d->hdcoef*(uleft  - TWO*ucenter + uright);
    hadv        = d->hacoef*(uright - uleft);
    udotdata[i] = hdiff + hadv;
  }
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  unsigned block = 256;
  unsigned grid  = (d->local_M + block - 1) / block;
  fKernel<<<grid,block>>>(udata, udotdata, d->bufs_dev, d->local_M, d->hdcoef, d->hacoef);
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
