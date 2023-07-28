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
 * Homogeneous Dirichlet boundary conditions are imposed, and the
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
 * Computational approach notes:
 *  (1) Two hypre nvectors are used: u and udot. Depending on the
 *      backend (SUNDIALS_HYPRE_BACKENDS), their data will be stored
 *      either on the host or on the device.
 *  (2) The interior M gridpoints are divided amongst the nprocs
 *      processes so that local_M is either k+1 or k for all procs.
 *  (3) Buffers: The four buffers used for MPI communications are
 *      stored together in the array 'bufs'. The entries correspond
 *      to [Lrecvbuf,Lsendbuf,Rsendbuf,Rrecvbuf]. If a device is in
 *      use, we similarly allocate 'bufs_dev', though only the Recv
 *      buffers are used by the device.
 *  (4) The diagram below shows the relationship between the nvector
 *      'u' (or 'udot') and the 4-element 'bufs' array used for MPI
 *      transfers. When using CUDA/HIP, there is an additional array
 *      'bufs_dev' on the device.
 *  proc           p-1                  p                  p+1
 *              0 1    k-1    ┊     0 1    k-1    ┊     0 1    k-1  
 *  'udata'   × ■ ■ ┄ ■ ■ ×   ┊   × ■ ■ ┄ ■ ■ ×   ┊   × ■ ■ ┄ ■ ■ × 
 *            ^ v       v ^   ┊   ^ v       v ^   ┊   ^ v       v ^ 
 *  'bufs'    ■ ■       ■ ■  <->  ■ ■       ■ ■  <->  ■ ■       ■ ■ 
 *            0 1       2 3  MPI  0 1       2 3  MPI  0 1       2 3 
 *  send/recv R S       S R   ┊   R S       S R   ┊   R S       S R 
 * -----------------------------------------------------------------
 * GPU usage notes:
 *  (1) As written, this example should work with SERIAL, CUDA, or
 *      HIP backends (determined by SUNDIALS_HYPRE_BACKENDS).
 *  (2) The macro NV_ADD_LANG_PREFIX_PH(...) adds the prefix "cuda"
 *      or "hip" depending on the backend. The user may use this
 *      macro for portability, or explicitly write cuda/hip calls.
 *      Example: "NV_ADD_LANG_PREFIX_PH(Free)" expands to "cudaFree"
 *      when using the CUDA hypre backend.
 *  (3) The MPI installation is assumed to *not* be CUDA/HIP-aware,
 *      so no device-to-device Send/Recv's are used here.
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

#include <mpi.h>                                  /* MPI constants and types */

/* Problem Constants */

#define ZERO  RCONST(0.0)
#define HALF  RCONST(0.5)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define XMAX  RCONST(2.0)    /* domain boundary           */
#define ATOL  RCONST(1.0e-5) /* scalar absolute tolerance */
#define T0    ZERO           /* initial time              */
#define T1    RCONST(0.5)    /* first output time         */
#define DTOUT RCONST(0.5)    /* output time increment     */
#define NOUT  10             /* number of output times    */

#define MPI_ASSERT(expr,msg,comm,myproc,code) \
  if(!expr) {                                 \
    if (myproc==0) printf(msg);               \
    MPI_Abort(comm,code);                     \
  }

/* Type : UserData
   contains grid constants, parhyp machine parameters, work array. */

typedef struct {
  MPI_Comm  comm;               // MPI communicator
  int       M, nprocs, myproc;  // global problem size, # of processes and my process id
  HYPRE_Int local_M, mybase;    // local problem size and global starting index
  realtype  dx, hdcoef, hacoef; // gridpoint spacing and problem coefficients
  realtype  *bufs;              // Send/Recv buffers: [Lrecvbuf,Lsendbuf,Rsendbuf,Rrecvbuf]
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  realtype  *bufs_dev;          // Buffer space (4 realtypes) to be allocated on device
#endif
} *UserData;

/* Private Helper Functions */

static void InitUserData(UserData data, MPI_Comm comm, HYPRE_Int M);

static void FreeUserData(UserData data);

static void SetIC(HYPRE_IJVector Uij, UserData data);

static void PrintIntro(int M, int nprocs);

static void PrintData(realtype t, realtype umax, long int nst);

static void PrintFinalStats(void *cvode_mem);

/* GPU Kernels */

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
__global__ void fKernel(realtype *u, realtype *udot,
                        HYPRE_Int local_M, realtype hdcoef, realtype hacoef);
#endif

/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  void*              cvode_mem;
  int                M, nprocs, myproc, iout, retval;
  long int           nst;
  realtype           reltol, abstol, t, tout, umax;
  
  UserData           data;
  N_Vector           u;
  SUNNonlinearSolver NLS;
  SUNContext         sunctx;

  MPI_Comm           comm;
  HYPRE_ParVector    Upar; /* Declare HYPRE parallel vector */
  HYPRE_IJVector     Uij;  /* Declare "IJ" interface to HYPRE vector */
  
/* Optional: create stream and choose policies to use instead of default */
// #if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
//   NV_ADD_LANG_PREFIX_PH(Stream_t) stream;
//   NV_ADD_LANG_PREFIX_PH(Error_t)  gpuerr;
// #endif

  data      = NULL;
  u         = NULL;
  cvode_mem = NULL;

  /* MPI Init */
  MPI_Init(&argc, &argv);

  /* Get MPI info */
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myproc);

  /* Parse inputs (required: M) */
  MPI_ASSERT(argc>=2,"ERROR: ONE (1) input required: M (# of interior gridpoints)\n",comm,myproc,1)
  M = (HYPRE_Int) atoi(argv[1]);
  MPI_ASSERT(M>=nprocs,"ERROR: M < nprocs (We require at least one interior gridpoint per process)\n",comm,myproc,1)

  /* Allocate UserData */
  data = (UserData) malloc(sizeof *data);
  if(check_retval((void *)data, "malloc", 2, myproc)) MPI_Abort(comm, 1);

  /* Initialize UserData */
  InitUserData(data, comm, M);

  /* Create SUNDIALS context */
  retval = SUNContext_Create(&comm, &sunctx);
  if (check_retval(&retval, "SUNContex_Create", 1, myproc)) MPI_Abort(comm, 1);

  /* Allocate hypre vector */
  HYPRE_IJVectorCreate(comm, data->mybase, data->mybase + data->local_M - 1, &Uij);
  HYPRE_IJVectorSetObjectType(Uij, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(Uij);

  /* Set tolerances */
  reltol = ZERO;
  abstol = ATOL;

  /* Initialize solution vector. */
  SetIC(Uij, data);
  HYPRE_IJVectorAssemble(Uij);
  HYPRE_IJVectorGetObject(Uij, (void**) &Upar);

  u = N_VMake_ParHyp(Upar, sunctx);  /* Create wrapper u around hypre vector */
  if(check_retval((void *)u, "N_VNew", 0, myproc)) MPI_Abort(comm, 1);

  /* Call CVodeCreate to create the solver memory and specify the
   * Adams-Moulton LMM */
  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0, myproc)) MPI_Abort(comm, 1);

  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1, myproc)) MPI_Abort(comm, 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1, myproc)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1, myproc)) return(1);

  /* Create fixed point nonlinear solver object */
  NLS = SUNNonlinSol_FixedPoint(u, 0, sunctx);
  if(check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0, myproc)) return(1);

  /* Attach nonlinear solver object to CVode */
  retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if(check_retval(&retval, "CVodeSetNonlinearSolver", 1, myproc)) return(1);

  if (myproc == 0) PrintIntro(M, nprocs);

  /* Print initial statistics */
  umax = N_VMaxNorm(u);
  if (myproc == 0) {
    t = T0;
    PrintData(t, umax, 0);
  }

  /* In loop over output points, call CVode, print results, test for error */
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    // printf("checkpoint 1\n");
    /* Advance to next output time */
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_retval(&retval, "CVode", 1, myproc)) break;

    // printf("checkpoint 2\n");
    /* Get number of steps taken */
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    if(check_retval(&retval, "CVodeGetNumSteps", 1, myproc)) break;

    // printf("checkpoint 3\n");
    /* Print time, current max norm, and # steps taken since last output */
    umax = N_VMaxNorm(u);
    if (myproc == 0) PrintData(t, umax, nst);
  }

  /* Print some final statistics */
  if (myproc == 0) PrintFinalStats(cvode_mem);  

  /* Free memory */
  FreeUserData(data);         /* Free UserData */
  N_VDestroy(u);              /* Free hypre vector wrapper */
  HYPRE_IJVectorDestroy(Uij); /* Free the underlying hypre vector */
  CVodeFree(&cvode_mem);      /* Free the integrator memory */
  SUNNonlinSolFree(NLS);      /* Free the nonlinear solver */
  SUNContext_Free(&sunctx);   /* Free context */

  /* MPI Finalize */
  MPI_Finalize();

  return(0);
}

/************************ Private Helper Functions ***********************/

/* Load constants in data */

static void InitUserData(UserData data, MPI_Comm comm, HYPRE_Int M)
{
  int nperproc, nrem;

  /* Attach MPI info */
  data->comm = comm;
  MPI_Comm_size(comm, &data->nprocs);
  MPI_Comm_rank(comm, &data->myproc);

  /* Set partitioning */
  data->M       = M;
  nperproc      = M / data->nprocs;
  nrem          = M - data->nprocs*nperproc;
  data->local_M = (data->myproc < nrem) ? nperproc+1 : nperproc;
  data->mybase  = (data->myproc < nrem) ? data->myproc*data->local_M : data->myproc*nperproc + nrem;

  /* Populate constants */
  data->dx      = XMAX/((realtype)(M+1));
  data->hdcoef  = ONE/(data->dx*data->dx);
  data->hacoef  = HALF/(TWO*data->dx);

  /* Initialize Send/Recv buffers */
  // Note: We have zero Dirichlet b.c.; Bdry procs will always have zero in their outer Recv buffer
  data->bufs = (realtype*) malloc(4*sizeof(realtype));
  memset(data->bufs, 0, 4*sizeof(realtype));
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Malloc)(&data->bufs_dev,    4*sizeof(realtype)); // e.g. cudaMalloc
  NV_ADD_LANG_PREFIX_PH(Memset)( data->bufs_dev, 0, 4*sizeof(realtype)); // e.g. cudaMemset
#endif
}

/* Free user data memory */

static void FreeUserData(UserData data)
{
  free(data->bufs);
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Free)(data->bufs_dev); // e.g. cudaFree
#endif
  free(data);
}

/* Set initial conditions in u vector */

static void SetIC(HYPRE_IJVector Uij, UserData data)
{
  int i;
  HYPRE_Int *iglobal;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u. */
  iglobal = (HYPRE_Int*) malloc(data->local_M*sizeof(HYPRE_Int));
  udata   = (realtype*)  malloc(data->local_M*sizeof(realtype));

  /* Load initial profile into u vector */
  for (i = 0; i < data->local_M; i++) {
    iglobal[i] = data->mybase + i;
    x = (iglobal[i] + 1)*data->dx;
    udata[i] = x*(XMAX - x)*SUNRexp(TWO*x);
  }

#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  /* For serial backend, HYPRE_IJVectorSetValues expects host pointers */
  HYPRE_IJVectorSetValues(Uij, data->local_M, iglobal, udata);

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  /* With GPU backend, HYPRE_IJVectorSetValues expects device pointers */
  HYPRE_Int* iglobal_dev;
  realtype*  udata_dev;

  /* Allocate device arrays */
  NV_ADD_LANG_PREFIX_PH(Malloc)(&iglobal_dev, data->local_M*sizeof(HYPRE_Int));
  NV_ADD_LANG_PREFIX_PH(Malloc)(&udata_dev, data->local_M*sizeof(realtype));
  
  /* Copy host data to device */
  NV_ADD_LANG_PREFIX_PH(Memcpy)(iglobal_dev, iglobal, data->local_M*sizeof(HYPRE_Int), NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
  NV_ADD_LANG_PREFIX_PH(Memcpy)(udata_dev, udata, data->local_M*sizeof(realtype), NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));

  /* Set hypre vector values from device arrays, then free the device arrays */
  HYPRE_IJVectorSetValues(Uij, data->local_M, iglobal_dev, udata_dev);
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
  printf("\n Number of MPI processes = %3d \n\n", nprocs);

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
__global__ void fKernel(realtype* udata, realtype* udotdata, realtype* bufs_dev,
                        HYPRE_Int local_M, realtype hdcoef, realtype hacoef)
{
  realtype uleft, uright, ucenter, hdiff, hadv;
  sunindextype tid;

  /* Get tid */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  /* Get u values (take from Recv buffers as needed) */
  if (tid < local_M) {
    uleft   = (tid==0        ) ? bufs_dev[0] : udata[tid-1];
    uright  = (tid==local_M-1) ? bufs_dev[3] : udata[tid+1];
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
  // printf("checkpoint 1a\n");
  int       nprocs, myproc, procfirst, proclast, procleft, procright;
  realtype *udata, *udotdata; // Both are length M (# of gridpoints excluding bdry)

  UserData        data;
  MPI_Comm        comm;
  MPI_Status      status;
  MPI_Request     request;
  // MPI_Status      status_leftward_flow, status_rightward_flow;
  HYPRE_ParVector uhyp;
  HYPRE_ParVector udothyp;

  // printf("checkpoint 1b\n");
  /* Get reference to UserData */
  data     = (UserData) user_data;

  /* Extract MPI info */
  comm     = data->comm;
  nprocs   = data->nprocs;
  myproc   = data->myproc;

  /* Extract hypre vectors */
  uhyp     = N_VGetVector_ParHyp(u);
  udothyp  = N_VGetVector_ParHyp(udot);

  /* Access hypre vectors local data */
  udata    = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));
  udotdata = hypre_VectorData(hypre_ParVectorLocalVector(udothyp));

  // printf("checkpoint 1c\n");
  /* Relevant process ids (utilize behavior of MPI_PROC_NULL in Send/Recv's) */
  // Note: Send/Recv's with MPI_PROC_NULL immediately return successful without modifying buffers.
  procfirst = 0;
  proclast  = nprocs-1;
  procleft  = (myproc==procfirst) ? MPI_PROC_NULL : myproc-1;
  procright = (myproc==proclast ) ? MPI_PROC_NULL : myproc+1;

  /* Populate Send buffers */
  // Note: bufs order is [Lrecvbuf,Lsendbuf,Rsendbuf,Rrecvbuf]
  // Note: MPI transfers are host-to-host (if MPI is CUDA/HIP-aware, could do device-to-device)
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  data->bufs[1] = udata[0];          //  leftmost local gridpoint is sent left
  data->bufs[2] = udata[data->local_M]; // rightmost local gridpoint is sent right
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Memcpy)(&data->bufs[1],&udata[0],              sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyDeviceToHost));
  NV_ADD_LANG_PREFIX_PH(Memcpy)(&data->bufs[2],&udata[data->local_M-1],sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyDeviceToHost));
#endif

  // printf("checkpoint 1d\n");
  /* Nonblocking leftward flow of data */
  MPI_Irecv(&data->bufs[3],1,MPI_SUNREALTYPE,procright,0,comm,&request); // Receive from procright
  MPI_Send( &data->bufs[1],1,MPI_SUNREALTYPE,procleft ,0,comm);          // Send to procleft
  MPI_Wait( &request,&status);

  /* Nonblocking rightward flow of data */
  MPI_Irecv(&data->bufs[0],1,MPI_SUNREALTYPE,procleft ,0,comm,&request); // Receive from procleft
  MPI_Send( &data->bufs[2],1,MPI_SUNREALTYPE,procright,0,comm);          // Send to procright
  MPI_Wait( &request,&status);

  // printf("checkpoint 1e\n");
  /* Move received data from host to GPU */
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Memcpy)(&data->bufs_dev[0],&data->bufs[0],sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
  NV_ADD_LANG_PREFIX_PH(Memcpy)(&data->bufs_dev[3],&data->bufs[3],sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
#endif

  // printf("checkpoint 1f\n");
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  int      i;
  realtype uleft, uright, ucenter, hdiff, hadv;
  /* Loop over all grid points in current process. */
  for (i=0; i<data->local_M; i++) {
    /* Extract u at x_i and two neighboring points */
    uleft   = (i==0              ) ? data->bufs[0] : udata[i-1];
    uright  = (i==data->local_M-1) ? data->bufs[3] : udata[i+1];
    ucenter = udata[i];

    /* Set diffusion and advection terms and load into udot */
    hdiff       = data->hdcoef*(uleft  - TWO*ucenter + uright);
    hadv        = data->hacoef*(uright - uleft);
    udotdata[i] = hdiff + hadv;
  }
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  unsigned block = 256;
  unsigned grid  = (data->local_M + block - 1) / block;
  fKernel<<<grid,block>>>(udata, udotdata, data->bufs_dev, data->local_M, data->hdcoef, data->hacoef);
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
