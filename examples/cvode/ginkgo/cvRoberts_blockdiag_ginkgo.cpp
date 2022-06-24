/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * The following is a simple example problem based off of cvRoberts_klu.c. We
 * simulate a scenario where a set of independent ODEs are grouped together to
 * form a larger system. For simplicity, each set of ODEs is the same problem.
 * The problem is from chemical kinetics, and consists of the following three
 * rate equations: dy1/dt = -.04*y1 + 1.e4*y2*y3 dy2/dt = .04*y1 - 1.e4*y2*y3 -
 * 3.e7*(y2)^2 dy3/dt = 3.e7*(y2)^2 on the interval from t = 0.0 to t = 4.e10,
 * with initial conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff. This
 * program solves the problem with the BDF method, Newton iteration, a
 * user-supplied Jacobian routine, and since the grouping of the independent
 * systems results in a block diagonal linear system, with the Gingko
 * SUNLinearSolver which supports batched iterative methods. It uses a scalar
 * relative tolerance and a vector absolute tolerance. Output is printed in
 * decades from t = .4 to t = 4.e10. Run statistics (optional outputsare
 * printed at the end.
 *
 * The program takes one optional argument, the number of groups
 * of independent ODE systems:
 *
 *    ./cvRoberts_blockdiag_ginkgo [number of groups]
 *
 * This problem is comparable to the cvRoberts_block_klu.c example.
 * --------------------------------------------------------------------------*/

#include <stdio.h>

#include <cvode/cvode.h>                              /* prototypes for CVODE fcts., consts.           */
#include <sunmatrix/sunmatrix_ginkgoblock.hpp>           /* access to the MAGMA dense SUNMatrix           */
#include <sunlinsol/sunlinsol_ginkgoblock.hpp>           /* access to MAGMA dense SUNLinearSolver         */

#if defined(SUNDIALS_GINKGO_BACKENDS_HIP)
#define HIP_OR_CUDA(a,b) a
#elif defined(SUNDIALS_GINKGO_BACKENDS_CUDA)
#define HIP_OR_CUDA(a,b) b
#else
#define HIP_OR_CUDA(a,b((void)0);
#endif

#if defined(SUNDIALS_GINKGO_BACKENDS_CUDA)
#include <nvector/nvector_cuda.h>
#include <sunmemory/sunmemory_cuda.h>
#elif defined(SUNDIALS_GINKGO_BACKENDS_HIP)
#include <nvector/nvector_hip.h>
#include <sunmemory/sunmemory_hip.h>
#endif


/* Problem Constants */

#define GROUPSIZE 3            /* number of equations per group */
#define Y1    RCONST(1.0)      /* initial y components */
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1.0e-4)   /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-8)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-14)
#define ATOL3 RCONST(1.0e-6)
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.4)      /* first output time      */
#define TMULT RCONST(10.0)     /* output time factor     */
#define NOUT  12               /* number of output times */

#define ZERO  RCONST(0.0)

/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

__global__
static void f_kernel(realtype t, realtype* y, realtype* ydot,
                     int neq, int ngroups);

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

__global__
static void j_kernel(int ngroups, int nnzper, realtype* ydata, realtype *Jdata);

/* Private function to initialize Jacobian */

static int JacInit(SUNMatrix J);

/* Private function to output results */

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem, SUNLinearSolver LS);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

/* user data structure */
typedef struct {
  int ngroups;
  int neq;
  sunindextype nnzper;
} UserData;

using GkoMatType = gko::matrix::BatchCsr<sunrealtype, sunindextype>;
using SUNMatrixView = sundials::ginkgo::BlockMatrix<GkoMatType>;
using SUNLinearSolverView = sundials::ginkgo::BlockLinearSolver<gko::solver::BatchGmres<sunrealtype>, SUNMatrixView>;

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(int argc, char *argv[])
{
  realtype reltol, t, tout;
  realtype *ydata, *abstol_data;
  N_Vector y, abstol;
  // SUNMatrix A;
  // SUNLinearSolver LS;
  void *cvode_mem;
  int retval, iout;
  int neq, ngroups, groupj;
  UserData udata;

  y = abstol = NULL;
  cvode_mem = NULL;

  /* Create the SUNDIALS context */
  sundials::Context sunctx;

  /* Create the Ginkgo executor */
  auto gko_exec =
      HIP_OR_CUDA(gko::HipExecutor::create(0, gko::OmpExecutor::create(), true),
                  gko::CudaExecutor::create(0, gko::OmpExecutor::create(), true));

  /* Parse command line arguments */
  if (argc > 1) {
    ngroups = atoi(argv[1]);
  } else {
    ngroups = 100;
  }
  neq = ngroups * GROUPSIZE;

  udata.ngroups = ngroups;
  udata.nnzper = GROUPSIZE*GROUPSIZE;
  udata.neq = neq;

  /* Create CUDA or HIP vector of length neq for I.C. and abstol */
  y = HIP_OR_CUDA( N_VNew_Hip(neq, sunctx);,
                   N_VNew_Cuda(neq, sunctx); )
  if (check_retval((void *)y, "N_VNew", 0)) { return(1); }
  abstol = N_VClone(y);
  if (check_retval((void *)abstol, "N_VClone", 0)) { return(1); }

  /* Initialize y */
  ydata = N_VGetArrayPointer(y);
  for (groupj = 0; groupj < neq; groupj += GROUPSIZE) {
    ydata[groupj]   = Y1;
    ydata[groupj+1] = Y2;
    ydata[groupj+2] = Y3;
  }
  HIP_OR_CUDA( N_VCopyToDevice_Hip(y);,
               N_VCopyToDevice_Cuda(y); )

  /* Set the scalar relative tolerance */
  reltol = RTOL;

  /* Set the vector absolute tolerance */
  abstol_data = N_VGetArrayPointer(abstol);
  for (groupj = 0; groupj < neq; groupj += GROUPSIZE) {
    abstol_data[groupj]   = ATOL1;
    abstol_data[groupj+1] = ATOL2;
    abstol_data[groupj+2] = ATOL3;
  }
  HIP_OR_CUDA( N_VCopyToDevice_Hip(abstol);,
               N_VCopyToDevice_Cuda(abstol); )

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) { return(1); }

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { return(1); }

  /* Call CVodeSetUserData to attach the user data structure */
  retval = CVodeSetUserData(cvode_mem, &udata);
  if (check_retval(&retval, "CVodeSetUserData", 1)) { return(1); }

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) { return(1); }

  /* Create SUNMatrix for use in linear solves */
  sunindextype block_rows = GROUPSIZE;
  sunindextype block_cols = GROUPSIZE;
  sunindextype nnzper = GROUPSIZE*GROUPSIZE;
  auto block_size = gko::batch_dim<2>(ngroups, gko::dim<2>(block_rows , block_cols));
  auto gko_matA = gko::share(GkoMatType::create(gko_exec, block_size, nnzper));
  SUNMatrixView A{gko_matA, sunctx};

  /* Initialiize the Jacobian with its fixed sparsity pattern */
  JacInit(A);

  /* Create the SUNLinearSolver object for use by CVode */
  auto precond_factory = gko::share(gko::preconditioner::BatchJacobi<sunrealtype>::build().on(gko_exec));
  SUNLinearSolverView LS{gko_exec, gko::stop::batch::ToleranceType::absolute,
                         nullptr, ngroups, sunctx};
  LS.setEnableScaling(true);

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS.get(), A.get());
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) { return(1); }

  /* Set the user-supplied Jacobian routine Jac */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
  
  //CVodeSetEpsLin(cvode_mem, 1e-10);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  printf(" \nGroup of independent 3-species kinetics problems\n\n");
  printf("number of groups = %d\n\n", ngroups);

  iout = 0;  tout = T1;
  while(1) {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    HIP_OR_CUDA( N_VCopyFromDevice_Hip(y);,
                 N_VCopyFromDevice_Cuda(y); )
    for (groupj = 0; groupj < ngroups; groupj += 10) {
      printf("group %d: ", groupj);
      PrintOutput(t, ydata[GROUPSIZE*groupj],
                     ydata[1+GROUPSIZE*groupj],
                     ydata[2+GROUPSIZE*groupj]);
    }

    if (check_retval(&retval, "CVode", 1)) break;
    if (retval == CV_SUCCESS) {
      iout++;
      tout *= TMULT;
    }

    if (iout == NOUT) break;
  }

  /* Print some final statistics */
  PrintFinalStats(cvode_mem, LS);

  /* Free y and abstol vectors */
  N_VDestroy(y);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  return(0);
}


/*
 * Jacobian initialization routine. This sets the sparisty pattern of
 * the blocks of the Jacobian J(t,y) = df/dy. This is performed on the CPU,
 * and only occurs at the beginning of the simulation.
 */

static int JacInit(SUNMatrix J)
{
  auto Jgko = static_cast<SUNMatrixView*>(J->content)->gkomtx();
  auto gko_exec = Jgko->get_executor();

  sunindextype* rowptrs = Jgko->get_row_ptrs();
  sunindextype* colvals = Jgko->get_col_idxs();

  /* there are 3 entries per row */
  rowptrs[0] = 0;
  rowptrs[1] = 3;
  rowptrs[2] = 6;
  rowptrs[3] = 9;

  /* first row of block */
  colvals[0] = 0;
  colvals[1] = 1;
  colvals[2] = 2;

  /* second row of block */
  colvals[3] = 0;
  colvals[4] = 1;
  colvals[5] = 2;

  /* third row of block */
  colvals[6] = 0;
  colvals[7] = 1;
  colvals[8] = 2;

  return(0);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* Right hand side function. This just launches the CUDA or HIP kernel
   to do the actual computation. At the very least, doing this
   saves moving the vector data in y and ydot to/from the device
   every evaluation of f. */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData *udata;
  realtype *ydata, *ydotdata;

  udata = (UserData*) user_data;
  ydata = N_VGetDeviceArrayPointer(y);
  ydotdata = N_VGetDeviceArrayPointer(ydot);

  unsigned block_size = HIP_OR_CUDA( 64, 32 );
  unsigned grid_size = (udata->neq + block_size - 1) / block_size;
  f_kernel<<<grid_size, block_size>>>(t, ydata, ydotdata, udata->neq, udata->ngroups);

  HIP_OR_CUDA( hipDeviceSynchronize();, cudaDeviceSynchronize(); )
  HIP_OR_CUDA( hipError_t cuerr = hipGetLastError();,
               cudaError_t cuerr = cudaGetLastError(); )
  if (cuerr != HIP_OR_CUDA( hipSuccess, cudaSuccess )) {
    fprintf(stderr,
            ">>> ERROR in f: cudaGetLastError returned %s\n",
            HIP_OR_CUDA( hipGetErrorName(cuerr), cudaGetErrorName(cuerr) ));
    return(-1);
  }

  return(0);
}

/* Right hand side function evalutation kernel. */
__global__
static void f_kernel(realtype t, realtype* ydata, realtype* ydotdata,
                     int neq, int ngroups)
{
  realtype y1, y2, y3, yd1, yd3;
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int groupj = i*GROUPSIZE;

  if (i < neq) {
    y1 = ydata[groupj]; y2 = ydata[groupj+1]; y3 = ydata[groupj+2];

    yd1 = ydotdata[groupj]   = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
    yd3 = ydotdata[groupj+2] = RCONST(3.0e7)*y2*y2;
          ydotdata[groupj+1] = -yd1 - yd3;
  }
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy.
 * This is done on the GPU.
 */

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData *udata = (UserData*) user_data;
  realtype *Jdata, *ydata;
  int *rowptrs, *colvals;
  unsigned block_size, grid_size;
  auto Jgko = static_cast<SUNMatrixView*>(J->content)->gkomtx();

  Jdata = Jgko->get_values();
  rowptrs = Jgko->get_row_ptrs();
  colvals = Jgko->get_col_idxs();
  ydata = N_VGetDeviceArrayPointer(y);

  block_size = HIP_OR_CUDA( 64, 32 );
  grid_size = (udata->neq + block_size - 1) / block_size;

  j_kernel<<<grid_size, block_size>>>(udata->ngroups, udata->nnzper, ydata, Jdata);

  HIP_OR_CUDA( hipDeviceSynchronize();, cudaDeviceSynchronize(); )
  HIP_OR_CUDA( hipError_t cuerr = hipGetLastError();,
               cudaError_t cuerr = cudaGetLastError(); )
  if (cuerr != HIP_OR_CUDA( hipSuccess, cudaSuccess )) {
    fprintf(stderr,
            ">>> ERROR in Jac: cudaGetLastError returned %s\n",
            HIP_OR_CUDA( hipGetErrorName(cuerr), cudaGetErrorName(cuerr) ));
    return(-1);
  }

  return(0);
}

/* Jacobian evaluation GPU kernel */
__global__
static void j_kernel(int ngroups, int nnzper, realtype* ydata, realtype *Jdata)
{
  int groupj;
  realtype y2, y3;

  for (groupj = blockIdx.x*blockDim.x + threadIdx.x;
       groupj < ngroups;
       groupj += blockDim.x * gridDim.x)
  {
    /* get y values */
    y2 = ydata[GROUPSIZE*groupj + 1];
    y3 = ydata[GROUPSIZE*groupj + 2];

    /* first row of block */
    Jdata[nnzper*groupj]       = RCONST(-0.04);
    Jdata[nnzper*groupj + 1]   = RCONST(1.0e4)*y3;
    Jdata[nnzper*groupj + 2]   = RCONST(1.0e4)*y2;

    /* second row of block */
    Jdata[nnzper*groupj + 3]   = RCONST(0.04);
    Jdata[nnzper*groupj + 4]   = (RCONST(-1.0e4)*y3) - (RCONST(6.0e7)*y2);
    Jdata[nnzper*groupj + 5]   = RCONST(-1.0e4)*y2;

    /* third row of block */
    Jdata[nnzper*groupj + 6]   = ZERO;
    Jdata[nnzper*groupj + 7]   = RCONST(6.0e7)*y2;
    Jdata[nnzper*groupj + 8]   = ZERO;

    // printf("%.6f %.6f\n", y2, y3);
    // printf("%.6f ", Jdata[nnzper*groupj]);
    // printf("%.6f ", Jdata[nnzper*groupj + 1]);
    // printf("%.6f\n", Jdata[nnzper*groupj + 2]);
    // printf("%.6f ", Jdata[nnzper*groupj + 3]);
    // printf("%.6f ", Jdata[nnzper*groupj + 4]);
    // printf("%.6f\n", Jdata[nnzper*groupj + 5]);
    // printf("%.6f ", Jdata[nnzper*groupj + 6]);
    // printf("%.6f ", Jdata[nnzper*groupj + 7]);
    // printf("%.6f\n\n", Jdata[nnzper*groupj + 8]);
  }
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}

/*
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem, SUNLinearSolver LS)
{
  long int nst, nfe, nsetups, nje, nni, ncfn, netf, nge;
  int retval;

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

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);

  retval = CVodeGetNumGEvals(cvode_mem, &nge);
  check_retval(&retval, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld    nge = %ld\n \n",
   nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s(failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

