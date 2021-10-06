/* -----------------------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ UIUC/LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Example problem:
 *
 * The following test implements the expectation-maximazation algorithm
 * for mixture densities.
 * Here, we consider a mixture density composed of 3 univariate normal
 * densities.
 *
 * .... MORE PROBLEM INFO ...
 *
 * Several command line options are available to change the problem parameters
 * and KINSOL settings. Use the flag --help for more information.
 * ---------------------------------------------------------------------------*/

#include "kin_nonlin_mpicuda.hpp"     // header file containing UserData
                                      //   and function declartions

//#define SERIALIZE

// -----------------------------------------------------------------------------
// Cuda Kernels
// -----------------------------------------------------------------------------

__global__
void PxKernel(realtype *mu, realtype *Px, realtype *x,
              realtype a1, realtype a2, realtype a3, realtype scale,
              sunindextype N)
{
  // Calculate all P(x_k) for each x value
  sunindextype tid;
  realtype val1, val2, val3;

  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < N) {
    val1 = x[tid] - mu[0];
    val2 = x[tid] - mu[1];
    val3 = x[tid] - mu[2];

    Px[tid] = a1 * scale * exp( -(val1 * val1)/TWO );
    Px[tid] += a2 * scale * exp( -(val2 * val2)/TWO );
    Px[tid] += a3 * scale * exp( -(val3 * val3)/TWO );
  }
}

__global__
void EMKernel(realtype *mu, realtype *mu_top, realtype *mu_bottom,
              realtype *x, realtype *Px,
              realtype a1, realtype a2, realtype a3, realtype scale,
              sunindextype N)
{
  sunindextype tid;
  realtype val1, val2, val3;
  realtype frac1, frac2, frac3;

  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < N) {
    val1 = x[tid] - mu[0];
    val2 = x[tid] - mu[1];
    val3 = x[tid] - mu[2];

    frac1 = a1 * scale * exp( -(val1 * val1)/TWO ) / Px[tid];
    frac2 = a2 * scale * exp( -(val2 * val2)/TWO ) / Px[tid];
    frac3 = a3 * scale * exp( -(val3 * val3)/TWO ) / Px[tid];

    // MAKE ALL OF THESE ATOMIC ADDS
    atomicAdd(mu_top,     x[tid] * frac1);
    atomicAdd(mu_top + 1, x[tid] * frac2);
    atomicAdd(mu_top + 2, x[tid] * frac3);

    atomicAdd(mu_bottom,     frac1);
    atomicAdd(mu_bottom + 1, frac2);
    atomicAdd(mu_bottom + 2, frac3);
  }
}

__global__
void EMKernelFin(realtype *mu, realtype *mu_top, realtype *mu_bottom,
                 sunindextype localn)
{
  sunindextype tid, max;

  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < localn) {
    mu[3*tid]   = mu_top[0] / mu_bottom[0];
    mu[3*tid+1] = mu_top[1] / mu_bottom[1];
    mu[3*tid+2] = mu_top[2] / mu_bottom[2];
  }
}

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int flag;                   // reusable error-checking flag
  UserData *udata    = NULL;  // user data structure
  N_Vector ulocal    = NULL;  // local vector for storing solution
  N_Vector u         = NULL;  // global vector for storing solution
  N_Vector scale     = NULL;  // vector for scaling initial value
  void *kin_mem      = NULL;  // KINSOL memory structure

  // Timing variables
  double t1 = 0.0;
  double t2 = 0.0;

  // MPI variables
  MPI_Comm comm_w = MPI_COMM_WORLD; // MPI communicator
  int myid;                         // MPI process ID

  // Initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;

  flag = MPI_Comm_rank(comm_w, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

  // Set output process flag
  bool outproc = (myid == 0);

  // ------------------------------------------
  // Setup UserData and parallel decomposition
  // ------------------------------------------

  // Allocate and initialize user data structure with default values. The
  // defaults may be overwritten by command line inputs in ReadInputs below.
  udata = new UserData;
  flag = InitUserData(udata);
  if (check_flag(&flag, "InitUserData", 1)) return 1;

  // Parse command line inputs
  flag = ReadInputs(&argc, &argv, udata, outproc);
  if (flag != 0) return 1;

  // Output problem setup/options
  if (outproc)
  {
    flag = PrintUserData(udata);
    if (check_flag(&flag, "PrintUserData", 1)) return 1;
  }

  // --------------------------
  // Create MPI + Cuda vectors
  // --------------------------

  // Create vector for solution
  ulocal = N_VNew_Cuda(3 * udata->nodes_loc);
  if (check_flag((void *) ulocal, "N_VNew_Cuda", 0)) return 1;

  u = N_VMake_MPIPlusX(udata->comm, ulocal);
  if (check_flag((void *) u, "N_VMake_MPIPlusX", 0)) return 1;

  // Enable fused vector operations for low sync orthogonalization routines
  // Turned off by default
  if (udata->fusedops)
  {
    flag = N_VEnableFusedOps_Cuda(ulocal, SUNTRUE);
    if (check_flag(&flag, "N_VEnableFusedOps_Cuda", 1)) return 1;

    flag = N_VEnableFusedOps_MPIPlusX(u, SUNTRUE);
    if (check_flag(&flag, "N_VEnableFusedOps_MPIPlusX", 1)) return 1;
  }
  // Disable linear combination because it requires a copy to device
  flag = N_VEnableLinearCombination_Cuda(ulocal, false);
  if (check_flag(&flag, "N_VEnableLinearCombination", 1)) return 1;

  // Create vector for scaling initial value
  scale = N_VClone(u);
  if (check_flag((void *) scale, "N_VClone", 0)) return 1;
  N_VConst(ONE, scale);

  // Set initial condition
  flag = RandomVec(u, udata);
  if (check_flag(&flag, "RandomVec", 1)) return 1;

  // Create vector true mu values
  udata->mu_true = N_VClone(u);
  if (check_flag((void *) (udata->mu_true), "N_VClone", 0)) return 1;

  // Create temporary vector for residual and error output
  udata->vtemp = N_VClone(u);
  if (check_flag((void *) (udata->vtemp), "N_VClone", 0)) return 1;

  // Clone ulocal for temporary vector
  udata->mu_bottom = N_VNew_Cuda(3);
  if (check_flag((void *) (udata->mu_bottom), "N_VNewCuda", 0)) return 1;

  udata->mu_top = N_VNew_Cuda(3);
  if (check_flag((void *) (udata->mu_top), "N_VNewCuda", 0)) return 1;

  // Create vector for samples
  udata->samples_local = N_VNew_Cuda(udata->num_samples);
  if (check_flag((void *) udata->samples_local, "N_VNew_Cuda", 0)) return 1;
  // Clone samples for temporary vector
  udata->px = N_VClone(udata->samples_local);
  if (check_flag((void *) (udata->px), "N_VClone", 0)) return 1;

  // --------------
  // Setup Mus
  // --------------
  flag = SetMus(udata);
  if (check_flag(&flag, "SetMus", 1)) return 1;

  // --------------
  // Setup Samples
  // --------------
  flag = SetupSamples(udata);
  if (check_flag(&flag, "SetupSamples", 1)) return 1;

  // --------------
  // Setup KINSOL
  // --------------

  // Initialize KINSOL memory
  kin_mem = KINCreate();
  if (check_flag((void *) kin_mem, "KINCreate", 0)) return 1;

  // Set number of prior residuals used in Anderson Accleration
  flag = KINSetMAA(kin_mem, udata->maa);
  if (check_flag(&flag, "KINSetMAA", 0)) return 1;

  // Set orthogonlization routine used in Anderson Accleration
  flag = KINSetOrthAA(kin_mem, udata->orthaa);
  if (check_flag(&flag, "KINSetOrthAA", 0)) return 1;

  // Set Fixed Point Function
  flag = KINInit(kin_mem, FPFunction, u);
  if (check_flag(&flag, "KINInit", 1)) return 1;

  // Specify tolerances
  flag = KINSetFuncNormTol(kin_mem, udata->rtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return 1;

  // Set maximum number of iterations
  flag = KINSetNumMaxIters(kin_mem, udata->maxits);
  if (check_flag(&flag, "KINSetMaxNumIters", 1)) return 1;

  // Set Anderson Acceleration damping parameter
  flag = KINSetDampingAA(kin_mem, udata->damping);
  if (check_flag(&flag, "KINSetDampingAA", 1)) return 1;

  // Attach user data
  flag = KINSetUserData(kin_mem, (void *) udata);
  if (check_flag(&flag, "KINSetUserData", 1)) return 1;

  // Set debugging output file
  FILE* debugfp;
  if (udata->debug)
  {
    char fname[MXSTR];
    snprintf(fname, MXSTR, "kinsol_output_%06d.txt", myid);
    debugfp = fopen(fname,"w");

    flag = KINSetDebugFile(kin_mem, debugfp);
    if (check_flag(&flag, "KINSetDebugFile", 1)) return 1;
  }

  // ----------------------------
  // Call KINSol to solve problem
  // ----------------------------

  // No scaling used
  N_VConst(ONE, scale);

  if (udata->output > 1)
  {
    flag = OpenOutput(udata);
    if (check_flag(&flag, "OpenOutput", 1)) return 1;
  }

  // Start timer
  t1 = MPI_Wtime();

  // Call main solver
  flag = KINSol(kin_mem,        // KINSol memory block
                u,              // inital guess on input; solution vector
                KIN_FP,         // global strategy choice
                scale,          // scaling vector, for the variable u
                scale);         // scaling vector for function values fval
  if (check_flag(&flag, "KINSol", 1)) return(1);

  // Stop timer
  t2 = MPI_Wtime();

  // Update timer
  udata->totaltime = t2 - t1;

  // -----------------------
  // Get solver statistics
  // -----------------------
  if (udata->output > 0 && outproc)
  {
    cout << "Final statistics:" << endl;
    flag = OutputStats(kin_mem, udata);
    if (check_flag(&flag, "OutputStats", 1)) return 1;
  }
  if (udata->output > 1)
  {
    flag = CloseOutput(udata);
    if (check_flag(&flag, "CloseOutput", 1)) return 1;
  }

  // ------------------------
  // Calculate solution error
  // ------------------------
  // Output final error
  flag = SolutionError(udata->mu_true, u, udata->vtemp, udata);
  if (check_flag(&flag, "SolutionError", 1)) return 1;

  realtype maxerr = N_VMaxNorm(udata->vtemp);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(numeric_limits<realtype>::digits10);
    cout << "  Max error = " << maxerr << endl;
    cout << endl;
  }

  // ------------------------------
  // Print timing
  // ------------------------------
  if (udata->timing)
  {
    flag = OutputTiming(udata);
    if (check_flag(&flag, "OutputTiming", 1)) return 1;
  }

  // ------------------------------
  // Free memory
  // ------------------------------

  if (udata->debug) fclose(debugfp);
  KINFree(&kin_mem);         // Free solver memory
  N_VDestroy(u);             // Free vectors
  N_VDestroy(scale);
  FreeUserData(udata);       // Free user data
  delete udata;
  flag = MPI_Finalize();     // Finalize MPI
  return 0;

}

// -----------------------------------------------------------------------------
// Functions called by the solver
// -----------------------------------------------------------------------------

// Fixed point function to compute G(u)
static int FPFunction(N_Vector u, N_Vector f, void *user_data)
{
  int          flag;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Start timer
  double t1 = MPI_Wtime();

  // Call EM Algorithm
  flag = EM(u, f, user_data);
  if (check_flag(&flag, "EM", 1)) return -1;

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->fevaltime += t2 - t1;

  // Calculate and output residual and error history
  if (udata->output > 1)
  {
    flag = WriteOutput(u, f, udata);
    if (check_flag(&flag, "WriteOutput", 1)) return 1;
  }

  // Return success
  return 0;
}

// Setup mean distribution samples
static int SetupSamples(UserData *udata)
{
  sunindextype i, j, start, k, end;
  realtype mean, val;

  // Access problem data
  realtype *samples_local = N_VGetHostArrayPointer_Cuda(udata->samples_local);
  if (check_flag((void *) samples_local, "N_VGetHostArrayPointer_Cuda", 0)) return 1;
  realtype *mu_host = N_VGetHostArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(udata->mu_true));
  if (check_flag((void *) mu_host, "N_VGetHostArrayPointer_Cuda", 0)) return 1;

  realtype std_dev = ONE;

  for (i = 0; i < 3; i++) {
    // Set number of samples with this mean
    if (i == 0 || i == 1) {
      end = udata->num_samples * PTTHREE;
      start = i * end;
      end += start;
    }
    else {
      end = udata->num_samples * PTFOUR;
      start = 2 * udata->num_samples * PTTHREE;
      end += start;
    }

    // Setup distribution parameters
    mean = mu_host[i];
    std::default_random_engine generator;
    std::normal_distribution<realtype> distribution(mean, std_dev);

    // Get samples
    for (j = start; j < end; j++) {
      val = distribution(generator);
      samples_local[j] = val;
    }
  }

  N_VCopyToDevice_Cuda(udata->samples_local);

  // Return success
  return 0;
}

// Fill the vector u with random data
static int SetMus(UserData *udata)
{
  sunindextype i;
  realtype increment1, increment2;

  realtype *mu_host = N_VGetHostArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(udata->mu_true));
  if (check_flag((void *) mu_host, "N_VGetHostArrayPointer_Cuda", 0)) return 1;

  increment1 = (FIVE - HALF)/udata->nodes_loc;
  increment2 = (TEN - ONE)/udata->nodes_loc;

  // Fill vectors with uniform random data in [-1,1]
  for (i = 0; i < udata->nodes_loc; i++) {
    /*mu_host[3*i]   = ZERO;
    mu_host[3*i+1] = HALF + i * increment1;
    mu_host[3*i+2] = ONE + i * increment2;*/

    /*mu_host[3*i]   = ZERO;
    mu_host[3*i+1] = HALF;
    mu_host[3*i+2] = ONE;*/

    mu_host[3*i]   = ZERO;
    mu_host[3*i+1] = FIVE;
    mu_host[3*i+2] = TEN;
  }

  N_VCopyToDevice_Cuda(N_VGetLocalVector_MPIPlusX(udata->mu_true));

  // Return success
  return 0;
}

static int SetStartGuess(N_Vector u, UserData* udata)
{
  sunindextype i;

  realtype *u_host = N_VGetHostArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(u));
  if (check_flag((void *) u_host, "N_VGetHostArrayPointer_Cuda", 0)) return 1;

  // Fill vectors with uniform random data in [-1,1]
  for (i = 0; i < udata->nodes_loc; i++) {
    u_host[3*i]   = 0.25;
    u_host[3*i+1] = 3.0;
    u_host[3*i+2] = 0.75;
  }

  N_VCopyToDevice_Cuda(N_VGetLocalVector_MPIPlusX(u));

  // Return success
  return 0;
}

// Fill the vector u with random data between [0,1] as starting guess
static int RandomVec(N_Vector u, void *user_data)
{
  int flag;
  sunindextype i;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  realtype *uarray_local = N_VGetHostArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(u));
  if (check_flag((void *) uarray_local, "N_VGetHostArrayPointer_Cuda", 0)) return -1;

  // Fill vectors with uniform random data in [0,1]
  realtype val;
  for (sunindextype j = 0; j < 3; j++) {
    val = ((realtype) rand() / (realtype) RAND_MAX);
    for (i = 0; i < udata->nodes_loc; i++) {
      uarray_local[3*i + j] = val;
    }
  }

  N_VCopyToDevice_Cuda(N_VGetLocalVector_MPIPlusX(u));

  // Return success
  return 0;
}

static int EM(N_Vector u, N_Vector f, void *user_data)
{
  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Set grid and block sizes for kernel launch
  unsigned block = 256;
  unsigned grid1 = (udata->num_samples + block - 1) / block;
  unsigned grid2 = (udata->nodes_loc + block - 1) / block;

  // ---------
  // PX KERNEL
  // ---------

  // Scale value for functions
  realtype scale = ONE / sqrt(TWO * PI);

  // Get input device pointers
  realtype *u_dev  = N_VGetDeviceArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(u));
  realtype *x_dev  = N_VGetDeviceArrayPointer_Cuda(udata->samples_local);

  // Get output device pointer
  realtype *Px_dev = N_VGetDeviceArrayPointer_Cuda(udata->px);

  // Compute Px
  PxKernel<<<grid1, block>>>(u_dev, Px_dev, x_dev,
                             udata->alpha1, udata->alpha2, udata->alpha3, scale,
                             udata->num_samples);

  // ---------
  // EM KERNEL
  // ---------

  // Get output device pointers
  realtype *mu_bottom_dev = N_VGetDeviceArrayPointer_Cuda(udata->mu_bottom);
  realtype *mu_top_dev    = N_VGetDeviceArrayPointer_Cuda(udata->mu_top);

#ifdef SERIALIZE
  // Get input host pointers
  N_VCopyFromDevice_Cuda(N_VGetLocalVector_MPIPlusX(u));
  N_VCopyFromDevice_Cuda(udata->samples_local);
  N_VCopyFromDevice_Cuda(udata->px);

  realtype *u_host  = N_VGetHostArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(u));
  realtype *x_host  = N_VGetHostArrayPointer_Cuda(udata->samples_local);
  realtype *px_host = N_VGetHostArrayPointer_Cuda(udata->px);

  // Initilaize output vectors to zero (for sum reduction)
  mub_host[0] = ZERO;
  mub_host[1] = ZERO;
  mub_host[2] = ZERO;

  mut_host[0] = ZERO;
  mut_host[1] = ZERO;
  mut_host[2] = ZERO;

  // Get output host pointers
  realtype *mub_host = N_VGetHostArrayPointer_Cuda(udata->mu_bottom);
  realtype *mut_host = N_VGetHostArrayPointer_Cuda(udata->mu_top);

  // Serial EM kernel calculation
  realtype frac1, frac2, frac3;
  realtype val1, val2, val3;
  realtype a1 = udata->alpha1;
  realtype a2 = udata->alpha2;
  realtype a3 = udata->alpha3;

  for (int i = 0; i < udata->num_samples; i++)
  {
    val1 = x_host[i] - u_host[0];
    val2 = x_host[i] - u_host[1];
    val3 = x_host[i] - u_host[2];

    frac1 = a1 * scale * exp( -(val1 * val1)/TWO ) / px_host[i];
    frac2 = a2 * scale * exp( -(val2 * val2)/TWO ) / px_host[i];
    frac3 = a3 * scale * exp( -(val3 * val3)/TWO ) / px_host[i];

    mut_host[0] += x_host[i] * frac1;
    mut_host[1] += x_host[i] * frac2;
    mut_host[2] += x_host[i] * frac3;

    mub_host[0] += frac1;
    mub_host[1] += frac2;
    mub_host[2] += frac3;
  }

  // Copy output vectors to the device
  N_VCopyToDevice_Cuda(udata->mu_bottom);
  N_VCopyToDevice_Cuda(udata->mu_top);
#else
  // Initilaize output vectors to zero (for sum reduction)
  N_VConst(ZERO, udata->mu_bottom);
  N_VConst(ZERO, udata->mu_top);

  EMKernel<<<grid1, block>>>(u_dev, mu_top_dev, mu_bottom_dev, x_dev, Px_dev,
                             udata->alpha1, udata->alpha2, udata->alpha3, scale,
                             udata->num_samples);
#endif

  // ------------------
  // EM FINALIZE KERNEL
  // ------------------

  realtype *f_dev = N_VGetDeviceArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(f));

  EMKernelFin<<<grid2, block>>>(f_dev, mu_top_dev, mu_bottom_dev,
                                udata->nodes_loc);

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Initialize memory allocated within Userdata
static int InitUserData(UserData *udata)
{
  int flag;

  // Sigmas
  udata->sigma = ONE;

  // Alphas - mixture proportions
  udata->alpha1 = PTTHREE;
  udata->alpha2 = PTTHREE;
  udata->alpha3 = PTFOUR;

  // MPI variables
  udata->comm = MPI_COMM_WORLD;

  // Get the number of processes
  flag = MPI_Comm_size(udata->comm, &(udata->nprocs_w));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << flag << endl;
    return -1;
  }

  // Get my rank
  flag = MPI_Comm_rank(udata->comm, &(udata->myid));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << flag << endl;
    return -1;
  }

  // Local number of nodes
  udata->nodes_loc = 1;

  // Global total number of nodes
  udata->nodes = udata->nodes_loc * udata->nprocs_w;

  // Integrator settings
  udata->rtol        = RCONST(1.e-8);   // relative tolerance
  udata->maa         = 3;               // 3 vectors in Anderson Acceleration space
  udata->damping     = ONE;             // no damping for Anderson Acceleration
  udata->orthaa      = 0;               // use MGS for Anderson Acceleration
  udata->maxits      = 200;             // max number of fixed point iterations

  // Vectors
  udata->samples_local = NULL;
  udata->px            = NULL;
  udata->mu_bottom     = NULL;
  udata->mu_top        = NULL;
  udata->mu_true       = NULL;

  // Number samples
  udata->num_samples = 100000;

  // Output variables
  udata->output = 1;   // 0 = no output, 1 = stats output, 2 = output to disk
  udata->vtemp  = NULL;

  // Timing variables
  udata->timing       = false;
  udata->totaltime    = 0.0;
  udata->fevaltime    = 0.0;

  // Fused operations
  udata->fusedops = false;

  udata->debug = false;

  // Return success
  return 0;
}

// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{

  // Free samples vectors
  if (udata->samples_local)
  {
    N_VDestroy(udata->samples_local);
    udata->samples_local = NULL;
  }

  // Free temporary vectors
  if (udata->px)
  {
    N_VDestroy(udata->px);
    udata->px = NULL;
  }
  if (udata->mu_bottom)
  {
    N_VDestroy(udata->mu_bottom);
    udata->mu_bottom = NULL;
  }
  if (udata->mu_top)
  {
    N_VDestroy(udata->mu_top);
    udata->mu_top = NULL;
  }
  if (udata->mu_true)
  {
    N_VDestroy(udata->mu_true);
    udata->mu_true = NULL;
  }

  // Free error vector
  if (udata->vtemp)
  {
    N_VDestroy(udata->vtemp);
    udata->vtemp = NULL;
  }

  // Free MPI communicator
  udata->comm = MPI_COMM_NULL;

  // Return success
  return 0;
}

// Read command line inputs
static int ReadInputs(int *argc, char ***argv, UserData *udata, bool outproc)
{
  // Check for input args
  int arg_idx = 1;

  while (arg_idx < (*argc))
  {
    string arg = (*argv)[arg_idx++];

    // Mesh points
    if (arg == "--nodes_loc")
    {
      udata->nodes_loc = stoi((*argv)[arg_idx++]);
    }
    // Fixed Point settings
    else if (arg == "--rtol")
    {
      udata->rtol = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--maa")
    {
      udata->maa = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--damping")
    {
      udata->damping = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--orthaa")
    {
      udata->orthaa = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--maxits")
    {
      udata->maxits = stoi((*argv)[arg_idx++]);
    }
    // Output settings
    else if (arg == "--output")
    {
      udata->output = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--timing")
    {
      udata->timing = true;
    }
    else if (arg == "--fusedops")
    {
      udata->fusedops = true;
    }
    else if (arg == "--debug")
    {
      udata->debug = true;
    }
    // Help
    else if (arg == "--help")
    {
      if (outproc) InputHelp();
      return -1;
    }
    // Unknown input
    else
    {
      if (outproc)
      {
        cerr << "ERROR: Invalid input " << arg << endl;
        InputHelp();
      }
      return -1;
    }
  }

  // Recompute local number of nodes
  udata->nodes = udata->nodes_loc * udata->nprocs_w;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the solution error
static int SolutionError(N_Vector u_true, N_Vector u, N_Vector err,
                         UserData *udata)
{
  // Put true solution in error vector
  SetMus(udata);

  // Compute absolute error
  N_VLinearSum(ONE, u_true, -ONE, u, err);
  N_VAbs(err, err);

  return 0;
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --nodes                 : global number of values in vector" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --maa                   : size of Anderson Acceleration subspace" << endl;
  cout << "  --damping               : damping for Anderson Acceleration" << endl;
  cout << "  --orthaa                : orthogonalization routined used in Anderson Acceleration" << endl;
  cout << "  --maxits <iterations>   : max fixed point iterations" << endl;
  cout << "  --output                : output nonlinear solver statistics" << endl;
  cout << "  --timing                : print timing data" << endl;
  cout << "  --help                  : print this message and exit" << endl;
}

// Print user data
static int PrintUserData(UserData *udata)
{
  cout << endl;
  cout << "Expectation-Maximizaton Alg. for Mixture Densities Terms:" << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  nprocs             = " << udata->nprocs_w                << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  nodes              = " << udata->nodes                   << endl;
  cout << "  nodes_loc (proc 0) = " << udata->nodes_loc               << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  sigma              = {"
             << udata->sigma << ", " << udata->sigma << ", "
             << udata->sigma << "}" << endl;
  cout << "  alpha              = {"
             << udata->alpha1 << ", " << udata->alpha2 << ", "
             << udata->alpha3 << "}" << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  rtol               = " << udata->rtol                        << endl;
  cout << "  maa                = " << udata->maa                         << endl;
  cout << "  damping            = " << udata->damping                     << endl;
  cout << "  orthaa             = " << udata->orthaa                      << endl;
  cout << "  maxits             = " << udata->maxits                      << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  output             = " << udata->output                      << endl;
  cout << "  fusedops           = " << udata->fusedops                    << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << endl;

  return 0;
}

// Print nonlinear solver statistics
static int OutputStats(void *kinsol_mem, UserData* udata)
{
  int flag;

  // Get solver stats
  long int nfe, nni;
  flag = KINGetNumNonlinSolvIters(kinsol_mem, &nni);
  if (check_flag(&flag, "KINGetNumNonlinSolvIters", 1)) return(1);
  flag = KINGetNumFuncEvals(kinsol_mem, &nfe);
  if (check_flag(&flag, "KINGetNumFuncEvals", 1)) return(1);

  cout << setprecision(6);

  cout << "  Func evals       = " << nfe     << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << endl;

  return 0;
}

static int OutputTiming(UserData *udata)
{
  bool outproc = (udata->myid == 0);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(6);
  }

  double maxtime = 0.0;

  MPI_Reduce(&(udata->totaltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm);
  if (outproc)
  {
    cout << "  Total time                = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->fevaltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm);
  if (outproc)
  {
    cout << "  Function evaluation time  = " << maxtime << " sec" << endl;
  }

  return 0;
}

// Open residual and error output
static int OpenOutput(UserData *udata)
{
  bool outproc = (udata->myid == 0);

  if (outproc)
  {
    stringstream fname;

    // Open output stream for residual
    fname.str("");
    fname.clear();
    fname << "EM_res_m" << udata->maa << "_orth" << udata->orthaa
          << "_len" << udata->nodes_loc << ".txt";
    udata->rout.open(fname.str());

    udata->rout << scientific;
    udata->rout << setprecision(numeric_limits<realtype>::digits10);

    // Open output stream for error
    fname.str("");
    fname.clear();
    fname << "EM_err_m" << udata->maa << "_orth" << udata->orthaa
          << "_len" << udata->nodes_loc << ".txt";
    udata->eout.open(fname.str());

    udata->eout << scientific;
    udata->eout << setprecision(numeric_limits<realtype>::digits10);
  }

  return 0;
}

// Write residual and error out to file
static int WriteOutput(N_Vector u, N_Vector f, UserData *udata)
{
  int flag;
  bool outproc = (udata->myid == 0);

  // r = \|G(u) - u\|_inf
  N_VLinearSum(ONE, f, -ONE, u, udata->vtemp);
  realtype res = N_VMaxNorm(udata->vtemp);

  // e = \|u_exact - u\|_inf
  flag = SolutionError(udata->mu_true, u, udata->vtemp, udata);
  if (check_flag(&flag, "SolutionError", 1)) return 1;
  realtype err = N_VMaxNorm(udata->vtemp);

  if (outproc)
  {
    // Output residual
    udata->rout << res;
    udata->rout << endl;

    // Output error
    udata->eout << err;
    udata->eout << endl;
  }

  return 0;
}

// Close residual and error output files
static int CloseOutput(UserData *udata)
{
  bool outproc = (udata->myid == 0);

  if (outproc)
  {
    // Close output streams
    udata->rout.close();
    udata->eout.close();
  }

  return 0;
}

// Check function return value
static int check_flag(void *flagvalue, const string funcname, int opt)
{
  // Check if the function returned a NULL pointer
  if (opt == 0)
  {
    if (flagvalue == NULL)
    {
      cerr << endl << "ERROR: " << funcname << " returned NULL pointer" << endl
           << endl;
      return 1;
    }
  }
  // Check the function return flag value
  else if (opt == 1 || opt == 2)
  {
    int errflag = *((int *) flagvalue);
    if  ((opt == 1 && errflag < 0) || (opt == 2 && errflag != 0))
    {
      cerr << endl << "ERROR: " << funcname << " returned with flag = "
           << errflag << endl << endl;
      return 1;
    }
  }
  else
  {
    cerr << endl << "ERROR: check_flag called with an invalid option value"
         << endl;
    return 1;
  }

  return 0;
}

//---- end of file ----
