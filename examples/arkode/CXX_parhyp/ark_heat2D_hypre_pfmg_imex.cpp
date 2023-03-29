/* -----------------------------------------------------------------------------
 * Programmer(s): Rujeko Chinomona @ SMU
 *                David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
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
 * The following test simulates a simple 2D advection diffusion reaction problem,
 *
 *   u_t = kx u_xx + ky u_yy + ax ux + ay uy + r*(u - u*u*u),
 *
 * for t in [0, 1] and (x,y) in [0, 1]^2, with initial conditions
 *
 *   u(0,x,y) = sin^2(pi x) sin^2(pi y),
 *
 * stationary boundary conditions
 *
 *   u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
 *
 *
 * Under this setup, the problem does not have a simple analytical solution
 *
 *
 * The diffusion component spatial derivatives are computed using
 * second-order centered differences, with the data distributed over
 * nx * ny points on a uniform spatial grid.
 * The advection is discretized using a first order upwind scheme.
 * The problem is advanced in time with an IMEX method using
 * using an inexact Newton method paired with GMRES linear solver using
 * hypre's PFMG preconditioner based on diffusion only.
 * If the requested accuracy is 5th order, the problem is treated implicitly.
 * This option is included for computing a reference solution.
 * Several command line options are available to
 * change the problem parameters and ARKStep settings. Use the flag --help for
 * more information.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>

#include "arkode/arkode_arkstep.h"     // access to ARKStep
#include "nvector/nvector_parallel.h"  // access to the MPI N_Vector
#include "sunlinsol/sunlinsol_pcg.h"   // access to PCG SUNLinearSolver
#include "sunlinsol/sunlinsol_spgmr.h" // access to GMRES SUNLinearSolver
#include "sundials/sundials_math.h"    // def. math fcns eg. SUNMAX
#include "HYPRE_struct_ls.h"           // HYPRE structured grid solver interface
#include "mpi.h"                       // MPI header file


// Macros for problem constants
#define PI    RCONST(3.141592653589793238462643383279502884197169)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define EIGHT RCONST(8.0)

// Macro to access (x,y) location in 1D NVector array
#define IDX(x,y,n) ((n)*(y)+(x))


using namespace std;

// -----------------------------------------------------------------------------
// User data structure
// -----------------------------------------------------------------------------

struct UserData
{
  // Diffusion coefficients in the x and y directions
  realtype kx;
  realtype ky;

  // Advection coefficients in the x and y directions
  realtype ax;
  realtype ay;

  // Reactions coefficient
  realtype r;

  // Final time
  realtype tf;

  // Upper bounds in x and y directions
  realtype xu;
  realtype yu;

  // Global number of nodes in the x and y directions
  sunindextype nx;
  sunindextype ny;

  // Global total number of nodes
  sunindextype nodes;

  // Mesh spacing in the x and y directions
  realtype dx;
  realtype dy;

  // Local number of nodes in the x and y directions
  sunindextype nx_loc;
  sunindextype ny_loc;

  // Local total number of nodes
  sunindextype nodes_loc;

  // Global x and y indices of this subdomain
  sunindextype is;  // x starting index
  sunindextype ie;  // x ending index
  sunindextype js;  // y starting index
  sunindextype je;  // y ending index

  // MPI variables
  MPI_Comm comm_c; // Cartesian communicator in space

  int nprocs_w; // total number of MPI processes in Comm world
  int npx;      // number of MPI processes in the x-direction
  int npy;      // number of MPI processes in the y-direction

  int myid_c; // process ID in Cartesian communicator

  // Flags denoting if this process has a neighbor
  bool HaveNbrW;
  bool HaveNbrE;
  bool HaveNbrS;
  bool HaveNbrN;

  // Neighbor IDs for exchange
  int ipW;
  int ipE;
  int ipS;
  int ipN;

  // Receive buffers for neighbor exchange
  realtype *Wrecv;
  realtype *Erecv;
  realtype *Srecv;
  realtype *Nrecv;

  // Receive requests for neighbor exchange
  MPI_Request reqRW;
  MPI_Request reqRE;
  MPI_Request reqRS;
  MPI_Request reqRN;

  // Send buffers for neighbor exchange
  realtype *Wsend;
  realtype *Esend;
  realtype *Ssend;
  realtype *Nsend;

  // Send requests for neighor exchange
  MPI_Request reqSW;
  MPI_Request reqSE;
  MPI_Request reqSS;
  MPI_Request reqSN;

  // Integrator settings
  realtype rtol;    // relative tolerance
  realtype atol;    // absolute tolerance
  bool     linear;  // enable/disable linearly implicit option
  int      order;   // ARKode method order

  realtype hf;         // fixed step size
  int      controller; // step size adaptivity method
  int      maxsteps;   // max number of steps between outputs

  // Linear solver and preconditioner settings
  int      liniters; // number of linear iterations
  realtype epslin;   // linear solver tolerance factor
  int      prectype; // preconditioner type (NONE or LEFT)
  int      msbp;     // max number of steps between preconditioner setups

  // hypre objects
  HYPRE_StructGrid    grid;
  HYPRE_StructStencil stencil;
  HYPRE_StructMatrix  Amatrix;
  HYPRE_StructVector  bvec;
  HYPRE_StructVector  xvec;
  HYPRE_StructSolver  precond;

  // hypre grid extents
  HYPRE_Int ilower[2];
  HYPRE_Int iupper[2];

  // hypre workspace
  HYPRE_Int   nwork;
  HYPRE_Real *work;

  // hypre counters
  HYPRE_Int pfmg_its;

  // hypre PFMG settings (hypre defaults)
  HYPRE_Int pfmg_relax;  // type of relaxation:
                         //   0 - Jacobi
                         //   1 - Weighted Jacobi
                         //   2 - symmetric R/B Gauss-Seidel (*)
                         //   3 - nonsymmetric R/B Gauss-Seidel
  HYPRE_Int pfmg_nrelax; // number of pre and post relaxation sweeps (2)

  // Ouput variables
  int      output; // output level
  int      nout;   // number of output times
  ofstream uout;   // output file streams


  // Timing variables
  bool   timing;     // print timings
  double evolvetime;
  double rhstime;
  double explicitrhstime;
  double implicitrhstime;
  double matfilltime;
  double psetuptime;
  double psolvetime;
  double exchangetime;
};

// -----------------------------------------------------------------------------
// Setup the parallel decomposition
// -----------------------------------------------------------------------------

static int SetupDecomp(MPI_Comm comm_w, UserData *udata);

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// ODE right hand side functions
static int f(realtype t, N_Vector u, N_Vector f, void *user_data);
static int fi(realtype t, N_Vector u, N_Vector f, void *user_data);
static int fe(realtype t, N_Vector u, N_Vector f, void *user_data);
static int ffeval(realtype *fval, realtype uval, realtype rval);

// Preconditioner setup and solve functions
static int PSetup(realtype t, N_Vector u, N_Vector f, booleantype jok,
                  booleantype *jcurPtr, realtype gamma, void *user_data);

static int PSolve(realtype t, N_Vector u, N_Vector f, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,
                  void *user_data);

// -----------------------------------------------------------------------------
// Preconditioner and RHS helper functions
// -----------------------------------------------------------------------------

// Create hypre objects
static int HyprePFMG(UserData *udata);

// Fill A = I - gamma * J matrix
static int FillMatrix(UserData *udata, realtype gamma);

// Perform neighbor exchange
static int PostRecv(UserData *udata);
static int SendData(N_Vector y, UserData *udata);
static int WaitRecv(UserData *udata);

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Set the default values in the UserData structure
static int InitUserData(UserData *udata);

// Free memory allocated within UserData
static int FreeUserData(UserData *udata);

// Read the command line inputs and set UserData values
static int ReadInputs(int *argc, char ***argv, UserData *udata, bool outproc);

// Compute the initial condition
static int InitSolution(realtype t, N_Vector u, UserData *udata);

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Print the command line options
static void InputHelp();

// Print some UserData information
static int PrintUserData(UserData *udata);

// Output solution and error
static int OpenOutput(UserData *udata);
static int WriteOutput(realtype t, N_Vector u, UserData *udata);
static int CloseOutput(UserData *udata);

// Print integration statistics
static int OutputStats(void *arkode_mem, UserData *udata);

// Print integration timing
static int OutputTiming(UserData *udata);

// Check function return values
static int check_flag(void *flagvalue, const string funcname, int opt);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int flag;                      // reusable error-checking flag
  N_Vector u         = NULL;     // vector for storing solution
  SUNLinearSolver LS = NULL;     // linear solver memory structure
  void *arkode_mem   = NULL;     // ARKode memory structure

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

  // Add scope so objects object are destroyed before MPIFinalize
  {
    // Set output process flag
    bool outproc = (myid == 0);

    // Create SUNDIALS context
    sundials::Context sunctx(&comm_w);

    // ------------------------------------------
    // Setup UserData and parallel decomposition
    // ------------------------------------------

    // Allocate and initialize user data structure
    UserData udata;
    flag = InitUserData(&udata);
    if (check_flag(&flag, "InitUserData", 1)) return 1;

    // Parse command line inputs
    flag = ReadInputs(&argc, &argv, &udata, outproc);
    if (flag != 0) return 1;

    // Setup parallel decomposition
    flag = SetupDecomp(comm_w, &udata);
    if (check_flag(&flag, "SetupDecomp", 1)) return 1;

    // Output problem setup/options
    if (outproc)
    {
      flag = PrintUserData(&udata);
      if (check_flag(&flag, "PrintUserData", 1)) return 1;
    }

    // ------------------------
    // Create parallel vectors
    // ------------------------

    // Create vector for solution
    u = N_VNew_Parallel(udata.comm_c, udata.nodes_loc, udata.nodes, sunctx);
    if (check_flag((void *) u, "N_VNew_Parallel", 0)) return 1;

    // Set initial condition
    flag = InitSolution(ZERO, u, &udata);
    if (check_flag(&flag, "InitSolution", 1)) return 1;

    // -----------------------------------------
    // Create linear solvers and preconditioners
    // -----------------------------------------

    // Create linear solver
    LS = SUNLinSol_SPGMR(u, udata.prectype, udata.liniters, sunctx);
    if (check_flag((void *) LS, "SUNLinSol_SPGMR", 0)) return 1;

    // Create hypre objects
    if (udata.prectype != SUN_PREC_NONE)
    {
      flag = HyprePFMG(&udata);
      if (check_flag(&flag, "HyprePFMG", 1)) return 1;
    }

    // ----------------------------------------------
    // Setup ARKStep integrator and set options
    // ----------------------------------------------

    // Create integrator
    if (udata.order == 5)
    {
      arkode_mem = ARKStepCreate(NULL, f, ZERO, u, sunctx);
      if (check_flag((void *) arkode_mem, "ARKStepCreate", 0)) return 1;
    }
    else
    {
      arkode_mem = ARKStepCreate(fe, fi, ZERO, u, sunctx);
      if (check_flag((void *) arkode_mem, "ARKStepCreate", 0)) return 1;
    }

    // Specify tolerances
    flag = ARKStepSStolerances(arkode_mem, udata.rtol, udata.atol);
    if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

    // Attach linear solver
    flag = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
    if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;

    if (udata.prectype != SUN_PREC_NONE)
    {
      // Attach preconditioner
      flag = ARKStepSetPreconditioner(arkode_mem, PSetup, PSolve);
      if (check_flag(&flag, "ARKStepSetPreconditioner", 1)) return 1;

      // Set max steps between linear solver (preconditioner) setup calls
      flag = ARKStepSetLSetupFrequency(arkode_mem, udata.msbp);
      if (check_flag(&flag, "ARKStepSetLSetupFrequency", 1)) return 1;
    }

    // Set linear solver tolerance factor
    flag = ARKStepSetEpsLin(arkode_mem, udata.epslin);
    if (check_flag(&flag, "ARKStepSetEpsLin", 1)) return 1;

    // Use an ARKode provided table
    flag = ARKStepSetOrder(arkode_mem, udata.order);
    if (check_flag(&flag, "ARKStepSetOrder", 1)) return 1;

    // Set fixed step size or adaptivity method
    if (udata.hf > ZERO)
    {
      flag = ARKStepSetFixedStep(arkode_mem, udata.hf);
      if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;
    }
    else
    {
      flag = ARKStepSetAdaptivityMethod(arkode_mem, udata.controller, SUNTRUE,
                                        SUNFALSE, NULL);
      if (check_flag(&flag, "ARKStepSetAdaptivityMethod", 1)) return 1;
    }

    // Specify linearly implicit non-time-dependent RHS
    if (udata.linear)
    {
      flag = ARKStepSetLinear(arkode_mem, 0);
      if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;
    }

    // Attach user data
    flag = ARKStepSetUserData(arkode_mem, &udata);
    if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

    // Set max steps between outputs
    flag = ARKStepSetMaxNumSteps(arkode_mem, udata.maxsteps);
    if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;

    // Set stopping time
    flag = ARKStepSetStopTime(arkode_mem, udata.tf);
    if (check_flag(&flag, "ARKStepSetStopTime", 1)) return 1;

    // -----------------------
    // Loop over output times
    // -----------------------

    realtype t     = ZERO;
    realtype dTout = udata.tf / udata.nout;
    realtype tout  = dTout;

    // Inital output
    flag = OpenOutput(&udata);
    if (check_flag(&flag, "OpenOutput", 1)) return 1;

    flag = WriteOutput(t, u, &udata);
    if (check_flag(&flag, "WriteOutput", 1)) return 1;

    for (int iout = 0; iout < udata.nout; iout++)
    {
      // Start timer
      t1 = MPI_Wtime();

      // Evolve in time
      flag = ARKStepEvolve(arkode_mem, tout, u, &t, ARK_NORMAL);
      if (check_flag(&flag, "ARKStepEvolve", 1)) break;

      // Stop timer
      t2 = MPI_Wtime();

      // Update timer
      udata.evolvetime += t2 - t1;

      // Output solution and error
      flag = WriteOutput(t, u, &udata);
      if (check_flag(&flag, "WriteOutput", 1)) return 1;

      // Update output time
      tout += dTout;
      tout = (tout > udata.tf) ? udata.tf : tout;
    }

    // Close output
    flag = CloseOutput(&udata);
    if (check_flag(&flag, "CloseOutput", 1)) return 1;

    // --------------
    // Final outputs
    // --------------

    // Print final integrator stats
    if (udata.output > 0 && outproc)
    {
      cout << "Final fast integrator statistics:" << endl;
      flag = OutputStats(arkode_mem, &udata);
      if (check_flag(&flag, "OutputStats", 1)) return 1;
    }

    // Print timing
    if (udata.timing)
    {
      flag = OutputTiming(&udata);
      if (check_flag(&flag, "OutputTiming", 1)) return 1;
    }

    // --------------------
    // Clean up and return
    // --------------------

    ARKStepFree(&arkode_mem);       // Free integrator memory
    SUNLinSolFree(LS);              // Free linear solver
    N_VDestroy(u);                  // Free vectors
    FreeUserData(&udata);           // Free user data
  }

  // Finalize MPI
  MPI_Finalize();
  return 0;
}

// -----------------------------------------------------------------------------
// Setup the parallel decomposition
// -----------------------------------------------------------------------------

static int SetupDecomp(MPI_Comm comm_w, UserData *udata)
{
  int flag;

  // Check that this has not been called before
  if (udata->Erecv != NULL || udata->Wrecv != NULL ||
      udata->Srecv != NULL || udata->Nrecv != NULL)
  {
    cerr << "SetupDecomp error: parallel decomposition already set up" << endl;
    return -1;
  }

  // Get the number of processes
  flag = MPI_Comm_size(comm_w, &(udata->nprocs_w));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << flag << endl;
    return -1;
  }

  // Check the processor grid
  if ((udata->npx * udata->npy) != udata->nprocs_w)
  {
    cerr << "Error: npx * npy != nproc" << endl;
    return -1;
  }

  // Set up 2D Cartesian communicator
  int dims[2];
  dims[0] = udata->npx;
  dims[1] = udata->npy;

  int periods[2];
  periods[0] = 0;
  periods[1] = 0;

  flag = MPI_Cart_create(comm_w, 2, dims, periods, 0, &(udata->comm_c));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_create = " << flag << endl;
    return -1;
  }

  // Get my rank in the new Cartesian communicator
  flag = MPI_Comm_rank(udata->comm_c, &(udata->myid_c));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << flag << endl;
    return -1;
  }

  // Get dimension of the Cartesian communicator and my coordinates
  int coords[2];
  flag = MPI_Cart_get(udata->comm_c, 2, dims, periods, coords);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_get = " << flag << endl;
    return -1;
  }

  // Determine local extents in x-direction
  int idx         = coords[0];
  sunindextype qx = udata->nx / dims[0];
  sunindextype rx = udata->nx % dims[0];

  udata->is = qx * idx + (idx < rx ? idx : rx);
  udata->ie = udata->is + qx - 1 + (idx < rx ? 1 : 0);

  // Sanity check
  if (udata->ie > (udata->nx - 1))
  {
    cerr << "Error ie > nx - 1" << endl;
    return -1;
  }

  // Determine local extents in y-direction
  int idy         = coords[1];
  sunindextype qy = udata->ny / dims[1];
  sunindextype ry = udata->ny % dims[1];

  udata->js = qy * idy + (idy < ry ? idy : ry);
  udata->je = udata->js + qy - 1 + (idy < ry ? 1 : 0);

  // Sanity check
  if (udata->je > (udata->ny - 1))
  {
    cerr << "Error je > ny - 1" << endl;
    return -1;
  }

  // Number of local nodes
  udata->nx_loc = (udata->ie) - (udata->is) + 1;
  udata->ny_loc = (udata->je) - (udata->js) + 1;

  // Initialize global and local vector lengths
  udata->nodes     = udata->nx * udata->ny;
  udata->nodes_loc = udata->nx_loc * udata->ny_loc;

  // Determine if this proc has neighbors
  udata->HaveNbrW = (udata->is != 0);
  udata->HaveNbrE = (udata->ie != udata->nx-1);
  udata->HaveNbrS = (udata->js != 0);
  udata->HaveNbrN = (udata->je != udata->ny-1);

  // Allocate exchange buffers if necessary
  if (udata->HaveNbrW)
  {
    udata->Wrecv = new realtype[udata->ny_loc];
    udata->Wsend = new realtype[udata->ny_loc];
  }
  if (udata->HaveNbrE)
  {
    udata->Erecv = new realtype[udata->ny_loc];
    udata->Esend = new realtype[udata->ny_loc];
  }
  if (udata->HaveNbrS)
  {
    udata->Srecv = new realtype[udata->nx_loc];
    udata->Ssend = new realtype[udata->nx_loc];
  }
  if (udata->HaveNbrN)
  {
    udata->Nrecv = new realtype[udata->nx_loc];
    udata->Nsend = new realtype[udata->nx_loc];
  }

  // MPI neighborhood information
  int nbcoords[2];

  // West neighbor
  if (udata->HaveNbrW)
  {
    nbcoords[0] = coords[0]-1;
    nbcoords[1] = coords[1];
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // East neighbor
  if (udata->HaveNbrE)
  {
    nbcoords[0] = coords[0]+1;
    nbcoords[1] = coords[1];
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // South neighbor
  if (udata->HaveNbrS)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]-1;
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // North neighbor
  if (udata->HaveNbrN)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]+1;
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------
// f routine to compute the full ODE RHS function f(t,y).
static int f(realtype t, N_Vector u, N_Vector f, void *user_data)
{
  int          flag;
  sunindextype i, j;
  realtype ftemp;

  // Start timer
  double t1 = MPI_Wtime();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  realtype *farray = N_VGetArrayPointer(f);
  if (check_flag((void *) farray, "N_VGetArrayPointer", 0)) return -1;

  // Open exchange receives
  flag = PostRecv(udata);
  if (check_flag(&flag, "PostRecv", 1)) return -1;

  // Send exchange data
  flag = SendData(u, udata);
  if (check_flag(&flag, "SendData", 1)) return -1;

  // Constants for computing diffusion
  realtype cx = udata->kx / (udata->dx * udata->dx);
  realtype cy = udata->ky / (udata->dy * udata->dy);
  realtype cc = -TWO * (cx + cy);

  // Constants for computing advection
  realtype cxp = -SUNMAX(udata->ax, ZERO) / (udata->dx);
  realtype cxm = -SUNMIN(udata->ax, ZERO) / (udata->dx);
  realtype cyp = -SUNMAX(udata->ay, ZERO) / (udata->dy);
  realtype cym = -SUNMIN(udata->ay, ZERO) / (udata->dy);

  // Shortcut to reaction coefficient
  realtype rval = udata->r;

  // Initialize rhs to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  // Iterate over subdomain interior setting rhs values
  for (j = 1; j < ny_loc - 1; j++)
  {
    for (i = 1; i < nx_loc - 1; i++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)])
        + cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)])
        + cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)])
        + cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)])
        + cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // Wait for exchange receives
  flag = WaitRecv(udata);
  if (check_flag(&flag, "WaitRecv", 1)) return -1;

  // Iterate over subdomain boundaries (if not at overall domain boundary)
  realtype *Warray = udata->Wrecv;
  realtype *Earray = udata->Erecv;
  realtype *Sarray = udata->Srecv;
  realtype *Narray = udata->Nrecv;

  // West face
  if (udata->HaveNbrW)
  {
    i = 0;

    if (udata->HaveNbrS)  // South-West corner
    {
      j = 0;

      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)])
        + cxp * (uarray[IDX(i,j,nx_loc)] - Warray[j])
        + cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)])
        + cyp * (uarray[IDX(i,j,ny_loc)] - Sarray[i])
        + cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)])
        + cxp * (uarray[IDX(i,j,nx_loc)] - Warray[j])
        + cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)])
        + cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)])
        + cyp * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }

    if (udata->HaveNbrN)  // North-West corner
    {
      j = ny_loc - 1;

      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i])
        + cxp * (uarray[IDX(i,j,nx_loc)] - Warray[j])
        + cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)])
        + cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)])
        + cym * (Narray[i] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // East face
  if (udata->HaveNbrE)
  {
    i = nx_loc - 1;

    if (udata->HaveNbrS)  // South-East corner
    {
      j = 0;

      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)])
        + cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)])
        + cxm * (Earray[j] - uarray[IDX(i,j,nx_loc)])
        + cyp * (uarray[IDX(i,j,ny_loc)] - Sarray[i])
        + cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)])
        + cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)])
        + cxm * (Earray[j] - uarray[IDX(i,j,nx_loc)])
        + cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)])
        + cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }

    if (udata->HaveNbrN)  // North-East corner
    {
      j = ny_loc - 1;

      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i])
        + cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)])
        + cxm * (Earray[j] - uarray[IDX(i,j,nx_loc)])
        + cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)])
        + cym * (Narray[i] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // South face
  if (udata->HaveNbrS)
  {
    j = 0;

    for (i = 1; i < nx_loc - 1; i++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)])
        + cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)])
        + cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)])
        + cyp * (uarray[IDX(i,j,ny_loc)] - Sarray[i])
        + cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // North face
  if (udata->HaveNbrN)
  {
    j = udata->ny_loc - 1;

    for (i = 1; i < nx_loc - 1; i++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i])
        + cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)])
        + cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)])
        + cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)])
        + cym * (Narray[i] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->rhstime += t2 - t1;

  // Return success
  return 0;
}

// fi routine to compute the implicit part of ODE RHS function fi(t,u).
static int fi(realtype t, N_Vector u, N_Vector f, void *user_data)
{
  int          flag;
  sunindextype i, j;
  realtype ftemp;

  // Start timer
  double t1 = MPI_Wtime();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  realtype *farray = N_VGetArrayPointer(f);
  if (check_flag((void *) farray, "N_VGetArrayPointer", 0)) return -1;

  // Open exchange receives
  flag = PostRecv(udata);
  if (check_flag(&flag, "PostRecv", 1)) return -1;

  // Send exchange data
  flag = SendData(u, udata);
  if (check_flag(&flag, "SendData", 1)) return -1;

  // Constants for computing diffusion
  realtype cx = udata->kx / (udata->dx * udata->dx);
  realtype cy = udata->ky / (udata->dy * udata->dy);
  realtype cc = -TWO * (cx + cy);

  // Shortcut to reaction coefficient
  realtype rval = udata->r;

  // Initialize rhs to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  // Iterate over subdomain interior setting rhs values
  for (j = 1; j < ny_loc - 1; j++)
  {
    for (i = 1; i < nx_loc - 1; i++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }
  }

  // Wait for exchange receives
  flag = WaitRecv(udata);
  if (check_flag(&flag, "WaitRecv", 1)) return -1;

  // Iterate over subdomain boundaries (if not at overall domain boundary)
  realtype *Warray = udata->Wrecv;
  realtype *Earray = udata->Erecv;
  realtype *Sarray = udata->Srecv;
  realtype *Narray = udata->Nrecv;

  // West face
  if (udata->HaveNbrW)
  {
    i = 0;

    if (udata->HaveNbrS)  // South-West corner
    {
      j = 0;

      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }

    if (udata->HaveNbrN)  // North-West corner
    {
      j = ny_loc - 1;

      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // East face
  if (udata->HaveNbrE)
  {
    i = nx_loc - 1;

    if (udata->HaveNbrS)  // South-East corner
    {
      j = 0;

      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }

    if (udata->HaveNbrN)  // North-East corner
    {
      j = ny_loc - 1;

      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // South face
  if (udata->HaveNbrS)
  {
    j = 0;

    for (i = 1; i < nx_loc - 1; i++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }
  }

  // North face
  if (udata->HaveNbrN)
  {
    j = udata->ny_loc - 1;

    for (i = 1; i < nx_loc - 1; i++)
    {
      ffeval(&ftemp, uarray[IDX(i,j,nx_loc)], rval);
      farray[IDX(i,j,nx_loc)] = ftemp +
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->implicitrhstime += t2 - t1;

  // Return success
  return 0;
}

// fe routine to compute the explicit part of ODE RHS function.
static int fe(realtype t, N_Vector u, N_Vector f, void *user_data)
{
  int          flag;
  sunindextype i, j;

  // Start timer
  double t1 = MPI_Wtime();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  realtype *farray = N_VGetArrayPointer(f);
  if (check_flag((void *) farray, "N_VGetArrayPointer", 0)) return -1;

  // Open exchange receives
  flag = PostRecv(udata);
  if (check_flag(&flag, "PostRecv", 1)) return -1;

  // Send exchange data
  flag = SendData(u, udata);
  if (check_flag(&flag, "SendData", 1)) return -1;

  // Constants for computing advection
  realtype cxp = -SUNMAX(udata->ax, ZERO) / (udata->dx);
  realtype cxm = -SUNMIN(udata->ax, ZERO) / (udata->dx);
  realtype cyp = -SUNMAX(udata->ay, ZERO) / (udata->dy);
  realtype cym = -SUNMIN(udata->ay, ZERO) / (udata->dy);

  // Initialize rhs to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  // Iterate over subdomain interior setting rhs values
  for (j = 1; j < ny_loc - 1; j++)
  {
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] =
        cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)]) +
        cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)]) +
        cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)]) +
        cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // Wait for exchange receives
  flag = WaitRecv(udata);
  if (check_flag(&flag, "WaitRecv", 1)) return -1;

  // Iterate over subdomain boundaries (if not at overall domain boundary)
  realtype *Warray = udata->Wrecv;
  realtype *Earray = udata->Erecv;
  realtype *Sarray = udata->Srecv;
  realtype *Narray = udata->Nrecv;

  // West face
  if (udata->HaveNbrW)
  {
    i = 0;

    if (udata->HaveNbrS)  // South-West corner
    {
      j = 0;

      farray[IDX(i,j,nx_loc)] =
        cxp * (uarray[IDX(i,j,nx_loc)] - Warray[j]) +
        cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)]) +
        cyp * (uarray[IDX(i,j,ny_loc)] - Sarray[i]) +
        cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i,j,nx_loc)] =
        cxp * (uarray[IDX(i,j,nx_loc)] - Warray[j]) +
        cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)]) +
        cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)]) +
        cyp * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }

    if (udata->HaveNbrN)  // North-West corner
    {
      j = ny_loc - 1;

      farray[IDX(i,j,nx_loc)] =
        cxp * (uarray[IDX(i,j,nx_loc)] - Warray[j]) +
        cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)]) +
        cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)]) +
        cym * (Narray[i] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // East face
  if (udata->HaveNbrE)
  {
    i = nx_loc - 1;

    if (udata->HaveNbrS)  // South-East corner
    {
      j = 0;

      farray[IDX(i,j,nx_loc)] =
        cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)]) +
        cxm * (Earray[j] - uarray[IDX(i,j,nx_loc)]) +
        cyp * (uarray[IDX(i,j,ny_loc)] - Sarray[i]) +
        cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i,j,nx_loc)] =
        cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)]) +
        cxm * (Earray[j] - uarray[IDX(i,j,nx_loc)]) +
        cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)]) +
        cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }

    if (udata->HaveNbrN)  // North-East corner
    {
      j = ny_loc - 1;

      farray[IDX(i,j,nx_loc)] =
        cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)]) +
        cxm * (Earray[j] - uarray[IDX(i,j,nx_loc)]) +
        cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)]) +
        cym * (Narray[i] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // South face
  if (udata->HaveNbrS)
  {
    j = 0;

    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] =
        cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)]) +
        cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)]) +
        cyp * (uarray[IDX(i,j,ny_loc)] - Sarray[i]) +
        cym * (uarray[IDX(i,j+1,ny_loc)] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // North face
  if (udata->HaveNbrN)
  {
    j = udata->ny_loc - 1;

    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] =
        cxp * (uarray[IDX(i,j,nx_loc)] - uarray[IDX(i-1,j,nx_loc)]) +
        cxm * (uarray[IDX(i+1,j,nx_loc)] - uarray[IDX(i,j,nx_loc)]) +
        cyp * (uarray[IDX(i,j,ny_loc)] - uarray[IDX(i,j-1,ny_loc)]) +
        cym * (Narray[i] - uarray[IDX(i,j,ny_loc)]);
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->explicitrhstime += t2 - t1;

  // Return success
  return 0;
}

// ff helper function
static int ffeval(realtype *fval, realtype uval, realtype rval) {

  *fval = rval * (uval - uval*uval*uval);
  return 0;
}

// Preconditioner setup routine
static int PSetup(realtype t, N_Vector u, N_Vector f, booleantype jok,
                  booleantype *jcurPtr, realtype gamma, void *user_data)
{
  int flag;

  // Start timer
  double t1 = MPI_Wtime();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Fill matrix A = I - gamma * J
  flag = FillMatrix(udata, gamma);
  if (flag != 0) return -1;

  // Indicate that the jacobian is current
  *jcurPtr = SUNTRUE;

  // -----------
  // Setup PFMG
  // -----------

  // Set rhs/solution vectors as all zero for now
  flag = HYPRE_StructVectorSetConstantValues(udata->bvec, ZERO);
  if (flag != 0) return -1;

  flag = HYPRE_StructVectorAssemble(udata->bvec);
  if (flag != 0) return -1;

  flag = HYPRE_StructVectorSetConstantValues(udata->xvec, ZERO);
  if (flag != 0) return -1;

  flag = HYPRE_StructVectorAssemble(udata->xvec);
  if (flag != 0) return -1;

  // Set up the solver
  flag = HYPRE_StructPFMGSetup(udata->precond, udata->Amatrix,
                               udata->bvec, udata->xvec);
  if (flag != 0) return -1;

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->psetuptime += t2 - t1;

  // Return success
  return 0;
}

// Preconditioner solve routine for Pz = r
static int PSolve(realtype t, N_Vector u, N_Vector f, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,
                  void *user_data)
{
  int flag;

  // Start timer
  double t1 = MPI_Wtime();

  // Access user_data structure
  UserData *udata = (UserData *) user_data;

  // Insert rhs N_Vector entries into HYPRE vector b and assemble
  flag = HYPRE_StructVectorSetBoxValues(udata->bvec,
                                        udata->ilower, udata->iupper,
                                        N_VGetArrayPointer(r));
  if (flag != 0) return -1;

  flag = HYPRE_StructVectorAssemble(udata->bvec);
  if (flag != 0) return -1;

  // Set the initial guess into HYPRE vector x and assemble
  flag = HYPRE_StructVectorSetConstantValues(udata->xvec, ZERO);
  if (flag != 0) return -1;

  flag = HYPRE_StructVectorAssemble(udata->xvec);
  if (flag != 0) return -1;

  // Solve the linear system
  flag = HYPRE_StructPFMGSolve(udata->precond, udata->Amatrix,
                               udata->bvec, udata->xvec);

  // If a convergence error occured, clear the error and continue. For any
  // other error return with a recoverable error.
  if (flag == HYPRE_ERROR_CONV) HYPRE_ClearError(HYPRE_ERROR_CONV);
  else if (flag != 0) return 1;

  // Update precond statistics
  HYPRE_Int itmp;
  flag = HYPRE_StructPFMGGetNumIterations(udata->precond, &itmp);
  if (flag != 0) return -1;

  udata->pfmg_its += itmp;

  // Extract solution values
  flag = HYPRE_StructVectorGetBoxValues(udata->xvec,
                                        udata->ilower, udata->iupper,
                                        N_VGetArrayPointer(z));
  if (flag != 0) return -1;

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->psolvetime += t2 - t1;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Preconditioner and RHS helper functions
// -----------------------------------------------------------------------------

// Create hypre objects
static int HyprePFMG(UserData *udata)
{
  int flag, result;

  // Check input
  if (udata == NULL) return -1;

  // Check if the grid or stencil have been created
  if ((udata->grid != NULL || udata->stencil != NULL))
  {
    cerr << "CreateHyprePFMG error: grid or stencil already exists" << endl;
    return -1;
  }

  // Check for valid 2D Cartesian MPI communicator
  flag = MPI_Topo_test(udata->comm_c, &result);
  if ((flag != MPI_SUCCESS) || (result != MPI_CART))
  {
    cerr << "CreateHyprePFMG error: communicator is not Cartesian" << endl;
    return -1;
  }

  flag = MPI_Cartdim_get(udata->comm_c, &result);
  if ((flag != MPI_SUCCESS) || (result != 2))
  {
    cerr << "CreateHyprePFMG error: communicator is not 2D" << endl;
    return -1;
  }

  // -----
  // Grid
  // -----

  // Create 2D grid object
  flag = HYPRE_StructGridCreate(udata->comm_c, 2, &(udata->grid));
  if (flag != 0) { FreeUserData(udata); return -1; }

  // Set grid extents (lower left and upper right corners)
  udata->ilower[0] = udata->is;
  udata->ilower[1] = udata->js;

  udata->iupper[0] = udata->ie;
  udata->iupper[1] = udata->je;

  flag = HYPRE_StructGridSetExtents(udata->grid, udata->ilower, udata->iupper);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // Assemble the grid
  flag = HYPRE_StructGridAssemble(udata->grid);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // --------
  // Stencil
  // --------

  // Create the 2D 5 point stencil object
  flag = HYPRE_StructStencilCreate(2, 5, &(udata->stencil));
  if (flag != 0) { FreeUserData(udata); return -1; }

  // Set the stencil entries (center, left, right, bottom, top)
  HYPRE_Int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

  for (int entry = 0; entry < 5; entry++)
  {
    flag = HYPRE_StructStencilSetElement(udata->stencil, entry, offsets[entry]);
    if (flag != 0) { FreeUserData(udata); return -1; }
  }

  // -----------
  // Work array
  // -----------

  udata->nwork = 5 * udata->nodes_loc;
  udata->work  = NULL;
  udata->work  = new HYPRE_Real[udata->nwork];
  if (udata->work == NULL) { FreeUserData(udata); return -1; }

  // ---------
  // x vector
  // ---------

  flag = HYPRE_StructVectorCreate(udata->comm_c, udata->grid, &(udata->xvec));
  if (flag != 0) { FreeUserData(udata); return -1; }

  flag = HYPRE_StructVectorInitialize(udata->xvec);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // ---------
  // b vector
  // ---------

  flag = HYPRE_StructVectorCreate(udata->comm_c, udata->grid, &(udata->bvec));
  if (flag != 0) { FreeUserData(udata); return -1; }

  flag = HYPRE_StructVectorInitialize(udata->bvec);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // ---------
  // A matrix
  // ---------

  flag = HYPRE_StructMatrixCreate(udata->comm_c, udata->grid, udata->stencil,
                                  &(udata->Amatrix));
  if (flag != 0) { FreeUserData(udata); return -1; }

  flag = HYPRE_StructMatrixInitialize(udata->Amatrix);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // --------------------
  // PFMG preconditioner
  // --------------------

  flag = HYPRE_StructPFMGCreate(udata->comm_c, &(udata->precond));
  if (flag != 0) { FreeUserData(udata); return -1; }

  // signal that the inital guess is zero
  flag = HYPRE_StructPFMGSetZeroGuess(udata->precond);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // tol <= 0.0 means do the max number of iterations
  flag = HYPRE_StructPFMGSetTol(udata->precond, ZERO);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // use one v-cycle
  flag = HYPRE_StructPFMGSetMaxIter(udata->precond, 1);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // use non-Galerkin corase grid operator
  flag = HYPRE_StructPFMGSetRAPType(udata->precond, 1);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // set the relaxation type
  flag = HYPRE_StructPFMGSetRelaxType(udata->precond, udata->pfmg_relax);
  if (flag != 0) { FreeUserData(udata); return -1; }

  // set the number of pre and post relaxation sweeps
  flag = HYPRE_StructPFMGSetNumPreRelax(udata->precond, udata->pfmg_nrelax);
  if (flag != 0) { FreeUserData(udata); return -1; }

  flag = HYPRE_StructPFMGSetNumPostRelax(udata->precond, udata->pfmg_nrelax);
  if (flag != 0) { FreeUserData(udata); return -1; }

  return 0;
}

// Fill A = I - gamma * J matrix
static int FillMatrix(UserData *udata, realtype gamma)
{
  // Variable shortcuts
  HYPRE_Int ilower[2];
  HYPRE_Int iupper[2];

  ilower[0] = udata->ilower[0];
  ilower[1] = udata->ilower[1];

  iupper[0] = udata->iupper[0];
  iupper[1] = udata->iupper[1];

  HYPRE_Int   nwork = udata->nwork;
  HYPRE_Real *work  = udata->work;

  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Matrix stencil
  HYPRE_Int entries[5] = {0, 1, 2, 3, 4};
  HYPRE_Int entry[1];

  // Grid extents for setting boundary entries
  HYPRE_Int bc_ilower[2];
  HYPRE_Int bc_iupper[2];

  // Loop counters
  HYPRE_Int idx, ix, iy;

  // hypre return flag
  int flag;

  // --------------------------
  // Compute A = I - gamma * J
  // --------------------------

  // Start timer
  double t1 = MPI_Wtime();

  // Only do work if the box is non-zero in size
  if ((ilower[0] <= iupper[0]) &&
      (ilower[1] <= iupper[1]))
  {
    // Jacobian values
    realtype c1 = udata->kx / (udata->dx * udata->dx);
    realtype c2 = udata->ky / (udata->dy * udata->dy);
    realtype c3 = -TWO * (c1 + c2);

    // Linear system values
    realtype cx = -gamma * c1;
    realtype cy = -gamma * c2;
    realtype cc = ONE - gamma * c3;

    // --------------------------------
    // Set matrix values for all nodes
    // --------------------------------

    // Set the matrix interior entries (center, left, right, bottom, top)
    idx = 0;
    for (iy = 0; iy < ny_loc; iy++)
    {
      for (ix = 0; ix < nx_loc; ix++)
      {
        work[idx]     = cc;
        work[idx + 1] = cx;
        work[idx + 2] = cx;
        work[idx + 3] = cy;
        work[idx + 4] = cy;
        idx += 5;
      }
    }

    // Modify the matrix
    flag = HYPRE_StructMatrixSetBoxValues(udata->Amatrix,
                                          ilower, iupper,
                                          5, entries, work);
    if (flag != 0) return -1;

    // ----------------------------------------
    // Correct matrix values at boundary nodes
    // ----------------------------------------

    // Set the matrix boundary entries (center, left, right, bottom, top)
    if (ilower[1] == 0 ||
        iupper[1] == (udata->ny - 1) ||
        ilower[0] == 0 ||
        iupper[0] == (udata->nx - 1))
    {
      idx = 0;
      for (iy = 0; iy < ny_loc; iy++)
      {
        for (ix = 0; ix < nx_loc; ix++)
        {
          work[idx]     = ONE;
          work[idx + 1] = ZERO;
          work[idx + 2] = ZERO;
          work[idx + 3] = ZERO;
          work[idx + 4] = ZERO;
          idx += 5;
        }
      }
    }

    // Set cells on western boundary
    if (ilower[0] == 0)
    {
      // Grid cell on south-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];

      // Grid cell on north-west corner
      bc_iupper[0] = ilower[0];
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(udata->Amatrix,
                                              bc_ilower, bc_iupper,
                                              5, entries, work);
        if (flag != 0) return -1;
      }
    }

    // Set cells on eastern boundary
    if (iupper[0] == (udata->nx - 1))
    {
      // Grid cell on south-east corner
      bc_ilower[0] = iupper[0];
      bc_ilower[1] = ilower[1];

      // Grid cell on north-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(udata->Amatrix,
                                              bc_ilower, bc_iupper,
                                              5, entries, work);
        if (flag != 0) return -1;
      }
    }

    // Correct cells on southern boundary
    if (ilower[1] == 0)
    {
      // Grid cell on south-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];

      // Grid cell on south-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = ilower[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(udata->Amatrix,
                                              bc_ilower, bc_iupper,
                                              5, entries, work);
        if (flag != 0) return -1;
      }
    }

    // Set cells on northern boundary
    if (iupper[1] == (udata->ny - 1))
    {
      // Grid cell on north-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = iupper[1];

      // Grid cell on north-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(udata->Amatrix,
                                              bc_ilower, bc_iupper,
                                              5, entries, work);
        if (flag != 0) return -1;
      }
    }

    // -----------------------------------------------------------
    // Remove connections between the interior and boundary nodes
    // -----------------------------------------------------------

    // Zero out work array
    for (ix = 0; ix < nwork; ix++)
    {
      work[ix] = ZERO;
    }

    // Second column of nodes (depends on western boundary)
    if ((ilower[0] <= 1) && (iupper[0] >= 1))
    {
      // Remove western dependency
      entry[0] = 1;

      // Grid cell on south-west corner
      bc_ilower[0] = 1;
      bc_ilower[1] = ilower[1];

      // Grid cell on north-west corner
      bc_iupper[0] = 1;
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(udata->Amatrix,
                                              bc_ilower, bc_iupper,
                                              1, entry, work);
        if (flag != 0) return -1;
      }
    }

    // Next to last column (depends on eastern boundary)
    if ((ilower[0] <= (udata->nx - 2)) &&
        (iupper[0] >= (udata->nx - 2)))
    {
      // Remove eastern dependency
      entry[0] = 2;

      // Grid cell on south-east corner
      bc_ilower[0] = udata->nx - 2;
      bc_ilower[1] = ilower[1];

      // Grid cell on north-east corner
      bc_iupper[0] = udata->nx - 2;
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(udata->Amatrix,
                                              bc_ilower, bc_iupper,
                                              1, entry, work);
        if (flag != 0) return -1;
      }
    }

    // Second row of nodes (depends on southern boundary)
    if ((udata->ilower[1] <= 1) && (udata->iupper[1] >= 1))
    {
      // Remove southern dependency
      entry[0] = 3;

      // Grid cell on south-west corner
      bc_ilower[0] = udata->ilower[0];
      bc_ilower[1] = 1;

      // Grid cell on south-east corner
      bc_iupper[0] = udata->iupper[0];
      bc_iupper[1] = 1;

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(udata->Amatrix,
                                              bc_ilower, bc_iupper,
                                              1, entry, work);
        if (flag != 0) return -1;
      }
    }

    // Next to last row of nodes (depends on northern boundary)
    if ((ilower[1] <= (udata->ny - 2)) &&
        (iupper[1] >= (udata->ny - 2)))
    {
      // Remove northern dependency
      entry[0] = 4;

      // Grid cell on north-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = udata->ny - 2;

      // Grid cell on north-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = udata->ny - 2;

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(udata->Amatrix,
                                              bc_ilower, bc_iupper,
                                              1, entry, work);
        if (flag != 0) return -1;
      }
    }
  }

  // Assemble matrix
  flag = HYPRE_StructMatrixAssemble(udata->Amatrix);
  if (flag != 0) return -1;

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->matfilltime += t2 - t1;

  // Return success
  return 0;
}

// Post exchange receives
static int PostRecv(UserData *udata)
{
  int flag;

  // Start timer
  double t1 = MPI_Wtime();

  // Open Irecv buffers
  if (udata->HaveNbrW)
  {
    flag = MPI_Irecv(udata->Wrecv, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipW, MPI_ANY_TAG, udata->comm_c, &(udata->reqRW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    flag = MPI_Irecv(udata->Erecv, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipE, MPI_ANY_TAG, udata->comm_c, &(udata->reqRE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    flag = MPI_Irecv(udata->Srecv, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipS, MPI_ANY_TAG, udata->comm_c, &(udata->reqRS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    flag = MPI_Irecv(udata->Nrecv, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipN, MPI_ANY_TAG, udata->comm_c, &(udata->reqRN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// Send exchange data
static int SendData(N_Vector y, UserData *udata)
{
  int flag, i;
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Start timer
  double t1 = MPI_Wtime();

  // Access data array
  realtype *Y = N_VGetArrayPointer(y);
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return -1;

  // Send data
  if (udata->HaveNbrW)
  {
    for (i = 0; i < ny_loc; i++) udata->Wsend[i] = Y[IDX(0,i,nx_loc)];
    flag = MPI_Isend(udata->Wsend, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipW, 0, udata->comm_c, &(udata->reqSW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    for (i = 0; i < ny_loc; i++) udata->Esend[i] = Y[IDX(nx_loc-1,i,nx_loc)];
    flag = MPI_Isend(udata->Esend, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipE, 1, udata->comm_c, &(udata->reqSE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    for (i = 0; i < nx_loc; i++) udata->Ssend[i] = Y[IDX(i,0,nx_loc)];
    flag = MPI_Isend(udata->Ssend, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipS, 2, udata->comm_c, &(udata->reqSS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    for (i = 0; i < nx_loc; i++) udata->Nsend[i] = Y[IDX(i,ny_loc-1,nx_loc)];
    flag = MPI_Isend(udata->Nsend, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipN, 3, udata->comm_c, &(udata->reqSN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// Wait for exchange data
static int WaitRecv(UserData *udata)
{
  // Local variables
  int flag;
  MPI_Status stat;

  // Start timer
  double t1 = MPI_Wtime();

  // Wait for messages to finish
  if (udata->HaveNbrW)
  {
    flag = MPI_Wait(&(udata->reqRW), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSW), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    flag = MPI_Wait(&(udata->reqRE), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSE), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    flag = MPI_Wait(&(udata->reqRS), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSS), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    flag = MPI_Wait(&(udata->reqRN), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSN), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Initialize memory allocated within Userdata
static int InitUserData(UserData *udata)
{
  // Diffusion coefficient
  udata->kx = ONE;
  udata->ky = ONE;

  // Advection coefficient
  udata->ax = ONE;
  udata->ay = ONE;

  // Reaction coefficient
  udata->r = ONE;

  // Final time
  udata->tf = ONE;

  // Upper bounds in x and y directions
  udata->xu = ONE;
  udata->yu = ONE;

  // Global number of nodes in the x and y directions
  udata->nx    = 64;
  udata->ny    = 64;
  udata->nodes = udata->nx * udata->ny;

  // Mesh spacing in the x and y directions
  udata->dx = udata->xu / (udata->nx - 1);
  udata->dy = udata->yu / (udata->ny - 1);

  // Locals number of nodes in the x and y directions (set in SetupDecomp)
  udata->nx_loc    = 0;
  udata->ny_loc    = 0;
  udata->nodes_loc = 0;

  // Global indices of this subdomain (set in SetupDecomp)
  udata->is = 0;
  udata->ie = 0;
  udata->js = 0;
  udata->je = 0;

  // MPI variables (set in SetupDecomp)
  udata->comm_c = MPI_COMM_NULL;

  udata->nprocs_w = 1;
  udata->npx      = 1;
  udata->npy      = 1;

  udata->myid_c = 0;

  // Flags denoting neighbors (set in SetupDecomp)
  udata->HaveNbrW = true;
  udata->HaveNbrE = true;
  udata->HaveNbrS = true;
  udata->HaveNbrN = true;

  // Exchange receive buffers (allocated in SetupDecomp)
  udata->Erecv = NULL;
  udata->Wrecv = NULL;
  udata->Nrecv = NULL;
  udata->Srecv = NULL;

  // Exchange send buffers (allocated in SetupDecomp)
  udata->Esend = NULL;
  udata->Wsend = NULL;
  udata->Nsend = NULL;
  udata->Ssend = NULL;

  // Neighbors IDs (set in SetupDecomp)
  udata->ipW = -1;
  udata->ipE = -1;
  udata->ipS = -1;
  udata->ipN = -1;

  // Integrator settings
  udata->rtol     = RCONST(1.e-5);   // relative tolerance
  udata->atol     = RCONST(1.e-10);  // absolute tolerance
  udata->linear   = false;           // linearly implicit problem
  udata->order    = 3;               // method order

  udata->hf         = ZERO;           // using adaptive step sizes at fast
  udata->controller = 0;              // PID controller
  udata->maxsteps   = 0;              // use ARKode default

  // Linear solver and preconditioner options
  udata->liniters = 10;         // max linear iterations
  udata->epslin   = -ONE;       // use ARKode default (0.05)
  udata->prectype = SUN_PREC_LEFT;  // no preconditioning
  udata->msbp     = 0;          // use ARKode default (20 steps)

  // hypre objects
  udata->grid    = NULL;
  udata->stencil = NULL;
  udata->Amatrix = NULL;
  udata->bvec    = NULL;
  udata->xvec    = NULL;
  udata->precond = NULL;

  // hypre grid extents
  udata->ilower[0] = 0;
  udata->ilower[1] = 0;

  udata->iupper[0] = 0;
  udata->iupper[1] = 0;

  // hypre workspace
  udata->nwork = 0;
  udata->work  = NULL;

  // hypre counters
  udata->pfmg_its = 0;

  // hypre PFMG settings
  udata->pfmg_relax  = 2;
  udata->pfmg_nrelax = 2;

  // Output variables
  udata->output = 1;   // 0 = no output, 1 = stats output, 2 = output to disk
  udata->nout   = 10;  // Number of output times

  // Timing variables
  udata->timing          = false;
  udata->evolvetime      = 0.0;
  udata->rhstime         = 0.0;
  udata->explicitrhstime = 0.0;
  udata->implicitrhstime = 0.0;
  udata->matfilltime     = 0.0;
  udata->psetuptime      = 0.0;
  udata->psolvetime      = 0.0;
  udata->exchangetime    = 0.0;

  // Return success
  return 0;
}

// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{
  // Free exchange buffers
  if (udata->Wrecv != NULL)  delete[] udata->Wrecv;
  if (udata->Wsend != NULL)  delete[] udata->Wsend;
  if (udata->Erecv != NULL)  delete[] udata->Erecv;
  if (udata->Esend != NULL)  delete[] udata->Esend;
  if (udata->Srecv != NULL)  delete[] udata->Srecv;
  if (udata->Ssend != NULL)  delete[] udata->Ssend;
  if (udata->Nrecv != NULL)  delete[] udata->Nrecv;
  if (udata->Nsend != NULL)  delete[] udata->Nsend;

  // Free hypre preconditioner data
  if (udata->grid    != NULL) HYPRE_StructGridDestroy(udata->grid);
  if (udata->stencil != NULL) HYPRE_StructStencilDestroy(udata->stencil);
  if (udata->Amatrix != NULL) HYPRE_StructMatrixDestroy(udata->Amatrix);
  if (udata->bvec    != NULL) HYPRE_StructVectorDestroy(udata->bvec);
  if (udata->xvec    != NULL) HYPRE_StructVectorDestroy(udata->xvec);
  if (udata->precond != NULL) HYPRE_StructPFMGDestroy(udata->precond);
  if (udata->work    != NULL) delete[] udata->work;

  // Free MPI Cartesian communicator
  if (udata->comm_c != MPI_COMM_NULL)
    MPI_Comm_free(&(udata->comm_c));

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
    if (arg == "--mesh")
    {
      udata->nx = stoi((*argv)[arg_idx++]);
      udata->ny = stoi((*argv)[arg_idx++]);
    }
    // MPI processes
    else if (arg == "--np")
    {
      udata->npx = stoi((*argv)[arg_idx++]);
      udata->npy = stoi((*argv)[arg_idx++]);
    }
    // Domain upper bounds
    else if (arg == "--domain")
    {
      udata->xu = stoi((*argv)[arg_idx++]);
      udata->yu = stoi((*argv)[arg_idx++]);
    }
    // Diffusion parameters
    else if (arg == "--k")
    {
      udata->kx = stod((*argv)[arg_idx++]);
      udata->ky = stod((*argv)[arg_idx++]);
    }
    // Advection parameters
    else if (arg == "--a")
    {
      udata->ax = stod((*argv)[arg_idx++]);
      udata->ay = stod((*argv)[arg_idx++]);
    }
    // Reaction parameters
    else if (arg == "--r")
    {
      udata->r = stod((*argv)[arg_idx++]);
    }
    // Temporal domain settings
    else if (arg == "--tf")
    {
      udata->tf = stod((*argv)[arg_idx++]);
    }
    // Integrator settings
    else if (arg == "--rtol")
    {
      udata->rtol = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--atol")
    {
      udata->atol = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--linear")
    {
      udata->linear = true;
    }
    else if (arg == "--order")
    {
      udata->order = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--hf")
    {
      udata->hf = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--controller")
    {
      udata->controller = stoi((*argv)[arg_idx++]);
    }
    // Linear solver settings
    else if (arg == "--liniters")
    {
      udata->liniters = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--epslin")
    {
      udata->epslin = stod((*argv)[arg_idx++]);
    }
    // Preconditioner settings
    else if (arg == "--noprec")
    {
      udata->prectype = SUN_PREC_NONE;
    }
    else if (arg == "--msbp")
    {
      udata->msbp = stoi((*argv)[arg_idx++]);
    }
    // PFMG settings
    else if (arg == "--pfmg_relax")
    {
      udata->pfmg_relax = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--pfmg_nrelax")
    {
      udata->pfmg_nrelax = stoi((*argv)[arg_idx++]);
    }
    // Output settings
    else if (arg == "--output")
    {
      udata->output = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--nout")
    {
      udata->nout = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--maxsteps")
    {
      udata->maxsteps = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--timing")
    {
      udata->timing = true;
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

  // Recompute total number of nodes
  udata->nodes = udata->nx * udata->ny;

  // Recompute x and y mesh spacing
  udata->dx = (udata->xu) / (udata->nx - 1);
  udata->dy = (udata->yu) / (udata->ny - 1);

  // If the method order is 1 fixed time stepping must be used
  if (udata->order == 1 && !(udata->hf > ZERO))
  {
    cerr << "ERROR: Method order 1 requires fixed time stepping" << endl;
    return -1;
  }

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the exact solution
static int InitSolution(realtype t, N_Vector u, UserData *udata)
{
  realtype x, y;
  realtype cos_sqr_t;
  realtype sin_sqr_x, sin_sqr_y;

  // Constants for computing solution
  cos_sqr_t = cos(PI * t) * cos(PI * t);

  // Initialize u to zero (handles boundary conditions)
  N_VConst(ZERO, u);

  // Iterative over domain interior
  sunindextype istart = (udata->HaveNbrW) ? 0 : 1;
  sunindextype iend   = (udata->HaveNbrE) ? udata->nx_loc : udata->nx_loc - 1;

  sunindextype jstart = (udata->HaveNbrS) ? 0 : 1;
  sunindextype jend   = (udata->HaveNbrN) ? udata->ny_loc : udata->ny_loc - 1;

  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  for (sunindextype j = jstart; j < jend; j++)
  {
    for (sunindextype i = istart; i < iend; i++)
    {
      x  = (udata->is + i) * udata->dx;
      y  = (udata->js + j) * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      uarray[IDX(i,j,udata->nx_loc)] = sin_sqr_x * sin_sqr_y * cos_sqr_t;
    }
  }

  return 0;
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --mesh <nx> <ny>        : mesh points in the x and y directions" << endl;
  cout << "  --np <npx> <npy>        : number of MPI processes in the x and y directions" << endl;
  cout << "  --domain <xu> <yu>      : domain upper bound in the x and y direction" << endl;
  cout << "  --k <kx> <ky>           : diffusion coefficients" << endl;
  cout << "  --a <ax> <ay>           : advection coefficients" << endl;
  cout << "  --r <r>                 : reaction coefficient" << endl;
  cout << "  --tf <time>             : final time" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --atol <atol>           : absoltue tolerance" << endl;
  cout << "  --linear                : enable linearly implicit flag" << endl;
  cout << "  --order <ord>           : method order" << endl;
  cout << "  --hf <step>             : used fixed step size" << endl;
  cout << "  --controller <ctr>      : time step adaptivity controller" << endl;
  cout << "  --liniters <iters>      : max number of iterations" << endl;
  cout << "  --epslin <factor>       : linear tolerance factor" << endl;
  cout << "  --noprec                : disable PFMG preconditioner" << endl;
  cout << "  --msbp <steps>          : max steps between prec setups" << endl;
  cout << "  --pfmg_relax <types>    : relaxtion type in PFMG" << endl;
  cout << "  --pfmg_nrelax <iters>   : pre/post relaxtion sweeps in PFMG" << endl;
  cout << "  --output <level>        : output level" << endl;
  cout << "  --nout <nout>           : number of outputs" << endl;
  cout << "  --maxsteps <steps>      : max steps between outputs" << endl;
  cout << "  --timing                : print timing data" << endl;
  cout << "  --help                  : print this message and exit" << endl;
}

// Print user data
static int PrintUserData(UserData *udata)
{
  cout << endl;
  cout << "2D Heat PDE test problem:"                     << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  nprocs         = " << udata->nprocs_w        << endl;
  cout << "  npx            = " << udata->npx             << endl;
  cout << "  npy            = " << udata->npy             << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  kx             = " << udata->kx              << endl;
  cout << "  ky             = " << udata->ky              << endl;
  cout << "  ax             = " << udata->ax              << endl;
  cout << "  ay             = " << udata->ay              << endl;
  cout << "  r              = " << udata->r               << endl;
  cout << "  tf             = " << udata->tf              << endl;
  cout << "  xu             = " << udata->xu              << endl;
  cout << "  yu             = " << udata->yu              << endl;
  cout << "  nx             = " << udata->nx              << endl;
  cout << "  ny             = " << udata->ny              << endl;
  cout << "  nxl (proc 0)   = " << udata->nx_loc          << endl;
  cout << "  nyl (proc 0)   = " << udata->ny_loc          << endl;
  cout << "  dx             = " << udata->dx              << endl;
  cout << "  dy             = " << udata->dy              << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  rtol           = " << udata->rtol            << endl;
  cout << "  atol           = " << udata->atol            << endl;
  cout << "  order          = " << udata->order           << endl;
  cout << "  hf             = " << udata->hf              << endl;
  cout << "  controller     = " << udata->controller      << endl;
  cout << "  linear         = " << udata->linear          << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  lin iters      = " << udata->liniters        << endl;
  cout << "  eps lins       = " << udata->epslin          << endl;
  cout << "  prectype       = " << udata->prectype        << endl;
  cout << "  msbp           = " << udata->msbp            << endl;
  cout << "  pfmg_relax     = " << udata->pfmg_relax      << endl;
  cout << "  pfmg_nrelax    = " << udata->pfmg_nrelax     << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  output         = " << udata->output          << endl;
  cout << " --------------------------------- "           << endl;
  cout << endl;

  return 0;
}

// Initialize output
static int OpenOutput(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  // Header for status output
  if (udata->output > 0 && outproc)
  {
    cout << scientific;
    cout << setprecision(numeric_limits<realtype>::digits10);
    cout << "          t           ";
    cout << "          ||u||_rms      " << endl;
    cout << " ---------------------";
    cout << "-------------------------" << endl;
  }

  // Output problem information and open output streams
  if (udata->output == 2)
  {
    // Each processor outputs subdomain information
    stringstream fname;
    fname << "heat2d_info." << setfill('0') << setw(5) << udata->myid_c
          << ".txt";

    ofstream dout;
    dout.open(fname.str());
    dout <<  "xu  " << udata->xu       << endl;
    dout <<  "yu  " << udata->yu       << endl;
    dout <<  "nx  " << udata->nx       << endl;
    dout <<  "ny  " << udata->ny       << endl;
    dout <<  "px  " << udata->npx      << endl;
    dout <<  "py  " << udata->npy      << endl;
    dout <<  "np  " << udata->nprocs_w << endl;
    dout <<  "is  " << udata->is       << endl;
    dout <<  "ie  " << udata->ie       << endl;
    dout <<  "js  " << udata->js       << endl;
    dout <<  "je  " << udata->je       << endl;
    dout <<  "nt  " << udata->nout + 1 << endl;
    dout.close();

    // Open output streams for solution
    fname.str("");
    fname.clear();
    fname << "heat2d_solution." << setfill('0') << setw(5) << udata->myid_c
          << ".txt";
    udata->uout.open(fname.str());

    udata->uout << scientific;
    udata->uout << setprecision(numeric_limits<realtype>::digits10);
  }

  return 0;
}

// Write output
static int WriteOutput(realtype t, N_Vector u, UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  if (udata->output > 0)
  {

    // Compute rms norm
    realtype urms = sqrt(N_VDotProd(u, u) / udata->nx / udata->ny);

    // Output current status
    if (outproc)
    {
      cout << setw(22) << t << setw(25) << urms << endl;
    }

    // Write solution to disk
    if (udata->output == 2)
    {
      realtype *uarray = N_VGetArrayPointer(u);
      if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

      udata->uout << t << " ";
      for (sunindextype i = 0; i < udata->nodes_loc; i++)
      {
        udata->uout << uarray[i] << " ";
      }
      udata->uout << endl;
    }
  }

  return 0;
}

// Finalize output
static int CloseOutput(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  // Footer for status output
  if (outproc && (udata->output > 0))
  {
    cout << " ---------------------";
    cout << "-------------------------" << endl;
    cout << endl;
  }

  if (udata->output == 2)
  {
    // Close output streams
    udata->uout.close();
  }

  return 0;
}

// Print fast integrator statistics
static int OutputStats(void *arkode_mem, UserData* udata)
{
  int flag;

  // Get integrator and solver stats
  long int nst, nst_a, netf, nfe, nfi, nni, ncfn, nli, nlcf, nsetups, nJv;
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return -1;
  flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_flag(&flag, "ARKStepGetNumStepAttempts", 1)) return -1;
  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  if (check_flag(&flag, "ARKStepGetNumErrTestFails", 1)) return -1;
  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  if (check_flag(&flag, "ARKStepGetNumRhsEvals", 1)) return -1;
  flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1)) return -1;
  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1)) return -1;
  flag = ARKStepGetNumLinIters(arkode_mem, &nli);
  if (check_flag(&flag, "ARKStepGetNumLinIters", 1)) return -1;
  flag = ARKStepGetNumLinConvFails(arkode_mem, &nlcf);
  if (check_flag(&flag, "ARKStepGetNumLinConvFails", 1)) return -1;
  flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  if (check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1)) return -1;
  flag = ARKStepGetNumJtimesEvals(arkode_mem, &nJv);
  if (check_flag(&flag, "ARKStepGetNumJtimesEvals", 1)) return -1;


  cout << fixed;
  cout << setprecision(6);

  cout << "  Steps            = " << nst     << endl;
  cout << "  Step attempts    = " << nst_a   << endl;
  cout << "  Error test fails = " << netf    << endl;
  cout << "  RHS evals        = " << nfi     << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << "  NLS fails        = " << ncfn    << endl;
  cout << "  LS iters         = " << nli     << endl;
  cout << "  LS fails         = " << nlcf    << endl;
  cout << "  LS setups        = " << nsetups << endl;
  cout << "  Jv products      = " << nJv     << endl;
  cout << endl;

  // Compute average nls iters per step attempt and ls iters per nls iter
  realtype avgnli = (realtype) nni / (realtype) nst_a;
  realtype avgli  = (realtype) nli / (realtype) nni;
  cout << "  Avg NLS iters per step attempt = " << avgnli << endl;
  cout << "  Avg LS iters per NLS iter      = " << avgli  << endl;
  cout << endl;

  // Get preconditioner stats
  if (udata->prectype != SUN_PREC_NONE)
  {
    long int npe, nps;
    flag = ARKStepGetNumPrecEvals(arkode_mem, &npe);
    if (check_flag(&flag, "ARKStepGetNumPrecEvals", 1)) return -1;
    flag = ARKStepGetNumPrecSolves(arkode_mem, &nps);
    if (check_flag(&flag, "ARKStepGetNumPrecSolves", 1)) return -1;

    cout << "  Preconditioner setups = " << npe << endl;
    cout << "  Preconditioner solves = " << nps << endl;
    cout << "  PFMG iters            = " << udata->pfmg_its << endl;
    cout << endl;
  }

  return 0;
}

static int OutputTiming(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(6);
  }

  double maxtime = 0.0;

  MPI_Reduce(&(udata->evolvetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Evolve time            = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->rhstime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Full RHS time          = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->explicitrhstime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Explicit RHS time      = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->implicitrhstime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Implicit RHS time      = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->exchangetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Exchange time          = " << maxtime << " sec" << endl;
    cout << endl;
  }

  if (udata->prectype != SUN_PREC_NONE)
  {
    MPI_Reduce(&(udata->matfilltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
               udata->comm_c);
    if (outproc)
    {
      cout << "  MatFill time          = " << maxtime << " sec" << endl;
    }

    MPI_Reduce(&(udata->psetuptime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
               udata->comm_c);
    if (outproc)
    {
      cout << "  PSetup time           = " << maxtime << " sec" << endl;
    }

    MPI_Reduce(&(udata->psolvetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
               udata->comm_c);
    if (outproc)
    {
      cout << "  PSolve time           = " << maxtime << " sec" << endl;
      cout << endl;
    }
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
