/* -----------------------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ LLNL
 *                David Gardner @ LLNL
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
 * ---------------------------------------------------------------------------*/

#ifndef KIN_HEAT2D_NONLIN_HYPRE_PFMG_HPP
#define KIN_HEAT2D_NONLIN_HYPRE_PFMG_HPP

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include "kinsol/kinsol.h"                  // access to KINSOL
#include "nvector/nvector_parallel.h"       // access to the MPI N_Vector
#include "sunlinsol/sunlinsol_pcg.h"        // access to PCG SUNLinearSolver
#include "sunlinsol/sunlinsol_spgmr.h"      // access to SPGMR SUNLinearSolver
#include "HYPRE_struct_ls.h"                // HYPRE structured grid solver interface
#include "mpi.h"                            // MPI header file

// Macros for problem constants
#define PI    RCONST(3.141592653589793238462643383279502884197169)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define EIGHT RCONST(8.0)

// Macro to access (x,y) location in 1D NVector array
#define IDX(x,y,n) ((n)*(y)+(x))

// Define c function type
typedef int (*cFn)(N_Vector u, N_Vector z, void *user_data);

using namespace std;

// -----------------------------------------------------------------------------
// User data structure
// -----------------------------------------------------------------------------

struct UserData
{
  // Diffusion coefficients in the x and y directions
  realtype kx;
  realtype ky;

  // Enable/disable forcing
  bool forcing;

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

  // Overall number of local nodes
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

  // Fixed Point Solver settings
  realtype rtol;        // relative tolerance
  int      maa;         // m for Anderson Acceleration
  realtype damping;     // daming for Anderson Acceleration
  int      orthaa;      // orthogonalization routine for AA
  int      maxits;      // max number of fixed point iterations

  // Linear solver and preconditioner settings
  bool     lsinfo;    // output residual history
  int      liniters;  // number of linear iterations
  int      msbp;      // max number of steps between preconditioner setups
  realtype epslin;    // linear solver tolerance factor
  
  // Linear solver object
  SUNLinearSolver LS;  // linear solver memory structure

  // c(u) Function and integer for help setting
  cFn c;
  int c_int;

  // Vectors to help with FPFunction definition and execution
  N_Vector b;      // defined using c(u_exact)
  N_Vector vtemp;  // temporary vector for function evaluation

  // hypre objects
  HYPRE_StructGrid    grid;
  HYPRE_StructStencil stencil;
  HYPRE_StructMatrix  Jmatrix;
  HYPRE_StructVector  bvec;
  HYPRE_StructVector  xvec;
  HYPRE_StructVector  vvec;
  HYPRE_StructVector  Jvvec;
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
  N_Vector e;      // error vector

  // Timing variables
  bool   timing;     // print timings
  double totaltime;
  double fevaltime;
  double matfilltime;
  double jvtime;
  double psetuptime;
  double psolvetime;

  // Fused operations variable
  bool fusedops;
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS iterator
// -----------------------------------------------------------------------------

// Nonlinear fixed point function
static int FPFunction(N_Vector u, N_Vector f, void *user_data);

// Jacobian-vector product function
static int JTimes(void *user_data, N_Vector v, N_Vector Jv);

// Preconditioner setup and solve functions
static int PSetup(void *user_data);

static int PSolve(void *user_data, N_Vector r, N_Vector z,
                  realtype tol, int lr);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Setup the parallel decomposition
static int SetupDecomp(MPI_Comm comm_w, UserData *udata);

// Create rhs = b - c(u) for FPFunction
static int SetupRHS(void *user_data);

// Set nonlinear function c(u)
static int SetC(UserData *udata);

// Create hypre objects
static int SetupHypre(UserData *udata);

// Create PCG Linear Solver and attach hypre
static int SetupLS(N_Vector u, void *user_data);

// Fill Jacobian and A = I - gamma * J
static int Jac(UserData *udata);

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Set the default values in the UserData structure
static int InitUserData(UserData *udata);

// Free memory allocated within UserData
static int FreeUserData(UserData *udata);

// Read the command line inputs and set UserData values
static int ReadInputs(int *argc, char ***argv, UserData *udata, bool outproc);

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the true solution
static int Solution(N_Vector u, UserData *udata);

// Compute the solution error solution
static int SolutionError(N_Vector u,  N_Vector e, UserData *udata);

// Print the command line options
static void InputHelp();

// Print some UserData information
static int PrintUserData(UserData *udata);

// Print solver statistics
static int OutputStats(void *kinsol_mem, UserData *udata);

// Print solver timing
static int OutputTiming(UserData *udata);

// Check function return values
static int check_flag(void *flagvalue, const string funcname, int opt);

// -----------------------------------------------------------------------------
// Multiple nonlinear functions for testing
// -----------------------------------------------------------------------------

// c(u) = u
static int c1(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];

      zarray[IDX(i,j,nx_loc)] = u_val;
    }
  }

  // Return success
  return 0;
}

// c(u) = u^3 - u
static int c2(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];

      zarray[IDX(i,j,nx_loc)] = u_val * u_val * u_val - u_val;
    }
  }

  // Return success
  return 0;
}

// c(u) = u - u^2
static int c3(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];

      zarray[IDX(i,j,nx_loc)] = u_val - u_val * u_val;
    }
  }

  // Return success
  return 0;
}

// c(u) = e^u
static int c4(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];

      zarray[IDX(i,j,nx_loc)] = exp(u_val);
    }
  }

  // Return success
  return 0;
}

// c(u) = u^4
static int c5(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];

      zarray[IDX(i,j,nx_loc)] = u_val * u_val * u_val * u_val;
    }
  }

  // Return success
  return 0;
}

// c(u) = cos^2(u) - sin^2(u)
static int c6(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];

      zarray[IDX(i,j,nx_loc)] = (cos(u_val) * cos(u_val)) - (sin(u_val) * sin(u_val));
    }
  }

  // Return success
  return 0;
}

// c(u) = cos^2(u) - sin^2(u) - e^u
static int c7(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];

      zarray[IDX(i,j,nx_loc)] = (cos(u_val) * cos(u_val)) - (sin(u_val) * sin(u_val)) - exp(u_val);
    }
  }

  // Return success
  return 0;
}

// c(u) = e^u * u^4 - u * e^{cos(u)}
static int c8(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, u2;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      u2 = u_val * u_val;

      zarray[IDX(i,j,nx_loc)] = exp(u_val) * u2 * u2 - u_val * exp(cos(u_val));
    }
  }

  // Return success
  return 0;
}

// c(u) = e^(cos^2(u))
static int c9(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, cos2u;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      cos2u = cos(u_val) * cos(u_val);

      zarray[IDX(i,j,nx_loc)] = exp(cos2u);
    }
  }

  // Return success
  return 0;
}

// c(u) = 10(u - u^2) 
static int c10(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, u2;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      u2 = u_val * u_val;

      zarray[IDX(i,j,nx_loc)] = 10.0 * (u_val - u2);
    }
  }

  // Return success
  return 0;
}

// c(u) = -13 + u + ((5-u)u - 2)u
static int c11(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, temp;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      temp = ((5.0 - u_val) * u_val) - 2.0;

      zarray[IDX(i,j,nx_loc)] = -13.0 + u_val + temp * u_val;
    }
  }

  // Return success
  return 0;
}

// c(u) = sqrt(5) * (u - u^2) 
static int c12(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, u2;
  realtype temp = sqrt(5);

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      u2    = u_val * u_val;

      zarray[IDX(i,j,nx_loc)] = temp * (u_val - u2);
    }
  }

  // Return success
  return 0;
}

// c(u) = (u - e^u)^2 + (u + u * sin(u) - cos(u))^2 
static int c13(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, eu, usin, temp;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      eu    = u_val - exp(u_val);
      usin  = u_val * sin(u_val);
      temp  = (u_val + usin - cos(u_val));

      zarray[IDX(i,j,nx_loc)] = eu * eu + temp * temp;
    }
  }

  // Return success
  return 0;
}

// c(u) = u + ue^u + ue^{-u} 
static int c14(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, ueu, ue_u;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      ueu   = u_val * exp(u_val);
      ue_u  = u_val * exp(-u_val);

      zarray[IDX(i,j,nx_loc)] = u_val + ueu + ue_u;
    }
  }

  // Return success
  return 0;
}

// c(u) = u + ue^u + ue^{-u} + (u - e^u)^2 
static int c15(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, ueu, ue_u, temp;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      ueu  = u_val * exp(u_val);
      ue_u = u_val * exp(-u_val);
      temp = u_val - exp(u_val);

      zarray[IDX(i,j,nx_loc)] = u_val + ueu + ue_u + (temp * temp);
    }
  }

  // Return success
  return 0;
}

// c(u) = u + ue^u + ue^{-u} + (u - e^u)^2 + (u + usin(u) - cos(u))^2 
static int c16(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, ueu, ue_u, temp, temp2;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      ueu    = u_val * exp(u_val);
      ue_u   = u_val * exp(-u_val);
      temp   = u_val - exp(u_val);
      temp2  = u_val + (u_val * sin(u_val)) - cos(u_val);

      zarray[IDX(i,j,nx_loc)] = u_val + ueu + ue_u + (temp * temp) + (temp2 * temp2);
    }
  }

  // Return success
  return 0;
}

// c(u) = u + ue^{-u} + e^u*(u + sin(u) - cos(u))^3 
static int c17(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_flag((void *) zarray, "N_VGetArrayPointer", 0)) return 1;
  
  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u) 
  realtype u_val, ue_u, eu, temp;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];
      ue_u  = u_val * exp(-u_val);
      eu    = exp(u_val);
      temp  = u_val + sin(u_val) - cos(u_val);

      zarray[IDX(i,j,nx_loc)] = u_val + ue_u + eu * (temp * temp * temp);
    }
  }

  // Return success
  return 0;
}

#endif
