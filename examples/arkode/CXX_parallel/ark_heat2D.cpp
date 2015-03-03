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
 * The following test simulates a simple anisotropic 2D heat 
 * equation,
 *    u_t = kx*u_xx + ky*u_yy + h,
 * for t in [0, 10], (x,y) in [0, 1]^2, with initial conditions
 *    u(0,x,y) =  0,
 * stationary boundary conditions, i.e. 
 *    u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
 * and a heat source of the form
 *    h(x,y) = sin(pi*x)*sin(2*pi*y).
 *
 * Under this setup, the problem has an analytical solution:
 *    u(t,x,y) = a(t)*sin(pi*x)*sin(2*pi*y), where
 *    a(t) = (1 - exp(-(kx+4*ky)*pi^2*t)) / ((kx+4*ky)*pi^2).
 * 
 * The spatial derivatives are computed using second-order 
 * centered differences, with the data distributed over nx*ny
 * points on a uniform spatial grid.
 *
 * This program solves the problem with a DIRK method.  This 
 * employs a Newton iteration with the PCG iterative linear solver, 
 * which itself uses a Jacobi preconditioner.  The example uses the 
 * built-in finite-difference Jacobian-vector product routine, but 
 * supplies both the RHS and preconditioner setup/solve functions.
 *
 * 20 outputs are printed at equal intervals, and run statistics 
 * are printed at the end.
 *---------------------------------------------------------------*/

// Header files 
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "arkode/arkode.h"            // prototypes for ARKode fcts., consts. 
#include "nvector/nvector_parallel.h" // parallel N_Vector types, fcts., macros 
#include "arkode/arkode_pcg.h"        // prototype for ARKPcg solver 
#include "sundials/sundials_types.h"  // def. of type 'realtype' 
#include "mpi.h"                      // MPI header file

using namespace std;

// accessor macros between (x,y) location and 1D NVector array
#define IDX(x,y,n) ((n)*(y)+(x))
#define PI RCONST(3.1415926535897932)
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)

// user data structure 
typedef struct {
  long int nx;          // global number of x grid points 
  long int ny;          // global number of y grid points
  long int is;          // global x indices of this subdomain
  long int ie;
  long int js;          // global y indices of this subdomain
  long int je;
  long int nxl;         // local number of x grid points 
  long int nyl;         // local number of y grid points 
  realtype dx;          // x-directional mesh spacing 
  realtype dy;          // y-directional mesh spacing 
  realtype kx;          // x-directional diffusion coefficient 
  realtype ky;          // y-directional diffusion coefficient 
  N_Vector h;           // heat source vector
  N_Vector d;           // inverse of Jacobian diagonal
  MPI_Comm comm;        // communicator object
  int myid;             // MPI process ID
  int nprocs;           // total number of MPI processes
  bool HaveBdry[2][2];  // flags denoting if on physical boundary
  realtype *Erecv;      // receive buffers for neighbor exchange
  realtype *Wrecv;
  realtype *Nrecv;
  realtype *Srecv;
  realtype *Esend;      // send buffers for neighbor exchange
  realtype *Wsend;
  realtype *Nsend;
  realtype *Ssend;
} UserData;

// User-supplied Functions Called by the Solver 
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int PSet(realtype t, N_Vector y, N_Vector fy, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int PSol(realtype t, N_Vector y, N_Vector fy, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, 
		void *user_data, N_Vector tmp);

// Private functions 
//    checks function return values 
static int check_flag(void *flagvalue, const string funcname, int opt);
//    sets default values into UserData structure
static int InitUserData(UserData *udata);
//    sets up parallel decomposition
static int SetupDecomp(UserData *udata);
//    performs neighbor exchange
static int Exchange(N_Vector y, UserData *udata);
//    frees memory allocated within UserData
static int FreeUserData(UserData *udata);

// Main Program 
int main(int argc, char* argv[]) {

  // general problem parameters 
  realtype T0 = RCONST(0.0);   // initial time 
  realtype Tf = RCONST(0.3);   // final time 
  int Nt = 20;                 // total number of output times 
  long int nx = 60;            // spatial mesh size
  long int ny = 120;
  realtype kx = 0.5;           // heat conductivity coefficients
  realtype ky = 0.75;
  realtype rtol = 1.e-5;       // relative and absolute tolerances
  realtype atol = 1.e-10;
  UserData *udata = NULL;
  realtype *data;
  long int N, Ntot, i, j;

  // general problem variables 
  int flag;                      // reusable error-checking flag 
  int myid;                      // MPI process ID
  N_Vector y = NULL;             // empty vector for storing solution 
  void *arkode_mem = NULL;       // empty ARKode memory structure 

  // initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;
  flag = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

  // allocate and fill udata structure 
  udata = new UserData;
  flag = InitUserData(udata);
  if (check_flag(&flag, "InitUserData", 1)) return 1;
  udata->nx = nx;
  udata->ny = ny;
  udata->kx = kx;
  udata->ky = ky;
  udata->dx = RCONST(1.0)/(1.0*nx-1.0);   // x mesh spacing 
  udata->dy = RCONST(1.0)/(1.0*ny-1.0);   // y mesh spacing 

  // Set up parallel decomposition
  flag = SetupDecomp(udata);
  if (check_flag(&flag, "SetupDecomp", 1)) return 1;

  // Initial problem output 
  bool outproc = (udata->myid == 0);
  if (outproc) {
    cout << "\n2D Heat PDE test problem:\n";
    cout << "   nprocs = " << udata->nprocs << "\n";
    cout << "   nx = " << udata->nx << "\n";
    cout << "   ny = " << udata->ny << "\n";
    cout << "   kx = " << udata->kx << "\n";
    cout << "   ky = " << udata->ky << "\n";
    cout << "   rtol = " << rtol << "\n";
    cout << "   atol = " << atol << "\n";
    cout << "   nxl (proc 0) = " << udata->nxl << "\n";
    cout << "   nyl (proc 0) = " << udata->nyl << "\n\n";
  }

  // Initialize data structures 
  N = (udata->nxl)*(udata->nyl);
  Ntot = nx*ny;
  y = N_VNew_Parallel(udata->comm, N, Ntot);         // Create parallel vector for solution 
  if (check_flag((void *) y, "N_VNew_Parallel", 0)) return 1;
  N_VConst(0.0, y);                                  // Set initial conditions 
  udata->h = N_VNew_Parallel(udata->comm, N, Ntot);  // Create vector for heat source
  if (check_flag((void *) udata->h, "N_VNew_Parallel", 0)) return 1;
  udata->d = N_VNew_Parallel(udata->comm, N, Ntot);  // Create vector for Jacobian diagonal
  if (check_flag((void *) udata->d, "N_VNew_Parallel", 0)) return 1;
  arkode_mem = ARKodeCreate();                       // Create the solver memory 
  if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;

  // fill in the heat source array
  data = N_VGetArrayPointer(udata->h);
  for (j=0; j<udata->nyl; j++)
    for (i=0; i<udata->nxl; i++)
      data[IDX(i,j,udata->nxl)] = sin(PI*(udata->is+i)*udata->dx) 
                                * sin(TWO*PI*(udata->js+j)*udata->dy);

  /* Call ARKodeInit to initialize the integrator memory and specify the
     right-hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  // Set routines 
  flag = ARKodeSetUserData(arkode_mem, (void *) udata);   // Pass udata to user functions 
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
  flag = ARKodeSetNonlinConvCoef(arkode_mem, 1.e-7);      // Update solver convergence coeff.
  if (check_flag(&flag, "ARKodeSetNonlinConvCoef", 1)) return 1;
  flag = ARKodeSStolerances(arkode_mem, rtol, atol);      // Specify tolerances 
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  // Linear solver specification 
  flag = ARKPcg(arkode_mem, 1, 20);                           // Specify the PCG solver 
  if (check_flag(&flag, "ARKPcg", 1)) return 1;
  flag = ARKSpilsSetPreconditioner(arkode_mem, PSet, PSol);   // Specify the Preconditoner
  if (check_flag(&flag, "ARKSpilsSetPreconditioner", 1)) return 1;

  // Specify linearly implicit RHS, with non-time-dependent preconditioner
  flag = ARKodeSetLinear(arkode_mem, 0);
  if (check_flag(&flag, "ARKodeSetLinear", 1)) return 1;

  // Each processor outputs subdomain information
  char outname[100];
  sprintf(outname, "heat2d_subdomain.%03i.txt", udata->myid);
  FILE *UFID = fopen(outname,"w");
  fprintf(UFID, "%li  %li  %li  %li  %li  %li\n", 
	  udata->nx, udata->ny, udata->is, udata->ie, udata->js, udata->je);
  fclose(UFID);

  // Open output streams for results, access data array 
  sprintf(outname, "heat2d.%03i.txt", udata->myid);
  UFID = fopen(outname,"w");
  data = N_VGetArrayPointer(y);

  // output initial condition to disk 
  for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
  fprintf(UFID,"\n");

  /* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  realtype t = T0;
  realtype dTout = (Tf-T0)/Nt;
  realtype tout = T0+dTout;
  realtype urms = sqrt(N_VDotProd(y,y)/nx/ny);
  if (outproc) {
    cout << "        t      ||u||_rms\n";
    cout << "   ----------------------\n";
    printf("  %10.6f  %10.6f\n", t, urms);
  }
  int iout;
  for (iout=0; iout<Nt; iout++) {

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);         // call integrator 
    if (check_flag(&flag, "ARKode", 1)) break;
    urms = sqrt(N_VDotProd(y,y)/nx/ny);
    if (outproc)  printf("  %10.6f  %10.6f\n", t, urms);        // print solution stats 
    if (flag >= 0) {                                            // successful solve: update output time
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                    // unsuccessful solve: break 
      if (outproc)  
	cerr << "Solver failure, stopping integration\n";
      break;
    }

    // output results to disk 
    for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
    fprintf(UFID,"\n");
  }
  if (outproc)  cout << "   ----------------------\n";
  fclose(UFID);

  // Print some final statistics 
  long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nlcf, nni, ncfn, netf, npe, nps;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  if (check_flag(&flag, "ARKodeGetNumSteps", 1)) return 1;
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_flag(&flag, "ARKodeGetNumStepAttempts", 1)) return 1;
  flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) return 1;
  flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  if (check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1)) return 1;
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  if (check_flag(&flag, "ARKodeGetNumErrTestFails", 1)) return 1;
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1)) return 1;
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1)) return 1;
  flag = ARKSpilsGetNumLinIters(arkode_mem, &nli);
  if (check_flag(&flag, "ARKSpilsGetNumLinIters", 1)) return 1;
  flag = ARKSpilsGetNumJtimesEvals(arkode_mem, &nJv);
  if (check_flag(&flag, "ARKSpilsGetNumJtimesEvals", 1)) return 1;
  flag = ARKSpilsGetNumConvFails(arkode_mem, &nlcf);
  if (check_flag(&flag, "ARKSpilsGetNumConvFails", 1)) return 1;
  flag = ARKSpilsGetNumPrecEvals(arkode_mem, &npe);
  if (check_flag(&flag, "ARKSpilsGetNumPrecEvals", 1)) return 1;
  flag = ARKSpilsGetNumPrecSolves(arkode_mem, &nps);

  if (outproc) {
    cout << "\nFinal Solver Statistics:\n";
    cout << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n";
    cout << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n";
    cout << "   Total linear solver setups = " << nsetups << "\n";
    cout << "   Total linear iterations = " << nli << "\n";
    cout << "   Total number of Jacobian-vector products = " << nJv << "\n";
    cout << "   Total number of Preconditioner setups = " << npe << "\n";
    cout << "   Total number of Preconditioner solves = " << nps << "\n";
    cout << "   Total number of linear solver convergence failures = " << nlcf << "\n";
    cout << "   Total number of Newton iterations = " << nni << "\n";
    cout << "   Total number of nonlinear solver convergence failures = " << ncfn << "\n";
    cout << "   Total number of error test failures = " << netf << "\n";
  }

  // Clean up and return with successful completion 
  N_VDestroy_Parallel(y);        // Free vectors 
  N_VDestroy_Parallel(udata->h);
  N_VDestroy_Parallel(udata->d);
  FreeUserData(udata);         // Free user data 
  delete udata;       
  ARKodeFree(&arkode_mem);     // Free integrator memory 
  flag = MPI_Finalize();       // Finalize MPI
  return 0;
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

// f routine to compute the ODE RHS function f(t,y). 
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(0.0, ydot);                           // Initialize ydot to zero 
  UserData *udata = (UserData *) user_data;      // access problem data 
  long int nxl = udata->nxl;                     // set variable shortcuts 
  long int nyl = udata->nyl;
  realtype kx = udata->kx;
  realtype ky = udata->ky;
  realtype dx = udata->dx;
  realtype dy = udata->dy;
  realtype *Y = N_VGetArrayPointer(y);           // access data arrays 
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return -1;
  realtype *Ydot = N_VGetArrayPointer(ydot);
  if (check_flag((void *) Ydot, "N_VGetArrayPointer", 0)) return -1;

  // Exchange boundary data with neighbors
  int ierr = Exchange(y, udata);
  if (check_flag(&ierr, "Exchange", 1)) return -1;

  // iterate over subdomain interior, computing approximation to RHS
  realtype c1 = kx/dx/dx;
  realtype c2 = ky/dy/dy;
  realtype c3 = -TWO*(c1 + c2);
  long int i, j;
  for (j=1; j<nyl-1; j++)                        // diffusive terms
    for (i=1; i<nxl-1; i++)
      Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + Y[IDX(i+1,j,nxl)])
                         + c2*(Y[IDX(i,j-1,nxl)] + Y[IDX(i,j+1,nxl)])
                         + c3*Y[IDX(i,j,nxl)];
  
  // iterate over subdomain boundaries (if not at overall domain boundary)
  if (!udata->HaveBdry[0][0]) {    // West face
    i=0;
    for (j=1; j<nyl-1; j++)  
      Ydot[IDX(i,j,nxl)] = c1*(udata->Wrecv[j]   + Y[IDX(i+1,j,nxl)])
                         + c2*(Y[IDX(i,j-1,nxl)] + Y[IDX(i,j+1,nxl)])
                         + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][1]) {    // East face
    i=nxl-1;
    for (j=1; j<nyl-1; j++)  
      Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + udata->Erecv[j])
                         + c2*(Y[IDX(i,j-1,nxl)] + Y[IDX(i,j+1,nxl)])
                         + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[1][0]) {    // South face
    j=0;
    for (i=1; i<nxl-1; i++)  
      Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + Y[IDX(i+1,j,nxl)])
                         + c2*(udata->Srecv[i]   + Y[IDX(i,j+1,nxl)])
                         + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[1][1]) {    // West face
    j=nyl-1;
    for (i=1; i<nxl-1; i++)  
      Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + Y[IDX(i+1,j,nxl)])
                         + c2*(Y[IDX(i,j-1,nxl)] + udata->Nrecv[i])
                         + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][0] && !udata->HaveBdry[1][0]) {  // South-West corner
    i = 0; 
    j = 0;
    Ydot[IDX(i,j,nxl)] = c1*(udata->Wrecv[j] + Y[IDX(i+1,j,nxl)])
                       + c2*(udata->Srecv[i] + Y[IDX(i,j+1,nxl)])
                       + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][0] && !udata->HaveBdry[1][1]) {  // North-West corner
    i = 0; 
    j = nyl-1;
    Ydot[IDX(i,j,nxl)] = c1*(udata->Wrecv[j]   + Y[IDX(i+1,j,nxl)])
                       + c2*(Y[IDX(i,j-1,nxl)] + udata->Nrecv[i])
                       + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][1] && !udata->HaveBdry[1][0]) {  // South-East corner
    i = nxl-1; 
    j = 0;
    Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + udata->Erecv[j])
                       + c2*(udata->Srecv[i]   + Y[IDX(i,j+1,nxl)])
                       + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][1] && !udata->HaveBdry[1][1]) {  // North-East corner
    i = nxl-1; 
    j = nyl-1;
    Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + udata->Erecv[j])
                       + c2*(Y[IDX(i,j-1,nxl)] + udata->Nrecv[i])
                       + c3*Y[IDX(i,j,nxl)];
  }

  // add in heat source
  N_VLinearSum(1.0, ydot, 1.0, udata->h, ydot);
  return 0;                                      // Return with success 
}

// Preconditioner setup routine (fills inverse of Jacobian diagonal)
static int PSet(realtype t, N_Vector y, N_Vector fy, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData *udata = (UserData *) user_data;      // variable shortcuts 
  realtype kx = udata->kx;
  realtype ky = udata->ky;
  realtype dx = udata->dx;
  realtype dy = udata->dy;
  realtype *diag = N_VGetArrayPointer(tmp1);  // access data arrays 
  if (check_flag((void *) diag, "N_VGetArrayPointer", 0)) return -1;

  // set all entries of tmp1 to the diagonal values of interior
  // (since boundary RHS is 0, set boundary diagonals to the same)
  realtype c = ONE + gamma*TWO*(kx/dx/dx + ky/dy/dy);
  N_VConst(c, tmp1);
  N_VInv(tmp1, udata->d);      // set d to inverse of diagonal
  return 0;                    // Return with success 
}

// Preconditioner solve routine
static int PSol(realtype t, N_Vector y, N_Vector fy, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, 
		void *user_data, N_Vector tmp)
{
  UserData *udata = (UserData *) user_data;  // access user_data structure
  N_VProd(r, udata->d, z);                   // perform Jacobi iteration
  return 0;                                  // Return with success 
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer  
*/
static int check_flag(void *flagvalue, const string funcname, int opt)
{
  int *errflag;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL) {
    cerr << "\nSUNDIALS_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return 1; }

  // Check if flag < 0
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      cerr << "\nSUNDIALS_ERROR: " << funcname << " failed with flag = " << *errflag << "\n\n";
      return 1; 
    }
  }
  
  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && flagvalue == NULL) {
    cerr << "\nMEMORY_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return 1; }

  return 0;
}

// Set up parallel decomposition
static int SetupDecomp(UserData *udata)
{
  // check that this has not been called before
  if (udata->Erecv != NULL || udata->Wrecv != NULL || 
      udata->Srecv != NULL || udata->Nrecv != NULL) {
    cerr << "SetupDecomp warning: parallel decomposition already set up\n";
    return 1;
  }

  // get suggested parallel decomposition
  int ierr, dims[] = {0, 0};
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &(udata->nprocs));
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Comm_size = " << ierr << "\n";
    return -1;
  }
  ierr = MPI_Dims_create(udata->nprocs, 2, dims);
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Dims_create = " << ierr << "\n";
    return -1;
  }

  // set up 2D Cartesian communicator
  int periods[] = {0, 0};
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &(udata->comm));
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Cart_create = " << ierr << "\n";
    return -1;
  }
  ierr = MPI_Comm_rank(udata->comm, &(udata->myid));
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Comm_rank = " << ierr << "\n";
    return -1;
  }

  // determine local extents
  int coords[2];
  ierr = MPI_Cart_get(udata->comm, 2, dims, periods, coords);
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Cart_get = " << ierr << "\n";
    return -1;
  }
  udata->is = (udata->nx)*(coords[0])/(dims[0]);
  udata->ie = (udata->nx)*(coords[0]+1)/(dims[0])-1;
  udata->js = (udata->ny)*(coords[1])/(dims[1]);
  udata->je = (udata->ny)*(coords[1]+1)/(dims[1])-1;
  udata->nxl = (udata->ie)-(udata->is)+1;
  udata->nyl = (udata->je)-(udata->js)+1;

  // determine if I have neighbors, and allocate exchange buffers
  udata->HaveBdry[0][0] = (udata->is == 0);
  udata->HaveBdry[0][1] = (udata->ie == udata->nx-1);
  udata->HaveBdry[1][0] = (udata->js == 0);
  udata->HaveBdry[1][1] = (udata->je == udata->ny-1);
  if (!udata->HaveBdry[0][0]) {
    udata->Wrecv = new realtype[udata->nyl];
    udata->Wsend = new realtype[udata->nyl];
  }
  if (!udata->HaveBdry[0][1]) {
    udata->Erecv = new realtype[udata->nyl];
    udata->Esend = new realtype[udata->nyl];
  }
  if (!udata->HaveBdry[1][0]) {
    udata->Srecv = new realtype[udata->nxl];
    udata->Ssend = new realtype[udata->nxl];
  }
  if (!udata->HaveBdry[1][1]) {
    udata->Nrecv = new realtype[udata->nxl];
    udata->Nsend = new realtype[udata->nxl];
  }

  return 0;     // return with success flag
}

// Perform neighbor exchange
static int Exchange(N_Vector y, UserData *udata)
{
  // local variables
  MPI_Request reqSW, reqSE, reqSS, reqSN, reqRW, reqRE, reqRS, reqRN;
  MPI_Status stat;
  int ierr, i, ipW=-1, ipE=-1, ipS=-1, ipN=-1;
  int coords[2], dims[2], periods[2], nbcoords[2];
  int nyl = udata->nyl;
  int nxl = udata->nxl;

  // access data array
  realtype *Y = N_VGetArrayPointer(y);
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return -1;

  // MPI equivalent of realtype type
  #if defined(SUNDIALS_SINGLE_PRECISION)
  #define REALTYPE_MPI_TYPE MPI_FLOAT
  #elif defined(SUNDIALS_DOUBLE_PRECISION)
  #define REALTYPE_MPI_TYPE MPI_DOUBLE
  #elif defined(SUNDIALS_EXTENDED_PRECISION)
  #define REALTYPE_MPI_TYPE MPI_LONG_DOUBLE
  #endif

  // MPI neighborhood information
  ierr = MPI_Cart_get(udata->comm, 2, dims, periods, coords);
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Cart_get = " << ierr << "\n";
    return -1;
  }
  if (!udata->HaveBdry[0][0]) {
    nbcoords[0] = coords[0]-1; 
    nbcoords[1] = coords[1];
    ierr = MPI_Cart_rank(udata->comm, nbcoords, &ipW);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    nbcoords[0] = coords[0]+1; 
    nbcoords[1] = coords[1];
    ierr = MPI_Cart_rank(udata->comm, nbcoords, &ipE);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    nbcoords[0] = coords[0]; 
    nbcoords[1] = coords[1]-1;
    ierr = MPI_Cart_rank(udata->comm, nbcoords, &ipS);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    nbcoords[0] = coords[0]; 
    nbcoords[1] = coords[1]+1;
    ierr = MPI_Cart_rank(udata->comm, nbcoords, &ipN);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  
  // open Irecv buffers
  if (!udata->HaveBdry[0][0]) {
    ierr = MPI_Irecv(udata->Wrecv, udata->nyl, REALTYPE_MPI_TYPE, ipW,
                   MPI_ANY_TAG, udata->comm, &reqRW);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    ierr = MPI_Irecv(udata->Erecv, udata->nyl, REALTYPE_MPI_TYPE, ipE,
                   MPI_ANY_TAG, udata->comm, &reqRE);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    ierr = MPI_Irecv(udata->Srecv, udata->nxl, REALTYPE_MPI_TYPE, ipS,
                   MPI_ANY_TAG, udata->comm, &reqRS);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    ierr = MPI_Irecv(udata->Nrecv, udata->nxl, REALTYPE_MPI_TYPE, ipN,
                   MPI_ANY_TAG, udata->comm, &reqRN);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }

  // send data
  if (!udata->HaveBdry[0][0]) {
    for (i=0; i<nyl; i++)  udata->Wsend[i] = Y[IDX(0,i,nxl)];
    ierr = MPI_Isend(udata->Wsend, udata->nyl, REALTYPE_MPI_TYPE, ipW, 0,
              udata->comm, &reqSW);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    for (i=0; i<nyl; i++)  udata->Esend[i] = Y[IDX(nxl-1,i,nxl)];
    ierr = MPI_Isend(udata->Esend, udata->nyl, REALTYPE_MPI_TYPE, ipE, 1,
              udata->comm, &reqSE);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    for (i=0; i<nxl; i++)  udata->Ssend[i] = Y[IDX(i,0,nxl)];
    ierr = MPI_Isend(udata->Ssend, udata->nxl, REALTYPE_MPI_TYPE, ipS, 2,
              udata->comm, &reqSS);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    for (i=0; i<nxl; i++)  udata->Nsend[i] = Y[IDX(i,nyl-1,nxl)];
    ierr = MPI_Isend(udata->Nsend, udata->nxl, REALTYPE_MPI_TYPE, ipN, 3,
              udata->comm, &reqSN);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }

  // wait for messages to finish
  if (!udata->HaveBdry[0][0]) {
    ierr = MPI_Wait(&reqRW, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSW, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    ierr = MPI_Wait(&reqRE, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSE, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    ierr = MPI_Wait(&reqRS, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSS, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    ierr = MPI_Wait(&reqRN, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSN, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }

  return 0;     // return with success flag
}

// Initialize memory allocated within Userdata
static int InitUserData(UserData *udata)
{
  udata->nx = 0;
  udata->ny = 0;
  udata->is = 0;
  udata->ie = 0;  
  udata->js = 0;
  udata->je = 0;  
  udata->nxl = 0;
  udata->nyl = 0;
  udata->dx = 0.0;
  udata->dy = 0.0;
  udata->kx = 0.0;
  udata->ky = 0.0;
  udata->h = NULL;
  udata->d = NULL;
  udata->comm = MPI_COMM_WORLD;
  udata->myid = 0;
  udata->nprocs = 0;
  udata->HaveBdry[0][0] = 1;
  udata->HaveBdry[0][1] = 1;
  udata->HaveBdry[1][0] = 1;
  udata->HaveBdry[1][1] = 1;
  udata->Erecv = NULL;
  udata->Wrecv = NULL;
  udata->Nrecv = NULL;
  udata->Srecv = NULL;
  udata->Esend = NULL;
  udata->Wsend = NULL;
  udata->Nsend = NULL;
  udata->Ssend = NULL;

  return 0;     // return with success flag
}

// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{
  // free exchange buffers
  if (udata->Wrecv != NULL)  delete[] udata->Wrecv;
  if (udata->Wsend != NULL)  delete[] udata->Wsend;
  if (udata->Erecv != NULL)  delete[] udata->Erecv;
  if (udata->Esend != NULL)  delete[] udata->Esend;
  if (udata->Srecv != NULL)  delete[] udata->Srecv;
  if (udata->Ssend != NULL)  delete[] udata->Ssend;
  if (udata->Nrecv != NULL)  delete[] udata->Nrecv;
  if (udata->Nsend != NULL)  delete[] udata->Nsend;

  return 0;     // return with success flag
}


//---- end of file ----
