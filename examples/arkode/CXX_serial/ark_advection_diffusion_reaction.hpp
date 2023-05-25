/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * Header file for ARKODE advection-diffusion-reaction equation example, see
 * ark_advection_diffusion_reaction.cpp for more details.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_mristep.h"
#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_band.h"
#include "sunlinsol/sunlinsol_band.h"

// Macros for problem constants
#define PI    RCONST(3.141592653589793238462643383279502884197169)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define NSPECIES 3

#define WIDTH (10 + numeric_limits<realtype>::digits10)

// Macro to access each species at an x location
#define UIDX(i) (NSPECIES * (i))
#define VIDX(i) (NSPECIES * (i) + 1)
#define WIDX(i) (NSPECIES * (i) + 2)

using namespace std;

// -----------------------------------------------------------------------------
// Problem parameters
// -----------------------------------------------------------------------------

struct UserData
{
  // RHS options
  bool diffusion = true;
  bool advection = true;

  // Splitting option
  int splitting = 3;

  // Advection and diffusion coefficients
  realtype c = RCONST(1.0e-3);
  realtype d = RCONST(1.0e-2);

  // Feed and reaction rates
  realtype A = RCONST(0.6);
  realtype B = RCONST(2.0);

  // Stiffness parameter
  realtype eps = RCONST(1.0e-2);

  // Final simulation time
  realtype tf = RCONST(3.0);

  // Domain boundaries
  realtype xl = ZERO;
  realtype xu = ONE;

  // Number of nodes
  sunindextype nx = 512;

  // Mesh spacing
  realtype dx = (xu - xl) / (nx - 1);

  // Number of equations
  sunindextype neq = NSPECIES * nx;

  // Temporary workspace vector and matrix
  N_Vector  temp_v = nullptr;
  SUNMatrix temp_J = nullptr;

  // Inner stepper memory
  MRIStepInnerStepper fast_mem = nullptr;

  ~UserData();
};


UserData::~UserData()
{
  if (temp_v)
  {
    N_VDestroy(temp_v);
    temp_v = nullptr;
  }

  if (temp_J)
  {
    SUNMatDestroy(temp_J);
    temp_J = nullptr;
  }
}

// -----------------------------------------------------------------------------
// Problem options
// -----------------------------------------------------------------------------

struct UserOptions
{
  // Integration method (0 = ERK, 1 = ARK, 2 = MRIARK, 3 = MRICVODE)
  int integrator = 1;

  // Method order
  int order      = 3;
  int order_fast = 3;
  bool ark_dirk  = false;

  // Relative and absolute tolerances
  realtype rtol      = RCONST(1.e-4);
  realtype atol      = RCONST(1.e-9);
  realtype rtol_fast = RCONST(1.e-4);
  realtype atol_fast = RCONST(1.e-9);

  // Step size selection (ZERO = adaptive steps)
  realtype fixed_h      = ZERO;
  realtype fixed_h_fast = ZERO;

  // First step growth factor
  realtype etamx1_fast = ZERO;

  int maxsteps      = 10000;  // max steps between outputs
  int controller    = -1;     // step size adaptivity method
  int predictor     = 0;      // predictor for nonlinear systems
  int ls_setup_freq = 0;      // linear solver setup frequency

  int controller_fast    = -1;  // fast step size adaptivity method
  int predictor_fast     = 0;   // predictor for fast nonlinear systems
  int ls_setup_freq_fast = 0;   // fast linear solver setup frequency

  bool linear = false;  // signal that the problem is linearly implicit

  // save and reuse prior fast time step sizes
  bool save_hinit = false;
  bool save_hcur  = false;

  realtype hcur_factor = RCONST(0.7);

  int output = 1;   // 0 = none, 1 = stats, 2 = disk, 3 = disk with tstop
  int nout   = 10;  // number of output times
  ofstream uout;    // output file stream
};

// -----------------------------------------------------------------------------
// Custom inner stepper content and functions
// -----------------------------------------------------------------------------

struct CVodeInnerStepperContent
{
  void* cvode_mem = nullptr;  // CVODE memory structure
  void* user_data = nullptr;  // user data pointer

  // saved step sizes
  bool save_hinit = false;
  bool save_hcur  = false;

  realtype hcur_factor = ONE;

  realtype hinit = ZERO;  // initial step size
  realtype hcur  = ZERO;  // current step size

  // saved integrator stats
  long int nst     = 0;  // time steps
  long int netf    = 0;  // error test fails
  long int nfe     = 0;  // rhs evals
  long int nni     = 0;  // nonlinear iterations
  long int nncf    = 0;  // nonlinear convergence failures
  long int nsetups = 0;  // linear solver setups
  long int nje     = 0;  // Jacobian evals
};

int CVodeInnerStepper_Evolve(MRIStepInnerStepper fast_mem, realtype t0,
                             realtype tout, N_Vector y);

int CVodeInnerStepper_FullRhs(MRIStepInnerStepper fast_mem, realtype t,
                              N_Vector y, N_Vector f, int mode);

int CVodeInnerStepper_Reset(MRIStepInnerStepper fast_mem, realtype tR,
                            N_Vector yR);

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrators
// -----------------------------------------------------------------------------

// ODE right hand side (RHS) functions
int f_advection(realtype t, N_Vector y, N_Vector f, void* user_data);
int f_diffusion(realtype t, N_Vector y, N_Vector f, void* user_data);
int f_reaction(realtype t, N_Vector y, N_Vector f, void* user_data);

int f_adv_diff(realtype t, N_Vector y, N_Vector f, void* user_data);
int f_adv_react(realtype t, N_Vector y, N_Vector f, void* user_data);
int f_diff_react(realtype t, N_Vector y, N_Vector f, void* user_data);
int f_adv_diff_react(realtype t, N_Vector y, N_Vector f, void* user_data);

int f_react_forcing(realtype t, N_Vector y, N_Vector f, void* user_data);

// Jacobian of RHS functions
int J_advection(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void* user_data, N_Vector tmp1, N_Vector tmp2,
                N_Vector tmp3);
int J_diffusion(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void* user_data, N_Vector tmp1, N_Vector tmp2,
                N_Vector tmp3);
int J_reaction(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2,
               N_Vector tmp3);
int J_adv_diff(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2,
               N_Vector tmp3);
int J_adv_react(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void* user_data, N_Vector tmp1, N_Vector tmp2,
                N_Vector tmp3);
int J_diff_react(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                 void* user_data, N_Vector tmp1, N_Vector tmp2,
                 N_Vector tmp3);
int J_adv_diff_react(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                     void* user_data, N_Vector tmp1, N_Vector tmp2,
                     N_Vector tmp3);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Integrator setup functions
int SetupERK(SUNContext ctx, UserData &udata, UserOptions &uopts, N_Vector y,
             void** arkode_mem);

int SetupARK(SUNContext ctx, UserData &udata, UserOptions &uopts, N_Vector y,
             SUNMatrix* A, SUNLinearSolver* LS, void** arkode_mem);

int SetupMRIARK(SUNContext ctx, UserData &udata, UserOptions &uopts,
                N_Vector y, SUNMatrix* A, SUNLinearSolver* LS,
                SUNMatrix* A_fast, SUNLinearSolver* LS_fast,
                MRIStepInnerStepper* fast_mem, void** arkode_mem);

int SetupMRICVODE(SUNContext ctx, UserData &udata, UserOptions &uopts,
                  N_Vector y, SUNMatrix* A, SUNLinearSolver* LS,
                  SUNMatrix* A_fast, SUNLinearSolver* LS_fast,
                  MRIStepInnerStepper* fast_mem, void** arkode_mem);

// Compute the initial condition
int SetIC(N_Vector y, UserData &udata);

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Check function return flag
int check_flag(int flag, const string funcname)
{
  if (!flag) return 0;
  if (flag < 0) cerr << "ERROR: ";
  cerr << funcname << " returned " << flag << endl;
  return 1;
}


// Check if a function returned a NULL pointer
int check_ptr(void *ptr, const string funcname)
{
  if (ptr) return 0;
  cerr << "ERROR: " << funcname << " returned NULL" << endl;
  return 1;
}


// Print ERK integrator statistics
int OutputStatsERK(void* arkode_mem, UserData &udata)
{
  int flag;

  // Get integrator and solver stats
  long int nst, nst_a, netf, nfe;
  flag = ERKStepGetNumSteps(arkode_mem, &nst);
  if (check_flag(flag, "ERKStepGetNumSteps")) return -1;
  flag = ERKStepGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_flag(flag, "ERKStepGetNumStepAttempts")) return -1;
  flag = ERKStepGetNumErrTestFails(arkode_mem, &netf);
  if (check_flag(flag, "ERKStepGetNumErrTestFails")) return -1;
  flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
  if (check_flag(flag, "ERKStepGetNumRhsEvals")) return -1;

  cout << "  Steps            = " << nst     << endl;
  cout << "  Step attempts    = " << nst_a   << endl;
  cout << "  Error test fails = " << netf    << endl;
  cout << "  RHS evals        = " << nfe     << endl;

  return 0;
}


// Print ARK integrator statistics
int OutputStatsARK(void* arkode_mem, UserData &udata)
{
  int flag;

  // Get integrator and solver stats
  long int nst, nst_a, netf, nfe, nfi;
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  if (check_flag(flag, "ARKStepGetNumSteps")) return -1;
  flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_flag(flag, "ARKStepGetNumStepAttempts")) return -1;
  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  if (check_flag(flag, "ARKStepGetNumErrTestFails")) return -1;
  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  if (check_flag(flag, "ARKStepGetNumRhsEvals")) return -1;

  cout << fixed << setprecision(6);
  cout << "  Steps              = " << nst     << endl;
  cout << "  Step attempts      = " << nst_a   << endl;
  cout << "  Error test fails   = " << netf    << endl;
  cout << "  Explicit RHS evals = " << nfe     << endl;
  cout << "  Implicit RHS evals = " << nfi     << endl;

  if (udata.splitting)
  {
    long int nni, ncfn;
    flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
    if (check_flag(flag, "ARKStepGetNumNonlinSolvIters")) return -1;
    flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
    if (check_flag(flag, "ARKStepGetNumNonlinSolvConvFails")) return -1;

    long int nsetups, nje;
    flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
    if (check_flag(flag, "ARKStepGetNumLinSolvSetups")) return -1;
    flag = ARKStepGetNumJacEvals(arkode_mem, &nje);
    if (check_flag(flag, "ARKStepGetNumJacEvals")) return -1;

    cout << "  NLS iters          = " << nni     << endl;
    cout << "  NLS fails          = " << ncfn    << endl;
    cout << "  LS setups          = " << nsetups << endl;
    cout << "  J evals            = " << nje     << endl;
    cout << endl;

    realtype avgnli = (realtype) nni / (realtype) nst_a;
    realtype avgls  = (realtype) nsetups / (realtype) nni;
    cout << "  Avg NLS iters per step attempt = " << avgnli << endl;
    cout << "  Avg LS setups per NLS iter     = " << avgls  << endl;
  }
  cout << endl;

  return 0;
}


// Print MRI integrator statistics
int OutputStatsMRIARK(void* arkode_mem, MRIStepInnerStepper fast_mem,
                      UserData &udata)
{
  int flag;

  // Get slow integrator and solver stats
  long int nst, nst_a, netf, nfe, nfi;
  flag = MRIStepGetNumSteps(arkode_mem, &nst);
  if (check_flag(flag, "MRIStepGetNumSteps")) return -1;
  flag = MRIStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  if (check_flag(flag, "MRIStepGetNumRhsEvals")) return -1;

  cout << fixed << setprecision(6);
  cout << endl << "Slow Integrator:" << endl;
  cout << "  Steps                   = " << nst << endl;
  cout << "  Slow explicit RHS evals = " << nfe << endl;
  cout << "  Slow implicit RHS evals = " << nfi << endl;

  if (udata.diffusion)
  {
    long int nni, ncfn;
    flag = MRIStepGetNumNonlinSolvIters(arkode_mem, &nni);
    if (check_flag(flag, "MRIStepGetNumNonlinSolvIters")) return -1;
    flag = MRIStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
    if (check_flag(flag, "MRIStepGetNumNonlinSolvConvFails")) return -1;

    long int nsetups, nje;
    flag = MRIStepGetNumLinSolvSetups(arkode_mem, &nsetups);
    if (check_flag(flag, "MRIStepGetNumLinSolvSetups")) return -1;
    flag = MRIStepGetNumJacEvals(arkode_mem, &nje);
    if (check_flag(flag, "MRIStepGetNumJacEvals")) return -1;

    cout << "  NLS iters               = " << nni << endl;
    cout << "  NLS fails               = " << ncfn << endl;
    cout << "  LS setups               = " << nsetups << endl;
    cout << "  J evals                 = " << nje << endl;
    cout << endl;

    // Compute average nls iters per step and ls iters per nls iter
    realtype avgnli = (realtype) nni / (realtype) nst;
    realtype avgls  = (realtype) nsetups / (realtype) nni;
    cout << "  Avg NLS iters per slow step = " << avgnli << endl;
    cout << "  Avg LS setups per NLS iter  = " << avgls  << endl;
  }
  cout << endl;

  // Get fast integrator stats and solver stats
  void* fast_arkode_mem;
  MRIStepInnerStepper_GetContent(fast_mem, &fast_arkode_mem);

  // Get fast integrator and solver stats
  flag = ARKStepGetNumSteps(fast_arkode_mem, &nst);
  if (check_flag(flag, "ARKStepGetNumSteps")) return -1;
  flag = ARKStepGetNumStepAttempts(fast_arkode_mem, &nst_a);
  if (check_flag(flag, "ARKStepGetNumStepAttempts")) return -1;
  flag = ARKStepGetNumErrTestFails(fast_arkode_mem, &netf);
  if (check_flag(flag, "ARKStepGetNumErrTestFails")) return -1;
  flag = ARKStepGetNumRhsEvals(fast_arkode_mem, &nfe, &nfi);
  if (check_flag(flag, "ARKStepGetNumRhsEvals")) return -1;

  cout << fixed << setprecision(6);
  cout << endl << "Fast Integrator:" << endl;
  cout << "  Steps                   = " << nst << endl;
  cout << "  Step attempts           = " << nst_a << endl;
  cout << "  Error test fails        = " << netf << endl;
  cout << "  Slow explicit RHS evals = " << nfe << endl;
  cout << "  Slow implicit RHS evals = " << nfi << endl;

  if (udata.splitting)
  {
    long int nni, ncfn;
    flag = ARKStepGetNumNonlinSolvIters(fast_arkode_mem, &nni);
    if (check_flag(flag, "ARKStepGetNumNonlinSolvIters")) return -1;
    flag = ARKStepGetNumNonlinSolvConvFails(fast_arkode_mem, &ncfn);
    if (check_flag(flag, "ARKStepGetNumNonlinSolvConvFails")) return -1;

    long int nsetups, nje;
    flag = ARKStepGetNumLinSolvSetups(fast_arkode_mem, &nsetups);
    if (check_flag(flag, "ARKStepGetNumLinSolvSetups")) return -1;
    flag = ARKStepGetNumJacEvals(fast_arkode_mem, &nje);
    if (check_flag(flag, "ARKStepGetNumJacEvals")) return -1;

    cout << "  NLS iters               = " << nni     << endl;
    cout << "  NLS fails               = " << ncfn    << endl;
    cout << "  LS setups               = " << nsetups << endl;
    cout << "  J evals                 = " << nje     << endl;
    cout << endl;

    realtype avgnli = (realtype) nni / (realtype) nst;
    realtype avgls  = (realtype) nsetups / (realtype) nni;
    cout << "  Avg NLS iters per fast step = " << avgnli << endl;
    cout << "  Avg LS setups per NLS iter  = " << avgls  << endl;
  }
  cout << endl;

  return 0;
}


// Save current stats
int UpdateCVodeStats(CVodeInnerStepperContent* content)
{
  int flag;
  long int nst, netf, nfe, nni, nncf, nsetups, nje;

  flag = CVodeGetNumSteps(content->cvode_mem, &nst);
  if (check_flag(flag, "CVodeGetNumSteps")) return -1;
  content->nst += nst;

  flag = CVodeGetNumErrTestFails(content->cvode_mem, &netf);
  if (check_flag(flag, "CVodeGetNumErrTestFails")) return -1;
  content->netf += netf;

  flag = CVodeGetNumRhsEvals(content->cvode_mem, &nfe);
  if (check_flag(flag, "CVodeGetNumRhsEvals")) return -1;
  content->nfe += nfe;

  flag = CVodeGetNumNonlinSolvIters(content->cvode_mem, &nni);
  if (check_flag(flag, "CVodeGetNumNonlinSolvIters")) return -1;
  content->nni += nni;

  flag = CVodeGetNumNonlinSolvConvFails(content->cvode_mem, &nncf);
  if (check_flag(flag, "CVodeGetNumNonlinSolvConvFails")) return -1;
  content->nncf += nncf;

  flag = CVodeGetNumLinSolvSetups(content->cvode_mem, &nsetups);
  if (check_flag(flag, "CVodeGetNumLinSolveSetups")) return -1;
  content->nsetups += nsetups;

  flag = CVodeGetNumJacEvals(content->cvode_mem, &nje);
  if (check_flag(flag, "CVodeGetNumJacEvals")) return -1;
  content->nje += nje;

  return 0;
}


// Print MRI integrator statistics
int OutputStatsMRICVODE(void* arkode_mem, MRIStepInnerStepper fast_mem,
                        UserData &udata)
{
  int flag;

  // Get slow integrator and solver stats
  long int nsts, nfse, nfsi;
  flag = MRIStepGetNumSteps(arkode_mem, &nsts);
  if (check_flag(flag, "MRIStepGetNumSteps")) return -1;
  flag = MRIStepGetNumRhsEvals(arkode_mem, &nfse, &nfsi);
  if (check_flag(flag, "MRIStepGetNumRhsEvals")) return -1;

  long int nni, ncfn;
  flag = MRIStepGetNumNonlinSolvIters(arkode_mem, &nni);
  if (check_flag(flag, "MRIStepGetNumNonlinSolvIters")) return -1;
  flag = MRIStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  if (check_flag(flag, "MRIStepGetNumNonlinSolvConvFails")) return -1;

  long int nsetups, nje;
  flag = MRIStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  if (check_flag(flag, "MRIStepGetNumLinSolvSetups")) return -1;
  flag = MRIStepGetNumJacEvals(arkode_mem, &nje);
  if (check_flag(flag, "MRIStepGetNumJacEvals")) return -1;

  cout << fixed << setprecision(6);
  cout << endl << "Slow Integrator:" << endl;
  cout << "  Steps                   = " << nsts    << endl;
  cout << "  Slow explicit RHS evals = " << nfse    << endl;
  cout << "  Slow implicit RHS evals = " << nfsi    << endl;
  cout << "  NLS iters               = " << nni     << endl;
  cout << "  NLS fails               = " << ncfn    << endl;
  cout << "  LS setups               = " << nsetups << endl;
  cout << "  J evals                 = " << nje     << endl;
  cout << endl;

  // Compute average nls iters per step and ls iters per nls iter
  realtype avgnli = (realtype) nni / (realtype) nsts;
  realtype avgls  = (realtype) nsetups / (realtype) nni;
  cout << "  Avg NLS iters per slow step = " << avgnli << endl;
  cout << "  Avg LS setups per NLS iter  = " << avgls  << endl;
  cout << endl;

  // Get fast integrator stats and solver stats
  void* inner_content;
  MRIStepInnerStepper_GetContent(fast_mem, &inner_content);
  CVodeInnerStepperContent* content = (CVodeInnerStepperContent*) inner_content;

  // Update CVODE stats
  flag = UpdateCVodeStats(content);
  if (check_flag(flag, "UpdateCVodeStats")) return -1;

  cout << fixed << setprecision(6);
  cout << "Fast Integrator:" << endl;
  cout << "  Steps            = " << content->nst     << endl;
  cout << "  Error test fails = " << content->netf    << endl;
  cout << "  Fast RHS evals   = " << content->nfe     << endl;
  cout << "  NLS iters        = " << content->nni     << endl;
  cout << "  NLS fails        = " << content->nncf    << endl;
  cout << "  LS setups        = " << content->nsetups << endl;
  cout << "  J evals          = " << content->nje     << endl;
  cout << endl;

  avgnli = (realtype) content->nni / (realtype) content->nst;
  avgls  = (realtype) content->nsetups / (realtype) content->nni;
  cout << "  Avg NLS iters per fast step = " << avgnli << endl;
  cout << "  Avg LS setups per NLS iter  = " << avgls  << endl;
  cout << endl;

  return 0;
}


// Print command line options
void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --no-advection           : disable advection\n";
  cout << "  --no-diffusion           : disable diffusion\n";
  cout << "  --splitting <int>        : RHS splitting\n";
  cout << "  --c <real>               : advection coefficient\n";
  cout << "  --d <real>               : diffusion coefficient\n";
  cout << "  --A <real>               : species A concentration\n";
  cout << "  --B <real>               : species B concentration\n";
  cout << "  --eps <real>             : stiffness parameter\n";
  cout << "  --tf <real>              : final time\n";
  cout << "  --xl <real>              : domain lower boundary\n";
  cout << "  --xu <real>              : domain upper boundary\n";
  cout << "  --nx <int>               : number of mesh points\n";
  cout << "  --integrator <int>       : integrator option\n";
  cout << "  --order <int>            : method order\n";
  cout << "  --order_fast <int>       : MRI fast method order\n";
  cout << "  --ark_dirk               : Use DIRK method from ARK method\n";
  cout << "  --rtol <real>            : relative tolerance\n";
  cout << "  --atol <real>            : absoltue tolerance\n";
  cout << "  --rtol_fast <real>       : MRI fast relative tolerance\n";
  cout << "  --atol_fast <real>       : MRI fast absolute tolerance\n";
  cout << "  --fixed_h <real>         : fixed step size\n";
  cout << "  --fixed_h_fast <real>    : MRI fast fixed step size\n";
  cout << "  --controller <int>       : time step adaptivity\n";
  cout << "  --predictor <int>        : nonlinear solver predictor\n";
  cout << "  --lssetupfreq <int>      : LS setup frequency\n";
  cout << "  --controller_fast <int>  : MRI fast time step adaptivity\n";
  cout << "  --predictor_fast <int>   : MRI fast nonlinear solver predictor\n";
  cout << "  --lssetupfreq_fast <int> : MRI fast LS setup frequency\n";
  cout << "  --maxsteps <int>         : max steps between outputs\n";
  cout << "  --etamx1_fast <real>     : max step size growth in first step\n";
  cout << "  --linear                 : linearly implicit\n";
  cout << "  --save_hinit             : reuse initial fast step\n";
  cout << "  --save_hcur              : reuse current fast step\n";
  cout << "  --hcur_factor            : current fast step safety factor\n";
  cout << "  --output <int>           : output level\n";
  cout << "  --nout <int>             : number of outputs\n";
  cout << "  --help                   : print options and exit\n";
}


inline void find_arg(vector<string> &args, const string key, realtype &dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    dest = stof(*(it + 1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    dest = stod(*(it + 1));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    dest = stold(*(it + 1));
#endif
    args.erase(it, it + 2);
  }
}


#if defined(SUNDIALS_INT64_T)
inline void find_arg(vector<string> &args, const string key, sunindextype &dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoll(*(it + 1));
    args.erase(it, it + 2);
  }
}
#endif


inline void find_arg(vector<string> &args, const string key, int &dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoi(*(it + 1));
    args.erase(it, it + 2);
  }
}


inline void find_arg(vector<string> &args, const string key, bool &dest,
                     bool store = true)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = store;
    args.erase(it);
  }
}


int ReadInputs(vector<string> &args, UserData &udata, UserOptions &uopts,
               SUNContext ctx)
{
  if (find(args.begin(), args.end(), "--help") != args.end())
  {
    InputHelp();
    return 1;
  }

  // Problem parameters
  find_arg(args, "--no-advection", udata.advection, false);
  find_arg(args, "--no-diffusion", udata.diffusion, false);
  find_arg(args, "--splitting", udata.splitting);
  find_arg(args, "--c", udata.c);
  find_arg(args, "--d", udata.d);
  find_arg(args, "--A", udata.A);
  find_arg(args, "--B", udata.B);
  find_arg(args, "--eps", udata.eps);
  find_arg(args, "--tf",  udata.tf);
  find_arg(args, "--xl",  udata.xl);
  find_arg(args, "--xu",  udata.xu);
  find_arg(args, "--nx",  udata.nx);

  // Integrator options
  find_arg(args, "--integrator", uopts.integrator);
  find_arg(args, "--order", uopts.order);
  find_arg(args, "--order_fast", uopts.order_fast);
  find_arg(args, "--ark_dirk", uopts.ark_dirk);
  find_arg(args, "--rtol", uopts.rtol);
  find_arg(args, "--atol", uopts.atol);
  find_arg(args, "--rtol_fast", uopts.rtol_fast);
  find_arg(args, "--atol_fast", uopts.atol_fast);
  find_arg(args, "--fixed_h", uopts.fixed_h);
  find_arg(args, "--fixed_h_fast", uopts.fixed_h_fast);
  find_arg(args, "--predictor", uopts.predictor);
  find_arg(args, "--controller", uopts.controller);
  find_arg(args, "--lssetupfreq", uopts.ls_setup_freq);
  find_arg(args, "--predictor_fast", uopts.predictor_fast);
  find_arg(args, "--controller_fast", uopts.controller_fast);
  find_arg(args, "--lssetupfreq_fast", uopts.ls_setup_freq_fast);
  find_arg(args, "--maxsteps", uopts.maxsteps);
  find_arg(args, "--etamx1_fast", uopts.etamx1_fast);
  find_arg(args, "--linear", uopts.linear);
  find_arg(args, "--save_hinit", uopts.save_hinit);
  find_arg(args, "--save_hcur", uopts.save_hcur);
  find_arg(args, "--hcur_factor", uopts.hcur_factor);
  find_arg(args, "--output", uopts.output);
  find_arg(args, "--nout", uopts.nout);

  // Recompute mesh spacing and total number of nodes
  udata.dx  = (udata.xu - udata.xl) / (udata.nx - 1);
  udata.neq = NSPECIES * udata.nx;

  // Create workspace
  switch(uopts.integrator)
  {
  case(0):
    // Create workspace vector
    udata.temp_v = N_VNew_Serial(udata.neq, ctx);
    if (check_ptr(udata.temp_v, "N_VNew_Serial")) return 1;
    N_VConst(ZERO, udata.temp_v);
    break;
  case(1):
    // Create workspace vector
    udata.temp_v = N_VNew_Serial(udata.neq, ctx);
    if (check_ptr(udata.temp_v, "N_VNew_Serial")) return 1;
    N_VConst(ZERO, udata.temp_v);
    // Create workspace matrix
    udata.temp_J = SUNBandMatrix(udata.neq, 3, 3, ctx);
    if (check_ptr(udata.temp_J, "SUNBandMatrix")) return 1;
    SUNMatZero(udata.temp_J);
    break;
  case(2):
    break;
  case(3):
    break;
  default:
    cerr << "Invalid integrator option" << endl;
    return 1;
  }

  // Input checks
  if (!udata.diffusion && !udata.advection)
  {
    cerr << "ERROR: Invalid problem configuration" << endl;
    return -1;
  }

  if (udata.diffusion && udata.advection)
  {
    if (udata.splitting < 0 || udata.splitting > 7)
    {
      cerr << "ERROR: Invalid splitting option" << endl;
      return -1;
    }
  }
  else
  {
    if (udata.splitting < 0 || udata.splitting > 4)
    {
      cerr << "ERROR: Invalid splitting option" << endl;
      return -1;
    }
  }

  if (uopts.integrator < 0 || uopts.integrator > 3)
  {
    cerr << "ERROR: Invalid integrator option" << endl;
    return -1;
  }

  // ERKStep setup requires splitting 0 (fully explicit)
  if (uopts.integrator == 0)
    udata.splitting = 0;

  // MRIStep + ARKStep requires splitting 0 or 1 (explicit v implicit reactions)
  if (uopts.integrator == 2 && udata.splitting)
    udata.splitting = 1;

  // MRIStep + CVODE requires splitting 1 (implicit reactions)
  if (uopts.integrator == 3)
    udata.splitting = 1;

  // MRIStep + CVODE requires adaptive fast step sizes
  if (uopts.integrator == 3)
    uopts.fixed_h_fast = ZERO;

  return 0;
}


// Print user data
int PrintSetup(UserData &udata, UserOptions &uopts)
{
  cout << endl;
  cout << "Problem parameters and options:"     << endl;
  cout << " --------------------------------- " << endl;
  cout << "  c                = " << udata.c << endl;
  cout << "  d                = " << udata.d << endl;
  cout << "  A                = " << udata.A << endl;
  cout << "  B                = " << udata.B << endl;
  cout << " --------------------------------- " << endl;
  cout << "  tf               = " << udata.tf << endl;
  cout << "  xl               = " << udata.xl << endl;
  cout << "  xu               = " << udata.xu << endl;
  cout << "  nx               = " << udata.nx << endl;
  cout << "  dx               = " << udata.dx << endl;
  cout << " --------------------------------- " << endl;

  if (uopts.integrator == 0)
  {
    cout << "  integrator       = ERK" << endl;
    if (udata.advection)
      cout << "  advection        = Explicit" << endl;
    else
      cout << "  advection        = OFF" << endl;
    if (udata.diffusion)
      cout << "  diffusion        = Explicit" << endl;
    else
      cout << "  diffusion        = OFF" << endl;
    cout << "  reaction         = Explicit" << endl;
  }
  else if (uopts.integrator == 1)
  {
    cout << "  integrator       = ARK" << endl;
    // advection-diffusion-reaction
    if (udata.diffusion && udata.advection)
    {
      switch(udata.splitting)
      {
      case(0):
        // ERK -- fully explicit
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case(1):
        // IMEX -- explicit advection-diffusion, implicit reaction
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case(2):
        // IMEX -- explicit advection-reaction, implicit diffusion
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case(3):
        // IMEX -- explicit advection, implicit diffusion-reaction
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case(4):
        // IMEX -- explicit diffusion-reaction, implicit advection
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case(5):
        // IMEX -- explicit diffusion, implicit advection-reaction
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case(6):
        // IMEX -- explicit reaction, implicit advection-diffusion
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case(7):
        // DIRK -- fully implicit
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      default:
        cerr << "ERROR: Invalid splitting option" << endl;
        return -1;
        break;
      }
    }
    // advection-reaction
    else if (!udata.diffusion && udata.advection)
    {
      switch(udata.splitting)
      {
      case(0):
        // ERK -- fully explicit
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = OFF" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case(1):
        // IMEX -- explicit advection, implicit reaction
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = OFF" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case(2):
        // IMEX -- explicit reaction, implicit advection
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = OFF" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case(3):
        // DIRK -- fully implicit
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = OFF" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      default:
        cerr << "ERROR: Invalid splitting option" << endl;
        return -1;
        break;
      }
    }
    // diffusion-reaction
    else if (udata.diffusion && !udata.advection)
    {
      switch(udata.splitting)
      {
      case(0):
        // ERK -- fully explicit
        cout << "  advection        = OFF" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case(1):
        // IMEX -- explicit diffusion, implicit reaction
        cout << "  advection        = OFF" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case(2):
        // IMEX -- explicit reaction, implicit diffusion
        cout << "  advection        = OFF" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case(4):
        // DIRK -- fully implicit
        cout << "  advection        = OFF" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      default:
        cerr << "ERROR: Invalid splitting option" << endl;
        return -1;
        break;
      }
    }
    else
    {
      cerr << "ERROR: Invalid problem configuration" << endl;
      return -1;
    }
  }
  else if (uopts.integrator == 2)
  {
    cout << "  integrator       = MRI + ARK" << endl;
    if (udata.diffusion && udata.advection)
    {
      // IMEX slow -- advection-diffusion
      cout << "  advection        = Slow-Explicit" << endl;
      cout << "  diffusion        = Slow-Implicit" << endl;
      if (udata.splitting)
        cout << "  reaction         = Fast-Implicit" << endl;
      else
        cout << "  reaction         = Fast-Explicit" << endl;
    }
    else if (!udata.diffusion && udata.advection)
    {
      // Explicit slow -- advection
      cout << "  advection        = Slow-Explicit" << endl;
      cout << "  diffusion        = OFF" << endl;
      if (udata.splitting)
        cout << "  reaction         = Fast-Implicit" << endl;
      else
        cout << "  reaction         = Fast-Explicit" << endl;
    }
    else if (udata.diffusion && !udata.advection)
    {
      cout << "  advection        = OFF" << endl;
      cout << "  diffusion        = Implicit" << endl;
      if (udata.splitting)
        cout << "  reaction         = Fast-Implicit" << endl;
      else
        cout << "  reaction         = Fast-Explicit" << endl;
    }
    else
    {
      // No slow time scale
      cerr << "ERROR: Invalid problem configuration" << endl;
      return -1;
    }
  }
  else if (uopts.integrator == 3)
  {
    cout << "  integrator       = MRI + CVODE" << endl;
    // Slow time scale
    if (udata.diffusion && udata.advection)
    {
      // IMEX slow -- advection-diffusion
      cout << "  advection        = Slow-Explicit" << endl;
      cout << "  diffusion        = Slow-Implicit" << endl;
      cout << "  reaction         = Fast-Implicit" << endl;
    }
    else if (!udata.diffusion && udata.advection)
    {
      // Explicit slow -- advection
      cout << "  advection        = Slow-Explicit" << endl;
      cout << "  diffusion        = OFF" << endl;
      cout << "  reaction         = Fast-Implicit" << endl;
    }
    else if (udata.diffusion && !udata.advection)
    {
      // Implicit slow -- diffusion
      cout << "  advection        = OFF" << endl;
      cout << "  diffusion        = Slow-Implicit" << endl;
      cout << "  reaction         = Fast-Implicit" << endl;
    }
    else
    {
      // No slow time scale
      cerr << "ERROR: Invalid problem configuration" << endl;
      return -1;
    }
  }
  else
  {
    cerr << "ERROR: Invalid integrator option" << endl;
    return -1;
  }
  if (uopts.ark_dirk)
    cout << "  order (ark_dirk) = " << uopts.order << endl;
  else
    cout << "  order            = " << uopts.order << endl;
  cout << "  rtol             = " << uopts.rtol << endl;
  cout << "  atol             = " << uopts.atol << endl;
  cout << "  fixed h          = " << uopts.fixed_h << endl;
  if (uopts.controller <= 0)
    cout << "  controller       = PID" << endl;
  else if (uopts.controller == 1)
    cout << "  controller       = PI" << endl;
  else if (uopts.controller == 2)
    cout << "  controller       = I" << endl;
  else if (uopts.controller == 3)
    cout << "  controller       = explicit Gustafsson" << endl;
  else if (uopts.controller == 4)
    cout << "  controller       = implicit Gustafsson" << endl;
  else if (uopts.controller == 5)
    cout << "  controller       = IMEX Gustafsson" << endl;
  else
    cout << "  controller       = " << uopts.controller << endl;
  if (uopts.integrator > 0)
  {
    if (uopts.predictor == 0)
      cout << "  predictor        = trivial" << endl;
    else if (uopts.predictor == 1)
      cout << "  predictor        = max order" << endl;
    else if (uopts.predictor == 2)
      cout << "  predictor        = variable order" << endl;
    else if (uopts.predictor == 3)
      cout << "  predictor        = cutoff order" << endl;
    else
      cout << "  predictor        = " << uopts.predictor << endl;
    cout << "  ls setup freq    = " << uopts.ls_setup_freq << endl;
    cout << "  linear           = " << uopts.linear << endl;
  }
  if (uopts.integrator == 2)
  {
    cout << " --------------------------------- " << endl;
    cout << "  fast integrator  = ARK" << endl;
    cout << "  rtol             = " << uopts.rtol_fast << endl;
    cout << "  atol             = " << uopts.atol_fast << endl;
    cout << "  order            = " << uopts.order_fast << endl;
    cout << "  fixed h          = " << uopts.fixed_h_fast << endl;
    if (uopts.controller_fast <= 0)
      cout << "  controller       = PID" << endl;
    else if (uopts.controller_fast == 1)
      cout << "  controller       = PI" << endl;
    else if (uopts.controller_fast == 2)
      cout << "  controller       = I" << endl;
    else if (uopts.controller_fast == 3)
      cout << "  controller       = explicit Gustafsson" << endl;
    else if (uopts.controller_fast == 4)
      cout << "  controller       = implicit Gustafsson" << endl;
    else if (uopts.controller_fast == 5)
      cout << "  controller       = IMEX Gustafsson" << endl;
    else
      cout << "  controller       = " << uopts.controller_fast << endl;
    if (uopts.predictor_fast == 0)
      cout << "  predictor        = trivial" << endl;
    else if (uopts.predictor_fast == 1)
      cout << "  predictor        = max order" << endl;
    else if (uopts.predictor_fast == 2)
      cout << "  predictor        = variable order" << endl;
    else if (uopts.predictor_fast == 3)
      cout << "  predictor        = cutoff order" << endl;
    else
      cout << "  predictor        = " << uopts.predictor_fast << endl;
    cout << "  ls setup freq    = " << uopts.ls_setup_freq_fast << endl;
  }
  else if (uopts.integrator == 3)
  {
    cout << " --------------------------------- " << endl;
    cout << "  fast integrator  = CVODE" << endl;
    cout << "  rtol             = " << uopts.rtol_fast << endl;
    cout << "  atol             = " << uopts.atol_fast << endl;
    cout << "  ls setup freq    = " << uopts.ls_setup_freq << endl;
    cout << "  etamx first step = " << uopts.etamx1_fast << endl;
    cout << "  reuse initial h  = " << uopts.save_hinit << endl;
    cout << "  reuse current h  = " << uopts.save_hcur << endl;
    cout << "  current h factor = " << uopts.hcur_factor << endl;
  }
  cout << " --------------------------------- " << endl;
  cout << "  output          = " << uopts.output << endl;
  cout << " --------------------------------- " << endl;
  cout << endl;

  return 0;
}


// Initialize output
int OpenOutput(UserData &udata, UserOptions &uopts)
{
  // Header for status output
  if (uopts.output)
  {
    cout << scientific;
    cout << setprecision(numeric_limits<realtype>::digits10);
    cout << "          t           ";
    cout << "          ||y||_rms      " << endl;
    cout << " ---------------------";
    cout << "-------------------------" << endl;
  }

  // Open output stream and output problem information
  if (uopts.output >= 2)
  {
    // Open output stream
    stringstream fname;
    fname << "advection_diffusion_reaction.out";
    uopts.uout.open(fname.str());

    uopts.uout << scientific;
    uopts.uout << setprecision(numeric_limits<realtype>::digits10);
    uopts.uout <<  "# title Advection-Diffusion-Reaction (Brusselator)" << endl;
    uopts.uout <<  "# nvar 3"     << endl;
    uopts.uout <<  "# vars u v w" << endl;
    uopts.uout <<  "# nt "        << uopts.nout + 1 << endl;
    uopts.uout <<  "# nx "        << udata.nx << endl;
    uopts.uout <<  "# xl "        << udata.xl << endl;
    uopts.uout <<  "# xu "        << udata.xu << endl;
  }

  return 0;
}


// Write output
int WriteOutput(realtype t, N_Vector y, UserData &udata, UserOptions &uopts)
{
  if (uopts.output)
  {
    // Compute rms norm of the state
    realtype urms = sqrt(N_VDotProd(y, y) / udata.nx);
    cout << setw(22) << t << setw(25) << urms << endl;

    // Write solution to disk
    if (uopts.output >= 2)
    {
      realtype* ydata = N_VGetArrayPointer(y);
      if (check_ptr(ydata, "N_VGetArrayPointer")) return -1;

      uopts.uout << t;
      for (sunindextype i = 0; i < udata.nx; i++)
      {
        uopts.uout << setw(WIDTH) << ydata[UIDX(i)];
        uopts.uout << setw(WIDTH) << ydata[VIDX(i)];
        uopts.uout << setw(WIDTH) << ydata[WIDX(i)];
      }
      uopts.uout << endl;
    }
  }

  return 0;
}


// Finalize output
int CloseOutput(UserOptions &uopts)
{
  // Footer for status output
  if (uopts.output)
  {
    cout << " ---------------------";
    cout << "-------------------------" << endl;
    cout << endl;
  }

  // Close output streams
  if (uopts.output >= 2)
  {
    uopts.uout.close();
  }

  return 0;
}


//---- end of file ----
