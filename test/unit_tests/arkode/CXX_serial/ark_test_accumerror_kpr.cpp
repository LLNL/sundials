/* ----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * Routine to test the accumulated temporal error estimation
 * approaches from ARKStep and ERKStep.  Uses a nonlinear
 * Kvaerno-Prothero-Robinson ODE test problem with analytical
 * solution,
 *
 *    [u]' =  [ G  e ] [(u^2-r-1)/(2u)] +  [ r'(t)/(2u) ]
 *    [v]     [ e -1 ] [(v^2-s-2)/(2v)]    [ s'(t)/(2v) ]
 *
 * where r(t) = 0.5*cos(t),  s(t) = sin(t),  -3 < t < 7. This
 * problem has analytical solution given by
 *    u(t) = sqrt(1+r(t)),  v(t) = sqrt(2+s(t)).
 *
 * We use the parameters: e = 0.5 and G = -100 [default]
 *
 * The stiffness of the problem is essentially determined
 * by G, for |G| > 50 it is 'stiff' and ideally suited to an
 * implicit method.
 *
 * We use either the ARKStep/DIRK/Newton/Dense solver (0) or
 * ERKStep (1).  Either defaults to using a 4th-order method.
 *
 * By default, all runs use temporal adaptivity; however, if the
 * requested 'ord' command-line input is negative, we run with
 * order |ord|, using fixed step sizes.
 *
 * The program should be run with arguments in the following order:
 *   $ a.out method ord G
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out method ord G
 *   $ a.out method ord
 *   $ a.out method
 *   $ a.out
 * are acceptable.  We require:
 *   * method = {0, 1}
 *   * G < 0.0
 *
 * For either temporally adaptive (ord >= 0) or fixed-step (ord < 0)
 * runs, we test a variety of tolerances/step sizes, and compare
 * the true error at the end of the run against the
 * integrator-reported accumulated error estimate.
 * ----------------------------------------------------------------*/

// Header files
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <vector>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_math.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)
#define PI4  RCONST(0.78539816339744830961566084581987572)

using namespace std;

// User data structure
struct UserData
{
  realtype G;
  realtype e;
};

// User-supplied functions called by the solver
static int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jn(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private utility functions
static int adaptive_run(void *arkode_mem, N_Vector y, realtype T0,
                        realtype Tf, int rk_type, UserData &udata);
static int fixed_run(void *arkode_mem, N_Vector y, realtype T0, realtype Tf,
                     int rk_type, UserData &udata);
static realtype r(realtype t);
static realtype s(realtype t);
static realtype rdot(realtype t);
static realtype sdot(realtype t);
static realtype utrue(realtype t);
static realtype vtrue(realtype t);
static int Ytrue(realtype t, N_Vector y);
static int check_retval(void *returnvalue, const char *funcname, int opt);

// Main Program
int main(int argc, char *argv[])
{
  // general problem parameters
  realtype T0 = RCONST(-3.0);       // initial time
  realtype Tf = RCONST(7.0);        // final time
  sunindextype NEQ = 2;             // number of dependent vars.
  int rk_type = 0;                  // type of RK method [DIRK=0, ERK=1]
  int order = 4;                    // order of accuracy for RK method
  booleantype adaptive = SUNTRUE;   // adaptive run vs convergence order

  // general problem variables
  int retval;                    // reusable error-checking flag
  N_Vector y = NULL;             // empty vector for the computed solution
  void *arkode_mem = NULL;       // empty ARKODE memory structure
  SUNMatrix A = NULL;            // empty system matrix
  SUNLinearSolver LS = NULL;     // empty system linear solver object
  UserData udata;                // user-data structure
  udata.G = RCONST(-100.0);      // stiffness parameter
  udata.e = RCONST(0.5);         // coupling strength

  //
  // Initialization
  //

  // Retrieve the command-line options:  method ord G
  if (argc > 1)  rk_type = (int) atoi(argv[1]);
  if (argc > 2)  order = (int) atoi(argv[2]);
  if (argc > 3)  udata.G = (realtype) atof(argv[3]);

  // Check arguments for validity
  //   0 <= rk_type <= 1
  //   G < 0.0
  if ((rk_type < 0) || (rk_type > 1)) {
    cerr << "ERROR: RK type be an integer in {0,1} \n";
    return(-1);
  }
  if (udata.G >= ZERO) {
    cerr << "ERROR: G must be a negative real number\n";
    return(-1);
  }

  // Handle adaptive run vs order-of-convergence run
  if (order < 0) {
    adaptive = SUNFALSE;
    order = abs(order);
  }
  if (order == 0)  order = 4;

  // Initial problem output (and set implicit solver tolerances as needed)
  cout << "\nAccumulated error estimation test (Nonlinear Kvaerno-Prothero-Robinson problem):\n";
  cout << "    time domain:  (" << T0 << "," << Tf << "]\n";
  cout << "    G = " << udata.G << endl;
  cout << "    e = " << udata.e << endl;
  if (rk_type == 0) {
    cout << "    DIRK solver, order = " << order << endl;
  } else if (rk_type == 1) {
    cout << "    ERK solver, order = " << order << endl;
  }

  //
  // Problem Setup
  //

  // Create SUNDIALS context
  sundials::Context ctx;

  // Create and initialize serial vector for the solution
  y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return 1;
  retval = Ytrue(T0, y);
  if (check_retval(&retval, "Ytrue", 1)) return 1;

  // Initialize ARKStep or ERKStep.
  if (rk_type == 0) {    // DIRK method

    arkode_mem = ARKStepCreate(NULL, fn, T0, y, ctx);
    if (check_retval((void *) arkode_mem, "ARKStepCreate", 0)) return 1;

    // Initialize/attach linear solvers (if required)
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void *)A, "SUNDenseMatrix", 0)) return 1;
    LS = SUNLinSol_Dense(y, A, ctx);
    if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) return 1;
    retval = ARKStepSetLinearSolver(arkode_mem, LS, A);
    if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return(1);
    retval = ARKStepSetJacFn(arkode_mem, Jn);
    if (check_retval(&retval, "ARKStepSetJacFn", 1)) return 1;

    // Set desired solver order
    retval = ARKStepSetOrder(arkode_mem, order);
    if (check_retval(&retval, "ARKStepSetOrder", 1)) return 1;

    // Set the user data pointer
    retval = ARKStepSetUserData(arkode_mem, (void *) &udata);
    if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;

  } else {               // ERK method

    arkode_mem = ERKStepCreate(fn, T0, y, ctx);
    if (check_retval((void *) arkode_mem, "ERKStepCreate", 0)) return 1;

    // Set maximum stepsize for ERK run
    retval = ERKStepSetMaxStep(arkode_mem, ONE/abs(udata.G));
    if (check_retval(&retval, "ERKStepSetMaxStep", 1)) return(1);

    // Set desired solver order
    retval = ERKStepSetOrder(arkode_mem, order);
    if (check_retval(&retval, "ERKStepSetOrder", 1)) return 1;

    // Set the user data pointer
    retval = ERKStepSetUserData(arkode_mem, (void *) &udata);
    if (check_retval(&retval, "ERKStepSetUserData", 1)) return 1;

  }

  // Integrate ODE, based on run type
  if (adaptive) {
    retval = adaptive_run(arkode_mem, y, T0, Tf, rk_type, udata);
    if (check_retval(&retval, "adaptive_run", 1)) return 1;
  } else {
    retval = fixed_run(arkode_mem, y, T0, Tf, rk_type, udata);
    if (check_retval(&retval, "fixed_run", 1)) return 1;
  }

  // Clean up and return
  if (rk_type == 0) {                 // Free integrator memory
    ARKStepFree(&arkode_mem);
  } else {
    ERKStepFree(&arkode_mem);
  }
  if (LS != NULL) SUNLinSolFree(LS);  // free system linear solver
  if (A  != NULL) SUNMatDestroy(A);   // free system matrix
  N_VDestroy(y);                      // Free y vector
  return 0;
}

//------------------------------
// Functions called by the solver
//------------------------------

static int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData *udata = (UserData *) user_data;
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype tmp1, tmp2;

  // fill in the RHS function:
  //   [G  e]*[(-1+u^2-r(t))/(2*u)] + [rdot(t)/(2u)]
  //   [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2v)]
  tmp1 = (-ONE+u*u-r(t))/(TWO*u);
  tmp2 = (-TWO+v*v-s(t))/(TWO*v);
  NV_Ith_S(ydot,0) = udata->G*tmp1 + udata->e*tmp2 + rdot(t)/(TWO*u);
  NV_Ith_S(ydot,1) = udata->e*tmp1 - tmp2 + sdot(t)/(TWO*v);

  // Return with success
  return 0;
}

static int Jn(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData *udata = (UserData *) user_data;
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype t11, t12, t21, t22;

  // fill in the Jacobian:
  //   [G  e]*[1-(u^2-r(t)-1)/(2*u^2),  0] + [-r'(t)/(2*u^2),  0]
  //   [e -1] [0,  1-(v^2-s(t)-2)/(2*v^2)]   [0,  -s'(t)/(2*v^2)]
  t11 = ONE-(u*u-r(t)-ONE)/(TWO*u*u);
  t12 = ZERO;
  t21 = ZERO;
  t22 = ONE-(v*v-s(t)-TWO)/(TWO*v*v);
  SM_ELEMENT_D(J,0,0) = udata->G*t11 + udata->e*t21 - rdot(t)/(TWO*u*u);
  SM_ELEMENT_D(J,0,1) = udata->G*t12 + udata->e*t22;
  SM_ELEMENT_D(J,1,0) = udata->e*t11 - t21;
  SM_ELEMENT_D(J,1,1) = udata->e*t12 - t22 - sdot(t)/(TWO*v*v);

  // Return with success
  return 0;
}



//------------------------------
// Private helper functions
//------------------------------

static int adaptive_run(void *arkode_mem, N_Vector y, realtype T0,
                        realtype Tf, int rk_type, UserData &udata)
{
  // Reused variables
  int retval;
  long int nsteps;
  realtype dsm_est;
  realtype t = T0;

  // Set testing tolerances
  realtype abstol = RCONST(1.e-12);
  vector<realtype> rtols = {RCONST(1.e-2), RCONST(1.e-4), RCONST(1.e-6)};
  vector<int> accum_types = {0, 1};

  // Loop over tolerances
  cout << "\n Adaptive-step runs:\n";
  for (size_t irtol=0; irtol<rtols.size(); irtol++) {

    cout << "   Rtol: " << rtols[irtol] << endl;

    // Loop over accumulation types
    for (size_t iaccum=0; iaccum<accum_types.size(); iaccum++) {

      // Reset integrator for this run, and evolve to Tf
      t = T0;
      retval = Ytrue(T0, y);
      if (check_retval(&retval, "Ytrue", 1)) return 1;
      if (rk_type == 0) {  // DIRK
        retval = ARKStepReInit(arkode_mem, NULL, fn, T0, y);
        if (check_retval(&retval, "ARKStepReInit", 1)) return 1;
        retval = ARKStepSetAccumulatedErrorType(arkode_mem, accum_types[iaccum]);
        if (check_retval(&retval, "ARKStepSetAccumulatedErrorType", 1)) return 1;
        retval = ARKStepResetAccumulatedError(arkode_mem);
        if (check_retval(&retval, "ARKStepResetAccumulatedError", 1)) return 1;
        retval = ARKStepSStolerances(arkode_mem, rtols[irtol], abstol);
        if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
        retval = ARKStepSetStopTime(arkode_mem, Tf);
        if (check_retval(&retval, "ARKStepSetStopTime", 1)) return 1;
        retval = ARKStepSetMaxNumSteps(arkode_mem, 1000000);
        if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return(1);
        retval = ARKStepEvolve(arkode_mem, Tf, y, &t, ARK_NORMAL);
        if (check_retval(&retval, "ARKStepEvolve", 1)) break;
        retval = ARKStepGetAccumulatedError(arkode_mem, &dsm_est);
        if (check_retval(&retval, "ARKStepGetAccumulatedError", 1)) break;
        retval = ARKStepGetNumSteps(arkode_mem, &nsteps);
        if (check_retval(&retval, "ARKStepGetNumSteps", 1)) break;
      } else {            // ERK
        retval = ERKStepReInit(arkode_mem, fn, T0, y);
        if (check_retval(&retval, "ERKStepReInit", 1)) return 1;
        retval = ERKStepSetAccumulatedErrorType(arkode_mem, accum_types[iaccum]);
        if (check_retval(&retval, "ERKStepSetAccumulatedErrorType", 1)) return 1;
        retval = ERKStepResetAccumulatedError(arkode_mem);
        if (check_retval(&retval, "ERKStepResetAccumulatedError", 1)) return 1;
        retval = ERKStepSStolerances(arkode_mem, rtols[irtol], abstol);
        if (check_retval(&retval, "ERKStepSStolerances", 1)) return 1;
        retval = ERKStepSetStopTime(arkode_mem, Tf);
        if (check_retval(&retval, "ERKStepSetStopTime", 1)) return 1;
        retval = ERKStepSetMaxNumSteps(arkode_mem, 1000000);
        if (check_retval(&retval, "ERKStepSetMaxNumSteps", 1)) return(1);
        retval = ERKStepEvolve(arkode_mem, Tf, y, &t, ARK_NORMAL);
        if (check_retval(&retval, "ERKStepEvolve", 1)) break;
        retval = ERKStepGetAccumulatedError(arkode_mem, &dsm_est);
        if (check_retval(&retval, "ERKStepGetAccumulatedError", 1)) break;
        retval = ERKStepGetNumSteps(arkode_mem, &nsteps);
        if (check_retval(&retval, "ERKStepGetNumSteps", 1)) break;
      }

      // Compute/print solution error
      realtype udsm = abs(NV_Ith_S(y,0)-utrue(t))/(abstol + rtols[irtol]*abs(utrue(t)));
      realtype vdsm = abs(NV_Ith_S(y,1)-vtrue(t))/(abstol + rtols[irtol]*abs(vtrue(t)));
      realtype dsm = rtols[irtol]*sqrt(0.5*(udsm*udsm + vdsm*vdsm));
      cout << "     acc type = " << accum_types[iaccum]
           << ",  dsm = " << dsm
           << ",  dsm_est = " << dsm_est
           << ",  dsm/dsm_est = " << dsm/dsm_est
           << ",  nsteps = " << nsteps
           << endl;
    }
  }
  cout << "   ------------------------------------------------------\n";

  return(0);
}

static int fixed_run(void *arkode_mem, N_Vector y, realtype T0, realtype Tf,
                     int rk_type, UserData &udata)
{
  // local variables
  int retval;
  long int nsteps, nsteps2;
  realtype dsm_est;
  realtype t, t2;
  N_Vector y2 = N_VClone(y);
  N_Vector ewt = N_VClone(y);

  // Set array of fixed step sizes to use, storage for corresponding errors/orders
  realtype hmax = (Tf - T0)/1000;
  if (rk_type == 1) hmax = min(hmax, ONE/abs(udata.G));
  vector<realtype> hvals = {hmax, hmax/4, hmax/16, hmax/64};
  vector<int> accum_types = {0, 1};

  // Loop over step sizes
  cout << "\n Fixed-step runs:\n";
  for (size_t ih=0; ih<hvals.size(); ih++) {

    cout << "   H: " << hvals[ih] << endl;

    // Loop over built-in accumulation types
    for (size_t iaccum=0; iaccum<accum_types.size(); iaccum++) {

      // Reset integrator for this run, and evolve to Tf
      t = T0;
      retval = Ytrue(T0, y);
      if (check_retval(&retval, "Ytrue", 1)) return 1;
      if (rk_type == 0) {  // DIRK
        retval = ARKStepReInit(arkode_mem, NULL, fn, T0, y);
        if (check_retval(&retval, "ARKStepReInit", 1)) return 1;
        retval = ARKStepSetAccumulatedErrorType(arkode_mem, accum_types[iaccum]);
        if (check_retval(&retval, "ARKStepSetAccumulatedErrorType", 1)) return 1;
        retval = ARKStepResetAccumulatedError(arkode_mem);
        if (check_retval(&retval, "ARKStepResetAccumulatedError", 1)) return 1;
        retval = ARKStepSetFixedStep(arkode_mem, hvals[ih]);
        if (check_retval(&retval, "ARKStepSetFixedStep", 1)) return 1;
        retval = ARKStepSetMaxNumSteps(arkode_mem, 1000000);
        if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return(1);
        retval = ARKStepSetStopTime(arkode_mem, Tf);
        if (check_retval(&retval, "ARKStepSetStopTime", 1)) return 1;
        retval = ARKStepSStolerances(arkode_mem, RCONST(1.e-9), RCONST(1.e-12));
        if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
        retval = ARKStepSetJacEvalFrequency(arkode_mem, 1);
        if (check_retval(&retval, "ARKStepSetJacEvalFrequency", 1)) return 1;
        retval = ARKStepSetLSetupFrequency(arkode_mem, 1);
        if (check_retval(&retval, "ARKStepSetLSetupFrequency", 1)) return 1;
        retval = ARKStepSetMaxNonlinIters(arkode_mem, 20);
        if (check_retval(&retval, "ARKStepSetMaxNonlinIters", 1)) return 1;
        retval = ARKStepEvolve(arkode_mem, Tf, y, &t, ARK_NORMAL);
        if (check_retval(&retval, "ARKStepEvolve", 1)) break;
        retval = ARKStepGetAccumulatedError(arkode_mem, &dsm_est);
        if (check_retval(&retval, "ARKStepGetAccumulatedError", 1)) break;
        retval = ARKStepGetNumSteps(arkode_mem, &nsteps);
        if (check_retval(&retval, "ARKStepGetNumSteps", 1)) break;
      } else {            // ERK
        retval = ERKStepReInit(arkode_mem, fn, T0, y);
        if (check_retval(&retval, "ERKStepReInit", 1)) return 1;
        retval = ERKStepSetAccumulatedErrorType(arkode_mem, accum_types[iaccum]);
        if (check_retval(&retval, "ERKStepSetAccumulatedErrorType", 1)) return 1;
        retval = ERKStepResetAccumulatedError(arkode_mem);
        if (check_retval(&retval, "ERKStepResetAccumulatedError", 1)) return 1;
        retval = ERKStepSetFixedStep(arkode_mem, hvals[ih]);
        if (check_retval(&retval, "ERKStepSetFixedStep", 1)) return 1;
        retval = ERKStepSetMaxNumSteps(arkode_mem, 1000000);
        if (check_retval(&retval, "ERKStepSetMaxNumSteps", 1)) return(1);
        retval = ERKStepSStolerances(arkode_mem, RCONST(1.e-9), RCONST(1.e-12));
        if (check_retval(&retval, "ERKStepSStolerances", 1)) return 1;
        retval = ERKStepSetStopTime(arkode_mem, Tf);
        if (check_retval(&retval, "ERKStepSetStopTime", 1)) return 1;
        retval = ERKStepEvolve(arkode_mem, Tf, y, &t, ARK_NORMAL);
        if (check_retval(&retval, "ERKStepEvolve", 1)) break;
        retval = ERKStepGetAccumulatedError(arkode_mem, &dsm_est);
        if (check_retval(&retval, "ERKStepGetAccumulatedError", 1)) break;
        retval = ERKStepGetNumSteps(arkode_mem, &nsteps);
        if (check_retval(&retval, "ERKStepGetNumSteps", 1)) break;
      }

      // Compute/print solution error
      realtype udsm = abs(NV_Ith_S(y,0)-utrue(t))/((1.e-12) + (1.e-9)*abs(utrue(t)));
      realtype vdsm = abs(NV_Ith_S(y,1)-vtrue(t))/((1.e-12) + (1.e-9)*abs(vtrue(t)));
      realtype dsm = (1.e-9)*sqrt(0.5*(udsm*udsm + vdsm*vdsm));
      cout << "     acc type = " << accum_types[iaccum]
           << ",  dsm = " << dsm
           << ",  dsm_est = " << dsm_est
           << ",  dsm/dsm_est = " << dsm/dsm_est
           << ",  nsteps = " << nsteps
           << endl;
    }

    // Test double fixed-step run error estimator
    t = T0;
    retval = Ytrue(T0, y);
    if (check_retval(&retval, "Ytrue", 1)) return 1;
    t2 = T0;
    retval = Ytrue(T0, y2);
    if (check_retval(&retval, "Ytrue", 1)) return 1;
    if (rk_type == 0) {  // DIRK
      retval = ARKStepReInit(arkode_mem, NULL, fn, T0, y);
      if (check_retval(&retval, "ARKStepReInit", 1)) return 1;
      retval = ARKStepSetAccumulatedErrorType(arkode_mem, 0);
      if (check_retval(&retval, "ARKStepSetAccumulatedErrorType", 1)) return 1;
      retval = ARKStepSetFixedStep(arkode_mem, hvals[ih]);
      if (check_retval(&retval, "ARKStepSetFixedStep", 1)) return 1;
      retval = ARKStepSetMaxNumSteps(arkode_mem, 1000000);
      if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return(1);
      retval = ARKStepSetStopTime(arkode_mem, Tf);
      if (check_retval(&retval, "ARKStepSetStopTime", 1)) return 1;
      retval = ARKStepSStolerances(arkode_mem, RCONST(1.e-9), RCONST(1.e-12));
      if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
      retval = ARKStepSetJacEvalFrequency(arkode_mem, 1);
      if (check_retval(&retval, "ARKStepSetJacEvalFrequency", 1)) return 1;
      retval = ARKStepSetLSetupFrequency(arkode_mem, 1);
      if (check_retval(&retval, "ARKStepSetLSetupFrequency", 1)) return 1;
      retval = ARKStepSetMaxNonlinIters(arkode_mem, 20);
      if (check_retval(&retval, "ARKStepSetMaxNonlinIters", 1)) return 1;
      retval = ARKStepEvolve(arkode_mem, Tf, y, &t, ARK_NORMAL);
      if (check_retval(&retval, "ARKStepEvolve", 1)) break;
      retval = ARKStepGetNumSteps(arkode_mem, &nsteps);
      if (check_retval(&retval, "ARKStepGetNumSteps", 1)) break;

      retval = ARKStepReInit(arkode_mem, NULL, fn, T0, y2);
      if (check_retval(&retval, "ARKStepReInit", 1)) return 1;
      retval = ARKStepSetFixedStep(arkode_mem, 2.0*hvals[ih]);
      if (check_retval(&retval, "ARKStepSetFixedStep", 1)) return 1;
      retval = ARKStepSetStopTime(arkode_mem, Tf);
      if (check_retval(&retval, "ARKStepSetStopTime", 1)) return 1;
      retval = ARKStepEvolve(arkode_mem, Tf, y2, &t2, ARK_NORMAL);
      if (check_retval(&retval, "ARKStepEvolve", 1)) break;
      retval = ARKStepGetNumSteps(arkode_mem, &nsteps2);
      if (check_retval(&retval, "ARKStepGetNumSteps", 1)) break;

      retval = ARKStepGetErrWeights(arkode_mem, ewt);
      if (check_retval(&retval, "ARKStepGetErrWeights", 1)) break;
      N_VLinearSum(1.0, y2, -1.0, y, y2);
      dsm_est = (1.e-9)*N_VWrmsNorm(y2, ewt);
      nsteps += nsteps2;
    } else {            // ERK
      retval = ERKStepReInit(arkode_mem, fn, T0, y);
      if (check_retval(&retval, "ERKStepReInit", 1)) return 1;
      retval = ERKStepSetAccumulatedErrorType(arkode_mem, 0);
      if (check_retval(&retval, "ERKStepSetAccumulatedErrorType", 1)) return 1;
      retval = ERKStepSetFixedStep(arkode_mem, hvals[ih]);
      if (check_retval(&retval, "ERKStepSetFixedStep", 1)) return 1;
      retval = ERKStepSetMaxNumSteps(arkode_mem, 1000000);
      if (check_retval(&retval, "ERKStepSetMaxNumSteps", 1)) return(1);
      retval = ERKStepSetStopTime(arkode_mem, Tf);
      if (check_retval(&retval, "ERKStepSetStopTime", 1)) return 1;
      retval = ERKStepSStolerances(arkode_mem, RCONST(1.e-9), RCONST(1.e-12));
      if (check_retval(&retval, "ERKStepSStolerances", 1)) return 1;
      retval = ERKStepEvolve(arkode_mem, Tf, y, &t, ARK_NORMAL);
      if (check_retval(&retval, "ERKStepEvolve", 1)) break;
      retval = ERKStepGetNumSteps(arkode_mem, &nsteps);
      if (check_retval(&retval, "ERKStepGetNumSteps", 1)) break;

      retval = ERKStepReInit(arkode_mem, fn, T0, y2);
      if (check_retval(&retval, "ERKStepReInit", 1)) return 1;
      retval = ERKStepSetFixedStep(arkode_mem, 2.0*hvals[ih]);
      if (check_retval(&retval, "ERKStepSetFixedStep", 1)) return 1;
      retval = ERKStepSetStopTime(arkode_mem, Tf);
      if (check_retval(&retval, "ERKStepSetStopTime", 1)) return 1;
      retval = ERKStepEvolve(arkode_mem, Tf, y2, &t2, ARK_NORMAL);
      if (check_retval(&retval, "ERKStepEvolve", 1)) break;
      retval = ERKStepGetNumSteps(arkode_mem, &nsteps2);
      if (check_retval(&retval, "ERKStepGetNumSteps", 1)) break;

      retval = ERKStepGetErrWeights(arkode_mem, ewt);
      if (check_retval(&retval, "ERKStepGetErrWeights", 1)) break;
      N_VLinearSum(1.0, y2, -1.0, y, y2);
      dsm_est = (1.e-9)*N_VWrmsNorm(y2, ewt);
      nsteps += nsteps2;
    }
    realtype udsm = abs(NV_Ith_S(y,0)-utrue(t))/((1.e-12) + (1.e-9)*abs(utrue(t)));
    realtype vdsm = abs(NV_Ith_S(y,1)-vtrue(t))/((1.e-12) + (1.e-9)*abs(vtrue(t)));
    realtype dsm = (1.e-9)*sqrt(0.5*(udsm*udsm + vdsm*vdsm));
    cout << "     acc type = " << 2
         << ",  dsm = " << dsm
         << ",  dsm_est = " << dsm_est
         << ",  dsm/dsm_est = " << dsm/dsm_est
         << ",  nsteps = " << nsteps
         << endl;


    // Test max-step run error estimator


  }
  cout << "   ------------------------------------------------------\n";

  N_VDestroy(y2);
  N_VDestroy(ewt);
  return(0);
}

static realtype r(realtype t)
{
  return( RCONST(0.5)*cos(t) );
}
static realtype s(realtype t)
{
  return( sin(t) );
}
static realtype rdot(realtype t)
{
  return( -RCONST(0.5)*sin(t) );
}
static realtype sdot(realtype t)
{
  return( cos(t) );
}
static realtype utrue(realtype t)
{
  return( SUNRsqrt(ONE+r(t)) );
}
static realtype vtrue(realtype t)
{
  return( SUNRsqrt(TWO+s(t)) );
}
static int Ytrue(realtype t, N_Vector y)
{
  NV_Ith_S(y,0) = utrue(t);
  NV_Ith_S(y,1) = vtrue(t);
  return(0);
}


/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval < 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  // Check if retval < 0
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1; }}

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


//---- end of file ----
