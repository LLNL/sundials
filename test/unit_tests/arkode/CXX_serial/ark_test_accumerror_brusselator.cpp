/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Routine to test the accumulated temporal error estimation
 * approaches from ARKStep.  Uses the "Brusselator" test problem,
 * an ODE system with 3 components, Y = [u,v,w], with IVP,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions
 * Y0 = [u0,v0,w0].
 *
 * We may run the problem in 3 different testing scenarios:
 *
 * Test 1:  u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5,  ep=1.0e-5
 *    Here, all three components exhibit a rapid transient change
 *    during the first 0.2 time units, followed by a slow and
 *    smooth evolution.
 *
 * Test 2:  u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5,  ep=5.0e-6
 *    Here, w experiences a fast initial transient, jumping 0.5
 *    within a few steps.  All values proceed smoothly until
 *    around t=6.5, when both u and v undergo a sharp transition,
 *    with u increaseing from around 0.5 to 5 and v decreasing
 *    from around 6 to 1 in less than 0.5 time units.  After this
 *    transition, both u and v continue to evolve somewhat
 *    rapidly for another 1.4 time units, and finish off smoothly.
 *
 * Test 3:  u0=3,  v0=3,  w0=3.5,  a=0.5,  b=3,  ep=5.0e-4
 *    Here, all components undergo very rapid initial transients
 *    during the first 0.3 time units, and all then proceed very
 *    smoothly for the remainder of the simulation.
 *
 * This program solves the problem with a DIRK/Newton/Dense
 * solver, and defaults to using a 4th-order method.
 *
 * By default, all runs use temporal adaptivity; however, if the
 * requested 'ord' command-line input is negative, we run with
 * order |ord|, using fixed step sizes.
 *
 * The program should be run with arguments in the following order:
 *   $ a.out ord test
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out ord test
 *   $ a.out ord
 *   $ a.out
 * are acceptable.  We require:
 *   * test = {1, 2, 3}
 *
 * For either temporally adaptive (ord >= 0) or fixed-step (ord < 0)
 * runs, we test a variety of tolerances/step sizes, and compare
 * the error at the end of the run (computed via a reference
 * solution) against the integrator-reported accumulated error
 * estimate.
 *-----------------------------------------------------------------*/

// Header files
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <vector>
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

using namespace std;

/* User-supplied Functions Called by the Solver */
static int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private utility functions
static int adaptive_run(void *arkode_mem, N_Vector y, realtype T0, realtype Tf, N_Vector yref);
static int fixed_run(void *arkode_mem, N_Vector y, realtype T0, realtype Tf, N_Vector yref);
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Main Program */
int main(int argc, char *argv[])
{
  // general problem parameters
  realtype T0 = RCONST(0.0);         // initial time
  realtype Tf = RCONST(10.0);        // final time
  sunindextype NEQ = 3;              // number of dependent vars.
  int order = 4;                     // order of accuracy for RK method
  booleantype adaptive = SUNTRUE;    // adaptive run vs convergence order
  int test = 2;                      // test problem to run
  realtype a, b, ep, u0, v0, w0;     // parameters

  // general problem variables
  int retval;                        // reusable error-checking flag
  N_Vector y = NULL;                 // empty vector for storing solution
  N_Vector yref = NULL;              // empty vector for storing reference solution
  void *arkode_mem = NULL;           // empty ARKStep memory structure
  SUNMatrix A = NULL;                // empty matrix for solver
  SUNLinearSolver LS = NULL;         // empty linear solver object
  realtype rdata[3];

  //
  // Initialization
  //

  // Retrieve the command-line options:  ord test
  if (argc > 1)  order = (int) atoi(argv[1]);
  if (argc > 2)  test = (int) atoi(argv[2]);

  // Check arguments for validity
  //   1 <= test <= 3
  if ((test < 1) || (test > 3)) {
    cerr << "ERROR: test type be an integer in {1,2,3} \n";
    return(-1);
  }

  // Handle adaptive run vs order-of-convergence run
  if (order < 0) {
    adaptive = SUNFALSE;
    order = abs(order);
  }
  if (order == 0)  order = 4;

  // set up the test problem according to the desired test
  if (test == 1) {
    u0 = RCONST(3.9);
    v0 = RCONST(1.1);
    w0 = RCONST(2.8);
    a  = RCONST(1.2);
    b  = RCONST(2.5);
    ep = RCONST(1.0e-5);
  } else if (test == 3) {
    u0 = RCONST(3.0);
    v0 = RCONST(3.0);
    w0 = RCONST(3.5);
    a  = RCONST(0.5);
    b  = RCONST(3.0);
    ep = RCONST(5.0e-4);
  } else {
    u0 = RCONST(1.2);
    v0 = RCONST(3.1);
    w0 = RCONST(3.0);
    a  = RCONST(1.0);
    b  = RCONST(3.5);
    ep = RCONST(5.0e-6);
  }
  rdata[0] = a;
  rdata[1] = b;
  rdata[2] = ep;

  // Initial problem output (and set implicit solver tolerances as needed)
  cout << "\nAccumulated error estimation test (Brusselator ODE problem):\n";
  cout << "    time domain:  (" << T0 << "," << Tf << "]\n";
  cout << "    initial conditions:  u0 = " << u0 << ",  v0 = " << v0 << ",  w0 = " << w0 << endl;
  cout << "    problem parameters:  a = " << a << ",  b = " << b << ",  ep = " << ep << endl;
  cout << "    DIRK solver, order = " << order << endl;

  //
  // Problem Setup
  //

  // Create SUNDIALS context
  sundials::Context ctx;

  // Create serial vectors for the solution and reference
  y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return 1;
  yref = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void *)yref, "N_VNew_Serial", 0)) return 1;

  // Generate reference solution.
  NV_Ith_S(yref,0) = u0;
  NV_Ith_S(yref,1) = v0;
  NV_Ith_S(yref,2) = w0;
  arkode_mem = ARKStepCreate(NULL, fn, T0, yref, ctx);
  if (check_retval((void *)arkode_mem, "ARKStepCreate", 0)) return 1;
  retval = ARKStepSetUserData(arkode_mem, (void *) rdata);
  if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_retval((void *)A, "SUNDenseMatrix", 0)) return 1;
  LS = SUNLinSol_Dense(yref, A, ctx);
  if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) return 1;
  retval = ARKStepSetLinearSolver(arkode_mem, LS, A);
  if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return 1;
  retval = ARKStepSetJacFn(arkode_mem, Jac);
  if (check_retval(&retval, "ARKStepSetJacFn", 1)) return 1;
  retval = ARKStepSetOrder(arkode_mem, 5);
  if (check_retval(&retval, "ARKStepSetOrder", 1)) return 1;
  retval = ARKStepSStolerances(arkode_mem, RCONST(1.e-9), RCONST(1.e-12));
  if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
  retval = ARKStepSetStopTime(arkode_mem, Tf);
  if (check_retval(&retval, "ARKStepSetStopTime", 1)) return 1;
  realtype t;
  retval = ARKStepEvolve(arkode_mem, Tf, yref, &t, ARK_NORMAL);
  if (check_retval(&retval, "ARKStepEvolve", 1)) return 1;

  // Reset ARKStep to specified order and default tolerances.
  NV_Ith_S(y,0) = u0;
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(y,2) = w0;
  retval = ARKStepSetOrder(arkode_mem, order);
  if (check_retval(&retval, "ARKStepSetOrder", 1)) return 1;
  retval = ARKStepSStolerances(arkode_mem, RCONST(1.e-4), RCONST(1.e-9));
  if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;

  // Integrate ODE, based on run type
  if (adaptive) {
    retval = adaptive_run(arkode_mem, y, T0, Tf, yref);
    if (check_retval(&retval, "adaptive_run", 1)) return 1;
  } else {
    retval = fixed_run(arkode_mem, y, T0, Tf, yref);
    if (check_retval(&retval, "fixed_run", 1)) return 1;
  }

  // Clean up and return
  ARKStepFree(&arkode_mem);
  SUNLinSolFree(LS);         // free system linear solver
  SUNMatDestroy(A);          // free system matrix
  N_VDestroy(y);             // Free y vector
  N_VDestroy(yref);          // Free yref vector
  return 0;

  return 0;
}

//------------------------------
// Functions called by the solver
//------------------------------

static int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype a  = rdata[0];                     // access data entries
  realtype b  = rdata[1];
  realtype ep = rdata[2];
  realtype u = NV_Ith_S(y,0);                 // access solution values
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  // fill in the RHS function
  NV_Ith_S(ydot,0) = a - (w+1.0)*u + v*u*u;
  NV_Ith_S(ydot,1) = w*u - v*u*u;
  NV_Ith_S(ydot,2) = (b-w)/ep - w*u;

  // Return with success
  return 0;
}

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype ep = rdata[2];                     // access data entries
  realtype u = NV_Ith_S(y,0);                 // access solution values
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  // fill in the Jacobian
  SM_ELEMENT_D(J,0,0) = -(w+1.0) + 2.0*u*v;
  SM_ELEMENT_D(J,0,1) = u*u;
  SM_ELEMENT_D(J,0,2) = -u;

  SM_ELEMENT_D(J,1,0) = w - 2.0*u*v;
  SM_ELEMENT_D(J,1,1) = -u*u;
  SM_ELEMENT_D(J,1,2) = u;

  SM_ELEMENT_D(J,2,0) = -w;
  SM_ELEMENT_D(J,2,1) = 0.0;
  SM_ELEMENT_D(J,2,2) = -1.0/ep - u;

  // Return with success
  return 0;
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

//------------------------------
// Private helper functions
//------------------------------

static int adaptive_run(void *arkode_mem, N_Vector y, realtype T0,
                        realtype Tf, N_Vector yref)
{
  // Reused variables
  int retval;
  long int nsteps;
  realtype dsm_est;
  realtype t = T0;
  N_Vector y0 = N_VClone(y);
  N_VScale(RCONST(1.0), y, y0);

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
      N_VScale(RCONST(1.0), y0, y);
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

      // Compute/print solution error
      realtype udsm = abs(NV_Ith_S(y,0)-NV_Ith_S(yref,0))/(abstol + rtols[irtol]*abs(NV_Ith_S(yref,0)));
      realtype vdsm = abs(NV_Ith_S(y,1)-NV_Ith_S(yref,1))/(abstol + rtols[irtol]*abs(NV_Ith_S(yref,1)));
      realtype wdsm = abs(NV_Ith_S(y,2)-NV_Ith_S(yref,2))/(abstol + rtols[irtol]*abs(NV_Ith_S(yref,2)));
      realtype dsm = rtols[irtol]*sqrt((udsm*udsm + vdsm*vdsm + wdsm*wdsm)/3);
      cout << "     acc type = " << accum_types[iaccum]
           << ",  dsm = " << dsm
           << ",  dsm_est = " << dsm_est
           << ",  dsm/dsm_est = " << dsm/dsm_est
           << ",  nsteps = " << nsteps
           << endl;
    }
  }
  cout << "   ------------------------------------------------------\n";

  N_VDestroy(y0);
  return(0);
}

static int fixed_run(void *arkode_mem, N_Vector y, realtype T0, realtype Tf,
                     N_Vector yref)
{
  // Reused variables
  int retval;
  long int nsteps, nsteps2;
  realtype dsm_est;
  realtype t, t2;
  N_Vector y0 = N_VClone(y);
  N_Vector y2 = N_VClone(y);
  N_Vector ewt = N_VClone(y);;
  N_VScale(RCONST(1.0), y, y0);

  // Set array of fixed step sizes to use, storage for corresponding errors/orders
  realtype hmax = (Tf - T0)/400;
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
      N_VScale(RCONST(1.0), y0, y);
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

      // Compute/print solution error
      realtype udsm = abs(NV_Ith_S(y,0)-NV_Ith_S(yref,0))/((1.e-12) + (1.e-9)*abs(NV_Ith_S(yref,0)));
      realtype vdsm = abs(NV_Ith_S(y,1)-NV_Ith_S(yref,1))/((1.e-12) + (1.e-9)*abs(NV_Ith_S(yref,1)));
      realtype wdsm = abs(NV_Ith_S(y,2)-NV_Ith_S(yref,2))/((1.e-12) + (1.e-9)*abs(NV_Ith_S(yref,2)));
      realtype dsm = (1.e-9)*sqrt((udsm*udsm + vdsm*vdsm + wdsm*wdsm)/3);
      cout << "     acc type = " << accum_types[iaccum]
           << ",  dsm = " << dsm
           << ",  dsm_est = " << dsm_est
           << ",  dsm/dsm_est = " << dsm/dsm_est
           << ",  nsteps = " << nsteps
           << endl;
    }

    // Test double fixed-step run error estimator
    t = T0;
    N_VScale(RCONST(1.0), y0, y);
    t2 = T0;
    N_VScale(RCONST(1.0), y0, y2);
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

    realtype udsm = abs(NV_Ith_S(y,0)-NV_Ith_S(yref,0))/((1.e-12) + (1.e-9)*abs(NV_Ith_S(yref,0)));
    realtype vdsm = abs(NV_Ith_S(y,1)-NV_Ith_S(yref,1))/((1.e-12) + (1.e-9)*abs(NV_Ith_S(yref,1)));
    realtype wdsm = abs(NV_Ith_S(y,2)-NV_Ith_S(yref,2))/((1.e-12) + (1.e-9)*abs(NV_Ith_S(yref,2)));
    realtype dsm = (1.e-9)*sqrt((udsm*udsm + vdsm*vdsm + wdsm*wdsm)/3);
    cout << "     acc type = " << 2
         << ",  dsm = " << dsm
         << ",  dsm_est = " << dsm_est
         << ",  dsm/dsm_est = " << dsm/dsm_est
         << ",  nsteps = " << nsteps
         << endl;

  }
  cout << "   ------------------------------------------------------\n";

  N_VDestroy(y0);
  N_VDestroy(y2);
  N_VDestroy(ewt);
  return(0);
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval >= 0
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


/*---- end of file ----*/
