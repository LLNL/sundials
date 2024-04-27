/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Routine to test an MRI method's embedding-based error
 * estimate.  Uses a simple polynomial test problem with 1
 * component,
 *    dy/dt = a*t^2 + b*t + c,  y(0) = 1.
 *
 * We run a single time step, using a variety of slow step sizes,
 * H = {1, 1/2, 1/4, 1/8, 1/16}.
 *
 * We place the entire ODE in the "slow" RHS partition.  All tests
 * use ARKODE's default fifth-order ERK method, with relative and
 * absolute tolerances set to 1e-10 and 1e-12, respectively.
 *
 * We select the slow integrator based on a command-line argument,
 * with the default being ARKODE_MRI_GARK_ERK22a.
 *
 * The program should be run with arguments in the following order:
 *   $ a.out method a b c
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out method a b
 *   $ a.out method a
 *   $ a.out method
 *   $ a.out
 * are acceptable.  We require that the "method" argument be a
 * string corresponding to a valid embedded ARKODE_MRITableID.
 *-----------------------------------------------------------------*/

// Header files
#include <algorithm>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>
#include <cmath>
#include <iostream>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <vector>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

using namespace std;

// User data structure
struct UserData
{
  sunrealtype a;
  sunrealtype b;
  sunrealtype c;
};

// User-supplied Functions Called by the Solver
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private utility functions
static int Ytrue(sunrealtype t, N_Vector y, UserData& udata);
static int computeErrorWeights(N_Vector ycur, N_Vector weight, sunrealtype rtol,
                               sunrealtype atol, N_Vector vtemp);
static int check_retval(void* returnvalue, const char* funcname, int opt);
static int run_test(void* mristep_mem, N_Vector y, sunrealtype T0,
                    vector<sunrealtype> Hvals, char* method, sunrealtype reltol,
                    sunrealtype abstol, UserData& udata);

// Main Program
int main(int argc, char* argv[])
{
  // general problem parameters
  sunrealtype T0   = ZERO;                 // initial time
  sunindextype NEQ = 1;                    // number of dependent vars.
  char* method;                            // MRI method name
  sunrealtype reltol = SUN_RCONST(1.e-10); // fast solver tolerances
  sunrealtype abstol = SUN_RCONST(1.e-12);

  // general problem variables
  int retval;     // reusable error-checking flag
  UserData udata; // user-data structure

  //
  // Initialization
  //

  // Retrieve the command-line options:  method Npart ep test
  if (argc > 1) { method = argv[1]; }
  else { method = "ARKODE_MRI_GARK_ERK22a"; }
  if (argc > 2) { udata.a = (sunrealtype)atof(argv[2]); }
  else { udata.a = ONE; }
  if (argc > 3) { udata.b = (sunrealtype)atof(argv[3]); }
  else { udata.b = ONE; }
  if (argc > 4) { udata.c = (sunrealtype)atof(argv[4]); }
  else { udata.c = ONE; }

  // Initial problem output (and set implicit solver tolerances as needed)
  cout << "\nSlow error estimation test (polynomial ODE problem):\n";
  cout << "    problem parameters:  a = " << udata.a << ",  b = " << udata.b
       << ",  c = " << udata.c << endl;
  cout << "    MRI method: " << method << endl;

  //
  // Problem Setup
  //

  // Create SUNDIALS context
  sundials::Context ctx;

  // Create serial vectors for the solution and reference
  N_Vector y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) return 1;
  retval = Ytrue(T0, y, udata);
  if (check_retval(&retval, "Ytrue", 1)) return 1;

  // Set up fast ARKStep integrator as fifth-order adaptive ERK
  void* inner_arkode_mem = ARKStepCreate(f0, NULL, T0, y, ctx);
  if (check_retval((void*)inner_arkode_mem, "ARKStepCreate", 0)) return 1;
  retval = ARKStepSetOrder(inner_arkode_mem, 5);
  if (check_retval(&retval, "ARKStepSetOrder", 1)) return 1;
  retval = ARKStepSStolerances(inner_arkode_mem, reltol, abstol);
  if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
  retval = ARKStepSetMaxNumSteps(inner_arkode_mem, 1000000);
  if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return (1);

  // Create inner stepper wrapper
  MRIStepInnerStepper inner_stepper = NULL; // inner stepper
  retval = ARKStepCreateMRIStepInnerStepper(inner_arkode_mem, &inner_stepper);
  if (check_retval(&retval, "ARKStepCreateMRIStepInnerStepper", 1)) return 1;

  // Set up slow MRIStep integrator
  sunbooleantype implicit = SUNFALSE;
  if ((strcmp(method, "ARKODE_MRI_GARK_IRK21a") == 0) ||
      (strcmp(method, "ARKODE_MRI_GARK_ESDIRK34a") == 0) ||
      (strcmp(method, "ARKODE_MRI_GARK_ESDIRK46a") == 0))
  {
    implicit = SUNTRUE;
  }
  void* mristep_mem = NULL;
  if (implicit)
  {
    mristep_mem = MRIStepCreate(NULL, fn, T0, y, inner_stepper, ctx);
  }
  else { mristep_mem = MRIStepCreate(fn, NULL, T0, y, inner_stepper, ctx); }
  if (check_retval((void*)mristep_mem, "MRIStepCreate", 0)) return 1;
  MRIStepCoupling C = MRIStepCoupling_LoadTableByName(method);
  if (check_retval((void*)C, "MRIStepCoupling_LoadTableByName", 0)) return 1;
  retval = MRIStepSetCoupling(mristep_mem, C);
  if (check_retval(&retval, "MRIStepSetCoupling", 1)) return 1;
  SUNMatrix A        = NULL; // matrix for slow solver
  SUNLinearSolver LS = NULL; // slow linear solver object
  if (implicit)
  {
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return 1;
    LS = SUNLinSol_Dense(y, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return 1;
    retval = MRIStepSetLinearSolver(mristep_mem, LS, A);
    if (check_retval(&retval, "MRIStepSetLinearSolver", 1)) return 1;
    retval = MRIStepSetJacFn(mristep_mem, Jn);
    if (check_retval(&retval, "MRIStepSetJacFn", 1)) return 1;
    retval = MRIStepSetJacEvalFrequency(mristep_mem, 1);
    if (check_retval(&retval, "MRIStepSetJacEvalFrequency", 1)) return 1;
    retval = MRIStepSetLSetupFrequency(mristep_mem, 1);
    if (check_retval(&retval, "MRIStepSetLSetupFrequency", 1)) return 1;
    retval = MRIStepSetMaxNonlinIters(mristep_mem, 20);
    if (check_retval(&retval, "MRIStepSetMaxNonlinIters", 1)) return 1;
  }
  retval = MRIStepSStolerances(mristep_mem, reltol, abstol);
  if (check_retval(&retval, "MRIStepSStolerances", 1)) return 1;
  retval = MRIStepSetUserData(mristep_mem, (void*)&udata);
  if (check_retval(&retval, "MRIStepSetUserData", 1)) return 1;
  retval = MRIStepSetAccumulatedErrorType(mristep_mem, 0);
  if (check_retval(&retval, "MRIStepSetAccumulatedErrorType", 1)) return 1;

  // Run test for various H values
  //  vector<sunrealtype> Hvals = {ONE, ONE/2.0, ONE/4.0, ONE/8.0, ONE/16.0};
  vector<sunrealtype> Hvals = {ONE / 100.0, ONE / 200.0, ONE / 400.0,
                               ONE / 800.0, ONE / 1600.0};
  retval = run_test(mristep_mem, y, T0, Hvals, method, reltol, abstol, udata);
  if (check_retval(&retval, "run_test", 1)) return 1;

  // Clean up and return
  ARKStepFree(&inner_arkode_mem);
  MRIStepInnerStepper_Free(&inner_stepper);
  MRIStepFree(&mristep_mem);
  if (LS) { SUNLinSolFree(LS); } // free system linear solver
  if (A) { SUNMatDestroy(A); }   // free system matrix
  N_VDestroy(y);                 // Free y and yref vectors
  return 0;
}

//------------------------------
// Functions called by the solver
//------------------------------

static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  // fill in the RHS function with zeros and return with success
  N_VConst(ZERO, ydot);
  return 0;
}

static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  // fill in the RHS function and return with success
  UserData* udata   = (UserData*)user_data;
  NV_Ith_S(ydot, 0) = udata->a * t * t + udata->b * t + udata->c;
  return 0;
}

static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // fill in the Jacobian and return with success
  SM_ELEMENT_D(J, 0, 0) = ZERO;
  return 0;
}

//------------------------------
// Private helper functions
//------------------------------

static int run_test(void* mristep_mem, N_Vector y, sunrealtype T0,
                    vector<sunrealtype> Hvals, char* method, sunrealtype reltol,
                    sunrealtype abstol, UserData& udata)
{
  // Reused variables
  int retval;
  sunrealtype t;
  N_Vector vtemp = N_VClone(y);
  N_Vector ele   = N_VClone(y);
  N_Vector ewt   = N_VClone(y);

  // Set storage for errors
  vector<sunrealtype> dsm(Hvals.size(), ZERO);
  vector<sunrealtype> dsm_est(Hvals.size(), ZERO);

  // Loop over step sizes
  for (size_t iH = 0; iH < Hvals.size(); iH++)
  {
    // Reset integrator for this run
    t      = T0;
    retval = Ytrue(t, y, udata);
    if (check_retval(&retval, "Ytrue", 1)) return 1;
    retval = MRIStepReset(mristep_mem, t, y);
    if (check_retval(&retval, "MRIStepReset", 1)) return 1;
    retval = MRIStepSetFixedStep(mristep_mem, Hvals[iH]);
    if (check_retval(&retval, "MRIStepSetFixedStep", 1)) return 1;
    retval = MRIStepResetAccumulatedError(mristep_mem);
    if (check_retval(&retval, "MRIStepResetAccumulatedError", 1)) return 1;

    // Run MRIStep to compute one step
    retval = MRIStepEvolve(mristep_mem, t + Hvals[iH], y, &t, ARK_ONE_STEP);
    if (check_retval(&retval, "MRIStepEvolve", 1)) return 1;
    retval = MRIStepGetEstLocalErrors(mristep_mem, ele);
    if (check_retval(&retval, "MRIStepGetEstLocalErrors", 1)) return 1;
    retval = computeErrorWeights(y, ewt, reltol, abstol, vtemp);
    if (check_retval(&retval, "computeErrorWeights", 1)) return 1;
    dsm_est[iH] = N_VWrmsNorm(ewt, ele);

    // Compute/print solution error
    retval = Ytrue(t, vtemp, udata);
    if (check_retval(&retval, "Ytrue", 1)) return 1;
    dsm[iH] = abs(NV_Ith_S(y, 0) - NV_Ith_S(vtemp, 0)) /
              (abstol + reltol * abs(NV_Ith_S(vtemp, 0)));
    if (strcmp(method, "ARKODE_MRI_GARK_ERK22a") == 0)
    {
      printf("       H  %.5f    dsm  %.2e    dsm_est  %.2e    dsm_anal  %.2e   "
             " dsm_est_anal  %.2e\n",
             Hvals[iH], dsm[iH], dsm_est[iH],
             Hvals[iH] * Hvals[iH] * Hvals[iH] * abs(udata.a / 12.0) /
               (abstol + reltol * abs(NV_Ith_S(vtemp, 0))),
             Hvals[iH] * Hvals[iH] *
               abs(udata.a * Hvals[iH] / 4.0 + udata.b / 2.0) /
               (abstol + reltol * abs(NV_Ith_S(vtemp, 0))));
    }
    else if (strcmp(method, "ARKODE_MRI_GARK_IRK21a") == 0)
    {
      printf("       H  %.5f    dsm  %.2e    dsm_est  %.2e    dsm_anal  %.2e   "
             " dsm_est_anal  %.2e\n",
             Hvals[iH], dsm[iH], dsm_est[iH],
             Hvals[iH] * Hvals[iH] * Hvals[iH] * abs(udata.a / 6.0) /
               (abstol + reltol * abs(NV_Ith_S(vtemp, 0))),
             Hvals[iH] * Hvals[iH] *
               abs(udata.a * Hvals[iH] / 2.0 + udata.b / 2.0) /
               (abstol + reltol * abs(NV_Ith_S(vtemp, 0))));
    }
    else
    {
      printf("       H  %.5f    dsm  %.2e    dsm_est  %.2e\n", Hvals[iH],
             dsm[iH], dsm_est[iH]);
    }
  }

  N_VDestroy(ele);
  N_VDestroy(ewt);
  N_VDestroy(vtemp);
  return (0);
}

static int Ytrue(sunrealtype t, N_Vector y, UserData& udata)
{
  NV_Ith_S(y, 0) = udata.a / SUN_RCONST(3.0) * t * t * t +
                   udata.b / SUN_RCONST(2.0) * t * t + udata.c * t + ONE;
  return (0);
}

/* Error weight calculation routine (mimics what's in ARKODE already) */
static int computeErrorWeights(N_Vector ycur, N_Vector weight, sunrealtype rtol,
                               sunrealtype atol, N_Vector vtemp)
{
  N_VAbs(ycur, vtemp);
  N_VScale(rtol, vtemp, vtemp);
  N_VAddConst(vtemp, atol, vtemp);
  N_VInv(vtemp, weight);
  return (0);
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  // Check if retval < 0
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

/*---- end of file ----*/
