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
 * Kvaerno-Prothero-Robinson ODE test problem:
 *
 *   [u]' = [ a  b ] [ (-1 + u^2 - r(t)) / (2u) ] + [ r'(t) / (2u) ]
 *   [v]    [ c  d ] [ (-2 + v^2 - s(t)) / (2v) ]   [ s'(t) / (2v) ]
 *
 * This problem has analytical solution given by
 *
 *   u(t) = sqrt(1 + r(t))
 *   v(t) = sqrt(2 + s(t))
 *
 * where, in this test, we use the functions
 *
 *   r(t) = 0.5 * cos(t)
 *   s(t) = cos(2t)
 * ---------------------------------------------------------------------------*/

#include "cvs_test_kpr.hpp"

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Read input options
  TestOptions opts;

  vector<string> args(argv + 1, argv + argc);

  int flag = ReadInputs(args, opts, sunctx);
  if (check_flag(flag, "ReadInputs")) return 1;

  // Create initial condition
  N_Vector y = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) return 1;

  realtype utrue, vtrue;
  flag = true_sol(ZERO, &utrue, &vtrue);
  if (check_flag(flag, "true_sol")) return 1;

  realtype* ydata = N_VGetArrayPointer(y);
  ydata[0] = utrue;
  ydata[1] = vtrue;

  // Create CVODES memory structure
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_ptr(cvode_mem, "CVodeCreate")) return 1;

  flag = CVodeInit(cvode_mem, f, ZERO, y);
  if (check_flag(flag, "CVodeInit")) return 1;

  flag = CVodeSStolerances(cvode_mem, opts.rtol, opts.atol);
  if (check_flag(flag, "CVodeSStolerances")) return 1;

  SUNMatrix A = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(A, "SUNDenseMatrix")) return 1;

  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_ptr(LS, "SUNLinSol_Dense")) return 1;

  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_flag(flag, "CVodeSetLinearSolver")) return 1;

  flag = CVodeSetJacFn(cvode_mem, J);
  if (check_flag(flag, "CVodeSetJacFn")) return 1;

  realtype udata[4] = {-TWO, HALF, HALF, -ONE};
  flag = CVodeSetUserData(cvode_mem, udata);
  if (check_flag(flag, "CVodeSetUserData")) return 1;

  flag = CVodeSetEtaFixedStepBounds(cvode_mem,
                                    opts.eta_min_fx, opts.eta_max_fx);
  if (check_flag(flag, "CVodeSetEtaFixeStepBounds")) return 1;

  flag = CVodeSetEtaMaxFirstStep(cvode_mem, opts.eta_max_fs);
  if (check_flag(flag, "CVodeSetEtaMaxFirstStep")) return 1;

  flag = CVodeSetEtaMaxEarlyStep(cvode_mem, opts.eta_max_es);
  if (check_flag(flag, "CVodeSetEtaMaxEarlyStep")) return 1;

  flag = CVodeSetNumStepsEtaMaxEarlyStep(cvode_mem, opts.small_nst);
  if (check_flag(flag, "CVodeSetNumStepsEtaMaxEarlyStep")) return 1;

  flag = CVodeSetEtaMax(cvode_mem, opts.eta_max_gs);
  if (check_flag(flag, "CVodeSetEtaMax")) return 1;

  flag = CVodeSetEtaMin(cvode_mem, opts.eta_min);
  if (check_flag(flag, "CVodeSetEtaMin")) return 1;

  flag = CVodeSetEtaMinErrFail(cvode_mem, opts.eta_min_ef);
  if (check_flag(flag, "CVodeSetEtaMinErrFail")) return 1;

  flag = CVodeSetEtaMaxErrFail(cvode_mem, opts.eta_max_ef);
  if (check_flag(flag, "CVodeSetEtaMaxErrFail")) return 1;

  flag = CVodeSetNumFailsEtaMaxErrFail(cvode_mem, opts.small_nef);
  if (check_flag(flag, "CVodeSetNumFailsEtaMaxErrFail")) return 1;

  flag = CVodeSetEtaConvFail(cvode_mem, opts.eta_cf);
  if (check_flag(flag, "CVodeSetEtaConvFail")) return 1;

  flag = CVodeSetDeltaGammaMaxLSetup(cvode_mem, opts.dgmax_lsetup);
  if (check_flag(flag, "CVodeSetDeltaGammaMaxLSetup")) return 1;

  flag = CVodeSetDeltaGammaMaxBadJac(cvode_mem, opts.dgmax_jbad);
  if (check_flag(flag, "CVodeSetDeltaGammaMaxBadJac")) return 1;

  // Initial time and fist output time
  realtype tret  = ZERO;
  realtype tout  = tret + opts.dtout;

  // Output initial contion
  cout << scientific;
  cout << setprecision(numeric_limits<realtype>::digits10);
  cout << "           t              ";
  cout << "          u              ";
  cout << "          v              ";
  cout << "        u err            ";
  cout << "        v err            " << endl;
  for (int i = 0; i < 9; i++)
    cout << "--------------";
  cout << endl;

  cout << setw(22) << tret
       << setw(25) << ydata[0]
       << setw(25) << ydata[1]
       << setw(25) << abs(ydata[0] - utrue)
       << setw(25) << abs(ydata[1] - vtrue) << endl;

  // Advance in time
  for (int i = 0; i < opts.nout; i++)
  {
    flag = CVode(cvode_mem, tout, y, &tret, CV_NORMAL);
    if (check_flag(flag, "CVode")) return 1;

    flag = true_sol(tret, &utrue, &vtrue);
    if (check_flag(flag, "true_sol")) return 1;

    cout << setw(22) << tret
         << setw(25) << ydata[0]
         << setw(25) << ydata[1]
         << setw(25) << abs(ydata[0] - utrue)
         << setw(25) << abs(ydata[1] - vtrue) << endl;

    // update output time
    tout += opts.dtout;
  }
  for (int i = 0; i < 9; i++)
    cout << "--------------";
  cout << endl;

  // Print some final statistics
  flag = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "CVodePrintAllStats")) return 1;

  // Clean up and return with successful completion
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  CVodeFree(&cvode_mem);

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the ODE RHS function:
 *   [a  b] * [ (-1 + u^2 - r(t)) / (2*u) ] + [ r'(t) / (2u) ]
 *   [c  d]   [ (-2 + v^2 - s(t)) / (2*v) ]   [ s'(t) / (2v) ]
 * ---------------------------------------------------------------------------*/
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype* udata = (realtype *) user_data;
  const realtype a = udata[0];
  const realtype b = udata[1];
  const realtype c = udata[2];
  const realtype d = udata[3];

  realtype* ydata = N_VGetArrayPointer(y);
  const realtype u = ydata[0];
  const realtype v = ydata[1];

  const realtype tmp1 = (-ONE + u * u - r(t)) / (TWO * u);
  const realtype tmp2 = (-TWO + v * v - s(t)) / (TWO * v);

  realtype* fdata = N_VGetArrayPointer(ydot);
  fdata[0] = a * tmp1 + b * tmp2 + rdot(t) / (TWO * u);
  fdata[1] = c * tmp1 + d * tmp2 + sdot(t) / (TWO * v);

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the ODE RHS Jacobin:
 *   [a/2 + (a(1+r(t))-rdot(t))/(2u^2)     b/2 + b*(2+s(t))/(2*v^2)         ]
 *   [c/2 + c(1+r(t))/(2u^2)               d/2 + (d(2+s(t))-sdot(t))/(2u^2) ]
 * ---------------------------------------------------------------------------*/
int J(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype* udata = (realtype *) user_data;
  const realtype a = udata[0];
  const realtype b = udata[1];
  const realtype c = udata[2];
  const realtype d = udata[3];

  realtype* ydata = N_VGetArrayPointer(y);
  realtype* Jdata = SUNDenseMatrix_Data(J);

  const realtype u = ydata[0];
  const realtype v = ydata[1];

  Jdata[0] = a / TWO + (a * (ONE + r(t)) - rdot(t)) / (TWO * u * u);
  Jdata[1] = c / TWO +  c * (ONE + r(t)) / (TWO * u * u);
  Jdata[2] = b / TWO +  b * (TWO + s(t)) / (TWO * v * v);
  Jdata[3] = d / TWO + (d * (TWO + s(t)) - sdot(t)) / (TWO * v * v);

  return 0;
}

/*---- end of file ----*/
