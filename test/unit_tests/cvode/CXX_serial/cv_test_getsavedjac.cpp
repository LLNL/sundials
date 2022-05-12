/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
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

#include "cv_test_getsavedjac.hpp"

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Create initial condition
  N_Vector y = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) return 1;

  int flag = ytrue(ZERO, y);
  if (check_flag(flag, "ytue")) return 1;

  // Create CVODE memory structure
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_ptr(cvode_mem, "CVodeCreate")) return 1;

  flag = CVodeInit(cvode_mem, f, ZERO, y);
  if (check_flag(flag, "CVodeInit")) return 1;

  flag = CVodeSStolerances(cvode_mem, RCONST(1.0e-6), RCONST(1.0e-10));
  if (check_flag(flag, "CVodeSStolerances")) return 1;

  SUNMatrix A = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(A, "SUNDenseMatrix")) return 1;

  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_ptr(LS, "SUNLinSol_Dense")) return 1;

  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_flag(flag, "CVodeSetLinearSolver")) return 1;

  realtype udata[4] = {-TWO, HALF, HALF, -ONE};
  flag = CVodeSetUserData(cvode_mem, udata);
  if (check_flag(flag, "CVodeSetUserData")) return 1;

  // Initial time and fist output time
  realtype tret  = ZERO;
  realtype tout  = tret + RCONST(0.1);

  // Advance one step in time
  flag = CVode(cvode_mem, tout, y, &tret, CV_ONE_STEP);
  if (check_flag(flag, "CVode")) return 1;

  // Get the saved internal finite difference Jacobian = J(t_pred, y_pred) which
  // approximates J(t_ret, y(t_ret))
  SUNMatrix Jsaved;
  flag = CVodeGetSavedJac(cvode_mem, &Jsaved);
  if (check_flag(flag, "CVodeGetSavedJac")) return 1;

  // Compute the true Jacobian at the returned time and true state
  SUNMatrix Jtrue = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(Jtrue, "SUNDenseMatrix")) return 1;

  flag = ytrue(tret, y);
  if (check_flag(flag, "ytue")) return 1;

  flag = J(tret, y, nullptr, Jtrue, &udata, nullptr, nullptr, nullptr);
  if (check_flag(flag, "J")) return 1;

  // Compare finite difference and true Jacobian
  realtype* Jsaved_data = SUNDenseMatrix_Data(Jsaved);
  if (check_ptr(Jsaved_data, "SUNDenseMatrix_Data")) return 1;

  realtype* Jtrue_data = SUNDenseMatrix_Data(Jtrue);
  if (check_ptr(Jtrue_data, "SUNDenseMatrix_Data")) return 1;

  // Output Jacobian
  cout << scientific;
  cout << setprecision(numeric_limits<realtype>::digits10);
  cout << setw(8)  << right << "Index"
       << setw(25) << right << "J saved"
       << setw(25) << right << "J true"
       << setw(25) << right << "absolute difference"
       << endl;
  for (int i = 0; i < 3 * 25 + 8; i++)
    cout << "-";
  cout << endl;

  sunindextype ldata = SUNDenseMatrix_LData(Jtrue);
  for (sunindextype i = 0; i < ldata; i++)
  {
    cout << setw(8)  << right << i
         << setw(25) << right << Jsaved_data[i]
         << setw(25) << right << Jtrue_data[i]
         << setw(25) << right << abs(Jsaved_data[i] - Jtrue_data[i])
         << endl;
  }

  // Clean up and return with successful completion
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNMatDestroy(Jtrue);
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
