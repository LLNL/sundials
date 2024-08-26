/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
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

#ifndef KPR_
#define KPR_

#include <cmath>

#include "sunmatrix/sunmatrix_dense.h"

// Macros for problem constants
#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define TWO    SUN_RCONST(2.0)
#define TWENTY SUN_RCONST(20.0)

sunrealtype kpr_udata[4] = {-TWO, HALF, HALF, -ONE};

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute r(t)
inline sunrealtype kpr_r(sunrealtype t) { return HALF * std::cos(t); }

// Compute the derivative of r(t)
inline sunrealtype kpr_rdot(sunrealtype t) { return -HALF * std::sin(t); }

// Compute s(t)
inline sunrealtype kpr_s(sunrealtype t) { return cos(TWENTY * t); }

// Compute the derivative of s(t)
inline sunrealtype kpr_sdot(sunrealtype t)
{
  return -TWENTY * std::sin(TWENTY * t);
}

// Compute the true solution
inline int kpr_true_sol(sunrealtype t, sunrealtype* u, sunrealtype* v)
{
  *u = std::sqrt(ONE + kpr_r(t));
  *v = std::sqrt(TWO + kpr_s(t));

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE RHS function:
 *   [a  b] * [ (-1 + u^2 - r(t)) / (2*u) ] + [ r'(t) / (2u) ]
 *   [c  d]   [ (-2 + v^2 - s(t)) / (2*v) ]   [ s'(t) / (2v) ]
 * ---------------------------------------------------------------------------*/
inline int kpr_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  const sunrealtype tmp1 = (-ONE + u * u - kpr_r(t)) / (TWO * u);
  const sunrealtype tmp2 = (-TWO + v * v - kpr_s(t)) / (TWO * v);

  fdata[0] = a * tmp1 + b * tmp2 + kpr_rdot(t) / (TWO * u);
  fdata[1] = c * tmp1 + d * tmp2 + kpr_sdot(t) / (TWO * v);

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE RHS Jacobin:
 *   [a/2 + (a(1+r(t))-rdot(t))/(2u^2)     b/2 + b*(2+s(t))/(2*v^2)         ]
 *   [c/2 + c(1+r(t))/(2u^2)               d/2 + (d(2+s(t))-sdot(t))/(2u^2) ]
 * ---------------------------------------------------------------------------*/
inline int kpr_rhs_jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                       void* user_data, N_Vector tmp1, N_Vector tmp2,
                       N_Vector tmp3)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  Jdata[0] = a / TWO + (a * (ONE + kpr_r(t)) - kpr_rdot(t)) / (TWO * u * u);
  Jdata[1] = c / TWO + c * (ONE + kpr_r(t)) / (TWO * u * u);
  Jdata[2] = b / TWO + b * (TWO + kpr_s(t)) / (TWO * v * v);
  Jdata[3] = d / TWO + (d * (TWO + kpr_s(t)) - kpr_sdot(t)) / (TWO * v * v);

  return 0;
}



#endif // KPR_
