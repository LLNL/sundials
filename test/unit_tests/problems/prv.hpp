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
 * Prothero-Robinson ODE test problem with a time-varying coefficient (PRV):
 *
 *   y' = L(t) (y - phi(t)) + phi'(t),
 *
 * where phi(t) = atan(t), phi'(t) = 1 / (1 + t^2), and the coefficient is
 * given by
 *
 *   L(t) = lambda - alpha * cos((10 - t) / 10 * pi).
 *
 * This problem has analytical solution y(t) = atan(t).
 *
 * The stiffness of the problem depends on the value of L(t) where the lambda
 * determines the center of the parameter and alpha the radius of the interval
 * in which the stiffness parameter lies. For a well-posed problem, the values
 * should be chosen such that L(t) is negative.
 * ---------------------------------------------------------------------------*/

#ifndef PRV_
#define PRV_

#include <cmath>
#include <sundials/sundials_core.hpp>

namespace problems {
namespace prv {

// Problem constants
static const sunrealtype pi = std::acos(SUN_RCONST(-1.0));

constexpr sunrealtype zero = SUN_RCONST(0.0);
constexpr sunrealtype one  = SUN_RCONST(1.0);
constexpr sunrealtype ten  = SUN_RCONST(10.0);

// lambda and alpha
sunrealtype prv_data[2] = {SUN_RCONST(-1000.0), SUN_RCONST(10.0)};

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute L(t)
inline sunrealtype l_coef(sunrealtype t, sunrealtype c[2])
{
  return c[0] - c[1] * std::cos((ten - t) / ten * pi);
}

// Compute phi(t)
inline sunrealtype phi(sunrealtype t) { return std::atan(t); }

// Compute phi'(t)
inline sunrealtype phi_prime(sunrealtype t) { return one / (one + t * t); }

// Compute the true solution
inline sunrealtype true_solution(sunrealtype t) { return phi(t); }

// ODE RHS function
inline int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* u_data  = static_cast<sunrealtype*>(user_data);
  sunrealtype* y_data  = N_VGetArrayPointer(y);
  sunrealtype* yd_data = N_VGetArrayPointer(ydot);

  yd_data[0] = l_coef(t, u_data) * (y_data[0] - phi(t)) + phi_prime(t);

  return 0;
}

// Dominant eigenvalue function
inline int ode_dom_eig(sunrealtype t, N_Vector y, N_Vector fn,
                       sunrealtype* lambdaR, sunrealtype* lambdaI,
                       void* user_data, N_Vector temp1, N_Vector temp2,
                       N_Vector temp3)
{
  sunrealtype* u_data = static_cast<sunrealtype*>(user_data);

  *lambdaR = l_coef(t, u_data);
  *lambdaI = zero;

  return 0;
}

} // namespace prv
} // namespace problems

#endif