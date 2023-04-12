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
 * Prothero-Robinson ODE test problem:
 *
 *   y' = lambda (y - r(t)) + r'(t)
 *
 * with the analytic solution
 *
 *   y(t) = r(t)
 *
 * In this test, we let
 *
 *   r(t) = atan(t)
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>

#include "sundials/sundials_nvector.h"
#include "sunmatrix/sunmatrix_dense.h"

// Macros for problem constants
#define ZERO             SUN_RCONST(0.0)
#define ONE              SUN_RCONST(1.0)
#define TWO              SUN_RCONST(2.0)
#define FIVE             SUN_RCONST(5.0)
#define SIX              SUN_RCONST(6.0)
#define TEN              SUN_RCONST(10.0)
#define TWENTY           SUN_RCONST(20.0)
#define TWENTYFOUR       SUN_RCONST(24.0)
#define ONEHUNDREDTWENTY SUN_RCONST(120.0)

// PR problem data
struct PRData
{
  sunrealtype lambda = -ONE;
};

// Compute y(t)
static sunrealtype r(sunrealtype t)
{
  return std::atan(t);
}

// Compute the first five derivatives of y(t)
static sunrealtype r_dot(sunrealtype t)
{
  return ONE / (t * t + ONE);
}

static sunrealtype r_dot2(sunrealtype t)
{
  return -(TWO * t) / std::pow((t * t + ONE), 2);
}

static sunrealtype r_dot3(sunrealtype t)
{
  return (SIX * t * t - TWO) / std::pow((t * t + ONE), 3);
}

static sunrealtype r_dot4(sunrealtype t)
{
  return -(TWENTYFOUR * t * (t * t - ONE)) / std::pow((t * t + 1), 4);
}

static sunrealtype r_dot5(sunrealtype t)
{
  return (TWENTYFOUR * (FIVE * std::pow(t, 4) - TEN * t * t + ONE)) /
         std::pow((t * t + ONE), 5);
}

static int PR_true(sunrealtype t, N_Vector y, PRData& pr_data)
{
  const sunrealtype r_t = r(t);
  const sunindextype N  = N_VGetLength(y);

  sunrealtype* ydata = N_VGetArrayPointer(y);

  for (sunindextype i = 0; i < N; i++)
  {
    ydata[i] = r_t;
  }

  return 0;
}

static int PR_true_dot(sunrealtype t, N_Vector ydot, int order, PRData& pr_data)
{
  const sunindextype N = N_VGetLength(ydot);
  sunrealtype* ydata   = N_VGetArrayPointer(ydot);

  switch(order)
  {
  case 0:
    for (sunindextype i = 0; i < N; i++)
    {
      ydata[i] = r(t);
    }
    break;
  case 1:
    for (sunindextype i = 0; i < N; i++)
    {
      ydata[i] = r_dot(t);
    }
    break;
  case 2:
    for (sunindextype i = 0; i < N; i++)
    {
      ydata[i] = r_dot2(t);
    }
    break;
  case 3:
    for (sunindextype i = 0; i < N; i++)
    {
      ydata[i] = r_dot3(t);
    }
    break;
  case 4:
    for (sunindextype i = 0; i < N; i++)
    {
      ydata[i] = r_dot4(t);
    }
    break;
  case 5:
    for (sunindextype i = 0; i < N; i++)
    {
      ydata[i] = r_dot5(t);
    }
    break;
  default:
    std::cerr << "Invalid derivative order: " << order << std::endl;
    return -1;
  }

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the ODE RHS function:  y' = lambda (y - r(t)) + r'(t)
 * ---------------------------------------------------------------------------*/
int PR_Rhs(sunrealtype t, N_Vector y, N_Vector ydot, PRData* pr_data)
{
  const sunrealtype lambda  = pr_data->lambda;
  const sunrealtype r_t     = r(t);
  const sunrealtype r_dot_t = r_dot(t);
  const sunindextype N      = N_VGetLength(y);

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  for (sunindextype i = 0; i < N; i++)
  {
    fdata[i] = lambda * (ydata[i] - r_t) + r_dot_t;
  }

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the ODE RHS Jacobin: J = lambda
 * ---------------------------------------------------------------------------*/
int PR_Jac(sunrealtype t, N_Vector y, SUNMatrix J, PRData* pr_data)
{
  const sunrealtype lambda = pr_data->lambda;
  const sunindextype N     = N_VGetLength(y);

  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  for (sunindextype i = 0; i < N; i++)
  {
    Jdata[i + i * N] = lambda;
  }

  return 0;
}

/*---- end of file ----*/
