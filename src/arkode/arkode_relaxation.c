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
 * This is the implementation file for ARKODE's relaxation (in time)
 * functionality
 *
 * Temporary vectors utilized in the functions below:
 *   tempv1 - holds delta_y, the update direction vector
 *   tempv2 - holds y_relax, the relaxed solution vector
 *   tempv3 - holds J_relax, the Jacobian of the relaxation function
 * ---------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_relaxation_impl.h"
#include "sundials/sundials_types.h"

/* =============================================================================
 * Private Functions
 * ===========================================================================*/

/* Access the ARKODE and relaxation memory structures */
static int arkRelaxAccessMem(void* arkode_mem, const char* fname,
                             ARKodeMem* ark_mem, ARKodeRelaxMem* relax_mem)
{
  if (!arkode_mem)
  {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", fname, MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  if (!((*ark_mem)->relax_mem))
  {
    arkProcessError(*ark_mem, ARK_RELAX_MEM_NULL, "ARKODE", fname,
                    MSG_RELAX_MEM_NULL);
    return ARK_RELAX_MEM_NULL;
  }
  *relax_mem = (ARKodeRelaxMem)((*ark_mem)->relax_mem);

  return ARK_SUCCESS;
}

/* Evaluates the relaxation residual function */
static int arkRelaxResidual(sunrealtype relax_param, N_Vector y_relax,
                            sunrealtype* relax_res, ARKodeMem ark_mem)
{
  int retval;
  sunrealtype e_old   = ark_mem->relax_mem->e_old;
  sunrealtype delta_e = ark_mem->relax_mem->delta_e;
  void* user_data     = ark_mem->user_data;

  /* Evaluate entropy function */
  retval = ark_mem->relax_mem->relax_fn(y_relax, relax_res, user_data);
  ark_mem->relax_mem->num_relax_fn_evals++;
  if (retval < 0) return ARK_RELAX_FUNC_FAIL;
  if (retval > 0) return ARK_RELAX_FUNC_RECV;

  /* Compute relaxation residual */
  *relax_res = *relax_res - e_old - relax_param * delta_e;

  return ARK_SUCCESS;
}

/* Evaluates the Jacobian of the relaxation residual function */
static int arkRelaxResidualJacobian(sunrealtype relax_param, N_Vector y_relax,
                                    sunrealtype* relax_jac, ARKodeMem ark_mem)
{
  int retval;
  N_Vector delta_y    = ark_mem->tempv1;
  N_Vector J_relax    = ark_mem->tempv3;
  sunrealtype delta_e = ark_mem->relax_mem->delta_e;
  void* user_data     = ark_mem->user_data;

  /* Evaluate Jacobian of entropy functions */
  retval = ark_mem->relax_mem->relax_jac_fn(y_relax, J_relax, user_data);
  ark_mem->relax_mem->num_relax_jac_evals++;
  if (retval < 0) return ARK_RELAX_JAC_FAIL;
  if (retval > 0) return ARK_RELAX_JAC_RECV;

  /* Compute relaxation residual Jacobian */
  *relax_jac = N_VDotProd(delta_y, J_relax);
  *relax_jac -= delta_e;

  return ARK_SUCCESS;
}

/* Solve the relaxation residual equation using Newton's method */
static int arkRelaxNewtonSolve(ARKodeMem ark_mem)
{
  int i, retval;
  ARKodeRelaxMem relax_mem = ark_mem->relax_mem;
  N_Vector delta_y = ark_mem->tempv1;
  N_Vector y_relax = ark_mem->tempv2;

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO,
                     "ARKODE::arkStep_TakeStep_Z", "NewtonSolve",
                     "tolerance = %g", relax_mem->res_tol);
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO,
                     "ARKODE::arkStep_TakeStep_Z", "NewtonSolve",
                     "0: relaxation param = %g", relax_mem->relax_param);
#endif

  for (i = 0; i < ark_mem->relax_mem->max_iters; i++)
  {
    /* y_relax = y_n + r * delta_y */
    N_VLinearSum(ONE, ark_mem->yn, relax_mem->relax_param, delta_y, y_relax);

    /* Compute the current residual */
    retval = arkRelaxResidual(relax_mem->relax_param, y_relax,
                              &(relax_mem->res), ark_mem);
    if (retval) return retval;

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO,
                       "ARKODE::arkStep_TakeStep_Z", "NewtonSolve",
                       "%d: relaxation resid = %g", i, relax_mem->res);
#endif

    /* Check for convergence */
    if (SUNRabs(relax_mem->res) < relax_mem->res_tol) return ARK_SUCCESS;

    /* Compute Jacobian and update */
    retval = arkRelaxResidualJacobian(relax_mem->relax_param, y_relax,
                                      &(relax_mem->jac),
                                      ark_mem);
    if (retval) return retval;

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO,
                       "ARKODE::arkStep_TakeStep_Z", "NewtonSolve",
                       "%d: relaxation jac = %g", i, relax_mem->jac);
#endif

    relax_mem->relax_param -= relax_mem->res / relax_mem->jac;

    /* Update cumulative iteration count */
    relax_mem->nls_iters++;

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO,
                       "ARKODE::arkStep_TakeStep_Z", "NewtonSolve",
                       "%d: relaxation param = %g", i+1, relax_mem->relax_param);
#endif
  }

  return ARK_RELAX_SOLVE_RECV;
}

/* Solve the relaxation residual equation using Newton's method */
static int arkRelaxBrentsSolve(ARKodeMem ark_mem)
{
  int i, j, retval;
  sunrealtype prev, f_prev; /* previous solution and function value */
  sunrealtype curr, f_curr; /* current solution and function value  */
  sunrealtype brac, f_brac; /* together brac and curr bracket zero  */
  sunrealtype mid;          /* midpoint between brac and curr       */
  ARKodeRelaxMem relax_mem = ark_mem->relax_mem;
  N_Vector delta_y = ark_mem->tempv1;
  N_Vector y_relax = ark_mem->tempv2;

  /* Compute interval that brackets the root */
  prev = SUN_RCONST(0.9) * relax_mem->relax_param;
  curr = SUN_RCONST(1.1) * relax_mem->relax_param;
  for (i = 0; i < 10; i++)
  {
    /* y_relax = y_n + r * delta_y */
    N_VLinearSum(ONE, ark_mem->yn, prev, delta_y, y_relax);

    /* Compute relaxation residual */
    retval = arkRelaxResidual(prev, y_relax, &f_prev, ark_mem);
    ark_mem->relax_mem->num_relax_fn_evals++;
    if (retval < 0) { return ARK_RELAX_FUNC_FAIL; }
    if (retval > 0) { return ARK_RELAX_FUNC_RECV; }

    /* Check if we got lucky */
    if (SUNRabs(f_prev) < relax_mem->tol)
    {
      relax_mem->res = f_prev;
      relax_mem->relax_param = prev;
      return ARK_SUCCESS;
    }

    if (relax_mem->res < ZERO) break;

    f_curr = f_prev;
    curr = prev;
    prev *= SUN_RCONST(0.9);
  }
  if (relax_mem->res > ZERO) { return ARK_RELAX_SOLVE_RECV; }

  for (i = 0; i < 10; i++)
  {
    /* y_relax = y_n + r * delta_y */
    N_VLinearSum(ONE, ark_mem->yn, curr, delta_y, y_relax);

    /* Compute relaxation residual */
    retval = arkRelaxResidual(curr, y_relax, &f_curr, ark_mem);
    ark_mem->relax_mem->num_relax_fn_evals++;
    if (retval < 0) { return ARK_RELAX_FUNC_FAIL; }
    if (retval > 0) { return ARK_RELAX_FUNC_RECV; }

    /* Check if we got lucky */
    if (SUNRabs(f_curr) < relax_mem->tol)
    {
      relax_mem->res = f_curr;
      relax_mem->relax_param = curr;
      return ARK_SUCCESS;
    }

    if (relax_mem->res > ZERO) break;

    f_prev = f_curr;
    prev = curr;
    curr *= SUN_RCONST(1.1);
  }
  if (relax_mem->res < ZERO) { return ARK_RELAX_SOLVE_RECV; }

  /* Initialize values (prev = lower bound and curr = upper bound) */
  brac = prev;
  f_brac = f_prev;

  /* Find root */
  for (i = 0; i < ark_mem->relax_mem->max_iters; i++)
  {
    /* Ensure brac and curr bracket zero */
    if (signbit(f_prev) != signbit(f_curr))
    {
      brac = prev;
      f_brac = f_prev;
      // spre = scur = xcur - xpre
    }

    /* Ensure the current solution better of bracketing values */
    if (SUNRabs(f_brac) < SUNRabs(f_curr))
    {
      prev = curr;
      curr = brac;
      brac = prev;

      f_prev = f_curr;
      f_curr = f_brac;
      f_brac = f_prev;
    }

    /* Compute midpoint */
    mid = SUN_RCONST(0.5) * (c - b);

    if (SUNRabs(mid) < tol1)
    {
      relax_mem->relax_param = curr;
      return ARK_SUCCESS;
    }


    if (SUNRabs(e) > tol1 && SUNRabs(f_perv) > SUNRabs(f_cur))
    {
      /* Inverse quadratic interpolation */

    }
    else
    {
      /* Bisection */
      d = mid;
      e = d;
    }

    /* Update previous best solution */
    prev = curr;
    f_prev = f_curr;

    /* Compute new current solution */
    if (SUNRabs(d) > tol1)
    {
      curr += d;
    }
    else
    {
      curr += copysign(tol1, mid);
    }

    /* y_relax = y_n + r * delta_y */
    N_VLinearSum(ONE, ark_mem->yn, curr, delta_y, y_relax);

    /* Compute relaxation residual */
    retval = arkRelaxResidual(curr, y_relax, &f_curr, ark_mem);
    ark_mem->relax_mem->num_relax_fn_evals++;
    if (retval < 0) { return ARK_RELAX_FUNC_FAIL; }
    if (retval > 0) { return ARK_RELAX_FUNC_RECV; }
  }

  return ARK_RELAX_SOLVE_RECV;
}

/* Solve the relaxation residual equation using Newton's method */
static int arkRelaxFixedPointSolve(ARKodeMem ark_mem)
{
  int i, retval;
  ARKodeRelaxMem relax_mem = ark_mem->relax_mem;
  N_Vector delta_y = ark_mem->tempv1;
  N_Vector y_relax = ark_mem->tempv2;

  for (i = 0; i < ark_mem->relax_mem->max_iters; i++)
  {
    /* y_relax = y_n + r * delta_y */
    N_VLinearSum(ONE, ark_mem->yn, relax_mem->relax_param, delta_y, y_relax);

    /* Compute the current residual */
    retval = arkRelaxResidual(relax_mem->relax_param, y_relax,
                              &(relax_mem->res), ark_mem);
    if (retval) return retval;

    /* Check for convergence */
    if (relax_mem->res < relax_mem->res_tol) return ARK_SUCCESS;

    relax_mem->relax_param -= relax_mem->res;

    /* Update iteration count */
    relax_mem->nls_iters++;
  }

  return ARK_RELAX_SOLVE_RECV;
}

/* Compute and apply relaxation parameter */
int arkRelaxSolve(ARKodeMem ark_mem, ARKodeRelaxMem relax_mem,
                  sunrealtype* relax_val_out)
{
  int retval;

  /* Get the change in entropy (uses temp vectors 1 and 2) */
  retval = relax_mem->delta_e_fn(ark_mem,
                                 relax_mem->relax_jac_fn,
                                 &(relax_mem->num_relax_jac_evals),
                                 &(relax_mem->delta_e));
  if (retval) return retval;

  /* Get the change in state (delta_y = tempv1) */
  retval = relax_mem->delta_y_fn(ark_mem, ark_mem->tempv1);
  if (retval) return retval;

  /* Store the current relaxation function value */
  retval = relax_mem->relax_fn(ark_mem->yn, &(relax_mem->e_old),
                               ark_mem->user_data);
  relax_mem->num_relax_fn_evals++;
  if (retval < 0) return ARK_RELAX_FUNC_FAIL;
  if (retval > 0) return ARK_RELAX_FUNC_RECV;

  /* Initial guess for relaxation parameter */
  /* ADD OPTION TO USE GAMMA FROM LAST STEP */
  relax_mem->relax_param = ONE;

  switch(relax_mem->solver)
  {
  case(ARK_RELAX_NEWTON):
    retval = arkRelaxNewtonSolve(ark_mem);
    break;
  case(ARK_RELAX_FIXEDPOINT):
    retval = arkRelaxFixedPointSolve(ark_mem);
    break;
  default:
    return ARK_ILL_INPUT;
    break;
  }

  /* Check for solver failure */
  if (retval)
  {
    relax_mem->nls_fails++;
    return retval;
  }

  /* Check for bad relaxation value */
  if (ark_mem->relax_mem->relax_param < relax_mem->lower_bound ||
      ark_mem->relax_mem->relax_param > relax_mem->upper_bound)
    return ARK_RELAX_SOLVE_RECV;

  /* Return relaxation value */
  *relax_val_out = ark_mem->relax_mem->relax_param;

  return ARK_SUCCESS;
}

/* =============================================================================
 * User Functions
 * ===========================================================================*/

/* -----------------------------------------------------------------------------
 * Set functions
 * ---------------------------------------------------------------------------*/

int arkRelaxSetEtaFail(void* arkode_mem, sunrealtype eta_fail)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxSetEtaFail", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  if (eta_fail < SUN_RCONST(1.0)) relax_mem->eta_fail = eta_fail;
  else relax_mem->eta_fail = ARK_RELAX_DEFAULT_ETA_FAIL;

  return ARK_SUCCESS;
}

int arkRelaxSetLowerBound(void* arkode_mem, sunrealtype lower)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxSetLowerBound", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  if (lower < SUN_RCONST(1.0)) relax_mem->lower_bound = lower;
  else relax_mem->lower_bound = ARK_RELAX_DEFAULT_LOWER_BOUND;

  return ARK_SUCCESS;
}

int arkRelaxSetMaxFails(void* arkode_mem, int max_fails)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxSetMaxFails", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  if (max_fails > 0) relax_mem->max_fails = max_fails;
  else relax_mem->max_fails = ARK_RELAX_DEFAULT_MAX_FAILS;

  return ARK_SUCCESS;
}

int arkRelaxSetMaxIters(void* arkode_mem, int max_iters)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxSetMaxIters", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  if (max_iters > 0) relax_mem->max_iters = max_iters;
  else relax_mem->max_iters = ARK_RELAX_DEFAULT_MAX_ITERS;

  return ARK_SUCCESS;
}

int arkRelaxSetSolver(void* arkode_mem, ARKRelaxSolver solver)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxSetSolver", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  relax_mem->solver = solver;

  return ARK_SUCCESS;
}

int arkRelaxSetResTol(void* arkode_mem, sunrealtype res_tol)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxSetResTol", &ark_mem, &relax_mem);
  if (retval) return retval;

  if (res_tol > SUN_RCONST(0.0)) relax_mem->res_tol = res_tol;
  else relax_mem->res_tol = ARK_RELAX_DEFAULT_RES_TOL;

  return ARK_SUCCESS;
}

int arkRelaxSetUpperBound(void* arkode_mem, sunrealtype upper)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxSetUpperBound", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  if (upper > SUN_RCONST(1.0)) relax_mem->upper_bound = upper;
  else relax_mem->upper_bound = ARK_RELAX_DEFAULT_UPPER_BOUND;

  return ARK_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Get functions
 * ---------------------------------------------------------------------------*/

int arkRelaxGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxGetNumRelaxFnEvals", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  *r_evals = relax_mem->num_relax_fn_evals;

  return ARK_SUCCESS;
}

int arkRelaxGetNumRelaxJacEvals(void* arkode_mem, long int* J_evals)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxGetNumRelaxJacEvals", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  *J_evals = relax_mem->num_relax_jac_evals;

  return ARK_SUCCESS;
}

int arkRelaxGetNumRelaxFails(void* arkode_mem, long int* relax_fails)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxGetNumRelaxFails", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  *relax_fails = relax_mem->num_fails;

  return ARK_SUCCESS;
}

int arkRelaxGetNumRelaxSolveFails(void* arkode_mem, long int* fails)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxGetNumRelaxSolveFails",
                             &ark_mem, &relax_mem);
  if (retval) return retval;

  *fails = relax_mem->nls_fails;

  return ARK_SUCCESS;
}

int arkRelaxGetNumRelaxSolveIters(void* arkode_mem, long int* iters)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxGetNumRelaxSolveIters",
                             &ark_mem, &relax_mem);
  if (retval) return retval;

  *iters = relax_mem->nls_iters;

  return ARK_SUCCESS;
}

int arkRelaxPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxPrintAllStats", &ark_mem,
                             &relax_mem);
  if (retval) return retval;

  switch(fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "Relax fn evals               = %ld\n",
            relax_mem->num_relax_fn_evals);
    fprintf(outfile, "Relax Jac evals              = %ld\n",
            relax_mem->num_relax_jac_evals);
    fprintf(outfile, "Relax fails                  = %ld\n",
            relax_mem->num_fails);
    fprintf(outfile, "Relax NLS iters              = %ld\n",
            relax_mem->nls_iters);
    fprintf(outfile, "Relax NLS fails              = %ld\n",
            relax_mem->nls_fails);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, ",Relax fn evals,%ld", relax_mem->num_relax_fn_evals);
    fprintf(outfile, ",Relax Jac evals,%ld", relax_mem->num_relax_jac_evals);
    fprintf(outfile, ",Relax fails,%ld", relax_mem->num_fails);
    fprintf(outfile, ",Relax NLS iters,%ld", relax_mem->nls_iters);
    fprintf(outfile, ",Relax NLS fails,%ld", relax_mem->nls_fails);
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRelaxPrintAllStats",
                    "Invalid formatting option.");
    return ARK_ILL_INPUT;
  }

  return ARK_SUCCESS;
}

/* =============================================================================
 * Driver and Stepper Functions
 * ===========================================================================*/

/* Constructor called by stepper */
int arkRelaxCreate(void* arkode_mem, ARKRelaxFn relax_fn,
                   ARKRelaxJacFn relax_jac_fn, ARKRelaxDeltaYFn delta_y_fn,
                   ARKRelaxDeltaEFn delta_e_fn, ARKRelaxGetOrderFn get_order_fn)
{
  ARKodeMem ark_mem;

  /* Check inputs */
  if (!arkode_mem)
  {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", "arkRelaxCreate",
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Disable relaxation if both user inputs are NULL */
  if (!relax_fn && !relax_jac_fn)
  {
    ark_mem->relax_enabled = SUNFALSE;
    return ARK_SUCCESS;
  }

  /* Ensure both the relaxation function and Jacobian are provided */
  if (!relax_fn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRelaxCreate",
                    "The relaxation function is NULL.");
    return ARK_ILL_INPUT;
  }

  if (!relax_jac_fn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRelaxCreate",
                    "The relaxation Jacobian function is NULL.");
    return ARK_ILL_INPUT;
  }

  /* Ensure stepper supplied inputs are provided */
  if (!delta_y_fn || !delta_e_fn || !get_order_fn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRelaxCreate",
                    "The Delta y, Delta e, or get order function is NULL.");
    return ARK_ILL_INPUT;
  }

  /* Allocate and initialize relaxation memory structure */
  if (!(ark_mem->relax_mem))
  {
    ark_mem->relax_mem = (ARKodeRelaxMem)malloc(sizeof(*(ark_mem->relax_mem)));
    if (!(ark_mem->relax_mem)) return ARK_MEM_FAIL;
    memset(ark_mem->relax_mem, 0, sizeof(struct ARKodeRelaxMemRec));

    /* Initialize other values */
    ark_mem->relax_mem->max_fails   = ARK_RELAX_DEFAULT_MAX_FAILS;
    ark_mem->relax_mem->lower_bound = ARK_RELAX_DEFAULT_LOWER_BOUND;
    ark_mem->relax_mem->upper_bound = ARK_RELAX_DEFAULT_UPPER_BOUND;
    ark_mem->relax_mem->eta_fail    = ARK_RELAX_DEFAULT_ETA_FAIL;
    ark_mem->relax_mem->solver      = ARK_RELAX_NEWTON;
    ark_mem->relax_mem->res_tol     = ARK_RELAX_DEFAULT_RES_TOL;
    ark_mem->relax_mem->max_iters   = ARK_RELAX_DEFAULT_MAX_ITERS;
  }

  /* Set function pointers */
  ark_mem->relax_mem->relax_fn     = relax_fn;
  ark_mem->relax_mem->relax_jac_fn = relax_jac_fn;
  ark_mem->relax_mem->delta_y_fn   = delta_y_fn;
  ark_mem->relax_mem->delta_e_fn   = delta_e_fn;
  ark_mem->relax_mem->get_order_fn = get_order_fn;

  /* Enable relaxation */
  ark_mem->relax_enabled = SUNTRUE;

  return ARK_SUCCESS;
}

/* Destructor called by driver */
int arkRelaxDestroy(ARKodeRelaxMem relax_mem)
{
  if (!relax_mem) return ARK_SUCCESS;

  /* Free structure */
  free(relax_mem);

  return ARK_SUCCESS;
}

/* Compute and apply relaxation, called by driver */
int arkRelax(ARKodeMem ark_mem, int* relax_fails, realtype* dsm_inout,
             int* nflag_out)
{
  int retval;
  sunrealtype relax_val;
  ARKodeRelaxMem relax_mem = ark_mem->relax_mem;

  /* Get the relaxation memory structure */
  if (!relax_mem)
  {
    arkProcessError(ark_mem, ARK_RELAX_MEM_NULL, "ARKODE",
                    "arkRelax", MSG_RELAX_MEM_NULL);
    return ARK_RELAX_MEM_NULL;
  }

  /* Compute the relaxation parameter */
  retval = arkRelaxSolve(ark_mem, relax_mem, &relax_val);
  if (retval < 0) return retval;
  if (retval > 0)
  {
    /* Update failure counts */
    relax_mem->num_fails++;
    (*relax_fails)++;

    /* Check for max fails in a step */
    if (*relax_fails == relax_mem->max_fails)
      return ARK_RELAX_FAIL;

    /* Return with an error if |h| == hmin */
    if (SUNRabs(ark_mem->h) <= ark_mem->hmin * ONEPSM)
      return ARK_RELAX_FAIL;

    /* Cut step size and try again */
    ark_mem->eta = relax_mem->eta_fail;

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO,
                       "ARKODE::arkStep_TakeStep_Z", "relaxation",
                       "relaxation failed");
#endif

    return TRY_AGAIN;
  }

  /* Update step size and error estimate */
  ark_mem->h *= relax_val;
  *dsm_inout *= SUNRpowerI(relax_val, relax_mem->get_order_fn(ark_mem));

  /* Relax solution */
  N_VLinearSum(relax_val, ark_mem->ycur, (ONE - relax_val), ark_mem->yn,
               ark_mem->ycur);

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO,
                     "ARKODE::arkStep_TakeStep_Z", "relaxation",
                     "relaxation parameter = %"RSYM", relaxed h = %"RSYM
                     ", relaxed error = %"RSYM,
                     relax_val, ark_mem->h, *dsm_inout);
#endif

  return ARK_SUCCESS;
}

/* =============================================================================
 * EOF
 * ===========================================================================*/
