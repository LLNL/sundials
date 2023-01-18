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

/* Evaluate the relaxation residual function */
static int arkRelaxResidual(sunrealtype* relax_vals, N_Vector* y_relax,
                            sunrealtype* res_vals, ARKodeMem ark_mem)
{
  int i, retval;
  int num_relax_fn            = ark_mem->relax_mem->num_relax_fn;
  long int num_relax_fn_evals = ark_mem->relax_mem->num_relax_fn_evals;
  sunrealtype* e_old          = ark_mem->relax_mem->e_old;
  sunrealtype* delta_e        = ark_mem->relax_mem->delta_e;
  void* user_data             = ark_mem->user_data;

  /* Evaluate entropy functions */
  retval = ark_mem->relax_mem->relax_fn(y_relax, res_vals, user_data);
  num_relax_fn_evals++;
  if (retval < 0) return ARK_RELAX_FUNC_FAIL;
  if (retval > 0) return ARK_RELAX_FUNC_RECV;

  /* Compute relaxation residual */
  for (i = 0; i < num_relax_fn; ++i)
    res_vals[i] = res_vals[i] - e_old[i] - relax_vals[i] * delta_e[i];

  return ARK_SUCCESS;
}

/* Evaluate the relaxation Jacobian function */
static int arkRelaxJacobian(sunrealtype* relax_vals, N_Vector* y_relax,
                            sunrealtype* jac_vals, ARKodeMem ark_mem)
{
  int i, retval;
  int num_relax_fn                = ark_mem->relax_mem->num_relax_fn;
  long int num_relax_jac_fn_evals = ark_mem->relax_mem->num_relax_jac_evals;
  N_Vector delta_y                = ark_mem->relax_mem->delta_y;
  N_Vector* J_vecs                = ark_mem->relax_mem->J_vecs;
  sunrealtype* delta_e            = ark_mem->relax_mem->delta_e;
  void* user_data                 = ark_mem->user_data;

  /* Evaluate Jacobian of entropy functions */
  retval = ark_mem->relax_mem->relax_jac_fn(y_relax, J_vecs, user_data);
  num_relax_jac_fn_evals++;
  if (retval < 0) return ARK_RELAX_JAC_FAIL;
  if (retval > 0) return ARK_RELAX_JAC_RECV;

  /* Compute relaxation residual Jacobian */
  retval = N_VDotProdMulti(num_relax_fn, delta_y, J_vecs, jac_vals);
  if (retval) return retval;

  for (i = 0; i < num_relax_fn; ++i) jac_vals[i] -= delta_e[i];

  return ARK_SUCCESS;
}

/* Solve the relaxation residual equation using Newton's method */
static int arkRelaxNewtonSolve(ARKodeMem ark_mem)
{
  int i, j, retval;
  sunrealtype max_res;
  ARKodeRelaxMem relax_mem = ark_mem->relax_mem;

  for (i = 0; i < ark_mem->relax_mem->max_iters; i++)
  {
    /* y_relax = y_n + r * delta_y */
    for (j = 0; j < ark_mem->relax_mem->num_relax_fn; ++j)
      N_VLinearSum(ONE, ark_mem->yn, relax_mem->relax_vals[j],
                   ark_mem->relax_mem->delta_y, ark_mem->relax_mem->y_relax[j]);

    /* Compute the current residual */
    retval = arkRelaxResidual(relax_mem->relax_vals, relax_mem->y_relax,
                              relax_mem->res_vals, ark_mem);
    if (retval) return retval;

    /* Check for convergence of all values */
    max_res = SUNRabs(relax_mem->res_vals[0]);
    for (j = 1; j < ark_mem->relax_mem->num_relax_fn; ++j)
      max_res = SUNMAX(max_res, relax_mem->res_vals[j]);

    if (max_res < relax_mem->tol) return ARK_SUCCESS;

    /* Update iteration count */
    relax_mem->nls_iters++;

    /* Compute Jacobian and update */
    retval = arkRelaxJacobian(relax_mem->relax_vals, relax_mem->y_relax,
                              relax_mem->jac_vals, ark_mem);
    if (retval) return retval;

    for (j = 0; j < ark_mem->relax_mem->num_relax_fn; ++j)
      relax_mem->relax_vals[j] -= relax_mem->res_vals[j] / relax_mem->jac_vals[j];
  }

  return ARK_RELAX_SOLVE_RECV;
}

/* Solve the relaxation residual equation using Newton's method */
static int arkRelaxFixedPointSolve(ARKodeMem ark_mem)
{
  int i, j, retval;
  sunrealtype max_res;
  ARKodeRelaxMem relax_mem = ark_mem->relax_mem;

  for (i = 0; i < ark_mem->relax_mem->max_iters; i++)
  {
    /* y_relax = y_n + r * delta_y */
    for (j = 0; j < ark_mem->relax_mem->num_relax_fn; ++j)
      N_VLinearSum(ONE, ark_mem->yn, relax_mem->relax_vals[j],
                   ark_mem->relax_mem->delta_y, ark_mem->relax_mem->y_relax[j]);

    /* Compute the current residual */
    retval = arkRelaxResidual(relax_mem->relax_vals, relax_mem->y_relax,
                              relax_mem->res_vals, ark_mem);
    if (retval) return retval;

    /* Check for convergence of all values */
    max_res = SUNRabs(relax_mem->res_vals[0]);
    for (j = 1; j < ark_mem->relax_mem->num_relax_fn; ++j)
      max_res = SUNMAX(max_res, relax_mem->res_vals[j]);

    if (max_res < relax_mem->tol) return ARK_SUCCESS;

    /* Update iteration count */
    relax_mem->nls_iters++;

    for (j = 0; j < ark_mem->relax_mem->num_relax_fn; ++j)
      relax_mem->relax_vals[j] -= relax_mem->res_vals[j];
  }

  return ARK_RELAX_SOLVE_RECV;
}

/* Compute and apply relaxation parameter */
int arkRelaxSolve(ARKodeMem ark_mem, ARKodeRelaxMem relax_mem,
                  sunrealtype* relax_val_out)
{
  int i, retval;
  sunrealtype relax_val_max, relax_val_min;

  /* Get the change in entropy (use y_relax and J_vecs as scratch) */
  retval = relax_mem->delta_e_fn(ark_mem, relax_mem->num_relax_fn,
                                 relax_mem->relax_jac_fn, relax_mem->y_relax,
                                 relax_mem->J_vecs,
                                 &(relax_mem->num_relax_jac_evals),
                                 relax_mem->delta_e);
  if (retval) return retval;

  /* Get the change in state */
  retval = relax_mem->delta_y_fn(ark_mem, &(relax_mem->delta_y));
  if (retval) return retval;

  /* Duplicate state to evaluate entropy functions at y_n */
  for (i = 0; i < ark_mem->relax_mem->num_relax_fn; ++i)
    N_VScale(ONE, ark_mem->yn, ark_mem->relax_mem->y_relax[i]);

  retval = relax_mem->relax_fn(ark_mem->relax_mem->y_relax, relax_mem->e_old,
                               ark_mem->user_data);
  relax_mem->num_relax_fn_evals++;
  if (retval < 0) return ARK_RELAX_FUNC_FAIL;
  if (retval > 0) return ARK_RELAX_FUNC_RECV;

  /* Initial guess for relaxation parameter */
  for (i = 0; i < ark_mem->relax_mem->num_relax_fn; ++i)
    relax_mem->relax_vals[i] = ONE;

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

  /* Compute max and min relaxation values */
  relax_val_max = relax_mem->relax_vals[0];
  relax_val_min = relax_mem->relax_vals[0];
  for (i = 1; i < ark_mem->relax_mem->num_relax_fn; ++i)
  {
    relax_val_min = SUNMIN(relax_val_min, relax_mem->relax_vals[i]);
    relax_val_max = SUNMAX(relax_val_max, relax_mem->relax_vals[i]);
  }

  /* Check for bad relaxation value */
  if (relax_val_min < relax_mem->lower_bound ||
      relax_val_max > relax_mem->upper_bound)
    return ARK_RELAX_SOLVE_RECV;

  /* Return min relaxation value */
  *relax_val_out = relax_val_min;

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

int arkRelaxSetTol(void* arkode_mem, sunrealtype tol)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, "arkRelaxSetTol", &ark_mem, &relax_mem);
  if (retval) return retval;

  if (tol > SUN_RCONST(0.0)) relax_mem->tol = tol;
  else relax_mem->tol = ARK_RELAX_DEFAULT_TOL;

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
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, ",Relax fn evals,%ld", relax_mem->num_relax_fn_evals);
    fprintf(outfile, ",Relax Jac evals,%ld", relax_mem->num_relax_jac_evals);
    fprintf(outfile, ",Relax fails,%ld", relax_mem->num_fails);
    fprintf(outfile, ",Relax NLS iters,%ld", relax_mem->nls_iters);
    fprintf(outfile, ",Relax NLS fails,%ld", relax_mem->nls_fails);
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
int arkRelaxCreate(void* arkode_mem, int num_relax_fn, ARKRelaxFn relax_fn,
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

  /* Check for valid inputs */
  if (num_relax_fn < 0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRelaxCreate",
                    "The number of relaxation functions must be non-negative");
    return ARK_ILL_INPUT;
  }

  if (num_relax_fn == 0 && (relax_fn || relax_jac_fn))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRelaxCreate",
                    "The number of relaxation functions is 0 but the relaxation"
                    " or Jacobian function was not NULL.");
    return ARK_ILL_INPUT;
  }

  if (num_relax_fn > 0 && !relax_fn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRelaxCreate",
                    "The relaxation function is NULL.");
    return ARK_ILL_INPUT;
  }

  if (num_relax_fn > 0 && !relax_jac_fn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRelaxCreate",
                    "The relaxation Jacobian function is NULL.");
    return ARK_ILL_INPUT;
  }

  /* Disable relaxation if user inputs are 0 and NULL */
  if (num_relax_fn == 0 && !relax_fn && !relax_jac_fn)
  {
    ark_mem->relax_enabled = SUNFALSE;
    return ARK_SUCCESS;
  }

  /* Check for stepper supplied inputs */
  if (!delta_y_fn || !delta_e_fn || !get_order_fn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRelaxCreate",
                    "The Delta y, Delta e, or get order function is NULL.");
    return ARK_ILL_INPUT;
  }

  /* Check if sufficient memory has already been allocated */
  if (ark_mem->relax_mem)
  {
    if (ark_mem->relax_mem->num_relax_fn_alloc < num_relax_fn)
    {
      arkRelaxDestroy(arkode_mem);
      ark_mem->relax_mem = NULL;
    }
  }

  /* Allocate the relaxation memory structure if necessary */
  if (!(ark_mem->relax_mem))
  {
    ark_mem->relax_mem = (ARKodeRelaxMem)malloc(sizeof(*(ark_mem->relax_mem)));
    if (!(ark_mem->relax_mem)) return ARK_MEM_FAIL;

    /* Zero out relax_mem */
    memset(ark_mem->relax_mem, 0, sizeof(struct ARKodeRelaxMemRec));

    /* Allocate vectors */
    ark_mem->relax_mem->y_relax = N_VCloneVectorArray(num_relax_fn, ark_mem->yn);
    if (!(ark_mem->relax_mem->y_relax))
    {
      arkRelaxDestroy(arkode_mem);
      return ARK_MEM_FAIL;
    }

    ark_mem->relax_mem->J_vecs = N_VCloneVectorArray(num_relax_fn, ark_mem->yn);
    if (!(ark_mem->relax_mem->J_vecs))
    {
      arkRelaxDestroy(arkode_mem);
      return ARK_MEM_FAIL;
    }

    /* Allocate arrays */
    ark_mem->relax_mem->delta_e =
      (sunrealtype*)malloc(num_relax_fn * sizeof(sunrealtype));
    if (!(ark_mem->relax_mem->delta_e))
    {
      arkRelaxDestroy(arkode_mem);
      return ARK_MEM_FAIL;
    }

    ark_mem->relax_mem->e_old =
      (sunrealtype*)malloc(num_relax_fn * sizeof(sunrealtype));
    if (!(ark_mem->relax_mem->e_old))
    {
      arkRelaxDestroy(arkode_mem);
      return ARK_MEM_FAIL;
    }

    ark_mem->relax_mem->res_vals =
      (sunrealtype*)malloc(num_relax_fn * sizeof(sunrealtype));
    if (!(ark_mem->relax_mem->res_vals))
    {
      arkRelaxDestroy(arkode_mem);
      return ARK_MEM_FAIL;
    }

    ark_mem->relax_mem->jac_vals =
      (sunrealtype*)malloc(num_relax_fn * sizeof(sunrealtype));
    if (!(ark_mem->relax_mem->jac_vals))
    {
      arkRelaxDestroy(arkode_mem);
      return ARK_MEM_FAIL;
    }

    ark_mem->relax_mem->relax_vals =
      (sunrealtype*)malloc(num_relax_fn * sizeof(sunrealtype));
    if (!(ark_mem->relax_mem->relax_vals))
    {
      arkRelaxDestroy(arkode_mem);
      return ARK_MEM_FAIL;
    }

    ark_mem->relax_mem->num_relax_fn_alloc = num_relax_fn;
  }

  /* Initialize other values */
  ark_mem->relax_mem->relax_fn     = relax_fn;
  ark_mem->relax_mem->relax_jac_fn = relax_jac_fn;
  ark_mem->relax_mem->delta_y_fn   = delta_y_fn;
  ark_mem->relax_mem->delta_e_fn   = delta_e_fn;
  ark_mem->relax_mem->get_order_fn = get_order_fn;
  ark_mem->relax_mem->num_relax_fn = num_relax_fn;
  ark_mem->relax_mem->max_fails    = ARK_RELAX_DEFAULT_MAX_FAILS;
  ark_mem->relax_mem->lower_bound  = ARK_RELAX_DEFAULT_LOWER_BOUND;
  ark_mem->relax_mem->upper_bound  = ARK_RELAX_DEFAULT_UPPER_BOUND;
  ark_mem->relax_mem->eta_fail     = ARK_RELAX_DEFAULT_ETA_FAIL;
  ark_mem->relax_mem->solver       = ARK_RELAX_NEWTON;
  ark_mem->relax_mem->tol          = ARK_RELAX_DEFAULT_TOL;
  ark_mem->relax_mem->max_iters    = ARK_RELAX_DEFAULT_MAX_ITERS;

  /* Enable relaxation */
  ark_mem->relax_enabled = SUNTRUE;

  return ARK_SUCCESS;
}

/* Destructor called by driver */
int arkRelaxDestroy(ARKodeRelaxMem relax_mem)
{
  if (!relax_mem) return ARK_SUCCESS;

  /* Free vectors */
  if (relax_mem->y_relax)
  {
    N_VDestroyVectorArray(relax_mem->y_relax, relax_mem->num_relax_fn_alloc);
    relax_mem->y_relax = NULL;
  }

  if (relax_mem->y_relax)
  {
    N_VDestroyVectorArray(relax_mem->J_vecs, relax_mem->num_relax_fn_alloc);
    relax_mem->J_vecs = NULL;
  }

  /* Free arrays */
  if (relax_mem->delta_e)
  {
    free(relax_mem->delta_e);
    relax_mem->delta_e = NULL;
  }

  if (relax_mem->e_old)
  {
    free(relax_mem->e_old);
    relax_mem->e_old = NULL;
  }

  if (relax_mem->res_vals)
  {
    free(relax_mem->res_vals);
    relax_mem->res_vals = NULL;
  }

  if (relax_mem->jac_vals)
  {
    free(relax_mem->jac_vals);
    relax_mem->jac_vals = NULL;
  }

  if (relax_mem->relax_vals)
  {
    free(relax_mem->relax_vals);
    relax_mem->relax_vals = NULL;
  }

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
