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
 * Implementation header file for ARKODE's relaxation (in time) functionality.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_RELAX_IMPL_H
#define _ARKODE_RELAX_IMPL_H

#include <stdarg.h>

#include <arkode/arkode.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode_types_impl.h"

/* -----------------------------------------------------------------------------
 * Relaxation Constants
 * ---------------------------------------------------------------------------*/

#define ARK_RELAX_DEFAULT_MAX_FAILS   10
#define ARK_RELAX_DEFAULT_RES_TOL     (4 * SUN_UNIT_ROUNDOFF)
#define ARK_RELAX_DEFAULT_REL_TOL     (4 * SUN_UNIT_ROUNDOFF)
#define ARK_RELAX_DEFAULT_ABS_TOL     SUN_RCONST(1.0e-14)
#define ARK_RELAX_DEFAULT_MAX_ITERS   10
#define ARK_RELAX_DEFAULT_LOWER_BOUND SUN_RCONST(0.8)
#define ARK_RELAX_DEFAULT_UPPER_BOUND SUN_RCONST(1.2)
#define ARK_RELAX_DEFAULT_ETA_FAIL    SUN_RCONST(0.25)

/* -----------------------------------------------------------------------------
 * Relaxation Private Return Values (see arkode/arkode.h for public values)
 * ---------------------------------------------------------------------------*/

#define ARK_RELAX_FUNC_RECV  1;
#define ARK_RELAX_JAC_RECV   2;
#define ARK_RELAX_SOLVE_RECV 3;

/* -----------------------------------------------------------------------------
 * Stepper Supplied Relaxation Functions
 * ---------------------------------------------------------------------------*/

/* Compute the change in state for the current step y_new = y_old + delta_y */
typedef int (*ARKRelaxDeltaYFn)(ARKodeMem ark_mem, N_Vector delta_y);

/* Compute the estimated change in entropy for this step delta_e */
typedef int (*ARKRelaxDeltaEFn)(ARKodeMem ark_mem, ARKRelaxJacFn relax_jac_fn,
                                long int* evals_out, sunrealtype* delta_e_out);

/* Get the method order */
typedef int (*ARKRelaxGetOrderFn)(ARKodeMem ark_mem);

/* -----------------------------------------------------------------------------
 * Relaxation Data Structure
 * ---------------------------------------------------------------------------*/

struct ARKodeRelaxMemRec
{
  /* user-supplied and stepper supplied functions */
  ARKRelaxFn relax_fn;             /* user relaxation function ("entropy") */
  ARKRelaxJacFn relax_jac_fn;      /* user relaxation Jacobian             */
  ARKRelaxDeltaYFn delta_y_fn;     /* get delta y from stepper             */
  ARKRelaxDeltaEFn delta_e_fn;     /* get delta entropy from stepper       */
  ARKRelaxGetOrderFn get_order_fn; /* get the method order                 */

  /* relaxation variables */
  int max_fails;                /* max allowed relax fails in a step   */
  long int num_relax_fn_evals;  /* counter for total function evals    */
  long int num_relax_jac_evals; /* counter for total jacobian evals    */
  long int num_fails;           /* counter for total relaxation fails  */
  sunrealtype e_old;            /* entropy at start of step y(t_{n-1}) */
  sunrealtype delta_e;          /* change in entropy                   */
  sunrealtype res;              /* relaxation residual value           */
  sunrealtype jac;              /* relaxation Jacobian value           */
  sunrealtype relax_param;      /* current relaxation parameter value  */
  sunrealtype relax_param_prev; /* previous relaxation parameter value */
  sunrealtype lower_bound;      /* smallest allowed relaxation value   */
  sunrealtype upper_bound;      /* largest allowed relaxation value    */
  sunrealtype eta_fail;         /* failed relaxation step size factor  */

  /* nonlinear solver settings */
  ARKRelaxSolver solver;      /* choice of relaxation solver          */
  sunrealtype res_tol;        /* nonlinear residual solve tolerance   */
  sunrealtype rel_tol;        /* nonlinear iterate relative tolerance */
  sunrealtype abs_tol;        /* nonlinear iterate absolute tolerance */
  int max_iters;              /* nonlinear solve max iterations       */
  long int nls_iters;         /* total nonlinear iterations           */
  long int nls_fails;         /* number of nonlinear solver fails     */
  long int bound_fails;       /* number of relax param bound fails    */
};

/* -----------------------------------------------------------------------------
 * Relaxation Functions
 * ---------------------------------------------------------------------------*/

/* Driver and Stepper Functions */
int arkRelaxCreate(void* arkode_mem, ARKRelaxFn relax_fn,
                   ARKRelaxJacFn relax_jac_fn, ARKRelaxDeltaYFn delta_y_fn,
                   ARKRelaxDeltaEFn delta_e_fn,
                   ARKRelaxGetOrderFn get_order_fn);
int arkRelaxDestroy(ARKodeRelaxMem relax_mem);
int arkRelax(ARKodeMem ark_mem, int* relax_fails, sunrealtype* dsm_inout,
             int* nflag_out);

/* User Functions */
int arkRelaxSetEtaFail(void* arkode_mem, sunrealtype eta_fail);
int arkRelaxSetLowerBound(void* arkode_mem, sunrealtype lower);
int arkRelaxSetMaxFails(void* arkode_mem, int max_fails);
int arkRelaxSetMaxIters(void* arkode_mem, int max_iters);
int arkRelaxSetSolver(void* arkode_mem, ARKRelaxSolver solver);
int arkRelaxSetResTol(void* arkode_mem, sunrealtype res_tol);
int arkRelaxSetTol(void* arkode_mem, sunrealtype rel_tol, sunrealtype abs_tol);
int arkRelaxSetUpperBound(void* arkode_mem, sunrealtype upper);

int arkRelaxGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals);
int arkRelaxGetNumRelaxJacEvals(void* arkode_mem, long int* j_evals);
int arkRelaxGetNumRelaxFails(void* arkode_mem, long int* relax_fails);
int arkRelaxGetNumRelaxBoundFails(void* arkode_mem, long int* fails);
int arkRelaxGetNumRelaxSolveFails(void* arkode_mem, long int* fails);
int arkRelaxGetNumRelaxSolveIters(void* arkode_mem, long int* iters);

int arkRelaxPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt);

/* -----------------------------------------------------------------------------
 * Error Messages
 * ---------------------------------------------------------------------------*/

#define MSG_RELAX_MEM_NULL "Relaxation memory is NULL."

#endif
