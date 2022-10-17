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

#define ARK_RELAX_DEFAULT_TOL         SUN_RCONST(1.0e-14)
#define ARK_RELAX_DEFAULT_MAX_ITERS   100
#define ARK_RELAX_DEFAULT_LOWER_BOUND SUN_RCONST(0.5)
#define ARK_RELAX_DEFAULT_UPPER_BOUND SUN_RCONST(1.5)
#define ARK_RELAX_DEFAULT_ETA_FAIL    SUN_RCONST(0.25)

/* -----------------------------------------------------------------------------
 * Relaxation Private Return Values (see arkode/arkode.h for public values)
 * ---------------------------------------------------------------------------*/

#define ARK_RELAX_FUNC_RECV  1;
#define ARK_RELAX_JAC_RECV   2;
#define ARK_RELAX_SOLVE_RECV -1; /* <<< should be recoverable */

/* -----------------------------------------------------------------------------
 * Stepper Supplied Relaxation Functions
 * ---------------------------------------------------------------------------*/

/* Compute the change in state for the current step y_new = y_old + delta_y */
typedef int (*ARKRelaxDeltaYFn)(ARKodeMem ark_mem, N_Vector* delta_y);

/* Compute the estimated change in entropy for this step delta_e */
typedef int (*ARKRelaxDeltaEFn)(ARKodeMem ark_mem, int num_relax_fn,
                                ARKRelaxJacFn relax_jac_fn,
                                N_Vector* work_space_1, N_Vector* work_space_2,
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

  /* stepper computed quantities used in the residual and Jacobian */
  N_Vector delta_y;     /* change in solution */
  sunrealtype* delta_e; /* change in entropy  */

  /* relaxation variables */
  int num_relax_fn;             /* number of entropy functions         */
  int num_relax_fn_alloc;       /* allocated workspce size             */
  long int num_relax_fn_evals;  /* counter for function evals          */
  long int num_relax_jac_evals; /* counter for jacobian evals          */
  N_Vector* y_relax;            /* relaxed state y_n + relax * delta_y */
  N_Vector* J_vecs;             /* relaxation Jacobian vectors         */
  sunrealtype* e_old;           /* entropy at start of step y(t_{n-1}) */
  sunrealtype* res_vals;        /* relaxation residual values          */
  sunrealtype* jac_vals;        /* relaxation Jacobian values          */
  sunrealtype* relax_vals;      /* relaxation parameter values         */
  sunrealtype lower_bound;      /* smallest allowed relaxation value   */
  sunrealtype upper_bound;      /* largest allowed relaxation value    */
  sunrealtype eta_fail;         /* failed relaxation step size factor  */

  /* nonlinear solver settings */
  ARKRelaxationSolver solver; /* choice of relaxation solver      */
  sunrealtype tol;            /* nonlinear solve tolerance        */
  int max_iters;              /* nonlinear solve max iterations   */
  long int total_iters;       /* total nonlinear iterations       */
  long int total_fails;       /* number of nonlinear solver fails */
};

/* -----------------------------------------------------------------------------
 * Relaxation Functions
 * ---------------------------------------------------------------------------*/

/* Driver and Stepper Functions */
int arkRelaxCreate(void* arkode_mem, int num_relax_fn, ARKRelaxFn relax_fn,
                   ARKRelaxJacFn relax_jac_fn, ARKRelaxDeltaYFn delta_y_fn,
                   ARKRelaxDeltaEFn delta_e_fn,
                   ARKRelaxGetOrderFn get_order_fn);
int arkRelaxDestroy(ARKodeRelaxMem relax_mem);
int arkRelax(ARKodeMem ark_mem, sunrealtype* dsm_inout, int* nflag_out);

/* User Functions */
int arkRelaxSetEtaFail(void* arkode_mem, sunrealtype eta_fail);
int arkRelaxSetLowerBound(void* arkode_mem, sunrealtype lower);
int arkRelaxSetMaxIters(void* arkode_mem, int max_iters);
int arkRelaxSetSolver(void* arkode_mem, ARKRelaxationSolver solver);
int arkRelaxSetTol(void* arkode_mem, sunrealtype tol);
int arkRelaxSetUpperBound(void* arkode_mem, sunrealtype upper);

int arkRelaxGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals);
int arkRelaxGetNumRelaxJacEvals(void* arkode_mem, long int* j_evals);
int arkRelaxGetNumSolveFails(void* arkode_mem, long int* fails);
int arkRelaxGetNumSolveIters(void* arkode_mem, long int* iters);

/* -----------------------------------------------------------------------------
 * Error Messages
 * ---------------------------------------------------------------------------*/

#define MSG_RELAX_MEM_NULL "Relaxation memory is NULL."

#endif
