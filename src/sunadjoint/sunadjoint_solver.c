/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include <sunadjoint/sunadjoint_solver.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sunadjoint/sunadjoint_checkpointscheme.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_stepper.h"
#include "sundials/sundials_types.h"

SUNErrCode SUNAdjointSolver_Create(SUNStepper fwd_stepper, SUNStepper adj_stepper,
                                   int64_t final_step_idx, N_Vector sf,
                                   sunrealtype tf,
                                   SUNAdjointCheckpointScheme checkpoint_scheme,
                                   SUNContext sunctx,
                                   SUNAdjointSolver* adj_solver_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointSolver adj_solver = malloc(sizeof(struct SUNAdjointSolver_));
  SUNAssert(adj_solver, SUN_ERR_MALLOC_FAIL);

  adj_solver->fwd_stepper       = fwd_stepper;
  adj_solver->adj_stepper       = adj_stepper;
  adj_solver->Jac               = NULL;
  adj_solver->JacP              = NULL;
  adj_solver->JacFn             = NULL;
  adj_solver->JacPFn            = NULL;
  adj_solver->JvpFn             = NULL;
  adj_solver->JPvpFn            = NULL;
  adj_solver->vJpFn             = NULL;
  adj_solver->vJPpFn            = NULL;
  adj_solver->checkpoint_scheme = checkpoint_scheme;
  adj_solver->tf                = tf;
  adj_solver->user_data         = NULL;
  adj_solver->sunctx            = sunctx;
  adj_solver->step_idx          = final_step_idx;

  *adj_solver_ptr = adj_solver;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_Solve(SUNAdjointSolver adj_solver, sunrealtype tout,
                                  N_Vector sens, sunrealtype* tret,
                                  int* stop_reason)

{
  SUNFunctionBegin(adj_solver->sunctx);
  SUNStepper adj_stepper = adj_solver->adj_stepper;

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = adj_solver->tf;
  *stop_reason       = 0;
  while (t > tout)
  {
    SUNCheckCall(SUNAdjointSolver_Step(adj_solver, tout, sens, tret, stop_reason));
    if (*stop_reason < 0)
    {
      retcode = SUN_ERR_ADJOINT_STEPPERFAILED;
      break;
    }
    else { t = *tret; }
  }

  return retcode;
}

SUNErrCode SUNAdjointSolver_Step(SUNAdjointSolver adj_solver, sunrealtype tout,
                                 N_Vector sens, sunrealtype* tret,
                                 int* stop_reason)

{
  SUNFunctionBegin(adj_solver->sunctx);
  SUNStepper adj_stepper = adj_solver->adj_stepper;

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = adj_solver->tf;
  *stop_reason       = 0;
  SUNCheckCall(
    SUNStepper_Step(adj_stepper, adj_solver->tf, tout, sens, &t, stop_reason));
  if (*stop_reason < 0) { retcode = SUN_ERR_ADJOINT_STEPPERFAILED; }
  else
  {
    // TODO(CJB): what reasons could this happen, and are they valid?
    // 1==TSTOP_RETURN
    // 2==ROOT_RETURN
    fprintf(stderr, ">>>> HERE, stop_reason = %d\n", *stop_reason);
  }
  adj_solver->step_idx--;
  *tret = t;

  return retcode;
}

SUNErrCode SUNAdjointSolver_SetRecompute(SUNAdjointSolver adj_solver,
                                         int64_t start_idx, int64_t stop_idx,
                                         sunrealtype t0, sunrealtype tf,
                                         N_Vector y0)
{
  SUNFunctionBegin(adj_solver->sunctx);
  SUNStepper adj_stepper = adj_solver->adj_stepper;

  adj_solver->recompute_start_step = start_idx;
  adj_solver->recompute_end_step   = stop_idx;
  adj_solver->recompute_t0         = t0;
  adj_solver->recompute_tf         = tf;
  adj_solver->recompute_y0         = y0;
  adj_solver->recompute_flag       = SUNTRUE;

  fprintf(stderr, ">>> here\n");
  int fwd_stop_reason = 0;
  sunrealtype fwd_t   = adj_solver->recompute_t0;
  SUNCheckCall(SUNStepper_Reset(adj_solver->fwd_stepper, adj_solver->recompute_t0,
                                adj_solver->recompute_y0));

  // TODO(CJB): need to get interval first then restore it after advance
  SUNCheckCall(
    SUNAdjointCheckpointScheme_SetInterval(adj_solver->checkpoint_scheme, 1));

  SUNCheckCall(
    SUNStepper_Advance(adj_solver->fwd_stepper, adj_solver->recompute_t0,
                       adj_solver->recompute_tf, adj_solver->recompute_y0,
                       &fwd_t, &fwd_stop_reason));

  adj_solver->recompute_flag = SUNFALSE;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_Destroy(SUNAdjointSolver* adj_solver_ptr)
{
  SUNAdjointSolver adj_solver = *adj_solver_ptr;
  // SUNAdjointCheckpointScheme_Destroy(adj_solver->checkpoint_scheme);
  SUNStepper_Destroy(&adj_solver->fwd_stepper);
  SUNStepper_Destroy(&adj_solver->adj_stepper);
  free(adj_solver);
  *adj_solver_ptr = NULL;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_SetJacFn(SUNAdjointSolver adj_solver,
                                     SUNJacFn JacFn, SUNMatrix Jac,
                                     SUNJacFn JacPFn, SUNMatrix JacP)
{
  SUNFunctionBegin(adj_solver->sunctx);

  SUNAssert(JacFn && Jac, SUN_ERR_ARG_CORRUPT);
  SUNAssert(!JacPFn || (JacPFn && JacP), SUN_ERR_ARG_CORRUPT);

  adj_solver->JacFn  = JacFn;
  adj_solver->Jac    = Jac;
  adj_solver->JacPFn = JacPFn;
  adj_solver->JacP   = JacP;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_SetJacTimesVecFn(SUNAdjointSolver adj_solver,
                                             SUNJacTimesFn Jvp,
                                             SUNJacTimesFn JPvp)
{
  SUNFunctionBegin(adj_solver->sunctx);

  SUNAssert(Jvp, SUN_ERR_ARG_CORRUPT);
  SUNAssert(!JPvp || (JPvp && Jvp), SUN_ERR_ARG_CORRUPT);

  adj_solver->JvpFn  = Jvp;
  adj_solver->JPvpFn = JPvp;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_SetVecTimesJacFn(SUNAdjointSolver adj_solver,
                                             SUNJacTimesFn vJp,
                                             SUNJacTimesFn vJPp)
{
  SUNFunctionBegin(adj_solver->sunctx);

  SUNAssert(vJp, SUN_ERR_ARG_CORRUPT);
  SUNAssert(!vJPp || (vJp && vJPp), SUN_ERR_ARG_CORRUPT);

  adj_solver->vJpFn  = vJp;
  adj_solver->vJPpFn = vJPp;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_SetUserData(SUNAdjointSolver adj_solver,
                                        void* user_data)
{
  SUNFunctionBegin(adj_solver->sunctx);

  adj_solver->user_data = user_data;

  return SUN_SUCCESS;
}
