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

#include <sunadjoint/sunadjoint_stepper.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sunadjoint/sunadjoint_checkpointscheme.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_stepper.h"
#include "sundials/sundials_types.h"

SUNErrCode SUNAdjointStepper_Create(
  SUNStepper fwd_stepper, SUNStepper adj_stepper, int64_t final_step_idx,
  N_Vector sf, sunrealtype tf, SUNAdjointCheckpointScheme checkpoint_scheme,
  SUNContext sunctx, SUNAdjointStepper* adj_solver_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointStepper adj_solver = malloc(sizeof(struct SUNAdjointStepper_));
  SUNAssert(adj_solver, SUN_ERR_MALLOC_FAIL);

  adj_solver->fwd_stepper       = fwd_stepper;
  adj_solver->adj_stepper       = adj_stepper;
  adj_solver->checkpoint_scheme = checkpoint_scheme;
  adj_solver->Jac               = NULL;
  adj_solver->JacP              = NULL;

  adj_solver->JacFn  = NULL;
  adj_solver->JacPFn = NULL;
  adj_solver->JvpFn  = NULL;
  adj_solver->JPvpFn = NULL;
  adj_solver->vJpFn  = NULL;
  adj_solver->vJPpFn = NULL;

  adj_solver->tf       = tf;
  adj_solver->step_idx = final_step_idx;

  adj_solver->njeval     = 0;
  adj_solver->njpeval    = 0;
  adj_solver->njtimesv   = 0;
  adj_solver->njptimesv  = 0;
  adj_solver->nvtimesj   = 0;
  adj_solver->nvtimesjp  = 0;
  adj_solver->nrecompute = 0;

  adj_solver->user_data = NULL;
  adj_solver->sunctx    = sunctx;

  *adj_solver_ptr = adj_solver;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_Solve(SUNAdjointStepper adj_solver,
                                   sunrealtype tout, N_Vector sens,
                                   sunrealtype* tret, int* stop_reason)

{
  SUNFunctionBegin(adj_solver->sunctx);
  SUNStepper adj_stepper = adj_solver->adj_stepper;

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = adj_solver->tf;
  *stop_reason       = 0;
  while (t > tout)
  {
    SUNCheckCall(
      SUNAdjointStepper_Step(adj_solver, tout, sens, tret, stop_reason));
    if (*stop_reason < 0)
    {
      retcode = SUN_ERR_ADJOINT_STEPPERFAILED;
      break;
    }
    else { t = *tret; }
  }

  return retcode;
}

SUNErrCode SUNAdjointStepper_Step(SUNAdjointStepper adj_solver,
                                  sunrealtype tout, N_Vector sens,
                                  sunrealtype* tret, int* stop_reason)

{
  SUNFunctionBegin(adj_solver->sunctx);
  SUNStepper adj_stepper = adj_solver->adj_stepper;

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = adj_solver->tf;
  *stop_reason       = 0;
  SUNCheckCall(SUNStepper_OneStep(adj_stepper, adj_solver->tf, tout, sens, &t,
                                  stop_reason));
  if (*stop_reason < 0) { retcode = SUN_ERR_ADJOINT_STEPPERFAILED; }
  else if (*stop_reason > 0) { retcode = SUN_ERR_ADJOINT_STEPPERINVALIDSTOP; }
  adj_solver->step_idx--;
  *tret = t;

  return retcode;
}

SUNErrCode SUNAdjointStepper_SetRecompute(SUNAdjointStepper adj_solver,
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

  int fwd_stop_reason = 0;
  sunrealtype fwd_t   = adj_solver->recompute_t0;
  SUNCheckCall(SUNStepper_Reset(adj_solver->fwd_stepper, adj_solver->recompute_t0,
                                adj_solver->recompute_y0));

  SUNCheckCall(
    SUNAdjointCheckpointScheme_EnableDense(adj_solver->checkpoint_scheme, 1));

  SUNCheckCall(
    SUNStepper_Advance(adj_solver->fwd_stepper, adj_solver->recompute_t0,
                       adj_solver->recompute_tf, adj_solver->recompute_y0,
                       &fwd_t, &fwd_stop_reason));
  adj_solver->nrecompute++;

  SUNCheckCall(
    SUNAdjointCheckpointScheme_EnableDense(adj_solver->checkpoint_scheme, 0));

  adj_solver->recompute_flag = SUNFALSE;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_Destroy(SUNAdjointStepper* adj_solver_ptr)
{
  SUNAdjointStepper adj_solver = *adj_solver_ptr;
  // SUNAdjointCheckpointScheme_Destroy(adj_solver->checkpoint_scheme);
  SUNStepper_Destroy(&adj_solver->fwd_stepper);
  SUNStepper_Destroy(&adj_solver->adj_stepper);
  free(adj_solver);
  *adj_solver_ptr = NULL;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetJacFn(SUNAdjointStepper adj_solver,
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

SUNErrCode SUNAdjointStepper_SetJacTimesVecFn(SUNAdjointStepper adj_solver,
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

SUNErrCode SUNAdjointStepper_SetVecTimesJacFn(SUNAdjointStepper adj_solver,
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

SUNErrCode SUNAdjointStepper_SetUserData(SUNAdjointStepper adj_solver,
                                         void* user_data)
{
  SUNFunctionBegin(adj_solver->sunctx);

  adj_solver->user_data = user_data;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_PrintAllStats(SUNAdjointStepper adj_solver,
                                           FILE* outfile, SUNOutputFormat fmt)
{
  switch (fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "Recomputed steps   = %llu\n", adj_solver->nrecompute);
    if (adj_solver->JacFn)
    {
      fprintf(outfile, "Jac fn evals       = %llu\n", adj_solver->njeval);
    }
    if (adj_solver->JacPFn)
    {
      fprintf(outfile, "JacP fn evals      = %llu\n", adj_solver->njpeval);
    }
    if (adj_solver->JvpFn)
    {
      fprintf(outfile, "Jac-times-v evals  = %llu\n", adj_solver->njtimesv);
    }
    if (adj_solver->JPvpFn)
    {
      fprintf(outfile, "JacP-times-v evals = %llu\n", adj_solver->njptimesv);
    }
    if (adj_solver->vJpFn)
    {
      fprintf(outfile, "v-times-Jac evals  = %llu\n", adj_solver->nvtimesj);
    }
    if (adj_solver->vJPpFn)
    {
      fprintf(outfile, "v-times-Jacp evals = %llu\n", adj_solver->nvtimesjp);
    }
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, "Recomputed steps,%llu", adj_solver->nrecompute);
    if (adj_solver->JacFn)
    {
      fprintf(outfile, ",Jac fn evals,%llu", adj_solver->njeval);
    }
    if (adj_solver->JacPFn)
    {
      fprintf(outfile, ",JacP fn evals,%llu", adj_solver->njpeval);
    }
    if (adj_solver->JvpFn)
    {
      fprintf(outfile, ",Jac-times-v evals,%llu", adj_solver->njtimesv);
    }
    if (adj_solver->JPvpFn)
    {
      fprintf(outfile, ",JacP-times-v evals,%llu", adj_solver->njptimesv);
    }
    if (adj_solver->vJpFn)
    {
      fprintf(outfile, ",v-times-Jac evals,%llu", adj_solver->nvtimesj);
    }
    if (adj_solver->vJPpFn)
    {
      fprintf(outfile, ",v-times-Jacp evals,%llu", adj_solver->nvtimesjp);
    }

    break;
  default: return SUN_ERR_ARG_INCOMPATIBLE;
  }

  return SUN_SUCCESS;
}
