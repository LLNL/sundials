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

#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sunadjoint/sunadjoint_stepper.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>
#include "sundials/sundials_types.h"
#include "sundials_macros.h"
#include "sundials_stepper_impl.h"

SUNErrCode SUNAdjointStepper_Create(
  SUNStepper fwd_sunstepper, SUNStepper adj_sunstepper, int64_t final_step_idx,
  SUNDIALS_MAYBE_UNUSED N_Vector sf, sunrealtype tf,
  SUNAdjointCheckpointScheme checkpoint_scheme, SUNContext sunctx,
  SUNAdjointStepper* adj_stepper_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointStepper adj_stepper = malloc(sizeof(struct SUNAdjointStepper_));
  SUNAssert(adj_stepper, SUN_ERR_MALLOC_FAIL);

  adj_stepper->fwd_sunstepper    = fwd_sunstepper;
  adj_stepper->adj_sunstepper    = adj_sunstepper;
  adj_stepper->checkpoint_scheme = checkpoint_scheme;
  adj_stepper->Jac               = NULL;
  adj_stepper->JacP              = NULL;

  adj_stepper->JacFn  = NULL;
  adj_stepper->JacPFn = NULL;
  adj_stepper->JvpFn  = NULL;
  adj_stepper->JPvpFn = NULL;
  adj_stepper->vJpFn  = NULL;
  adj_stepper->vJPpFn = NULL;

  adj_stepper->tf             = tf;
  adj_stepper->step_idx       = final_step_idx;
  adj_stepper->final_step_idx = final_step_idx;
  adj_stepper->nst            = 0;

  adj_stepper->njeval     = 0;
  adj_stepper->njpeval    = 0;
  adj_stepper->njtimesv   = 0;
  adj_stepper->njptimesv  = 0;
  adj_stepper->nvtimesj   = 0;
  adj_stepper->nvtimesjp  = 0;
  adj_stepper->nrecompute = 0;

  adj_stepper->user_data = NULL;
  adj_stepper->sunctx    = sunctx;

  *adj_stepper_ptr = adj_stepper;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper adj, N_Vector y0,
                                    sunrealtype t0, N_Vector sf, sunrealtype tf)
{
  SUNFunctionBegin(adj->sunctx);
  adj->tf         = tf;
  adj->step_idx   = adj->final_step_idx;
  adj->njeval     = 0;
  adj->njpeval    = 0;
  adj->njtimesv   = 0;
  adj->njptimesv  = 0;
  adj->nvtimesj   = 0;
  adj->nvtimesjp  = 0;
  adj->nrecompute = 0;
  adj->nst        = 0;
  SUNStepper_Reset(adj->adj_sunstepper, tf, sf, 0);
  SUNStepper_Reset(adj->fwd_sunstepper, t0, y0, 0);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper adj_stepper,
                                    sunrealtype tout, N_Vector sens,
                                    sunrealtype* tret, int* stop_reason)

{
  SUNFunctionBegin(adj_stepper->sunctx);

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = adj_stepper->tf;
  *stop_reason       = 0;
  while (t > tout)
  {
    SUNCheckCall(
      SUNAdjointStepper_OneStep(adj_stepper, tout, sens, tret, stop_reason));
    if (*stop_reason < 0)
    {
      retcode = SUN_ERR_ADJOINT_STEPPERFAILED;
      break;
    }
    else { t = *tret; }
  }

  return retcode;
}

SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper adj_stepper,
                                     sunrealtype tout, N_Vector sens,
                                     sunrealtype* tret, int* stop_reason)

{
  SUNFunctionBegin(adj_stepper->sunctx);
  SUNStepper adj_sunstepper = adj_stepper->adj_sunstepper;

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = adj_stepper->tf;
  SUNCheckCall(
    SUNStepper_OneStep(adj_sunstepper, adj_stepper->tf, tout, sens, &t));
  *stop_reason = adj_sunstepper->last_flag;

  adj_stepper->step_idx--;
  adj_stepper->nst++;

  if (*stop_reason < 0) { retcode = SUN_ERR_ADJOINT_STEPPERFAILED; }
  else if (*stop_reason > 0) { retcode = SUN_ERR_ADJOINT_STEPPERINVALIDSTOP; }

  *tret = t;

  return retcode;
}

SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper adj_stepper,
                                          int64_t start_idx, sunrealtype t0,
                                          sunrealtype tf, N_Vector y0)
{
  SUNFunctionBegin(adj_stepper->sunctx);

  SUNErrCode retcode = SUN_SUCCESS;

  sunrealtype fwd_t      = t0;
  SUNStepper fwd_stepper = adj_stepper->fwd_sunstepper;
  SUNCheckCall(SUNStepper_Reset(fwd_stepper, t0, y0, start_idx));

  SUNCheckCall(
    SUNAdjointCheckpointScheme_EnableDense(adj_stepper->checkpoint_scheme, 1));

  SUNCheckCall(SUNStepper_SetStopTime(fwd_stepper, tf));

  SUNCheckCall(SUNStepper_Evolve(fwd_stepper, t0, tf, y0, &fwd_t));
  adj_stepper->nrecompute++;

  if (fwd_stepper->last_flag < 0) { retcode = SUN_ERR_ADJOINT_STEPPERFAILED; }
  else if (fwd_stepper->last_flag > 1)
  {
    /* if last_flags is not a successful (0) or tstop (1) return,
        we do not have a way to handle it */
    retcode = SUN_ERR_ADJOINT_STEPPERINVALIDSTOP;
  }

  SUNCheckCall(
    SUNAdjointCheckpointScheme_EnableDense(adj_stepper->checkpoint_scheme, 0));

  return retcode;
}

SUNErrCode SUNAdjointStepper_Destroy(SUNAdjointStepper* adj_stepper_ptr)
{
  SUNAdjointStepper adj_stepper = *adj_stepper_ptr;
  SUNStepper_Destroy(&adj_stepper->fwd_sunstepper);
  SUNStepper_Destroy(&adj_stepper->adj_sunstepper);
  free(adj_stepper);
  *adj_stepper_ptr = NULL;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetJacFn(SUNAdjointStepper adj_stepper,
                                      SUNJacFn JacFn, SUNMatrix Jac,
                                      SUNJacFn JacPFn, SUNMatrix JacP)
{
  SUNFunctionBegin(adj_stepper->sunctx);

  adj_stepper->JacFn  = JacFn;
  adj_stepper->Jac    = Jac;
  adj_stepper->JacPFn = JacPFn;
  adj_stepper->JacP   = JacP;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetJacTimesVecFn(SUNAdjointStepper adj_stepper,
                                              SUNJacTimesFn Jvp,
                                              SUNJacTimesFn JPvp)
{
  SUNFunctionBegin(adj_stepper->sunctx);

  adj_stepper->JvpFn  = Jvp;
  adj_stepper->JPvpFn = JPvp;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetVecTimesJacFn(SUNAdjointStepper adj_stepper,
                                              SUNJacTimesFn vJp,
                                              SUNJacTimesFn vJPp)
{
  SUNFunctionBegin(adj_stepper->sunctx);

  adj_stepper->vJpFn  = vJp;
  adj_stepper->vJPpFn = vJPp;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetUserData(SUNAdjointStepper adj_stepper,
                                         void* user_data)
{
  SUNFunctionBegin(adj_stepper->sunctx);

  adj_stepper->user_data = user_data;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_PrintAllStats(SUNAdjointStepper adj_stepper,
                                           FILE* outfile, SUNOutputFormat fmt)
{
  switch (fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "Num backwards steps    = %llu\n",
            (unsigned long long)adj_stepper->nst);
    fprintf(outfile, "Num recompute passes   = %llu\n",
            (unsigned long long)adj_stepper->nrecompute);
    if (adj_stepper->JacFn)
    {
      fprintf(outfile, "Jac fn evals           = %llu\n",
              (unsigned long long)adj_stepper->njeval);
    }
    if (adj_stepper->JacPFn)
    {
      fprintf(outfile, "JacP fn evals          = %llu\n",
              (unsigned long long)adj_stepper->njpeval);
    }
    if (adj_stepper->JvpFn)
    {
      fprintf(outfile, "Jac-times-v evals      = %llu\n",
              (unsigned long long)adj_stepper->njtimesv);
    }
    if (adj_stepper->JPvpFn)
    {
      fprintf(outfile, "JacP-times-v evals     = %llu\n",
              (unsigned long long)adj_stepper->njptimesv);
    }
    if (adj_stepper->vJpFn)
    {
      fprintf(outfile, "v-times-Jac evals      = %llu\n",
              (unsigned long long)adj_stepper->nvtimesj);
    }
    if (adj_stepper->vJPpFn)
    {
      fprintf(outfile, "v-times-Jacp evals     = %llu\n",
              (unsigned long long)adj_stepper->nvtimesjp);
    }
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, "Num backwards steps,%llu",
            (unsigned long long)adj_stepper->nst);
    fprintf(outfile, "Num recompute passes,%llu",
            (unsigned long long)adj_stepper->nrecompute);
    if (adj_stepper->JacFn)
    {
      fprintf(outfile, ",Jac fn evals,%llu",
              (unsigned long long)adj_stepper->njeval);
    }
    if (adj_stepper->JacPFn)
    {
      fprintf(outfile, ",JacP fn evals,%llu",
              (unsigned long long)adj_stepper->njpeval);
    }
    if (adj_stepper->JvpFn)
    {
      fprintf(outfile, ",Jac-times-v evals,%llu",
              (unsigned long long)adj_stepper->njtimesv);
    }
    if (adj_stepper->JPvpFn)
    {
      fprintf(outfile, ",JacP-times-v evals,%llu",
              (unsigned long long)adj_stepper->njptimesv);
    }
    if (adj_stepper->vJpFn)
    {
      fprintf(outfile, ",v-times-Jac evals,%llu",
              (unsigned long long)adj_stepper->nvtimesj);
    }
    if (adj_stepper->vJPpFn)
    {
      fprintf(outfile, ",v-times-Jacp evals,%llu",
              (unsigned long long)adj_stepper->nvtimesjp);
    }

    break;
  default: return SUN_ERR_ARG_INCOMPATIBLE;
  }

  return SUN_SUCCESS;
}
