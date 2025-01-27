/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>
#include "sundials/sundials_types.h"
#include "sundials_macros.h"
#include "sundials_stepper_impl.h"
#include "sundials_utils.h"

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

SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper self, N_Vector y0,
                                    sunrealtype t0, N_Vector sf, sunrealtype tf)
{
  SUNFunctionBegin(self->sunctx);
  self->tf         = tf;
  self->step_idx   = self->final_step_idx;
  self->njeval     = 0;
  self->njpeval    = 0;
  self->njtimesv   = 0;
  self->njptimesv  = 0;
  self->nvtimesj   = 0;
  self->nvtimesjp  = 0;
  self->nrecompute = 0;
  self->nst        = 0;
  SUNStepper_Reset(self->adj_sunstepper, tf, sf, 0);
  SUNStepper_Reset(self->fwd_sunstepper, t0, y0, 0);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper self, sunrealtype tout,
                                    N_Vector sens, sunrealtype* tret)

{
  SUNFunctionBegin(self->sunctx);

  SUNErrCode retcode     = SUN_SUCCESS;
  const sunrealtype zero = SUN_RCONST(0.0);
  const sunrealtype one  = SUN_RCONST(1.0);
  sunrealtype t          = self->tf;
  sunrealtype direction  = (t - tout) > zero ? -one : one;

  self->last_flag = 0;

  while ((direction == -one && t > tout) || (direction == one && t < tout))
  {
    SUNCheckCall(SUNAdjointStepper_OneStep(self, tout, sens, tret));
    if (self->last_flag < 0)
    {
      retcode = SUN_ERR_ADJOINT_STEPPERFAILED;
      break;
    }
    else { t = *tret; }
  }

  return retcode;
}

SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper self, sunrealtype tout,
                                     N_Vector sens, sunrealtype* tret)

{
  SUNFunctionBegin(self->sunctx);
  SUNStepper adj_sunstepper = self->adj_sunstepper;

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = self->tf;
  SUNCheckCall(SUNStepper_OneStep(adj_sunstepper, tout, sens, &t));
  self->last_flag = adj_sunstepper->last_flag;

  self->step_idx--;
  self->nst++;

  if (self->last_flag < 0) { retcode = SUN_ERR_ADJOINT_STEPPERFAILED; }
  else if (self->last_flag > 0)
  {
    retcode = SUN_ERR_ADJOINT_STEPPERINVALIDSTOP;
  }

  *tret = t;

  return retcode;
}

SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper self,
                                          int64_t start_idx, sunrealtype t0,
                                          sunrealtype tf, N_Vector y0)
{
  SUNFunctionBegin(self->sunctx);

  SUNErrCode retcode = SUN_SUCCESS;

  sunrealtype fwd_t      = t0;
  SUNStepper fwd_stepper = self->fwd_sunstepper;
  SUNCheckCall(SUNStepper_Reset(fwd_stepper, t0, y0, start_idx));

  SUNCheckCall(SUNAdjointCheckpointScheme_EnableDense(self->checkpoint_scheme, 1));

  SUNCheckCall(SUNStepper_SetStopTime(fwd_stepper, tf));

  SUNCheckCall(SUNStepper_Evolve(fwd_stepper, tf, y0, &fwd_t));
  self->nrecompute++;

  if (fwd_stepper->last_flag < 0) { retcode = SUN_ERR_ADJOINT_STEPPERFAILED; }
  else if (fwd_stepper->last_flag > 1)
  {
    /* if last_flags is not a successful (0) or tstop (1) return,
        we do not have a way to handle it */
    retcode = SUN_ERR_ADJOINT_STEPPERINVALIDSTOP;
  }

  SUNCheckCall(SUNAdjointCheckpointScheme_EnableDense(self->checkpoint_scheme, 0));

  return retcode;
}

SUNErrCode SUNAdjointStepper_Destroy(SUNAdjointStepper* self_ptr)
{
  SUNAdjointStepper self = *self_ptr;
  SUNStepper_Destroy(&self->fwd_sunstepper);
  SUNStepper_Destroy(&self->adj_sunstepper);
  free(self);
  *self_ptr = NULL;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetJacFn(SUNAdjointStepper self, SUNRhsJacFn JacFn,
                                      SUNMatrix Jac, SUNRhsJacFn JacPFn,
                                      SUNMatrix JacP)
{
  SUNFunctionBegin(self->sunctx);

  self->JacFn  = JacFn;
  self->Jac    = Jac;
  self->JacPFn = JacPFn;
  self->JacP   = JacP;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetJacTimesVecFn(SUNAdjointStepper self,
                                              SUNRhsJacTimesFn Jvp,
                                              SUNRhsJacTimesFn JPvp)
{
  SUNFunctionBegin(self->sunctx);

  self->JvpFn  = Jvp;
  self->JPvpFn = JPvp;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetVecTimesJacFn(SUNAdjointStepper self,
                                              SUNRhsJacTimesFn vJp,
                                              SUNRhsJacTimesFn vJPp)
{
  SUNFunctionBegin(self->sunctx);

  self->vJpFn  = vJp;
  self->vJPpFn = vJPp;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetUserData(SUNAdjointStepper self, void* user_data)
{
  SUNFunctionBegin(self->sunctx);

  self->user_data = user_data;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_PrintAllStats(SUNAdjointStepper self,
                                           FILE* outfile, SUNOutputFormat fmt)
{
  sunfprintf_long(outfile, fmt, SUNFALSE, "Num backwards steps", self->nst);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Num recompute passes",
                  self->nrecompute);
  if (self->JacFn)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac fn evals", self->njeval);
  }
  if (self->JacPFn)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "JacP fn evals", self->njpeval);
  }
  if (self->JvpFn)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times-v evals", self->njtimesv);
  }
  if (self->JPvpFn)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "JacP-times-v evals",
                    self->njptimesv);
  }
  if (self->vJpFn)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "v-times-Jac evals", self->nvtimesj);
  }
  if (self->vJPpFn)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "v-times-Jacp evals",
                    self->nvtimesjp);
  }
  return SUN_SUCCESS;
}
