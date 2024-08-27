

#include <stdlib.h>
#include <string.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sundials/sundials_errors.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"
#include "sundials_stepper_impl.h"

SUNErrCode SUNStepper_Create(SUNContext sunctx, SUNStepper* stepper_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNStepper stepper = NULL;
  stepper            = malloc(sizeof(*stepper));
  SUNAssert(stepper, SUN_ERR_MALLOC_FAIL);

  stepper->content            = NULL;
  stepper->sunctx             = sunctx;
  stepper->last_flag          = SUN_SUCCESS;
  stepper->forcing            = NULL;
  stepper->nforcing           = 0;
  stepper->nforcing_allocated = 0;
  stepper->tshift             = SUN_RCONST(0.0);
  stepper->tscale             = SUN_RCONST(0.0);
  stepper->fused_scalars      = NULL;
  stepper->fused_vectors      = NULL;

  stepper->ops = malloc(sizeof(*(stepper->ops)));
  SUNAssert(stepper->ops, SUN_ERR_MALLOC_FAIL);

  stepper->ops->evolve      = NULL;
  stepper->ops->onestep     = NULL;
  stepper->ops->trystep     = NULL;
  stepper->ops->fullrhs     = NULL;
  stepper->ops->reset       = NULL;
  stepper->ops->setstoptime = NULL;

  *stepper_ptr = stepper;

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_Destroy(SUNStepper* stepper_ptr)
{
  SUNFunctionBegin((*stepper_ptr)->sunctx);

  SUNStepper stepper = *stepper_ptr;

  /* free the inner forcing and fused op workspace vector */
  if (stepper->forcing)
  {
    N_VDestroyVectorArray(stepper->forcing, stepper->nforcing);
    N_VDestroyVectorArray(stepper->fused_vectors, stepper->nforcing);
    free(stepper->fused_scalars);
  }

  free(stepper->ops);
  free(stepper);
  *stepper_ptr = NULL;

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_Evolve(SUNStepper stepper, sunrealtype t0,
                             sunrealtype tout, N_Vector y, N_Vector yp,
                             sunrealtype* tret, int* stop_reason)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->evolve)
  {
    return stepper->ops->evolve(stepper, t0, tout, y, yp, tret, stop_reason);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_OneStep(SUNStepper stepper, sunrealtype t0,
                              sunrealtype tout, N_Vector y, N_Vector yp,
                              sunrealtype* tret, int* stop_reason)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->onestep)
  {
    return stepper->ops->onestep(stepper, t0, tout, y, yp, tret, stop_reason);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_TryStep(SUNStepper stepper, sunrealtype t0,
                              sunrealtype tout, N_Vector y, N_Vector yp,
                              sunrealtype* tret, int* stop_reason)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->trystep)
  {
    return stepper->ops->trystep(stepper, t0, tout, y, yp, tret, stop_reason);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_Reset(SUNStepper stepper, sunrealtype tR, N_Vector yR,
                            N_Vector ypR, int64_t ckptIdxR)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->reset)
  {
    return stepper->ops->reset(stepper, tR, yR, ypR, ckptIdxR);
  }
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNStepper_SetContent(SUNStepper stepper, void* content)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->content = content;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_GetContent(SUNStepper stepper, void** content)
{
  SUNFunctionBegin(stepper->sunctx);
  *content = stepper->content;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetEvolveFn(SUNStepper stepper, SUNStepperEvolveFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->evolve = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetOneStepFn(SUNStepper stepper, SUNStepperOneStepFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->onestep = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetTryStepFn(SUNStepper stepper, SUNStepperTryStepFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->trystep = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetFullRhsFn(SUNStepper stepper, SUNStepperFullRhsFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->fullrhs = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetResetFn(SUNStepper stepper, SUNStepperResetFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->reset = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetStopTimeFn(SUNStepper stepper, SUNStepperSetStopTimeFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->setstoptime = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetForcing(SUNStepper stepper, int count, N_Vector tmpl)
{
  SUNFunctionBegin(stepper->sunctx);

  stepper->nforcing = count;

  if (stepper->nforcing_allocated < stepper->nforcing)
  {
    if (stepper->nforcing_allocated)
    {
      N_VDestroyVectorArray(stepper->forcing, stepper->nforcing_allocated);
      SUNCheckLastErr();
    }
    stepper->forcing = N_VCloneVectorArray(stepper->nforcing, tmpl);
    SUNCheckLastErr();
    stepper->nforcing_allocated = stepper->nforcing;
  }

  if (!stepper->fused_vectors)
  {
    stepper->fused_vectors = N_VNewVectorArray(stepper->nforcing, SUNCTX_);
    SUNCheckLastErr();
  }

  if (!stepper->fused_scalars)
  {
    stepper->fused_scalars =
      (sunrealtype*)calloc(count + 1, sizeof(*stepper->fused_scalars));
    SUNAssert(stepper->fused_scalars, SUN_ERR_MEM_FAIL);
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_AddForcing(SUNStepper stepper, sunrealtype t, N_Vector f)
{
  SUNFunctionBegin(stepper->sunctx);

  /* if nforcing is 0, then forcing is 'disabled' */
  if (!stepper->nforcing) { return SUN_SUCCESS; }

  /* always append the constant forcing term */
  stepper->fused_scalars[0] = SUN_RCONST(1.0);
  stepper->fused_vectors[0] = f;

  /* compute normalized time tau and initialize tau^i */
  sunrealtype tau  = (t - stepper->tshift) / (stepper->tscale);
  sunrealtype taui = SUN_RCONST(1.0);

  for (int i = 0; i < stepper->nforcing; i++)
  {
    stepper->fused_scalars[i + 1] = taui;
    stepper->fused_vectors[i + 1] = stepper->forcing[i];
    taui *= tau;
  }

  N_VLinearCombination(stepper->nforcing + 1, stepper->fused_scalars,
                       stepper->fused_vectors, f);

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_GetForcingData(SUNStepper stepper, sunrealtype* tshift,
                                     sunrealtype* tscale, N_Vector** forcing,
                                     int* nforcing)
{
  *tshift   = stepper->tshift;
  *tscale   = stepper->tscale;
  *forcing  = stepper->forcing;
  *nforcing = stepper->nforcing;

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_SetStopTime(SUNStepper stepper, sunrealtype tstop)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->ops->setstoptime)
  {
    return stepper->ops->setstoptime(stepper, tstop);
  }
  else { return SUN_ERR_NOT_IMPLEMENTED; }
}