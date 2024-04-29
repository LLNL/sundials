

#include <stdlib.h>
#include <string.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sundials_stepper_impl.h"

SUNErrCode SUNStepper_Create(SUNContext sunctx, SUNStepper* stepper)
{
  SUNFunctionBegin(sunctx);

  *stepper = NULL;
  *stepper = (SUNStepper)malloc(sizeof(**stepper));
  SUNAssert(stepper, SUN_ERR_MALLOC_FAIL);

  memset(*stepper, 0, sizeof(**stepper));

  (*stepper)->ops = (SUNStepper_Ops)malloc(sizeof(*((*stepper)->ops)));
  SUNAssert((*stepper)->ops, SUN_ERR_MALLOC_FAIL);

  memset((*stepper)->ops, 0, sizeof(*((*stepper)->ops)));

  (*stepper)->priv_ops =
    (SUNStepper_PrivOps)malloc(sizeof(*((*stepper)->priv_ops)));
  SUNAssert((*stepper)->priv_ops, SUN_ERR_MALLOC_FAIL);

  memset((*stepper)->priv_ops, 0, sizeof(*((*stepper)->priv_ops)));

  /* initialize stepper data */
  (*stepper)->last_flag = SUN_SUCCESS;
  (*stepper)->sunctx    = sunctx;

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_Free(SUNStepper* stepper)
{
  SUNFunctionBegin((*stepper)->sunctx);

  /* free the inner forcing and fused op workspace vector */
  if ((*stepper)->forcing) { sunStepper_freeForcing(*stepper); }

  /* free operations structure */
  free((*stepper)->ops);

  free((*stepper)->priv_ops);

  /* free inner stepper mem */
  free(*stepper);
  *stepper = NULL;

  return SUN_SUCCESS;
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

SUNErrCode SUNStepper_SetAdvanceFn(SUNStepper stepper, SUNStepperAdvanceFn fn)
{
  SUNFunctionBegin(stepper->sunctx);
  stepper->ops->advance = fn;
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

SUNErrCode SUNStepper_AddForcing(SUNStepper stepper, sunrealtype t, N_Vector f)
{
  SUNFunctionBegin(stepper->sunctx);

  sunrealtype tau, taui;
  int i;

  /* always append the constant forcing term */
  stepper->vals[0] = SUN_RCONST(1.0);
  stepper->vecs[0] = f;

  /* compute normalized time tau and initialize tau^i */
  tau  = (t - stepper->tshift) / (stepper->tscale);
  taui = SUN_RCONST(1.0);

  for (i = 0; i < stepper->nforcing; i++)
  {
    stepper->vals[i + 1] = taui;
    stepper->vecs[i + 1] = stepper->forcing[i];
    taui *= tau;
  }

  N_VLinearCombination(stepper->nforcing + 1, stepper->vals, stepper->vecs, f);

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_GetForcingData(SUNStepper stepper, sunrealtype* tshift,
                                     sunrealtype* tscale, N_Vector** forcing,
                                     int* nforcing)
{
  sunrealtype tau, taui;
  int i;

  *tshift   = stepper->tshift;
  *tscale   = stepper->tscale;
  *forcing  = stepper->forcing;
  *nforcing = stepper->nforcing;

  return SUN_SUCCESS;
}

SUNErrCode sunStepper_allocForcing(SUNStepper stepper, int count, N_Vector tmpl)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->priv_ops->allocForcing)
  {
    return stepper->priv_ops->allocForcing(stepper, count, tmpl);
  }
  else { return SUN_ERR_NOT_IMPLEMENTED; }
}

SUNErrCode sunStepper_freeForcing(SUNStepper stepper)
{
  SUNFunctionBegin(stepper->sunctx);
  if (stepper->priv_ops->freeForcing)
  {
    return stepper->priv_ops->freeForcing(stepper);
  }
  else { return SUN_ERR_NOT_IMPLEMENTED; }
}
