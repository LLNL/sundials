

#include <stdlib.h>
#include <string.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sundials/sundials_errors.h"
#include "sundials/sundials_nvector.h"
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

  /* initialize stepper data */
  (*stepper)->last_flag = SUN_SUCCESS;
  (*stepper)->sunctx    = sunctx;

  return SUN_SUCCESS;
}

SUNErrCode SUNStepper_Free(SUNStepper* stepper_ptr)
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
  sunrealtype tau, taui;
  int i;

  *tshift   = stepper->tshift;
  *tscale   = stepper->tscale;
  *forcing  = stepper->forcing;
  *nforcing = stepper->nforcing;

  return SUN_SUCCESS;
}
