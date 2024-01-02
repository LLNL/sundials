/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a generic SUNAdaptController
 * package. It contains the implementation of the SUNAdaptController
 * operations listed in sundials_adaptcontroller.h
 * -----------------------------------------------------------------*/

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_adaptcontroller.h>

#include "sundials/sundials_errors.h"

/* -----------------------------------------------------------------
 * Create a new empty SUNAdaptController object
 * ----------------------------------------------------------------- */

SUNAdaptController SUNAdaptController_NewEmpty(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptController_Ops ops;

  /* a context is required */
  if (sunctx == NULL) { return (NULL); }

  SUNFunctionBegin(sunctx);

  /* create controller object */
  C = NULL;
  C = (SUNAdaptController)malloc(sizeof *C);
  SUNAssertNull(C, SUN_ERR_MALLOC_FAIL);

  /* create matrix ops structure */
  ops = NULL;
  ops = (SUNAdaptController_Ops)malloc(sizeof *ops);
  SUNAssertNull(ops, SUN_ERR_MALLOC_FAIL);

  /* initialize operations to NULL */
  ops->gettype      = NULL;
  ops->destroy      = NULL;
  ops->reset        = NULL;
  ops->estimatestep = NULL;
  ops->setdefaults  = NULL;
  ops->write        = NULL;
  ops->seterrorbias = NULL;
  ops->updateh      = NULL;
  ops->space        = NULL;

  /* attach ops and initialize content to NULL */
  C->ops     = ops;
  C->content = NULL;
  C->sunctx  = sunctx;

  return (C);
}

/* -----------------------------------------------------------------
 * Free a generic SUNAdaptController (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNAdaptController_DestroyEmpty(SUNAdaptController C)
{
  if (C == NULL) { return; }

  /* free non-NULL ops structure */
  if (C->ops) { free(C->ops); }
  C->ops = NULL;

  /* free overall SUNAdaptController object and return */
  free(C);
  return;
}

/* -----------------------------------------------------------------
 * Required functions in the 'ops' structure for non-NULL controller
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType(SUNAdaptController C)
{
  if (C == NULL) { return SUN_ADAPTCONTROLLER_NONE; }
  if (C->ops->gettype) { return C->ops->gettype(C); }
  return SUN_ADAPTCONTROLLER_NONE;
}

/* -----------------------------------------------------------------
 * Optional functions in the 'ops' structure
 * ----------------------------------------------------------------- */

SUNErrCode SUNAdaptController_Destroy(SUNAdaptController C)
{
  if (C == NULL) { return (SUN_SUCCESS); }

  /* if the destroy operation exists use it */
  if (C->ops)
  {
    if (C->ops->destroy) { return (C->ops->destroy(C)); }
  }

  /* if we reach this point, either ops == NULL or destroy == NULL,
     try to cleanup by freeing the content, ops, and matrix */
  if (C->content)
  {
    free(C->content);
    C->content = NULL;
  }
  if (C->ops)
  {
    free(C->ops);
    C->ops = NULL;
  }
  free(C);
  C = NULL;

  return (SUN_SUCCESS);
}

SUNErrCode SUNAdaptController_EstimateStep(SUNAdaptController C, sunrealtype h,
                                           int p, sunrealtype dsm,
                                           sunrealtype* hnew)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  SUNAssert(hnew, SUN_ERR_ARG_CORRUPT);
  *hnew = h; /* initialize output with identity */
  if (C->ops->estimatestep) { ier = C->ops->estimatestep(C, h, p, dsm, hnew); }
  return (ier);
}

SUNErrCode SUNAdaptController_Reset(SUNAdaptController C)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  if (C->ops->reset) { ier = C->ops->reset(C); }
  return (ier);
}

SUNErrCode SUNAdaptController_SetDefaults(SUNAdaptController C)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  if (C->ops->setdefaults) { ier = C->ops->setdefaults(C); }
  return (ier);
}

SUNErrCode SUNAdaptController_Write(SUNAdaptController C, FILE* fptr)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  if (C->ops->write) { ier = C->ops->write(C, fptr); }
  return (ier);
}

SUNErrCode SUNAdaptController_SetErrorBias(SUNAdaptController C, sunrealtype bias)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  if (C->ops->seterrorbias) { ier = C->ops->seterrorbias(C, bias); }
  return (ier);
}

SUNErrCode SUNAdaptController_UpdateH(SUNAdaptController C, sunrealtype h,
                                      sunrealtype dsm)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  if (C->ops->updateh) { ier = C->ops->updateh(C, h, dsm); }
  return (ier);
}

SUNErrCode SUNAdaptController_Space(SUNAdaptController C, long int* lenrw,
                                    long int* leniw)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = 0; /* initialize outputs with identity */
  *leniw = 0;
  if (C->ops->space) { ier = C->ops->space(C, lenrw, leniw); }
  return (ier);
}
