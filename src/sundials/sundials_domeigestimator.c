/* -----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a generic SUNDomEigEst
 * package. It contains the implementation of the SUNDomEigEstimator
 * operations listed in sundials_domeigestimator.h
 * -----------------------------------------------------------------*/

#include <stdlib.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include <sundials/priv/sundials_domeigestimator_impl.h>
#include <sundials/sundials_domeigestimator.h>

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(SUNDomEigEstimator DEE)
{
  return (DEE->sunctx->profiler);
}
#endif

/* -----------------------------------------------------------------
 * Create a new empty SUNDomEigEstimator object
 * ----------------------------------------------------------------- */

SUNDomEigEstimator SUNDomEigEstNewEmpty(SUNContext sunctx)
{
  SUNDomEigEstimator DEE;
  SUNDomEigEstimator_Ops ops;

  if (sunctx == NULL) { return NULL; }

  SUNFunctionBegin(sunctx);

  /* create dominant eigenvalue estimator object */
  DEE = NULL;
  DEE = (SUNDomEigEstimator)malloc(sizeof *DEE);
  SUNAssertNull(DEE, SUN_ERR_MALLOC_FAIL);

  /* create dominant eigenvalue estimator ops structure */
  ops = NULL;
  ops = (SUNDomEigEstimator_Ops)malloc(sizeof *ops);
  SUNAssertNull(ops, SUN_ERR_MALLOC_FAIL);

  /* initialize operations to NULL */
  ops->gettype            = NULL;
  ops->setatimes          = NULL;
  ops->setnumofperprocess = NULL;
  ops->initialize         = NULL;
  ops->preprocess         = NULL;
  ops->computehess        = NULL;
  ops->initialize         = NULL;
  ops->estimate           = NULL;
  ops->free               = NULL;

  /* attach ops and initialize content and context to NULL */
  DEE->ops     = ops;
  DEE->content = NULL;
  DEE->sunctx  = sunctx;

  return (DEE);
}

/* -----------------------------------------------------------------
 * Free a generic SUNDomEigEstimator (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNDomEigEstFreeEmpty(SUNDomEigEstimator DEE)
{
  if (DEE == NULL) { return; }

  /* free non-NULL ops structure */
  if (DEE->ops) { free(DEE->ops); }
  DEE->ops = NULL;

  /* free overall N_Vector object and return */
  free(DEE);
  return;
}

/* -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * -----------------------------------------------------------------*/

SUNDomEigEstimator_Type SUNDomEigEstGetType(SUNDomEigEstimator DEE)
{
  return (DEE->ops->gettype(DEE));
}

SUNErrCode SUNDomEigEstSetATimes(SUNDomEigEstimator DEE, void* A_data,
                                 SUNATimesFn ATimes)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->setatimes) { ier = DEE->ops->setatimes(DEE, A_data, ATimes); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEstInitialize(SUNDomEigEstimator DEE)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->initialize) { ier = DEE->ops->initialize(DEE); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEstPreProcess(SUNDomEigEstimator DEE)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->preprocess) { ier = DEE->ops->preprocess(DEE); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEstComputeHess(SUNDomEigEstimator DEE)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->computehess) { ier = DEE->ops->computehess(DEE); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEstimate(SUNDomEigEstimator DEE, suncomplextype* dom_eig)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->estimate) { ier = DEE->ops->estimate(DEE, dom_eig); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}