/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a generic SUNControl
 * package. It contains the implementation of the SUNControl
 * operations listed in sundials_control.h
 * -----------------------------------------------------------------*/

#include <sundials/sundials_control.h>

/* -----------------------------------------------------------------
 * Create a new empty SUNControl object
 * ----------------------------------------------------------------- */

SUNControl SUNControlNewEmpty(SUNContext sunctx)
{
  SUNControl C;
  SUNControl_Ops ops;

  /* a context is required */
  if (sunctx == NULL) return(NULL);

  /* create controller object */
  C = NULL;
  C = (SUNControl) malloc(sizeof *C);
  if (C == NULL) return(NULL);

  /* create matrix ops structure */
  ops = NULL;
  ops = (SUNControl_Ops) malloc(sizeof *ops);
  if (ops == NULL) { free(C); return(NULL); }

  /* initialize operations to NULL */
  ops->getid                = NULL;
  ops->destroy              = NULL;
  ops->reset                = NULL;
  ops->estimatestep         = NULL;
  ops->estimatestepandorder = NULL;
  ops->estimatemristeps     = NULL;
  ops->estimatesteptol      = NULL;
  ops->setdefaults          = NULL;
  ops->write                = NULL;
  ops->setmethodorder       = NULL;
  ops->setembeddingorder    = NULL;
  ops->seterrorbias         = NULL;
  ops->update               = NULL;
  ops->updatemrih           = NULL;
  ops->updatemritol         = NULL;
  ops->space                = NULL;

  /* attach ops and initialize content to NULL */
  C->ops     = ops;
  C->content = NULL;
  C->sunctx  = sunctx;

  return(C);
}


/* -----------------------------------------------------------------
 * Free a generic SUNControl (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNControlFreeEmpty(SUNControl C)
{
  if (C == NULL)  return;

  /* free non-NULL ops structure */
  if (C->ops)  free(C->ops);
  C->ops = NULL;

  /* free overall SUNControl object and return */
  free(C);
  return;
}


/* -----------------------------------------------------------------
 * Required functions in the 'ops' structure for non-NULL controller
 * ----------------------------------------------------------------- */

SUNControl_Type SUNControlGetType(SUNControl C)
{
  if (C == NULL) { return SUNDIALS_CONTROL_NONE; }
  return(C->ops->getid(C));
}

/* -----------------------------------------------------------------
 * Optional functions in the 'ops' structure
 * ----------------------------------------------------------------- */

void SUNControlDestroy(SUNControl C)
{
  if (C == NULL) return;

  /* if the destroy operation exists use it */
  if (C->ops)
    if (C->ops->destroy) { C->ops->destroy(C); return; }

  /* if we reach this point, either ops == NULL or destroy == NULL,
     try to cleanup by freeing the content, ops, and matrix */
  if (C->content) { free(C->content); C->content = NULL; }
  if (C->ops) { free(C->ops); C->ops = NULL; }
  free(C); C = NULL;

  return;
}

int SUNControlEstimateStep(SUNControl C, realtype h, realtype dsm,
                           realtype* hnew)
{
  int ier = 0;
  *hnew = h;   /* initialize output with identity */
  if (C == NULL) { return ier; }
  if (C->ops->estimatestep)
  {
    ier = C->ops->estimatestep(C, h, dsm, hnew);
  }
  return(ier);
}


int SUNControlEstimateStepAndOrder(SUNControl C, realtype h, int q,
                                   realtype dsm, realtype* hnew,
                                   int *qnew)
{
  int ier = 0;
  *hnew = h;   /* initialize outputs with identity */
  *qnew = q;
  if (C == NULL) { return ier; }
  if (C->ops->estimatestepandorder)
  {
    ier = C->ops->estimatestepandorder(C, h, q, dsm, hnew, qnew);
  }
  return(ier);
}

int SUNControlEstimateMRISteps(SUNControl C, realtype H, realtype h,
                               realtype DSM, realtype dsm,
                               realtype* Hnew, realtype *hnew)
{
  int ier = 0;
  *Hnew = H;   /* initialize outputs with identity */
  *hnew = h;
  if (C == NULL) { return ier; }
  if (C->ops->estimatemristeps)
  {
    ier = C->ops->estimatemristeps(C, H, DSM, dsm, h, Hnew, hnew);
  }
  return(ier);
}

int SUNControlEstimateStepTol(SUNControl C, realtype H,
                              realtype tolfac, realtype DSM,
                              realtype dsm, realtype *Hnew,
                              realtype* tolfacnew)
{
  int ier = 0;
  *Hnew = H;   /* initialize outputs with identity */
  *tolfacnew = tolfac;
  if (C == NULL) { return ier; }
  if (C->ops->estimatesteptol)
  {
    ier = C->ops->estimatesteptol(C, H, tolfac, DSM, dsm,
                                  Hnew, tolfacnew);
  }
  return(ier);
}

int SUNControlReset(SUNControl C)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->reset) { ier = C->ops->reset(C); }
  return(ier);
}

int SUNControlSetDefaults(SUNControl C)
{
  int ier = 0;
  if (C->ops->setdefaults) { ier = C->ops->setdefaults(C); }
  return(ier);
}

int SUNControlWrite(SUNControl C, FILE* fptr)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->write) { ier = C->ops->write(C, fptr); }
  return(ier);
}

int SUNControlSetMethodOrder(SUNControl C, int q)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->setmethodorder) { ier = C->ops->setmethodorder(C, q); }
  return(ier);
}

int SUNControlSetEmbeddingOrder(SUNControl C, int p)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->setembeddingorder) { ier = C->ops->setembeddingorder(C, p); }
  return(ier);
}

int SUNControlSetErrorBias(SUNControl C, realtype bias)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->seterrorbias) { ier = C->ops->seterrorbias(C, bias); }
  return(ier);
}

int SUNControlUpdate(SUNControl C, realtype h, realtype dsm)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->update) { ier = C->ops->update(C, h, dsm); }
  return(ier);
}

int SUNControlUpdateMRIH(SUNControl C, realtype H, realtype h,
                         realtype DSM, realtype dsm)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->updatemrih) { ier = C->ops->updatemrih(C, H, h, DSM, dsm); }
  return(ier);
}

int SUNControlUpdateMRITol(SUNControl C, realtype H, realtype tolfac,
                           realtype DSM, realtype dsm)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->updatemritol) { ier = C->ops->updatemritol(C, H, tolfac,
                                                         DSM, dsm); }
  return(ier);
}

int SUNControlSpace(SUNControl C, long int *lenrw, long int *leniw)
{
  int ier = 0;
  *lenrw = 0;   /* initialize outputs with identity */
  *leniw = 0;
  if (C == NULL) { return ier; }
  if (C->ops->space) { ier = C->ops->space(C, lenrw, leniw); }
  return(ier);
}
