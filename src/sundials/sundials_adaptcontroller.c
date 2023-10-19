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
 * This is the implementation file for a generic SUNAdaptController
 * package. It contains the implementation of the SUNAdaptController
 * operations listed in sundials_adaptcontroller.h
 * -----------------------------------------------------------------*/

#include <sundials/sundials_adaptcontroller.h>

/* -----------------------------------------------------------------
 * Create a new empty SUNAdaptController object
 * ----------------------------------------------------------------- */

SUNAdaptController SUNAdaptController_NewEmpty(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptController_Ops ops;

  /* a context is required */
  if (sunctx == NULL) return(NULL);

  /* create controller object */
  C = NULL;
  C = (SUNAdaptController) malloc(sizeof *C);
  if (C == NULL) return(NULL);

  /* create matrix ops structure */
  ops = NULL;
  ops = (SUNAdaptController_Ops) malloc(sizeof *ops);
  if (ops == NULL) { free(C); return(NULL); }

  /* initialize operations to NULL */
  ops->gettype               = NULL;
  ops->destroy               = NULL;
  ops->reset                 = NULL;
  ops->estimatestep          = NULL;
  ops->estimatestepandorder  = NULL;
  ops->estimatemristeps      = NULL;
  ops->estimatesteptol       = NULL;
  ops->setdefaults           = NULL;
  ops->write                 = NULL;
  ops->setmethodorder        = NULL;
  ops->adjustcontrollerorder = NULL;
  ops->seterrorbias          = NULL;
  ops->update                = NULL;
  ops->updatemrih            = NULL;
  ops->updatemritol          = NULL;
  ops->space                 = NULL;

  /* attach ops and initialize content to NULL */
  C->ops     = ops;
  C->content = NULL;
  C->sunctx  = sunctx;

  return(C);
}


/* -----------------------------------------------------------------
 * Free a generic SUNAdaptController (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNAdaptController_FreeEmpty(SUNAdaptController C)
{
  if (C == NULL)  return;

  /* free non-NULL ops structure */
  if (C->ops)  free(C->ops);
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
  return(C->ops->gettype(C));
}

/* -----------------------------------------------------------------
 * Optional functions in the 'ops' structure
 * ----------------------------------------------------------------- */

int SUNAdaptController_Destroy(SUNAdaptController C)
{
  if (C == NULL) return(SUNADAPTCONTROLLER_SUCCESS);

  /* if the destroy operation exists use it */
  if (C->ops)
    if (C->ops->destroy) { return(C->ops->destroy(C)); }

  /* if we reach this point, either ops == NULL or destroy == NULL,
     try to cleanup by freeing the content, ops, and matrix */
  if (C->content) { free(C->content); C->content = NULL; }
  if (C->ops) { free(C->ops); C->ops = NULL; }
  free(C); C = NULL;

  return(SUNADAPTCONTROLLER_SUCCESS);
}

int SUNAdaptController_EstimateStep(SUNAdaptController C, sunrealtype h, sunrealtype dsm,
                                    sunrealtype* hnew)
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


int SUNAdaptController_EstimateStepAndOrder(SUNAdaptController C, sunrealtype h, int q,
                                            sunrealtype dsm, sunrealtype* hnew,
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

int SUNAdaptController_EstimateMRISteps(SUNAdaptController C, sunrealtype H, sunrealtype h,
                                        sunrealtype DSM, sunrealtype dsm,
                                        sunrealtype* Hnew, sunrealtype *hnew)
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

int SUNAdaptController_EstimateStepTol(SUNAdaptController C, sunrealtype H,
                                       sunrealtype tolfac, sunrealtype DSM,
                                       sunrealtype dsm, sunrealtype *Hnew,
                                       sunrealtype* tolfacnew)
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

int SUNAdaptController_Reset(SUNAdaptController C)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->reset) { ier = C->ops->reset(C); }
  return(ier);
}

int SUNAdaptController_SetDefaults(SUNAdaptController C)
{
  int ier = 0;
  if (C->ops->setdefaults) { ier = C->ops->setdefaults(C); }
  return(ier);
}

int SUNAdaptController_Write(SUNAdaptController C, FILE* fptr)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->write) { ier = C->ops->write(C, fptr); }
  return(ier);
}

int SUNAdaptController_SetMethodOrder(SUNAdaptController C, int p)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->setmethodorder) { ier = C->ops->setmethodorder(C, p); }
  return(ier);
}

int SUNAdaptController_AdjustControllerOrder(SUNAdaptController C, int adj)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->adjustcontrollerorder) { ier = C->ops->adjustcontrollerorder(C, adj); }
  return(ier);
}

int SUNAdaptController_SetErrorBias(SUNAdaptController C, sunrealtype bias)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->seterrorbias) { ier = C->ops->seterrorbias(C, bias); }
  return(ier);
}

int SUNAdaptController_Update(SUNAdaptController C, sunrealtype h, sunrealtype dsm)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->update) { ier = C->ops->update(C, h, dsm); }
  return(ier);
}

int SUNAdaptController_UpdateMRIH(SUNAdaptController C, sunrealtype H, sunrealtype h,
                                  sunrealtype DSM, sunrealtype dsm)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->updatemrih) { ier = C->ops->updatemrih(C, H, h, DSM, dsm); }
  return(ier);
}

int SUNAdaptController_UpdateMRITol(SUNAdaptController C, sunrealtype H, sunrealtype tolfac,
                                    sunrealtype DSM, sunrealtype dsm)
{
  int ier = 0;
  if (C == NULL) { return ier; }
  if (C->ops->updatemritol) { ier = C->ops->updatemritol(C, H, tolfac, DSM, dsm); }
  return(ier);
}

int SUNAdaptController_Space(SUNAdaptController C, long int *lenrw, long int *leniw)
{
  int ier = 0;
  *lenrw = 0;   /* initialize outputs with identity */
  *leniw = 0;
  if (C == NULL) { return ier; }
  if (C->ops->space) { ier = C->ops->space(C, lenrw, leniw); }
  return(ier);
}
