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
 * This is the implementation file for a generic
 * SUNTimestepHeuristics package. It contains the implementation of
 * the SUNTimestepHeuristics operations listed in
 * sundials_timestepheuristics.h
 * -----------------------------------------------------------------*/

#include <sundials/sundials_timestepheuristics.h>
#include <sundials/sundials_nvector.h>

/* -----------------------------------------------------------------
 * Create a new empty SUNTimestepHeuristics object
 * ----------------------------------------------------------------- */

SUNTimestepHeuristics SUNTimestepHeuristics_NewEmpty(SUNContext sunctx)
{
  SUNTimestepHeuristics H;
  SUNTimestepHeuristics_Ops ops;

  /* a context is required */
  if (sunctx == NULL) return(NULL);

  /* create heuristics object */
  H = NULL;
  H = (SUNTimestepHeuristics) malloc(sizeof *H);
  if (H == NULL) return(NULL);

  /* create matrix ops structure */
  ops = NULL;
  ops = (SUNTimestepHeuristics_Ops) malloc(sizeof *ops);
  if (ops == NULL) { free(H); return(NULL); }

  /* initialize operations to NULL */
  ops->getid              = NULL;
  ops->destroy            = NULL;
  ops->constrainstep      = NULL;
  ops->etestfail          = NULL;
  ops->convfail           = NULL;
  ops->boundreduction     = NULL;
  ops->boundfirststep     = NULL;
  ops->reset              = NULL;
  ops->update             = NULL;
  ops->setdefaults        = NULL;
  ops->write              = NULL;
  ops->setmaxstep         = NULL;
  ops->setminstep         = NULL;
  ops->setexpstabfn       = NULL;
  ops->setcflfraction     = NULL;
  ops->setsafetyfactor    = NULL;
  ops->setmaxgrowth       = NULL;
  ops->setminreduction    = NULL;
  ops->setfixedstepbounds = NULL;
  ops->setmaxfirstgrowth  = NULL;
  ops->setmaxefailgrowth  = NULL;
  ops->setsmallnumefails  = NULL;
  ops->setmaxcfailgrowth  = NULL;
  ops->getnumexpsteps     = NULL;
  ops->getnumaccsteps     = NULL;
  ops->space              = NULL;

  /* attach ops and initialize content to NULL */
  H->ops     = ops;
  H->content = NULL;
  H->sunctx  = sunctx;

  return(H);
}


/* -----------------------------------------------------------------
 * Required functions in the 'ops' structure
 * ----------------------------------------------------------------- */

SUNTimestepHeuristics_ID SUNTimestepHeuristics_GetID(SUNTimestepHeuristics H)
{
  if (H == NULL) { return SUN_TIMESTEPHEURISTICS_NULL; }
  return(H->ops->getid(H));
}

/* -----------------------------------------------------------------
 * Optional functions in the 'ops' structure
 * ----------------------------------------------------------------- */

int SUNTimestepHeuristics_Destroy(SUNTimestepHeuristics H)
{
  if (H == NULL) return(SUNTIMESTEPHEURISTICS_SUCCESS);

  /* if the destroy operation exists use it */
  if (H->ops)
    if (H->ops->destroy) { return(H->ops->destroy(H)); }

  /* if we reach this point, either ops == NULL or destroy == NULL,
     try to cleanup by freeing the content, ops, and matrix */
  if (H->content) { free(H->content); H->content = NULL; }
  if (H->ops) { free(H->ops); H->ops = NULL; }
  free(H); H = NULL;

  return(SUNTIMESTEPHEURISTICS_SUCCESS);
}

int SUNTimestepHeuristics_ConstrainStep(SUNTimestepHeuristics H,
                                        sunrealtype hcur,
                                        sunrealtype hnew, sunrealtype *hconstr)
{
  int ier = 0;
  *hconstr = hnew;   /* initialize output with identity */
  if (H->ops->constrainstep)
  {
    ier = H->ops->constrainstep(H, hcur, hnew, hconstr);
  }
  return(ier);
}

int SUNTimestepHeuristics_ETestFail(SUNTimestepHeuristics H,
                                    sunrealtype hcur,
                                    sunrealtype hnew, int nef,
                                    sunrealtype *hconstr)
{
  int ier = 0;
  *hconstr = hnew;   /* initialize output with identity */
  if (H->ops->etestfail)
  {
    ier = H->ops->etestfail(H, hcur, hnew, nef, hconstr);
  }
  return(ier);
}

int SUNTimestepHeuristics_ConvFail(SUNTimestepHeuristics H,
                                   sunrealtype hcur,
                                   sunrealtype *hconstr)
{
  int ier = 0;
  *hconstr = hcur;   /* initialize output with identity */
  if (H->ops->convfail)
  {
    ier = H->ops->convfail(H, hcur, hconstr);
  }
  return(ier);
}

int SUNTimestepHeuristics_BoundReduction(SUNTimestepHeuristics H,
                                         sunrealtype hcur,
                                         sunrealtype hnew,
                                         sunrealtype *hconstr)
{
  int ier = 0;
  *hconstr = hnew;   /* initialize output with identity */
  if (H->ops->boundreduction)
  {
    ier = H->ops->boundreduction(H, hcur, hnew, hconstr);
  }
  return(ier);
}


int SUNTimestepHeuristics_BoundFirstStep(SUNTimestepHeuristics H,
                                         sunrealtype h0,
                                         sunrealtype *h0constr)
{
  int ier = 0;
  *h0constr = h0;   /* initialize output with identity */
  if (H->ops->boundfirststep)
  {
    ier = H->ops->boundfirststep(H, h0, h0constr);
  }
  return(ier);
}


int SUNTimestepHeuristics_Reset(SUNTimestepHeuristics H)
{
  int ier = 0;
  if (H->ops->reset) { ier = H->ops->reset(H); }
  return(ier);
}

int SUNTimestepHeuristics_Update(SUNTimestepHeuristics H)
{
  int ier = 0;
  if (H->ops->update) { ier = H->ops->update(H); }
  return(ier);
}

int SUNTimestepHeuristics_SetDefaults(SUNTimestepHeuristics H)
{
  int ier = 0;
  if (H->ops->setdefaults) { ier = H->ops->setdefaults(H); }
  return(ier);
}

int SUNTimestepHeuristics_Write(SUNTimestepHeuristics H, FILE* fptr)
{
  int ier = 0;
  if (H->ops->write) { ier = H->ops->write(H, fptr); }
  return(ier);
}

int SUNTimestepHeuristics_SetMaxStep(SUNTimestepHeuristics H,
                                     sunrealtype hmax)
{
  int ier = 0;
  if (H->ops->setmaxstep)
  {
    ier = H->ops->setmaxstep(H, hmax);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetMinStep(SUNTimestepHeuristics H,
                                     sunrealtype hmin)
{
  int ier = 0;
  if (H->ops->setminstep)
  {
    ier = H->ops->setminstep(H, hmin);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetExpStabFn(SUNTimestepHeuristics H,
                                       SUNExpStabFn EStab,
                                       void* estab_data)
{
  int ier = 0;
  if (H->ops->setexpstabfn)
  {
    ier = H->ops->setexpstabfn(H, EStab, estab_data);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetCFLFraction(SUNTimestepHeuristics H,
                                         sunrealtype cfl_frac)
{
  int ier = 0;
  if (H->ops->setcflfraction)
  {
    ier = H->ops->setcflfraction(H, cfl_frac);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetSafetyFactor(SUNTimestepHeuristics H,
                                          sunrealtype safety)
{
  int ier = 0;
  if (H->ops->setsafetyfactor)
  {
    ier = H->ops->setsafetyfactor(H, safety);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetMaxGrowth(SUNTimestepHeuristics H,
                                       sunrealtype mx_growth)
{
  int ier = 0;
  if (H->ops->setmaxgrowth)
  {
    ier = H->ops->setmaxgrowth(H, mx_growth);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetMinReduction(SUNTimestepHeuristics H,
                                          sunrealtype eta_min)
{
  int ier = 0;
  if (H->ops->setminreduction)
  {
    ier = H->ops->setminreduction(H, eta_min);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetFixedStepBounds(SUNTimestepHeuristics H,
                                             sunrealtype lb, sunrealtype ub)
{
  int ier = 0;
  if (H->ops->setfixedstepbounds)
  {
    ier = H->ops->setfixedstepbounds(H, lb, ub);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetMaxFirstGrowth(SUNTimestepHeuristics H,
                                            sunrealtype etamx1)
{
  int ier = 0;
  if (H->ops->setmaxfirstgrowth)
  {
    ier = H->ops->setmaxfirstgrowth(H, etamx1);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetMaxEFailGrowth(SUNTimestepHeuristics H,
                                            sunrealtype etamxf)
{
  int ier = 0;
  if (H->ops->setmaxefailgrowth)
  {
    ier = H->ops->setmaxefailgrowth(H, etamxf);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetSmallNumEFails(SUNTimestepHeuristics H,
                                            int small_nef)
{
  int ier = 0;
  if (H->ops->setsmallnumefails)
  {
    ier = H->ops->setsmallnumefails(H, small_nef);
  }
  return(ier);
}

int SUNTimestepHeuristics_SetMaxCFailGrowth(SUNTimestepHeuristics H,
                                            sunrealtype etacf)
{
  int ier = 0;
  if (H->ops->setmaxcfailgrowth)
  {
    ier = H->ops->setmaxcfailgrowth(H, etacf);
  }
  return(ier);
}

int SUNTimestepHeuristics_GetNumExpSteps(SUNTimestepHeuristics H,
                                         long int* expsteps)
{
  int ier = 0;
  *expsteps = 0;   /* initialize output with identity */
  if (H->ops->getnumexpsteps)
  {
    ier = H->ops->getnumexpsteps(H, expsteps);
  }
  return(ier);
}

int SUNTimestepHeuristics_GetNumAccSteps(SUNTimestepHeuristics H,
                                         long int* accsteps)
{
  int ier = 0;
  *accsteps = 0;   /* initialize output with identity */
  if (H->ops->getnumaccsteps)
  {
    ier = H->ops->getnumaccsteps(H, accsteps);
  }
  return(ier);
}

int SUNTimestepHeuristics_Space(SUNTimestepHeuristics H,
                                long int *lenrw,
                                long int *leniw)
{
  int ier = 0;
  *lenrw = 0;   /* initialize outputs with identity */
  *leniw = 0;
  if (H->ops->space) { ier = H->ops->space(H, lenrw, leniw); }
  return(ier);
}
