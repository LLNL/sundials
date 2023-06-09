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
 * This is the implementation file for a generic SUNHeuristics
 * package. It contains the implementation of the SUNHeuristics
 * operations listed in sundials_heuristics.h
 * -----------------------------------------------------------------*/

#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

/* -----------------------------------------------------------------
 * Create a new empty SUNHeuristics object
 * ----------------------------------------------------------------- */

SUNHeuristics SUNHeuristicsNewEmpty(SUNContext sunctx)
{
  SUNHeuristics H;
  SUNHeuristics_Ops ops;

  /* a context is required */
  if (sunctx == NULL) return(NULL);

  /* create heuristics object */
  H = NULL;
  H = (SUNHeuristics) malloc(sizeof *H);
  if (H == NULL) return(NULL);

  /* create matrix ops structure */
  ops = NULL;
  ops = (SUNHeuristics_Ops) malloc(sizeof *ops);
  if (ops == NULL) { free(H); return(NULL); }

  /* initialize operations to NULL */
  ops->getid              = NULL;
  ops->destroy            = NULL;
  ops->constrainstep      = NULL;
  ops->reset              = NULL;
  ops->setdefaults        = NULL;
  ops->write              = NULL;
  ops->setexpstabfn       = NULL;
  ops->setcflfraction     = NULL;
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
 * Free a generic SUNHeuristics (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNHeuristicsFreeEmpty(SUNHeuristics H)
{
  if (H == NULL)  return;

  /* free non-NULL ops structure */
  if (H->ops)  free(H->ops);
  H->ops = NULL;

  /* free overall SUNHeuristics object and return */
  free(H);
  return;
}


/* -----------------------------------------------------------------
 * Required functions in the 'ops' structure
 * ----------------------------------------------------------------- */

SUNHeuristics_ID SUNHeuristicsGetID(SUNHeuristics H)
{
  return(H->ops->getid(H));
}

void SUNHeuristicsDestroy(SUNHeuristics H)
{
  if (H == NULL) return;

  /* if the destroy operation exists use it */
  if (H->ops)
    if (H->ops->destroy) { H->ops->destroy(H); return; }

  /* if we reach this point, either ops == NULL or destroy == NULL,
     try to cleanup by freeing the content, ops, and matrix */
  if (H->content) { free(H->content); H->content = NULL; }
  if (H->ops) { free(H->ops); H->ops = NULL; }
  free(H); H = NULL;

  return;
}

/* -----------------------------------------------------------------
 * Optional functions in the 'ops' structure
 * ----------------------------------------------------------------- */

int SUNHeuristicsConstrainStep(SUNHeuristics H, realtype hcur,
                               realtype hnew, booleantype efail,
                               booleantype cfail, realtype *hconstr)
{
  int ier = 0;
  *hconstr = hnew;   /* initialize output with identity */
  if (H->ops->constrainstep)
  {
    ier = H->ops->constrainstep(H, hcur, hnew,
                                efail, cfail, hconstr);
  }
  return(ier);
}

int SUNHeuristicsReset(SUNHeuristics H)
{
  int ier = 0;
  if (H->ops->reset) { ier = H->ops->reset(H); }
  return(ier);
}

int SUNHeuristicsSetDefaults(SUNHeuristics H)
{
  int ier = 0;
  if (H->ops->setdefaults) { ier = H->ops->setdefaults(H); }
  return(ier);
}

int SUNHeuristicsWrite(SUNHeuristics H, FILE* fptr)
{
  int ier = 0;
  if (H->ops->write) { ier = H->ops->write(H, fptr); }
  return(ier);
}

int SUNHeuristicsSetExpStabFn(SUNHeuristics H,
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

int SUNHeuristicsSetCFLFraction(SUNHeuristics H,
                                realtype cfl_frac)
{
  int ier = 0;
  if (H->ops->setcflfraction)
  {
    ier = H->ops->setcflfraction(H, cfl_frac);
  }
  return(ier);
}

int SUNHeuristicsSetMaxGrowth(SUNHeuristics H,
                              realtype mx_growth)
{
  int ier = 0;
  if (H->ops->setmaxgrowth)
  {
    ier = H->ops->setmaxgrowth(H, mx_growth);
  }
  return(ier);
}

int SUNHeuristicsSetMinReduction(SUNHeuristics H,
                                 realtype eta_min)
{
  int ier = 0;
  if (H->ops->setminreduction)
  {
    ier = H->ops->setminreduction(H, eta_min);
  }
  return(ier);
}

int SUNHeuristicsSetFixedStepBounds(SUNHeuristics H,
                                    realtype lb, realtype ub)
{
  int ier = 0;
  if (H->ops->setfixedstepbounds)
  {
    ier = H->ops->setfixedstepbounds(H, lb, ub);
  }
  return(ier);
}

int SUNHeuristicsSetMaxFirstGrowth(SUNHeuristics H,
                                   realtype etamx1)
{
  int ier = 0;
  if (H->ops->setmaxfirstgrowth)
  {
    ier = H->ops->setmaxfirstgrowth(H, etamx1);
  }
  return(ier);
}

int SUNHeuristicsSetMaxEFailGrowth(SUNHeuristics H,
                                   realtype etamxf)
{
  int ier = 0;
  if (H->ops->setmaxefailgrowth)
  {
    ier = H->ops->setmaxefailgrowth(H, etamxf);
  }
  return(ier);
}

int SUNHeuristicsSetSmallNumEFails(SUNHeuristics H,
                                   int small_nef)
{
  int ier = 0;
  if (H->ops->setsmallnumefails)
  {
    ier = H->ops->setsmallnumefails(H, small_nef);
  }
  return(ier);
}

int SUNHeuristicsSetMaxCFailGrowth(SUNHeuristics H,
                                   realtype etacf)
{
  int ier = 0;
  if (H->ops->setmaxcfailgrowth)
  {
    ier = H->ops->setmaxcfailgrowth(H, etacf);
  }
  return(ier);
}

int SUNHeuristicsGetNumExpSteps(SUNHeuristics H,
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

int SUNHeuristicsGetNumAccSteps(SUNHeuristics H,
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

int SUNHeuristicsSpace(SUNHeuristics H, long int *lenrw,
                       long int *leniw)
{
  int ier = 0;
  *lenrw = 0;   /* initialize outputs with identity */
  *leniw = 0;
  if (H->ops->space) { ier = H->ops->space(H, lenrw, leniw); }
  return(ier);
}
