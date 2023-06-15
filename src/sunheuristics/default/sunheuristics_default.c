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
 * This is the implementation file for the SUNHeuristics_Default
 * module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sunheuristics/sunheuristics_default.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SH_CONTENT(H)     ( (SUNHeuristicsContent_Default)(H->content) )
#define SH_HMAX_INV(H)    ( SH_CONTENT(H)->hmax_inv )
#define SH_HMIN(H)        ( SH_CONTENT(H)->hmin )
#define SH_ETAMAX(H)      ( SH_CONTENT(H)->etamax )
#define SH_ETAMX1(H)      ( SH_CONTENT(H)->etamx1 )
#define SH_ETAMXF(H)      ( SH_CONTENT(H)->etamxf )
#define SH_ETAMIN(H)      ( SH_CONTENT(H)->etamin )
#define SH_SMALL_NEF(H)   ( SH_CONTENT(H)->small_nef )
#define SH_ETACF(H)       ( SH_CONTENT(H)->etacf )
#define SH_EXPSTAB(H)     ( SH_CONTENT(H)->expstab )
#define SH_ESTAB_DATA(H)  ( SH_CONTENT(H)->estab_data )
#define SH_CFL(H)         ( SH_CONTENT(H)->cfl )
#define SH_GROWTH(H)      ( SH_CONTENT(H)->growth )
#define SH_LBOUND(H)      ( SH_CONTENT(H)->lbound )
#define SH_UBOUND(H)      ( SH_CONTENT(H)->ubound )
#define SH_NST_ACC(H)     ( SH_CONTENT(H)->nst_acc )
#define SH_NST_EXP(H)     ( SH_CONTENT(H)->nst_exp )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_CFLFAC     RCONST(0.5)
#define DEFAULT_GROWTH     RCONST(20.0)
#define DEFAULT_HFIXED_LB  RCONST(1.0)
#define DEFAULT_HFIXED_UB  RCONST(1.5)
#define DEFAULT_ETAMX1     RCONST(10000.0)
#define DEFAULT_ETAMXF     RCONST(0.3)
#define DEFAULT_ETAMIN     RCONST(0.1)
#define DEFAULT_ETACF      RCONST(0.25)
#define DEFAULT_SMALL_NEF  2
#define ONEPSM             RCONST(1.000001)
#define ONEMSM             RCONST(0.999999)


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new default heuristics module
 */

SUNHeuristics SUNHeuristicsDefault(SUNContext sunctx)
{
  SUNHeuristics H;
  SUNHeuristicsContent_Default content;

  /* Create an empty heuristics object */
  H = NULL;
  H = SUNHeuristicsNewEmpty(sunctx);
  if (H == NULL) { return (NULL); }

  /* Attach operations */
  H->ops->getid              = SUNHeuristicsGetID_Default;
  H->ops->constrainstep      = SUNHeuristicsConstrainStep_Default;
  H->ops->etestfail          = SUNHeuristicsETestFail_Default;
  H->ops->convfail           = SUNHeuristicsConvFail_Default;
  H->ops->boundreduction     = SUNHeuristicsBoundReduction_Default;
  H->ops->reset              = SUNHeuristicsReset_Default;
  H->ops->update             = SUNHeuristicsUpdate_Default;
  H->ops->setdefaults        = SUNHeuristicsSetDefaults_Default;
  H->ops->write              = SUNHeuristicsWrite_Default;
  H->ops->setmaxstep         = SUNHeuristicsSetMaxStep_Default;
  H->ops->setminstep         = SUNHeuristicsSetMinStep_Default;
  H->ops->setexpstabfn       = SUNHeuristicsSetExpStabFn_Default;
  H->ops->setcflfraction     = SUNHeuristicsSetCFLFraction_Default;
  H->ops->setmaxgrowth       = SUNHeuristicsSetMaxGrowth_Default;
  H->ops->setminreduction    = SUNHeuristicsSetMinReduction_Default;
  H->ops->setfixedstepbounds = SUNHeuristicsSetFixedStepBounds_Default;
  H->ops->setmaxfirstgrowth  = SUNHeuristicsSetMaxFirstGrowth_Default;
  H->ops->setmaxefailgrowth  = SUNHeuristicsSetMaxEFailGrowth_Default;
  H->ops->setsmallnumefails  = SUNHeuristicsSetSmallNumEFails_Default;
  H->ops->setmaxcfailgrowth  = SUNHeuristicsSetMaxCFailGrowth_Default;
  H->ops->getnumexpsteps     = SUNHeuristicsGetNumExpSteps_Default;
  H->ops->getnumaccsteps     = SUNHeuristicsGetNumAccSteps_Default;
  H->ops->space              = SUNHeuristicsSpace_Default;

  /* Create content */
  content = NULL;
  content = (SUNHeuristicsContent_Default)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNHeuristicsDestroy(H);
    return (NULL);
  }

  /* Attach content */
  H->content = content;

  /* Fill content with default/reset values */
  SUNHeuristicsSetDefaults_Default(H);
  SUNHeuristicsReset_Default(H);

  /* Initialize explicit stability function and data */
  content->expstab = NULL;
  content->estab_data = NULL;

  return (H);
}


/* -----------------------------------------------------------------
 * implementation of heuristic operations
 * ----------------------------------------------------------------- */

SUNHeuristics_ID SUNHeuristicsGetID_Default(SUNHeuristics H)
{ return SUNDIALS_HEURISTICS_STD; }

int SUNHeuristicsConstrainStep_Default(SUNHeuristics H, realtype hcur,
                                       realtype h_acc, realtype* hconstr)
{
  /* Determine direction of integration, int_dir */
  realtype int_dir = hcur / SUNRabs(hcur);

  /* Call explicit stability function if present, multiply result by CFL factor * int_dir */
  int retval = 0;
  realtype h_cfl = RCONST(1.0e3) * SUNRabs(hcur);
  if (SH_EXPSTAB(H) != NULL)
  {
    retval = SH_EXPSTAB(H)(&h_cfl, SH_ESTAB_DATA(H));
    if (retval != 0) { return SUNHEURISTICS_USER_FCN_FAIL; }
    h_cfl *= int_dir * SH_CFL(H);
  }

  /* Enforce maximum bound on h_acc growth (etamax) */
  h_acc = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(SH_ETAMAX(H)*hcur));

  /* Enforce minimum bound on time step reduction (etamin) */
  h_acc = int_dir * SUNMAX(SUNRabs(h_acc), SUNRabs(SH_ETAMIN(H)*hcur));

  /* Increment nst_acc vs nst_exp counter, and set desired step */
  if (SUNRabs(h_acc) < SUNRabs(h_cfl)) { SH_NST_ACC(H)++; }
  else { SH_NST_EXP(H)++; }
  realtype hnew = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(h_cfl));

  /* Enforce lbound*ONEMSM and ubound*ONEPSM bounds if hnew>hcur */
  if (SUNRabs(hnew) > SUNRabs(hcur)) {
    if ( (SUNRabs(hnew) > SUNRabs(hcur*SH_LBOUND(H)*ONEMSM)) &&
         (SUNRabs(hnew) < SUNRabs(hcur*SH_UBOUND(H)*ONEPSM)) )
      hnew = hcur;
  }

  /* Enforce max step size */
  hnew /= SUNMAX(RCONST(1.0), SUNRabs(hnew) * SH_HMAX_INV(H));

  /* Bound any requested stepsize reduction, and return */
  return (SUNHeuristicsBoundReduction_Default(H, hcur, hnew, hconstr));
}

int SUNHeuristicsETestFail_Default(SUNHeuristics H, realtype hcur,
                                   realtype hnew, int nef, realtype* hconstr)
{
  /* Set etamax to 1.0 for next step attempt */
  SH_ETAMAX(H) = RCONST(1.0);

  /* Enforce failure bounds on step */
  if (nef >= SH_SMALL_NEF(H)) { hnew = SUNMIN(hnew, hcur*SH_ETAMXF(H)); }

  /* Bound any requested stepsize reduction, and return */
  return (SUNHeuristicsBoundReduction_Default(H, hcur, hnew, hconstr));
}

int SUNHeuristicsConvFail_Default(SUNHeuristics H, realtype hcur,
                                  realtype* hconstr)
{
  /* Set etamax to 1.0 for next step attempt */
  SH_ETAMAX(H) = RCONST(1.0);

  /* Enforce failure bounds on next step */
  realtype hnew = hcur*SH_ETACF(H);

  /* Bound any requested stepsize reduction, and return */
  return (SUNHeuristicsBoundReduction_Default(H, hcur, hnew, hconstr));
}

int SUNHeuristicsBoundReduction_Default(SUNHeuristics H, realtype hcur,
                                        realtype hnew, realtype *hconstr)
{
  /* If reduction requested first ensure that it its possible. */
  if ((SUNRabs(hcur) <= SH_HMIN(H)*ONEPSM) && (SUNRabs(hnew) <= SUNRabs(hcur)))
  { return SUNHEURISTICS_CANNOT_DECREASE; }

  /* Enforce minimum step size and return. */
  *hconstr = SUNMAX(SUNRabs(hnew), SH_HMIN(H));
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsReset_Default(SUNHeuristics H)
{
  SH_ETAMAX(H) = SH_ETAMX1(H);
  SH_NST_ACC(H) = 0;
  SH_NST_EXP(H) = 0;
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsUpdate_Default(SUNHeuristics H)
{
  SH_ETAMAX(H) = SH_GROWTH(H);
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetDefaults_Default(SUNHeuristics H)
{
  SH_HMAX_INV(H) = RCONST(0.0);
  SH_HMIN(H) = RCONST(0.0);
  SH_ETAMX1(H) = DEFAULT_ETAMX1;
  SH_ETAMXF(H) = DEFAULT_ETAMXF;
  SH_ETAMIN(H) = DEFAULT_ETAMIN;
  SH_SMALL_NEF(H) = DEFAULT_SMALL_NEF;
  SH_ETACF(H) = DEFAULT_ETACF;
  SH_CFL(H) = DEFAULT_CFLFAC;
  SH_GROWTH(H) = DEFAULT_GROWTH;
  SH_LBOUND(H) = DEFAULT_HFIXED_LB;
  SH_UBOUND(H) = DEFAULT_HFIXED_UB;
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsWrite_Default(SUNHeuristics H, FILE *fptr)
{
  fprintf(fptr, "SUNHeuristics_Default module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  hmin = %12Lg\n", SH_HMIN(H));
  fprintf(fptr, "  hmax_inv = %12Lg\n", SH_HMAX_INV(H));
  fprintf(fptr, "  etamax = %12Lg\n", SH_ETAMAX(H));
  fprintf(fptr, "  etamx1 = %12Lg\n", SH_ETAMX1(H));
  fprintf(fptr, "  etamxf = %12Lg\n", SH_ETAMXF(H));
  fprintf(fptr, "  etamin = %12Lg\n", SH_ETAMIN(H));
  fprintf(fptr, "  etacf = %12Lg\n", SH_ETACF(H));
  fprintf(fptr, "  cfl = %12Lg\n", SH_CFL(H));
  fprintf(fptr, "  growth = %12Lg\n", SH_GROWTH(H));
  fprintf(fptr, "  lbound = %12Lg\n", SH_LBOUND(H));
  fprintf(fptr, "  ubound = %12Lg\n", SH_UBOUND(H));
#else
  fprintf(fptr, "  hmin = %12g\n", SH_HMIN(H));
  fprintf(fptr, "  hmax_inv = %12g\n", SH_HMAX_INV(H));
  fprintf(fptr, "  etamax = %12g\n", SH_ETAMAX(H));
  fprintf(fptr, "  etamx1 = %12g\n", SH_ETAMX1(H));
  fprintf(fptr, "  etamxf = %12g\n", SH_ETAMXF(H));
  fprintf(fptr, "  etamin = %12g\n", SH_ETAMIN(H));
  fprintf(fptr, "  etacf = %12g\n", SH_ETACF(H));
  fprintf(fptr, "  cfl = %12g\n", SH_CFL(H));
  fprintf(fptr, "  growth = %12g\n", SH_GROWTH(H));
  fprintf(fptr, "  lbound = %12g\n", SH_LBOUND(H));
  fprintf(fptr, "  ubound = %12g\n", SH_UBOUND(H));
#endif
  fprintf(fptr, "  small_nef = %i\n", SH_SMALL_NEF(H));
  if (SH_EXPSTAB(H) != NULL) {
    fprintf(fptr, "  user-provided explicit stability function\n");
    fprintf(fptr, "  stability function data pointer = %p\n", SH_ESTAB_DATA(H));
  }
  fprintf(fptr, "  nst_acc = %li\n", SH_NST_ACC(H));
  fprintf(fptr, "  nst_exp = %li\n", SH_NST_EXP(H));
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetMaxStep_Default(SUNHeuristics H, realtype hmax)
{
  realtype hmax_inv;

  /* Passing a value <= 0 sets hmax = infinity */
  if (hmax <= RCONST(0.0))
  {
    SH_HMAX_INV(H) = RCONST(0.0);
    return SUNHEURISTICS_SUCCESS;
  }

  /* check that hmax and hmin are agreeable */
  hmax_inv = RCONST(1.0)/hmax;
  if (hmax_inv * SH_HMIN(H) > RCONST(1.0))
  {
    return SUNHEURISTICS_ILL_INPUT;
  }

  /* set the value and return */
  SH_HMAX_INV(H) = hmax_inv;

  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetMinStep_Default(SUNHeuristics H, realtype hmin)
{
  /* Passing a value <= 0 sets hmin = 0 */
  if (hmin <= RCONST(0.0))
  {
    SH_HMIN(H) = RCONST(0.0);
    return SUNHEURISTICS_SUCCESS;
  }

  /* check that hmin and hmax are agreeable */
  if (hmin * SH_HMAX_INV(H) > RCONST(1.0))
  {
    return SUNHEURISTICS_ILL_INPUT;
  }

  /* set the value and return */
  SH_HMIN(H) = hmin;
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetExpStabFn_Default(SUNHeuristics H, SUNExpStabFn EStab,
                                      void* estab_data)
{
  SH_EXPSTAB(H) = EStab;
  SH_ESTAB_DATA(H) = estab_data;
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetCFLFraction_Default(SUNHeuristics H, realtype cfl_frac)
{
  /* check for allowable parameters */
  if (cfl_frac >= RCONST(1.0))
  {
    return SUNHEURISTICS_ILL_INPUT;
  }

  /* set positive-valued parameters, otherwise set default */
  if (cfl_frac <= RCONST(0.0))
  {
    SH_CFL(H) = DEFAULT_CFLFAC;
  }
  else
  {
    SH_CFL(H) = cfl_frac;
  }
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetMaxGrowth_Default(SUNHeuristics H, realtype mx_growth)
{
  /* set allowed value, otherwise set default */
  if (mx_growth <= RCONST(1.0))
  {
    SH_GROWTH(H) = DEFAULT_GROWTH;
  }
  else
  {
    SH_GROWTH(H) = mx_growth;
  }
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetMinReduction_Default(SUNHeuristics H, realtype eta_min)
{
  /* set allowed value, otherwise set default */
  if (eta_min >= RCONST(1.0) || eta_min <= RCONST(0.0))
  {
    SH_ETAMIN(H) = DEFAULT_ETAMIN;
  } else {
    SH_ETAMIN(H) = eta_min;
  }
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetFixedStepBounds_Default(SUNHeuristics H, realtype lb, realtype ub)
{
  /* set allowable interval, otherwise set defaults */
  if ((lb <= RCONST(1.0)) && (ub >= RCONST(1.0))) {
    SH_LBOUND(H) = lb;
    SH_UBOUND(H) = ub;
  } else {
    SH_LBOUND(H) = DEFAULT_HFIXED_LB;
    SH_UBOUND(H) = DEFAULT_HFIXED_UB;
  }
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetMaxFirstGrowth_Default(SUNHeuristics H, realtype etamx1)
{
  /* if argument legal set it, otherwise set default */
  if (etamx1 <= RCONST(1.0))
  {
    SH_ETAMX1(H) = DEFAULT_ETAMX1;
  }
  else
  {
    SH_ETAMX1(H) = etamx1;
  }
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetMaxEFailGrowth_Default(SUNHeuristics H, realtype etamxf)
{
  /* if argument legal set it, otherwise set default */
  if ((etamxf <= RCONST(0.0)) || (etamxf > RCONST(1.0)))
  {
    SH_ETAMXF(H) = DEFAULT_ETAMXF;
  }
  else
  {
    SH_ETAMXF(H) = etamxf;
  }
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetSmallNumEFails_Default(SUNHeuristics H, int small_nef)
{
  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0)
  {
    SH_SMALL_NEF(H) = DEFAULT_SMALL_NEF;
  }
  else
  {
    SH_SMALL_NEF(H) = small_nef;
  }
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSetMaxCFailGrowth_Default(SUNHeuristics H, realtype etacf)
{
  /* if argument legal set it, otherwise set default */
  if ((etacf <= RCONST(0.0)) || (etacf > RCONST(1.0)))
  {
    SH_ETACF(H) = DEFAULT_ETACF;
  }
  else
  {
    SH_ETACF(H) = etacf;
  }
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsGetNumExpSteps_Default(SUNHeuristics H, long int* expsteps)
{
  *expsteps = SH_NST_EXP(H);
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsGetNumAccSteps_Default(SUNHeuristics H, long int* accsteps)
{
  *accsteps = SH_NST_ACC(H);
  return SUNHEURISTICS_SUCCESS;
}

int SUNHeuristicsSpace_Default(SUNHeuristics H, long int* lenrw, long int* leniw)
{
  *lenrw = 11;
  *leniw = 5;
  return SUNHEURISTICS_SUCCESS;
}
