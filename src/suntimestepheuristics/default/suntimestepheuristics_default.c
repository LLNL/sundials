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
 * This is the implementation file for the
 * SUNTimestepHeuristics_Default module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suntimestepheuristics/suntimestepheuristics_default.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SH_CONTENT(H)     ( (SUNTimestepHeuristicsContent_Default)(H->content) )
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
#define SH_SAFETY(H)      ( SH_CONTENT(H)->safety )
#define SH_GROWTH(H)      ( SH_CONTENT(H)->growth )
#define SH_LBOUND(H)      ( SH_CONTENT(H)->lbound )
#define SH_UBOUND(H)      ( SH_CONTENT(H)->ubound )
#define SH_NST_ACC(H)     ( SH_CONTENT(H)->nst_acc )
#define SH_NST_EXP(H)     ( SH_CONTENT(H)->nst_exp )

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_CFLFAC     RCONST(0.5)
#define DEFAULT_SAFETY     RCONST(0.96)
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

SUNTimestepHeuristics SUNTimestepHeuristics_Default(SUNContext sunctx)
{
  SUNTimestepHeuristics H;
  SUNTimestepHeuristicsContent_Default content;

  /* Create an empty heuristics object */
  H = NULL;
  H = SUNTimestepHeuristics_NewEmpty(sunctx);
  if (H == NULL) { return (NULL); }

  /* Attach operations */
  H->ops->getid              = SUNTimestepHeuristics_GetID_Default;
  H->ops->constrainstep      = SUNTimestepHeuristics_ConstrainStep_Default;
  H->ops->etestfail          = SUNTimestepHeuristics_ETestFail_Default;
  H->ops->convfail           = SUNTimestepHeuristics_ConvFail_Default;
  H->ops->boundreduction     = SUNTimestepHeuristics_BoundReduction_Default;
  H->ops->boundfirststep     = SUNTimestepHeuristics_BoundFirstStep_Default;
  H->ops->reset              = SUNTimestepHeuristics_Reset_Default;
  H->ops->update             = SUNTimestepHeuristics_Update_Default;
  H->ops->setdefaults        = SUNTimestepHeuristics_SetDefaults_Default;
  H->ops->write              = SUNTimestepHeuristics_Write_Default;
  H->ops->setmaxstep         = SUNTimestepHeuristics_SetMaxStep_Default;
  H->ops->setminstep         = SUNTimestepHeuristics_SetMinStep_Default;
  H->ops->setexpstabfn       = SUNTimestepHeuristics_SetExpStabFn_Default;
  H->ops->setcflfraction     = SUNTimestepHeuristics_SetCFLFraction_Default;
  H->ops->setsafetyfactor    = SUNTimestepHeuristics_SetSafetyFactor_Default;
  H->ops->setmaxgrowth       = SUNTimestepHeuristics_SetMaxGrowth_Default;
  H->ops->setminreduction    = SUNTimestepHeuristics_SetMinReduction_Default;
  H->ops->setfixedstepbounds = SUNTimestepHeuristics_SetFixedStepBounds_Default;
  H->ops->setmaxfirstgrowth  = SUNTimestepHeuristics_SetMaxFirstGrowth_Default;
  H->ops->setmaxefailgrowth  = SUNTimestepHeuristics_SetMaxEFailGrowth_Default;
  H->ops->setsmallnumefails  = SUNTimestepHeuristics_SetSmallNumEFails_Default;
  H->ops->setmaxcfailgrowth  = SUNTimestepHeuristics_SetMaxCFailGrowth_Default;
  H->ops->getnumexpsteps     = SUNTimestepHeuristics_GetNumExpSteps_Default;
  H->ops->getnumaccsteps     = SUNTimestepHeuristics_GetNumAccSteps_Default;
  H->ops->space              = SUNTimestepHeuristics_Space_Default;

  /* Create content */
  content = NULL;
  content = (SUNTimestepHeuristicsContent_Default)malloc(sizeof *content);
  if (content == NULL)
  {
    (void) SUNTimestepHeuristics_Destroy(H);
    return (NULL);
  }

  /* Attach content */
  H->content = content;

  /* Fill content with default/reset values */
  SUNTimestepHeuristics_SetDefaults_Default(H);
  SUNTimestepHeuristics_Reset_Default(H);

  /* Initialize explicit stability function and data */
  content->expstab = NULL;
  content->estab_data = NULL;

  return (H);
}


/* -----------------------------------------------------------------
 * implementation of heuristic operations
 * ----------------------------------------------------------------- */

SUNTimestepHeuristics_ID SUNTimestepHeuristics_GetID_Default(SUNTimestepHeuristics H)
{ return SUN_TIMESTEPHEURISTICS_STD; }

int SUNTimestepHeuristics_ConstrainStep_Default(SUNTimestepHeuristics H,
                                                sunrealtype hcur,
                                                sunrealtype h_acc,
                                                sunrealtype* hconstr)
{
  /* Determine direction of integration, int_dir */
  sunrealtype int_dir = hcur / SUNRabs(hcur);

  /* Call explicit stability function if present, multiply result by CFL factor * int_dir */
  int retval = 0;
  sunrealtype h_cfl = RCONST(1.0e30) * SUNRabs(hcur);
  if (SH_EXPSTAB(H) != NULL)
  {
    retval = SH_EXPSTAB(H)(&h_cfl, SH_ESTAB_DATA(H));
    if (retval != 0) { return SUNTIMESTEPHEURISTICS_USER_FCN_FAIL; }
    h_cfl *= int_dir * SH_CFL(H);
  }

  /* Enforce safety factor on h_acc prediction. */
  h_acc *= SH_SAFETY(H);

  /* Enforce maximum bound on h_acc growth (etamax) */
  h_acc = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(SH_ETAMAX(H)*hcur));

  /* Enforce minimum bound on time step reduction (etamin) */
  h_acc = int_dir * SUNMAX(SUNRabs(h_acc), SUNRabs(SH_ETAMIN(H)*hcur));

  /* Increment nst_acc vs nst_exp counter, and set desired step */
  if (SUNRabs(h_acc) < SUNRabs(h_cfl)) { SH_NST_ACC(H)++; }
  else { SH_NST_EXP(H)++; }
  sunrealtype hnew = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(h_cfl));

  /* Enforce lbound*ONEMSM and ubound*ONEPSM bounds if hnew>hcur */
  if (SUNRabs(hnew) > SUNRabs(hcur)) {
    if ( (SUNRabs(hnew) > SUNRabs(hcur*SH_LBOUND(H)*ONEMSM)) &&
         (SUNRabs(hnew) < SUNRabs(hcur*SH_UBOUND(H)*ONEPSM)) )
      hnew = hcur;
  }

  /* Enforce max step size */
  hnew /= SUNMAX(RCONST(1.0), SUNRabs(hnew) * SH_HMAX_INV(H));

  /* Bound any requested stepsize reduction, and return */
  return (SUNTimestepHeuristics_BoundReduction_Default(H, hcur, hnew, hconstr));
}

int SUNTimestepHeuristics_ETestFail_Default(SUNTimestepHeuristics H,
                                            sunrealtype hcur,
                                            sunrealtype hnew, int nef,
                                            sunrealtype* hconstr)
{
  /* Set etamax to 1.0 for next step attempt */
  SH_ETAMAX(H) = RCONST(1.0);

  /* Enforce failure bounds on step */
  if (nef >= SH_SMALL_NEF(H)) { hnew = SUNMIN(hnew, hcur*SH_ETAMXF(H)); }

  /* Bound any requested stepsize reduction, and return */
  return (SUNTimestepHeuristics_BoundReduction_Default(H, hcur, hnew, hconstr));
}

int SUNTimestepHeuristics_ConvFail_Default(SUNTimestepHeuristics H,
                                           sunrealtype hcur,
                                           sunrealtype* hconstr)
{
  /* Set etamax to 1.0 for next step attempt */
  SH_ETAMAX(H) = RCONST(1.0);

  /* Enforce failure bounds on next step */
  sunrealtype hnew = hcur*SH_ETACF(H);

  /* Bound any requested stepsize reduction, and return */
  return (SUNTimestepHeuristics_BoundReduction_Default(H, hcur, hnew, hconstr));
}

int SUNTimestepHeuristics_BoundReduction_Default(SUNTimestepHeuristics H,
                                                 sunrealtype hcur,
                                                 sunrealtype hnew,
                                                 sunrealtype *hconstr)
{
  /* If reduction requested first ensure that it its possible. */
  if ((SUNRabs(hcur) <= SH_HMIN(H)*ONEPSM) && (SUNRabs(hnew) <= SUNRabs(hcur)))
  { return SUNTIMESTEPHEURISTICS_CANNOT_DECREASE; }

  /* Enforce minimum step size and return. */
  *hconstr = SUNMAX(SUNRabs(hnew), SH_HMIN(H));
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_BoundFirstStep_Default(SUNTimestepHeuristics H,
                                                 sunrealtype h0,
                                                 sunrealtype *h0constr)
{
  /* Enforce max step size */
  h0 /= SUNMAX(RCONST(1.0), SUNRabs(h0) * SH_HMAX_INV(H));

  /* Enforce minimum step size and return. */
  *h0constr = SUNMAX(SUNRabs(h0), SH_HMIN(H));
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_Reset_Default(SUNTimestepHeuristics H)
{
  SH_ETAMAX(H) = SH_ETAMX1(H);
  SH_NST_ACC(H) = 0;
  SH_NST_EXP(H) = 0;
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_Update_Default(SUNTimestepHeuristics H)
{
  SH_ETAMAX(H) = SH_GROWTH(H);
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetDefaults_Default(SUNTimestepHeuristics H)
{
  SH_HMAX_INV(H) = RCONST(0.0);
  SH_HMIN(H) = RCONST(0.0);
  SH_ETAMX1(H) = DEFAULT_ETAMX1;
  SH_ETAMXF(H) = DEFAULT_ETAMXF;
  SH_ETAMIN(H) = DEFAULT_ETAMIN;
  SH_SMALL_NEF(H) = DEFAULT_SMALL_NEF;
  SH_ETACF(H) = DEFAULT_ETACF;
  SH_CFL(H) = DEFAULT_CFLFAC;
  SH_SAFETY(H) = DEFAULT_SAFETY;
  SH_GROWTH(H) = DEFAULT_GROWTH;
  SH_LBOUND(H) = DEFAULT_HFIXED_LB;
  SH_UBOUND(H) = DEFAULT_HFIXED_UB;
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_Write_Default(SUNTimestepHeuristics H, FILE *fptr)
{
  fprintf(fptr, "Default SUNTimestepHeuristics module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  Minimum step size = %32Lg\n", SH_HMIN(H));
  fprintf(fptr, "  Maximum step size = %32Lg\n", RCONST(1.0)/SH_HMAX_INV(H));
  fprintf(fptr, "  Maximum step increase factor (first step) = %32Lg\n", SH_ETAMX1(H));
  fprintf(fptr, "  Maximum step increase factor (general steps) = %32Lg\n", SH_GROWTH(H));
  fprintf(fptr, "  Current maximum step increase factor = %32Lg\n", SH_ETAMAX(H));
  fprintf(fptr, "  Minimum step reduction factor = %32Lg\n", SH_ETAMIN(H));
  fprintf(fptr, "  Step reduction factor on algebraic solver convergence failure = %32Lg\n", SH_ETACF(H));
  fprintf(fptr, "  Step reduction factor on multiple error fails = %32Lg\n", SH_ETAMXF(H));
  fprintf(fptr, "  Minimum number of error fails before above factor is used = %i\n", SH_SMALL_NEF(H));
  fprintf(fptr, "  Explicit stability safety factor = %32Lg\n", SH_CFL(H));
  fprintf(fptr, "  Step adaptivity safety factor = %32Lg\n", SH_SAFETY(H));
  fprintf(fptr, "  Unchanged step growth interval = [%32Lg,  %32Lg]\n", SH_LBOUND(H), SH_UBOUND(H));
#else
  fprintf(fptr, "  Minimum step size = %16g\n", SH_HMIN(H));
  fprintf(fptr, "  Maximum step size = %16g\n", RCONST(1.0)/SH_HMAX_INV(H));
  fprintf(fptr, "  Maximum step increase factor (first step) = %16g\n", SH_ETAMX1(H));
  fprintf(fptr, "  Maximum step increase factor (general steps) = %16g\n", SH_GROWTH(H));
  fprintf(fptr, "  Current maximum step increase factor = %16g\n", SH_ETAMAX(H));
  fprintf(fptr, "  Minimum step reduction factor = %16g\n", SH_ETAMIN(H));
  fprintf(fptr, "  Step reduction factor on algebraic solver convergence failure = %16g\n", SH_ETACF(H));
  fprintf(fptr, "  Step reduction factor on multiple error fails = %16g\n", SH_ETAMXF(H));
  fprintf(fptr, "  Minimum number of error fails before above factor is used = %i\n", SH_SMALL_NEF(H));
  fprintf(fptr, "  Explicit stability safety factor = %16g\n", SH_CFL(H));
  fprintf(fptr, "  Step adaptivity safety factor = %16g\n", SH_SAFETY(H));
  fprintf(fptr, "  Unchanged step growth interval = [%16g,  %16g]\n", SH_LBOUND(H), SH_UBOUND(H));
#endif
  if (SH_EXPSTAB(H) != NULL) {
    fprintf(fptr, "  User-provided explicit stability function\n");
    fprintf(fptr, "  Stability function data pointer = %p\n", SH_ESTAB_DATA(H));
  }
  fprintf(fptr, "  Current number of accuracy-limited steps = %li\n", SH_NST_ACC(H));
  fprintf(fptr, "  Current number of stability-limited steps = %li\n", SH_NST_EXP(H));
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetMaxStep_Default(SUNTimestepHeuristics H,
                                             sunrealtype hmax)
{
  /* Passing a value <= 0 sets hmax = infinity */
  if (hmax <= RCONST(0.0))
  {
    SH_HMAX_INV(H) = RCONST(0.0);
    return SUNTIMESTEPHEURISTICS_SUCCESS;
  }

  /* set the value and return */
  SH_HMAX_INV(H) = RCONST(1.0)/hmax;

  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetMinStep_Default(SUNTimestepHeuristics H,
                                             sunrealtype hmin)
{
  /* Passing a value <= 0 sets hmin = 0 */
  if (hmin <= RCONST(0.0))
  {
    SH_HMIN(H) = RCONST(0.0);
    return SUNTIMESTEPHEURISTICS_SUCCESS;
  }

  /* set the value and return */
  SH_HMIN(H) = hmin;
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetExpStabFn_Default(SUNTimestepHeuristics H,
                                               SUNExpStabFn EStab,
                                               void* estab_data)
{
  SH_EXPSTAB(H) = EStab;
  SH_ESTAB_DATA(H) = estab_data;
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetCFLFraction_Default(SUNTimestepHeuristics H,
                                                 sunrealtype cfl_frac)
{
  /* check for allowable parameters */
  if (cfl_frac >= RCONST(1.0))
  {
    return SUNTIMESTEPHEURISTICS_ILL_INPUT;
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
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetSafetyFactor_Default(SUNTimestepHeuristics H,
                                                  sunrealtype safety)
{
  /* check for legal input */
  if (safety >= RCONST(1.0)) { return SUNTIMESTEPHEURISTICS_ILL_INPUT; }

  /* set positive-valued parameters, otherwise set default */
  if (safety <= RCONST(0.0)) {
    SH_SAFETY(H) = DEFAULT_SAFETY;
  } else {
    SH_SAFETY(H) = safety;
  }

  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetMaxGrowth_Default(SUNTimestepHeuristics H,
                                               sunrealtype mx_growth)
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
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetMinReduction_Default(SUNTimestepHeuristics H,
                                                  sunrealtype eta_min)
{
  /* set allowed value, otherwise set default */
  if (eta_min >= RCONST(1.0) || eta_min <= RCONST(0.0))
  {
    SH_ETAMIN(H) = DEFAULT_ETAMIN;
  } else {
    SH_ETAMIN(H) = eta_min;
  }
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetFixedStepBounds_Default(SUNTimestepHeuristics H,
                                                     sunrealtype lb, sunrealtype ub)
{
  /* set allowable interval, otherwise set defaults */
  if ((lb <= RCONST(1.0)) && (ub >= RCONST(1.0))) {
    SH_LBOUND(H) = lb;
    SH_UBOUND(H) = ub;
  } else {
    SH_LBOUND(H) = DEFAULT_HFIXED_LB;
    SH_UBOUND(H) = DEFAULT_HFIXED_UB;
  }
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetMaxFirstGrowth_Default(SUNTimestepHeuristics H,
                                                    sunrealtype etamx1)
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
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetMaxEFailGrowth_Default(SUNTimestepHeuristics H,
                                                    sunrealtype etamxf)
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
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetSmallNumEFails_Default(SUNTimestepHeuristics H,
                                                    int small_nef)
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
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_SetMaxCFailGrowth_Default(SUNTimestepHeuristics H,
                                                    sunrealtype etacf)
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
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_GetNumExpSteps_Default(SUNTimestepHeuristics H,
                                                 long int* expsteps)
{
  *expsteps = SH_NST_EXP(H);
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_GetNumAccSteps_Default(SUNTimestepHeuristics H,
                                                 long int* accsteps)
{
  *accsteps = SH_NST_ACC(H);
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}

int SUNTimestepHeuristics_Space_Default(SUNTimestepHeuristics H,
                                        long int* lenrw, long int* leniw)
{
  *lenrw = 11;
  *leniw = 5;
  return SUNTIMESTEPHEURISTICS_SUCCESS;
}
