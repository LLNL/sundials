/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for ARKode's temporal
 * interpolation utility.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

#define NO_DEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
#include <nvector/nvector_serial.h>
#endif


/*---------------------------------------------------------------
 arkInterpCreate:

 This routine creates an ARKodeInterpMem structure, through
 cloning an input template N_Vector.  This returns a non-NULL
 structure if no errors occurred, or a NULL value otherwise.
---------------------------------------------------------------*/
ARKodeInterpMem arkInterpCreate(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeInterpMem interp_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* allocate structure */
  interp_mem = (ARKodeInterpMem) malloc(sizeof(struct ARKodeInterpMemRec));
  if (interp_mem == NULL)  return(NULL);
  memset(interp_mem, 0, sizeof(struct ARKodeInterpMemRec));

  /* set interpolation order based on user request (if possible) */
  /* if ((ark_mem->dense_q < 0) || (ark_mem->dense_q > 3)) { */
  if ((ark_mem->dense_q < 0) || (ark_mem->dense_q > 5)) {
    interp_mem->order = QDENSE_DEF;
  } else {
    interp_mem->order = ark_mem->dense_q;
  }

  /* vector allocation */
  if (!arkAllocVec(ark_mem, ark_mem->yn, &interp_mem->fold)) {
    arkInterpFree(&interp_mem); return(NULL);
  }
  if (!arkAllocVec(ark_mem, ark_mem->yn, &interp_mem->fnew)) {
    arkInterpFree(&interp_mem); return(NULL);
  }
  if (!arkAllocVec(ark_mem, ark_mem->yn, &interp_mem->yold)) {
    arkInterpFree(&interp_mem); return(NULL);
  }
  if (!arkAllocVec(ark_mem, ark_mem->yn, &interp_mem->fa)) {
    arkInterpFree(&interp_mem); return(NULL);
  }
  if (!arkAllocVec(ark_mem, ark_mem->yn, &interp_mem->fb)) {
    arkInterpFree(&interp_mem); return(NULL);
  }

  /* set ynew pointer to ark_mem->yn */
  interp_mem->ynew = ark_mem->yn;

  /* update workspace sizes */
  ark_mem->lrw += ARK_INTERP_LRW;
  ark_mem->liw += ARK_INTERP_LIW;

  /* copy ark_mem->yn into yold */
  N_VScale(ONE, ark_mem->yn, interp_mem->yold);

  /* initialize time values */
  interp_mem->told = ark_mem->tcur;
  interp_mem->tnew = ark_mem->tcur;
  interp_mem->t_fa = RCONST(0.0);
  interp_mem->t_fb = RCONST(0.0);
  interp_mem->h    = RCONST(0.0);

  return(interp_mem);
}


/*---------------------------------------------------------------
 arkInterpResize:

 This routine resizes the internal vectors in an ARKodeInterpMem
 structure.
---------------------------------------------------------------*/
int arkInterpResize(void* arkode_mem, ARKodeInterpMem interp_mem,
                    ARKVecResizeFn resize, void *resize_data,
                    sunindextype lrw_diff, sunindextype liw_diff,
                    N_Vector y0)
{
  int ier;
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* resize vectors */
  if (interp_mem == NULL)  return(ARK_SUCCESS);
  if (interp_mem->fold != NULL) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &interp_mem->fold);
    if (ier != ARK_SUCCESS)  return(ier);
  }
  if (interp_mem->fnew != NULL) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &interp_mem->fnew);
    if (ier != ARK_SUCCESS)  return(ier);
  }
  if (interp_mem->yold != NULL) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &interp_mem->yold);
    if (ier != ARK_SUCCESS)  return(ier);
  }
  if (interp_mem->fa != NULL) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &interp_mem->fa);
    if (ier != ARK_SUCCESS)  return(ier);
  }
  if (interp_mem->fb != NULL) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &interp_mem->fb);
    if (ier != ARK_SUCCESS)  return(ier);
  }

  /* update yold with current solution */
  N_VScale(ONE, y0, interp_mem->yold);

  /* update ynew pointer to point to current ark_mem->yn */
  interp_mem->ynew = ark_mem->yn;

  /* reinitialize time values */
  interp_mem->told = ark_mem->tcur;
  interp_mem->tnew = ark_mem->tcur;
  interp_mem->t_fa = RCONST(0.0);
  interp_mem->t_fb = RCONST(0.0);
  interp_mem->h    = RCONST(0.0);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkInterpFree:

 This routine frees an ARKodeInterpMem structure.
---------------------------------------------------------------*/
void arkInterpFree(ARKodeInterpMem *interp_mem)
{
  if (*interp_mem != NULL) {
    if ((*interp_mem)->fold != NULL) N_VDestroy((*interp_mem)->fold);
    if ((*interp_mem)->fnew != NULL) N_VDestroy((*interp_mem)->fnew);
    if ((*interp_mem)->yold != NULL) N_VDestroy((*interp_mem)->yold);
    if ((*interp_mem)->fa   != NULL) N_VDestroy((*interp_mem)->fa);
    if ((*interp_mem)->fb   != NULL) N_VDestroy((*interp_mem)->fb);
    free(*interp_mem);
  }
}


/*---------------------------------------------------------------
 arkPrintInterpMem

 This routine outputs the temporal interpolation memory structure
 to a specified file pointer.
---------------------------------------------------------------*/
void arkPrintInterpMem(ARKodeInterpMem interp_mem, FILE *outfile)
{
  if (interp_mem != NULL) {
    fprintf(outfile, "ark_interp: order = %d\n", interp_mem->order);
    fprintf(outfile, "ark_interp: told = %"RSYM"\n", interp_mem->told);
    fprintf(outfile, "ark_interp: tnew = %"RSYM"\n", interp_mem->tnew);
    fprintf(outfile, "ark_interp: t_fa = %"RSYM"\n", interp_mem->t_fa);
    fprintf(outfile, "ark_interp: t_fb = %"RSYM"\n", interp_mem->t_fb);
    fprintf(outfile, "ark_interp: h = %"RSYM"\n", interp_mem->h);
#ifdef DEBUG_OUTPUT
    if (interp_mem->fold != NULL) {
      fprintf(outfile, "ark_interp: fold:\n");
      N_VPrint_Serial(interp_mem->fold);
    }
    if (interp_mem->fnew != NULL) {
      fprintf(outfile, "ark_interp: fnew:\n");
      N_VPrint_Serial(interp_mem->fnew);
    }
    if (interp_mem->yold != NULL) {
      fprintf(outfile, "ark_interp: yold:\n");
      N_VPrint_Serial(interp_mem->yold);
    }
    if (interp_mem->ynew != NULL) {
      fprintf(outfile, "ark_interp: ynew:\n");
      N_VPrint_Serial(interp_mem->ynew);
    }
    if (interp_mem->fa != NULL) {
      fprintf(outfile, "ark_interp: fa:\n");
      N_VPrint_Serial(interp_mem->fa);
    }
    if (interp_mem->fb != NULL) {
      fprintf(outfile, "ark_interp: fb:\n");
      N_VPrint_Serial(interp_mem->fb);
    }
#endif
  }
}


/*---------------------------------------------------------------
 arkInterpInit

 This routine performs the following steps:
 1. Sets tnew and told to the input time
 1. Copies ark_mem->yn into yold
 2. Calls the full RHS routine to fill fnew
 3. Copies fnew into fold
---------------------------------------------------------------*/
int arkInterpInit(void* arkode_mem, ARKodeInterpMem interp,
                  realtype tnew)
{
  int ier;
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* return with success if no interpolation structure is allocated */
  if (interp == NULL)  return(ARK_SUCCESS);

  /* initialize time values */
  interp->told = tnew;
  interp->tnew = tnew;
  interp->h    = RCONST(0.0);

  /* copy current solution into yold */
  N_VScale(ONE, ark_mem->yn, interp->yold);

  /* fill fnew */
  ier = ark_mem->step_fullrhs(ark_mem, tnew, interp->ynew,
                              interp->fnew, 0);
  if (ier != 0)  return(ARK_RHSFUNC_FAIL);

  /* copy fnew into fold */
  N_VScale(ONE, interp->fnew, interp->fold);

  /* return with success */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkInterpUpdate

 This routine performs the following steps:
 1. Copies ynew into yold, and swaps the fnew <-> fold pointers,
    so that yold and fold contain the previous values
 2. Calls the full RHS routine to fill fnew, using ark_mem->ycur
    for the time-evolved solution (since ynew==ark_mem->yn
    has not been updated yet).

 Note: if forceRHS==SUNTRUE, then any previously-stored RHS
 function data in the time step module is suspect, and all RHS
 function(s) require recomputation; we therefore signal the
 fullrhs function with a corresponding flag.
---------------------------------------------------------------*/
int arkInterpUpdate(void* arkode_mem, ARKodeInterpMem interp,
                    realtype tnew, booleantype forceRHS)
{
  int ier, mode;
  N_Vector tempvec;
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* return with success if no interpolation structure is allocated */
  if (interp == NULL)  return(ARK_SUCCESS);

  /* copy ynew into yold */
  N_VScale(ONE, interp->ynew, interp->yold);

  /* swap fold & fnew N_Vector pointers */
  tempvec = interp->fold;
  interp->fold = interp->fnew;
  interp->fnew = tempvec;

  /* update time values */
  interp->told = interp->tnew;
  interp->tnew = tnew;
  interp->h    = ark_mem->h;

  /* determine mode for calling fullrhs */
  mode = (forceRHS) ? 0 : 1;

  /* fill fnew */
  ier = ark_mem->step_fullrhs(ark_mem, tnew, ark_mem->ycur,
                              interp->fnew, mode);
  if (ier != 0)  return(ARK_RHSFUNC_FAIL);

  /* return with success */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkInterpEvaluate

  This routine evaluates a temporal interpolation/extrapolation
  based on the data in the interpolation structure:
     yold = y(told)
     ynew = y(tnew)
     fold = f(told, yold)
     fnew = f(told, ynew)
  This typically consists of using a cubic Hermite interpolating
  formula with this data.  If greater polynomial order than 3 is
  requested, then we can bootstrap up to a 5th-order accurate
  interpolant.  For lower order interpolants than cubic, we use:
     {yold,ynew,fnew} for quadratic
     {yold,ynew} for linear
     {0.5*(yold+ynew)} for constant.

  Derivatives have lower accuracy than the interpolant
  itself, losing one order per derivative.  We will provide
  derivatives up to d = min(5,q).

  The input 'tau' specifies the time at which to return derivative
  information, the formula is
               t = told + tau*(tnew-told),
  where h = tnew-told, i.e. values 0<tau<1 provide interpolation,
  other values result in extrapolation.
---------------------------------------------------------------*/
int arkInterpEvaluate(void* arkode_mem, ARKodeInterpMem interp,
                      realtype tau, int d, int order, N_Vector yout)
{
  /* local variables */
  int q, retval;
  realtype tval, a0, a1, tau2, tau3, tau4, tau5, tau6;
  realtype h, h2, h3, h4, h5, h6;
  realtype a[7];
  N_Vector X[7];
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* set constants */
  tau2 = tau*tau;
  tau3 = tau*tau2;
  tau4 = tau*tau3;
  tau5 = tau*tau4;
  tau6 = tau*tau5;
  h = interp->h;
  h2 = h*h;
  h3 = h*h2;
  h4 = h*h3;
  h5 = h*h4;
  h6 = h*h5;

  /* determine polynomial order q */
  q = SUNMAX(order, 0);        /* respect lower bound  */
  /* q = SUNMIN(q, 3);            /\* respect max possible *\/ */
  q = SUNMIN(q, 5);            /* respect max possible */

  /* error on illegal d */
  if (d < 0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInterpEvaluate", "Requested illegal derivative.");
    return (ARK_ILL_INPUT);
  }

  /* if d is too high, just return zeros */
  if (d > q) {
    N_VConst(ZERO, yout);
    return(ARK_SUCCESS);
  }

  /* build polynomial based on order */
  switch (q) {

  case(0):    /* constant interpolant, yout = 0.5*(yn+yp) */
    N_VLinearSum(HALF, interp->yold, HALF, interp->ynew, yout);
    break;

  case(1):    /* linear interpolant */
    if (d == 0) {
      a0 = -tau;
      a1 = ONE+tau;
    } else {  /* d=1 */
      a0 = -ONE/h;
      a1 =  ONE/h;
    }
    N_VLinearSum(a0, interp->yold, a1, interp->ynew, yout);
    break;

  case(2):    /* quadratic interpolant */
    if (d == 0) {
      a[0] = tau2;
      a[1] = ONE - tau2;
      a[2] = h*(tau2 + tau);
    } else if (d == 1) {
      a[0] = TWO*tau/h;
      a[1] = -TWO*tau/h;
      a[2] = (ONE + TWO*tau);
    } else {  /* d == 2 */
      a[0] = TWO/h/h;
      a[1] = -TWO/h/h;
      a[2] = TWO/h;
    }
    X[0] = interp->yold;
    X[1] = interp->ynew;
    X[2] = interp->fnew;
    retval = N_VLinearCombination(3, a, X, yout);
    if (retval != 0)  return(ARK_VECTOROP_ERR);
    break;

  case(3):    /* cubic interpolant */
    if (d == 0) {
      a[0] = THREE*tau2 + TWO*tau3;
      a[1] = ONE - THREE*tau2 - TWO*tau3;
      a[2] = h*(tau2 + tau3);
      a[3] = h*(tau + TWO*tau2 + tau3);
    } else if (d == 1) {
      a[0] = SIX*(tau + tau2)/h;
      a[1] = -SIX*(tau + tau2)/h;
      a[2] = TWO*tau + THREE*tau2;
      a[3] = ONE + FOUR*tau + THREE*tau2;
    } else if (d == 2) {
      a[0] = SIX*(ONE + TWO*tau)/h2;
      a[1] = -SIX*(ONE + TWO*tau)/h2;
      a[2] = (TWO + SIX*tau)/h;
      a[3] = (FOUR + SIX*tau)/h;
    } else {  /* d == 3 */
      a[0] = TWELVE/h3;
      a[1] = -TWELVE/h3;
      a[2] = SIX/h2;
      a[3] = SIX/h2;
    }
    X[0] = interp->yold;
    X[1] = interp->ynew;
    X[2] = interp->fold;
    X[3] = interp->fnew;
    retval = N_VLinearCombination(4, a, X, yout);
    if (retval != 0) return(ARK_VECTOROP_ERR);
   break;

  case(4):    /* quartic interpolant */

    /* first, evaluate cubic interpolant at tau=-1/3 */
    tval = -ONE/THREE;
    retval = arkInterpEvaluate(arkode_mem, interp, tval, 0, 3, yout);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* second, evaluate RHS at tau=-1/3, storing the result in fa */
    tval = interp->tnew - h/THREE;
    /* if (SUNRabs(tval - interp->t_fa) > 10.0*UNIT_ROUNDOFF) { */
      retval = ark_mem->step_fullrhs(ark_mem, tval, yout, interp->fa, 2);
      if (retval != 0)  return(ARK_RHSFUNC_FAIL);
    /*   interp->t_fa = tval; */
    /* } */

    /* evaluate desired function */
    if (d == 0) {
      a[0] = -SIX*tau2 - RCONST(16.0)*tau3 - RCONST(9.0)*tau4;
      a[1] = ONE + SIX*tau2 + RCONST(16.0)*tau3 + RCONST(9.0)*tau4;
      a[2] = h*FOURTH*(-FIVE*tau2 - RCONST(14.0)*tau3 - RCONST(9.0)*tau4);
      a[3] = h*(tau + TWO*tau2 + tau3);
      a[4] = h*RCONST(27.0)*FOURTH*(-tau4 - TWO*tau3 - tau2);
    } else if (d == 1) {
      a[0] = (-TWELVE*tau - RCONST(48.0)*tau2 - RCONST(36.0)*tau3)/h;
      a[1] = (TWELVE*tau + RCONST(48.0)*tau2 + RCONST(36.0)*tau3)/h;
      a[2] = HALF*(-FIVE*tau - RCONST(21.0)*tau2 - RCONST(18.0)*tau3);
      a[3] = (ONE + FOUR*tau + THREE*tau2);
      a[4] = -RCONST(27.0)*HALF*(TWO*tau3 + THREE*tau2 + tau);
    } else if (d == 2) {
      a[0] = (-TWELVE - RCONST(96.0)*tau - RCONST(108.0)*tau2)/h2;
      a[1] = (TWELVE + RCONST(96.0)*tau + RCONST(108.0)*tau2)/h2;
      a[2] = (-FIVE*HALF - RCONST(21.0)*tau - RCONST(27.0)*tau2)/h;
      a[3] = (FOUR + SIX*tau)/h;
      a[4] = (-RCONST(27.0)*HALF - RCONST(81.0)*tau - RCONST(81.0)*tau2)/h;
    } else if (d == 3) {
      a[0] = (-RCONST(96.0) - RCONST(216.0)*tau)/h3;
      a[1] = (RCONST(96.0) + RCONST(216.0)*tau)/h3;
      a[2] = (-RCONST(21.0) - RCONST(54.0)*tau)/h2;
      a[3] = SIX/h2;
      a[4] = (-RCONST(81.0) - RCONST(162.0)*tau)/h2;
    } else {  /* d == 4 */
      a[0] = -RCONST(216.0)/h4;
      a[1] = RCONST(216.0)/h4;
      a[2] = -RCONST(54.0)/h3;
      a[3] = ZERO;
      a[4] = -RCONST(162.0)/h3;
    }
    X[0] = interp->yold;
    X[1] = interp->ynew;
    X[2] = interp->fold;
    X[3] = interp->fnew;
    X[4] = interp->fa;
    retval = N_VLinearCombination(5, a, X, yout);
    if (retval != 0) return(ARK_VECTOROP_ERR);
    break;

  case(5):    /* quintic interpolant */

    /* first, evaluate quartic interpolant at tau=-1/3 */
    tval = -ONE/THREE;
    retval = arkInterpEvaluate(arkode_mem, interp, tval, 0, 4, yout);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* second, evaluate RHS at tau=-1/3, storing the result in fa */
    tval = interp->tnew - h/THREE;
    /* if (SUNRabs(tval - interp->t_fa) > 10.0*UNIT_ROUNDOFF) { */
      retval = ark_mem->step_fullrhs(ark_mem, tval, yout, interp->fa, 2);
      if (retval != 0)  return(ARK_RHSFUNC_FAIL);
    /*   interp->t_fa = tval; */
    /* } */

    /* third, evaluate quartic interpolant at tau=-2/3 */
    tval = -TWO/THREE;
    retval = arkInterpEvaluate(arkode_mem, interp, tval, 0, 4, yout);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* fourth, evaluate RHS at tau=-2/3, storing the result in fb */
    tval = interp->tnew - h*TWO/THREE;
    /* if (SUNRabs(tval - interp->t_fb) > 10.0*UNIT_ROUNDOFF) { */
      retval = ark_mem->step_fullrhs(ark_mem, tval, yout, interp->fb, 2);
      if (retval != 0)  return(ARK_RHSFUNC_FAIL);
    /*   interp->t_fb = tval; */
    /* } */

    /* evaluate desired function */
    if (d == 0) {
      a[0] = RCONST(54.0)*tau5 + RCONST(135.0)*tau4 + RCONST(110.0)*tau3 + RCONST(30.0)*tau2;
      a[1] = ONE - a[0];
      a[2] = h/FOUR*(RCONST(27.0)*tau5 + RCONST(63.0)*tau4 + RCONST(49.0)*tau3 + RCONST(13.0)*tau2);
      a[3] = h/FOUR*(RCONST(27.0)*tau5 + RCONST(72.0)*tau4 + RCONST(67.0)*tau3 + RCONST(26.0)*tau2 + FOUR*tau);
      a[4] = h/FOUR*(RCONST(81.0)*tau5 + RCONST(189.0)*tau4 + RCONST(135.0)*tau3 + RCONST(27.0)*tau2);
      a[5] = h/FOUR*(RCONST(81.0)*tau5 + RCONST(216.0)*tau4 + RCONST(189.0)*tau3 + RCONST(54.0)*tau2);
    } else if (d == 1) {
      a[0] = (RCONST(270.0)*tau4 + RCONST(540.0)*tau3 + RCONST(330.0)*tau2 + RCONST(60.0)*tau)/h;
      a[1] = -a[0];
      a[2] = (RCONST(135.0)*tau4 + RCONST(252.0)*tau3 + RCONST(147.0)*tau2 + RCONST(26.0)*tau)/FOUR;
      a[3] = (RCONST(135.0)*tau4 + RCONST(288.0)*tau3 + RCONST(201.0)*tau2 + RCONST(52.0)*tau + FOUR)/FOUR;
      a[4] = (RCONST(405.0)*tau4 + RCONST(4.0)*189*tau3 + RCONST(405.0)*tau2 + RCONST(54.0)*tau)/FOUR;
      a[5] = (RCONST(405.0)*tau4 + RCONST(864.0)*tau3 + RCONST(567.0)*tau2 + RCONST(108.0)*tau)/FOUR;
    } else if (d == 2) {
      a[0] = (RCONST(1080.0)*tau3 + RCONST(1620.0)*tau2 + RCONST(660.0)*tau + RCONST(60.0))/h2;
      a[1] = -a[0];
      a[2] = (RCONST(270.0)*tau3 + RCONST(378.0)*tau2 + RCONST(147.0)*tau + RCONST(13.0))/(TWO*h);
      a[3] = (RCONST(270.0)*tau3 + RCONST(432.0)*tau2 + RCONST(201.0)*tau + RCONST(26.0))/(TWO*h);
      a[4] = (RCONST(810.0)*tau3 + RCONST(1134.0)*tau2 + RCONST(405.0)*tau + RCONST(27.0))/(TWO*h);
      a[5] = (RCONST(810.0)*tau3 + RCONST(1296.0)*tau2 + RCONST(567.0)*tau + RCONST(54.0))/(TWO*h);
    } else if (d == 3) {
      a[0] = (RCONST(3240.0)*tau2 + RCONST(3240.0)*tau + RCONST(660.0))/h3;
      a[1] = -a[0];
      a[2] = (RCONST(810.0)*tau2 + RCONST(756.0)*tau + RCONST(147.0))/(TWO*h2);
      a[3] = (RCONST(810.0)*tau2 + RCONST(864.0)*tau + RCONST(201.0))/(TWO*h2);
      a[4] = (RCONST(2430.0)*tau2 + RCONST(2268.0)*tau + RCONST(405.0))/(TWO*h2);
      a[5] = (RCONST(2430.0)*tau2 + RCONST(2592.0)*tau + RCONST(567.0))/(TWO*h2);
    } else if (d == 4) {
      a[0] = (RCONST(6480.0)*tau + RCONST(3240.0))/h4;
      a[1] = -a[0];
      a[2] = (RCONST(810.0)*tau + RCONST(378.0))/h3;
      a[3] = (RCONST(810.0)*tau + RCONST(432.0))/h3;
      a[4] = (RCONST(2430.0)*tau + RCONST(1134.0))/h3;
      a[5] = (RCONST(2430.0)*tau + RCONST(1296.0))/h3;
    } else {  /* d == 5 */
      a[0] = RCONST(6480.0)/h5;
      a[1] = -a[0];
      a[2] = RCONST(810.0)/h4;
      a[3] = a[2];
      a[4] = RCONST(2430.0)/h4;
      a[5] = a[4];
    }
    X[0] = interp->yold;
    X[1] = interp->ynew;
    X[2] = interp->fold;
    X[3] = interp->fnew;
    X[4] = interp->fa;
    X[5] = interp->fb;
    retval = N_VLinearCombination(6, a, X, yout);
    if (retval != 0) return(ARK_VECTOROP_ERR);
    break;

  /* case(6):    /\* hextic interpolant *\/ */

  /*   /\* first, evaluate quartic interpolant at tau=-1/3 *\/ */
  /*   tval = -ONE/THREE; */
  /*   retval = arkInterpEvaluate(arkode_mem, interp, tval, 0, 5, yout); */
  /*   if (retval != 0)  return(ARK_RHSFUNC_FAIL); */

  /*   /\* second, evaluate RHS at tau=-1/3, storing the result in fa *\/ */
  /*   tval = interp->tnew - h/THREE; */
  /*   if (SUNRabs(tval - interp->t_fa) > 10.0*UNIT_ROUNDOFF) { */
  /*     retval = ark_mem->step_fullrhs(ark_mem, tval, yout, interp->fa, 2); */
  /*     if (retval != 0)  return(ARK_RHSFUNC_FAIL); */
  /*     interp->t_fa = tval; */
  /*   } */

  /*   /\* third, evaluate quartic interpolant at tau=-2/3 *\/ */
  /*   tval = -TWO/THREE; */
  /*   retval = arkInterpEvaluate(arkode_mem, interp, tval, 0, 5, yout); */
  /*   if (retval != 0)  return(ARK_RHSFUNC_FAIL); */

  /*   /\* fourth, evaluate RHS at tau=-2/3, storing the result in fb *\/ */
  /*   tval = interp->tnew - h*TWO/THREE; */
  /*   if (SUNRabs(tval - interp->t_fb) > 10.0*UNIT_ROUNDOFF) { */
  /*     retval = ark_mem->step_fullrhs(ark_mem, tval, yout, interp->fb, 2); */
  /*     if (retval != 0)  return(ARK_RHSFUNC_FAIL); */
  /*     interp->t_fb = tval; */
  /*   } */

  /*   /\* fifth, evaluate quartic interpolant at tau=-1/4 *\/ */
  /*   tval = -ONE/FOUR; */
  /*   retval = arkInterpEvaluate(arkode_mem, interp, tval, 0, 5, yout); */
  /*   if (retval != 0)  return(ARK_RHSFUNC_FAIL); */

  /*   /\* fourth, evaluate RHS at tau=-2/3, storing the result in fc *\/ */
  /*   tval = interp->tnew - h*ONE/FOUR; */
  /*   if (SUNRabs(tval - interp->t_fc) > 10.0*UNIT_ROUNDOFF) { */
  /*     retval = ark_mem->step_fullrhs(ark_mem, tval, yout, interp->fc, 2); */
  /*     if (retval != 0)  return(ARK_RHSFUNC_FAIL); */
  /*     interp->t_fc = tval; */
  /*   } */

  /*   /\* evaluate desired function *\/ */
  /*   if (d == 0) { */
  /*     a[0] = (-RCONST(6924.0)*tau6 - RCONST(16506.0)*tau5 - RCONST(12415.0)*tau4 */
  /*             - RCONST(2850.0)*tau3 + RCONST(62.0)*tau2)/RCONST(79.0); */
  /*     a[1] = (RCONST(6924.0)*tau6 + RCONST(16506.0)*tau5 + RCONST(12415.0)*tau4 */
  /*             + RCONST(2850.0)*tau3 - RCONST(62.0)*tau2 + RCONST(79.0))/RCONST(79.0); */
  /*     a[2] = h*(-RCONST(918.0)*tau6 - RCONST(621.0)*tau5 + RCONST(1917.0)*tau4 */
  /*               + RCONST(2341.0)*tau3 + RCONST(721.0)*tau2)/RCONST(316.0); */
  /*     a[3] = h*(-RCONST(2166.0)*tau6 - RCONST(4365.0)*tau5 - RCONST(1532.0)*tau4 */
  /*               + RCONST(1683.0)*tau3 + RCONST(1332.0)*tau2 + RCONST(316.0)*tau)/RCONST(316.0); */
  /*     a[4] = h*(-RCONST(2754.0)*tau6 - RCONST(1863.0)*tau5 + RCONST(5751.0)*tau4 */
  /*               + RCONST(6075.0)*tau3 + RCONST(1215.0)*tau2)/RCONST(316.0); */
  /*     a[5] = h*(-RCONST(15714.0)*tau6 - RCONST(40743.0)*tau5 - RCONST(35316.0)*tau4 */
  /*               - RCONST(11259.0)*tau3 - RCONST(9072.0)*tau2)/RCONST(316.0); */
  /*     a[6] = h*(-RCONST(1536.0)*tau6 - RCONST(4608.0)*tau5 - RCONST(5120.0)*tau4 */
  /*               - RCONST(2560.0)*tau3 - RCONST(512.0)*tau2)/RCONST(79.0); */
  /*   } else if (d == 1) { */
  /*     a[0] = (-RCONST(41544.0)*tau5 - RCONST(82530.0)*tau4 - RCONST(49660.0)*tau3 */
  /*             - RCONST(8550.0)*tau2 + RCONST(124.0)*tau)/(RCONST(79.0)*h); */
  /*     a[1] = (RCONST(41544.0)*tau5 + RCONST(82530.0)*tau4 + RCONST(49660.0)*tau3 */
  /*             + RCONST(8550.0)*tau2 - RCONST(124.0)*tau)/(RCONST(79.0)*h); */
  /*     a[2] = (-RCONST(5508.0)*tau5 - RCONST(3105.0)*tau4 + RCONST(7668.0)*tau3 */
  /*             + RCONST(7023.0)*tau2 + RCONST(1442.0)*tau)/RCONST(316.0); */
  /*     a[3] = (-RCONST(12996.0)*tau5 - RCONST(21825.0)*tau4 - RCONST(6128.0)*tau3 */
  /*             + RCONST(5049.0)*tau2 + RCONST(2664.0)*tau + RCONST(316.0))/RCONST(316.0); */
  /*     a[4] = (-RCONST(16524.0)*tau5 - RCONST(9315.0)*tau4 + RCONST(22604.0)*tau3 */
  /*             + RCONST(18225.0)*tau2 + RCONST(2430.0)*tau)/RCONST(316.0); */
  /*     a[5] = (-RCONST(94284.0)*tau5 - RCONST(203715.0)*tau4 - RCONST(141264.0)*tau3 */
  /*             - RCONST(33777.0)*tau2 - RCONST(1944.0)*tau)/RCONST(316.0); */
  /*     a[6] = (-RCONST(9216.0)*tau5 - RCONST(23040.0)*tau4 - RCONST(20480.0)*tau3 */
  /*             - RCONST(7680.0)*tau2 - RCONST(1024.0)*tau)/RCONST(79.0); */
  /*   } else if (d == 2) { */
  /*     a[0] = (-RCONST(207720.0)*tau4 - RCONST(330120.0)*tau3 - RCONST(148980.0)*tau2 */
  /*             - RCONST(17100.0)*tau + RCONST(124.0))/(RCONST(79.0)*h2); */
  /*     a[1] = (RCONST(207720.0)*tau4 + RCONST(330120.0)*tau3 + RCONST(148980.0)*tau2 */
  /*             + RCONST(17100.0)*tau - RCONST(124.0))/(RCONST(79.0)*h2); */
  /*     a[2] = (-RCONST(13770.0)*tau4 - RCONST(6210.0)*tau3 + RCONST(11502.0)*tau2 */
  /*             + RCONST(7023.0)*tau + RCONST(721.0))/(RCONST(158.0)*h); */
  /*     a[3] = (-RCONST(32490.0)*tau4 - RCONST(43650.0)*tau3 - RCONST(9192.0)*tau2 */
  /*             + RCONST(5049.0)*tau + RCONST(1332.0))/(RCONST(158.0)*h); */
  /*     a[4] = (-RCONST(41310.0)*tau4 - RCONST(18630.0)*tau3 + RCONST(34506.0)*tau2 */
  /*             + RCONST(18225.0)*tau + RCONST(1215.0))/(RCONST(158.0)*h); */
  /*     a[5] = (-RCONST(235710.0)*tau4 - RCONST(407430.0)*tau3 - RCONST(211896.0)*tau2 */
  /*             - RCONST(33777.0)*tau - RCONST(972.0))/(RCONST(158.0)*h); */
  /*     a[6] = (-RCONST(46080.0)*tau4 - RCONST(92160.0)*tau3 - RCONST(61440.0)*tau2 */
  /*             - RCONST(15360.0)*tau - RCONST(1024.0))/(RCONST(79.0)*h); */
  /*   } else if (d == 3) { */
  /*     a[0] = (-RCONST(830880.0)*tau3 - RCONST(990360.0)*tau2 - RCONST(297960.0)*tau */
  /*             - RCONST(17100.0))/(RCONST(79.0)*h3); */
  /*     a[1] = (RCONST(830880.0)*tau3 + RCONST(990360.0)*tau2 + RCONST(297960.0)*tau */
  /*             + RCONST(17100.0))/(RCONST(79.0)*h3); */
  /*     a[2] = (-RCONST(55080.0)*tau3 - RCONST(18630.0)*tau2 + RCONST(23004.0)*tau */
  /*             + RCONST(7023.0))/(RCONST(158.0)*h2); */
  /*     a[3] = (-RCONST(129960.0)*tau3 - RCONST(130950.0)*tau2 - RCONST(18384.0)*tau */
  /*             + RCONST(5049.0))/(RCONST(158.0)*h2); */
  /*     a[4] = (-RCONST(165240.0)*tau3 - RCONST(55890.0)*tau2 + RCONST(69012.0)*tau */
  /*             + RCONST(18225.0))/(RCONST(158.0)*h2); */
  /*     a[5] = (-RCONST(942840.0)*tau3 - RCONST(1222290.0)*tau2 - RCONST(423792.0)*tau */
  /*             - RCONST(33777.0))/(RCONST(158.0)*h2); */
  /*     a[6] = (-RCONST(184320.0)*tau3 - RCONST(276480.0)*tau2 - RCONST(122880.0)*tau */
  /*             - RCONST(15360.0))/(RCONST(79.0)*h2); */
  /*   } else if (d == 4) { */
  /*     a[0] = (-RCONST(2492640.0)*tau2 - RCONST(1980720.0)*tau - RCONST(297960.0))/(RCONST(79.0)*h4); */
  /*     a[1] = (RCONST(2492640.0)*tau2 + RCONST(1980720.0)*tau + RCONST(297960.0))/(RCONST(79.0)*h4); */
  /*     a[2] = (-RCONST(82620.0)*tau2 - RCONST(18630.0)*tau + RCONST(11502.0))/(RCONST(79.0)*h3); */
  /*     a[3] = (-RCONST(194940.0)*tau2 - RCONST(130950.0)*tau - RCONST(9192.0))/(RCONST(79.0)*h3); */
  /*     a[4] = (-RCONST(247860.0)*tau2 - RCONST(55890.0)*tau + RCONST(34506.0))/(RCONST(79.0)*h3); */
  /*     a[5] = (-RCONST(1414260.0)*tau2 - RCONST(1222290.0)*tau - RCONST(211896.0))/(RCONST(79.0)*h3); */
  /*     a[6] = (-RCONST(552960.0)*tau2 - RCONST(552960.0)*tau - RCONST(122880.0))/(RCONST(79.0)*h3); */
  /*   } else if (d == 5) { */
  /*     a[0] = (-RCONST(4985280.0)*tau - RCONST(1980720.0))/(RCONST(79.0)*h5); */
  /*     a[1] = (RCONST(4985280.0)*tau + RCONST(1980720.0))/(RCONST(79.0)*h5); */
  /*     a[2] = (-RCONST(165240.0)*tau - RCONST(18630.0))/(RCONST(79.0)*h4); */
  /*     a[3] = (-RCONST(389880.0)*tau - RCONST(130950.0))/(RCONST(79.0)*h4); */
  /*     a[4] = (-RCONST(495720.0)*tau - RCONST(55890.0))/(RCONST(79.0)*h4); */
  /*     a[5] = (-RCONST(2828520.0)*tau - RCONST(1222290.0))/(RCONST(79.0)*h4); */
  /*     a[6] = (-RCONST(1105920.0)*tau - RCONST(552960.0))/(RCONST(79.0)*h4); */
  /*   } else {  /\* d == 6 *\/ */
  /*     a[0] = -RCONST(4985280.0)/(RCONST(79.0)*h6); */
  /*     a[1] = RCONST(4985280.0)/(RCONST(79.0)*h6); */
  /*     a[2] = -RCONST(165240.0)/(RCONST(79.0)*h5); */
  /*     a[3] = -RCONST(389880.0)/(RCONST(79.0)*h5); */
  /*     a[4] = -RCONST(495720.0)/(RCONST(79.0)*h5); */
  /*     a[5] = -RCONST(2828520.0)/(RCONST(79.0)*h5); */
  /*     a[6] = -RCONST(1105920.0)/(RCONST(79.0)*h5); */
  /*   } */
  /*   X[0] = interp->yold; */
  /*   X[1] = interp->ynew; */
  /*   X[2] = interp->fold; */
  /*   X[3] = interp->fnew; */
  /*   X[4] = interp->fa; */
  /*   X[5] = interp->fb; */
  /*   X[6] = interp->fc; */
  /*   retval = N_VLinearCombination(7, a, X, yout); */
  /*   if (retval != 0) return(ARK_VECTOROP_ERR); */
  /*   break; */

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInterpEvaluate",
                    "Illegal polynomial order");
    return (ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*===============================================================
   EOF
===============================================================*/
