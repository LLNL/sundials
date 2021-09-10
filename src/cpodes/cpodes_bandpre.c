/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file contains implementations of the banded difference
 * quotient Jacobian-based preconditioner and solver routines for
 * use with the CPSPILS linear solvers.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpodes/cpodes_sptfqmr.h>
#include <cpodes/cpodes_spbcgs.h>
#include <cpodes/cpodes_spgmr.h>

#include <sundials/sundials_math.h>

#include "cpodes_bandpre_impl.h"
#include "cpodes_private.h"
#include "cpodes_spils_impl.h"

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

#define MIN_INC_MULT RCONST(1000.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* cpBandPrecSetupExpl and cpBandPrecSetupImpl */
  
static int cpBandPrecSetupExpl(realtype t, N_Vector y, N_Vector fy, 
                               booleantype jok, booleantype *jcurPtr, 
                               realtype gamma, void *bp_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int cpBandPrecSetupImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                               realtype gamma, void *bp_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* cpBandPrecSolveExpl and cpBandPrecSolveImpl */

static int cpBandPrecSolveExpl(realtype t, N_Vector y, N_Vector fy, 
                               N_Vector b, N_Vector x, 
                               realtype gamma, realtype delta,
                               int lr, void *bp_data, N_Vector tmp);

static int cpBandPrecSolveImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                               N_Vector b, N_Vector x,
                               realtype gamma, realtype delta, 
                               void *bp_data, N_Vector tmp);

/* Prototype for cpBandPrecFree */

static int cpBandPrecFree(CPodeMem cp_mem);

/* Difference quotient Jacobian calculation routines */

static int cpBandPDQJacExpl(CPBandPrecData pdata, 
                            realtype t, N_Vector y, N_Vector fy, 
                            N_Vector ftemp, N_Vector ytemp);

static int cpBandPDQJacImpl(CPBandPrecData pdata, 
                            realtype t, realtype gamma,
                            N_Vector y, N_Vector yp, N_Vector r, 
                            N_Vector ftemp, N_Vector ytemp, N_Vector yptemp);

/* Redability replacements */

#define ode_type (cp_mem->cp_ode_type)
#define vec_tmpl (cp_mem->cp_tempv)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

int CPBandPrecInit(void *cpode_mem, int N, int mu, int ml)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  CPBandPrecData pdata;
  int mup, mlp, storagemu;
  int flag=0;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPBANDPRE", "CPBandPrecInit", MSGBP_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (cp_mem->cp_lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPBANDPRE", "CPBandPrecInit", MSGBP_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) cp_mem->cp_lmem;

  /* Test if the NVECTOR package is compatible with the BAND preconditioner */
    if(vec_tmpl->ops->nvgetarraypointer == NULL) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPBANDPRE", "CPBandPrecInit", MSGBP_BAD_NVECTOR);
    return(CPSPILS_ILL_INPUT);
  }

  pdata = NULL;
  pdata = (CPBandPrecData) malloc(sizeof *pdata);  /* Allocate data memory */
  if (pdata == NULL) {
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPBANDPRE", "CPBandPrecInit", MSGBP_MEM_FAIL);
    return(CPSPILS_MEM_FAIL);
  }

  /* Load pointers and bandwidths into pdata block. */
  pdata->cpode_mem = cpode_mem;
  pdata->N = N;
  pdata->mu = mup = SUNMIN(N-1, SUNMAX(0,mu));
  pdata->ml = mlp = SUNMIN(N-1, SUNMAX(0,ml));

  /* Initialize nfeBP counter */
  pdata->nfeBP = 0;

  /* Allocate memory for saved banded Jacobian approximation. */
  pdata->savedJ = NULL;
  pdata->savedJ = NewBandMat(N, mup, mlp, mup);
  if (pdata->savedJ == NULL) {
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPBANDPRE", "CPBandPrecInit", MSGBP_MEM_FAIL);
    return(CPSPILS_MEM_FAIL);
  }

  /* Allocate memory for banded preconditioner. */
  storagemu = SUNMIN(N-1, mup+mlp);
  pdata->savedP = NULL;
  pdata->savedP = NewBandMat(N, mup, mlp, storagemu);
  if (pdata->savedP == NULL) {
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPBANDPRE", "CPBandPrecInit", MSGBP_MEM_FAIL);
    return(CPSPILS_MEM_FAIL);
  }

  /* Allocate memory for pivot array. */
  pdata->pivots = NULL;
  pdata->pivots = NewIndexArray(N);
  if (pdata->savedJ == NULL) {
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPBANDPRE", "CPBandPrecInit", MSGBP_MEM_FAIL);
    return(CPSPILS_MEM_FAIL);
  }

  /* Overwrite the P_data field in the SPILS memory */
  cpspils_mem->s_P_data = pdata;

  /* Attach the pfree function */
  cpspils_mem->s_pfree = cpBandPrecFree;

  /* Attach preconditioner solve and setup functions */
  switch (ode_type) {
  case CP_EXPL:
    flag = CPSpilsSetPrecFnExpl(cpode_mem, cpBandPrecSetupExpl, cpBandPrecSolveExpl);
    break;
  case CP_IMPL:
    flag = CPSpilsSetPrecFnImpl(cpode_mem, cpBandPrecSetupImpl, cpBandPrecSolveImpl);
    break;
  }

  return(flag);
}

int CPBandPrecGetWorkSpace(void *cpode_mem, long int *lenrwBP, long int *leniwBP)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  CPBandPrecData pdata;
  int N, ml, mu, smu;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPBANDPRE", "CPBandPrecGetWorkSpace", MSGBP_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (cp_mem->cp_lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPBANDPRE", "CPBandPrecGetWorkSpace", MSGBP_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) cp_mem->cp_lmem;

  if (cpspils_mem->s_P_data == NULL) {
    cpProcessError(cp_mem, CPSPILS_PMEM_NULL, "CPBANDPRE", "CPBandPrecGetWorkSpace", MSGBP_PMEM_NULL);
    return(CPSPILS_PMEM_NULL);
  } 
  pdata = (CPBandPrecData) cpspils_mem->s_P_data;

  N   = pdata->N;
  mu  = pdata->mu;
  ml  = pdata->ml;
  smu = SUNMIN( N-1, mu + ml);

  *leniwBP = pdata->N;
  *lenrwBP = N * ( 2*ml + smu + mu + 2 );

  return(CPSPILS_SUCCESS);
}

int CPBandPrecGetNumRhsEvals(void *cpode_mem, long int *nfevalsBP)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  CPBandPrecData pdata;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPBANDPRE", "CPBandPrecGetNumRhsEvals", MSGBP_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (cp_mem->cp_lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPBANDPRE", "CPBandPrecGetNumRhsEvals", MSGBP_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) cp_mem->cp_lmem;

  if (cpspils_mem->s_P_data == NULL) {
    cpProcessError(cp_mem, CPSPILS_PMEM_NULL, "CPBANDPRE", "CPBandPrecGetNumRhsEvals", MSGBP_PMEM_NULL);
    return(CPSPILS_PMEM_NULL);
  } 
  pdata = (CPBandPrecData) cpspils_mem->s_P_data;

  *nfevalsBP = pdata->nfeBP;

  return(CPSPILS_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS
 * =================================================================
 */

/* Readability Replacements */

#define N      (pdata->N)
#define mu     (pdata->mu)
#define ml     (pdata->ml)
#define pivots (pdata->pivots)
#define savedJ (pdata->savedJ)
#define savedP (pdata->savedP)
#define nfeBP  (pdata->nfeBP)

/*
 * -----------------------------------------------------------------
 * cpBandPrecSetupExpl
 * -----------------------------------------------------------------
 * Together cpBandPrecSetupExpl and cpBandPrecSolveExpl use a banded
 * difference quotient Jacobian to create a preconditioner.
 * cpBandPrecSetupExpl calculates a new J, if necessary, then
 * calculates P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of cpBandPrecSetupExpl are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *         namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == SUNFALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == SUNTRUE means that Jacobian data from the
 *                  previous cpBandPrecSetupExpl call will be reused
 *                  (with the current value of gamma).
 *         A cpBandPrecSetupExpl call with jok == SUNTRUE should only
 *         occur after a call with jok == SUNFALSE.
 *
 * *jcurPtr is a pointer to an output integer flag which is
 *          set by cpBandPrecsetupExpl as follows:
 *            *jcurPtr = SUNTRUE if Jacobian data was recomputed.
 *            *jcurPtr = SUNFALSE if Jacobian data was not recomputed,
 *                       but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * bp_data is a pointer to preconditoner data - the same as the
 *         bp_data parameter passed to CPSp*.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *           for vectors of length N for work space. This
 *           routine uses only tmp1 and tmp2.
 *
 * The value to be returned by the cpBandPrecSetupExpl function is
 *   0  if successful, or
 *   1  if the band factorization failed.
 * -----------------------------------------------------------------
 */

static int cpBandPrecSetupExpl(realtype t, N_Vector y, N_Vector fy, 
                               booleantype jok, booleantype *jcurPtr, 
                               realtype gamma, void *bp_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CPBandPrecData pdata;
  CPodeMem cp_mem;
  int ier, retval;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  if (jok) {

    /* If jok = SUNTRUE, use saved copy of J. */
    *jcurPtr = SUNFALSE;
    BandCopy(savedJ, savedP, mu, ml);

  } else {

    /* If jok = SUNFALSE, call cpBandPDQJac for new J value. */
    *jcurPtr = SUNTRUE;
    retval = cpBandPDQJacExpl(pdata, t, y, fy, tmp1, tmp2);
    if (retval < 0) {
      cpProcessError(cp_mem, -1, "CPBANDPRE", "cpBandPrecSetupExpl", MSGBP_FUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    BandCopy(savedJ, savedP, mu, ml);

  }
  
  /* Scale and add I to get savedP = I - gamma*J. */
  BandScale(-gamma, savedP);
  AddIdentity(savedP);
 
  /* Do LU factorization of matrix. */
  ier = BandGBTRF(savedP, pivots);
 
  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandPrecSolveExpl
 * -----------------------------------------------------------------
 * cpBandPrecSolveExpl solves a linear system P x = b, where P is the
 * matrix computed by cpBandPrecSetupExpl.
 *
 * The parameters of cpBandPrecSolveExpl used here are as follows:
 *
 * b       is the right-hand side vector of the linear system.
 *
 * bp_data is a pointer to preconditioner data - the same as the
 *         bp_data parameter passed to CPSp*.
 *
 * x       is the output vector computed by cpBandPrecSolveExpl.
 *
 * The value returned by the cpBandPrecSolveExpl function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */ 

static int cpBandPrecSolveExpl(realtype t, N_Vector y, N_Vector fy, 
                               N_Vector b, N_Vector x, 
                               realtype gamma, realtype delta,
                               int lr, void *bp_data, N_Vector tmp)
{
  CPBandPrecData pdata;
  realtype *xd;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  /* Copy b to x. */
  N_VScale(ONE, b, x);

  /* Do band backsolve on the vector z. */
  xd = N_VGetArrayPointer(x);

  BandGBTRS(savedP, pivots, xd);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandPrecSetupImpl
 * -----------------------------------------------------------------
 * Together cpBandPrecSetupImpl and cpBandPrecSolveImpl use a banded
 * difference quotient Jacobian to create a preconditioner.
 * cpBandPrecSetupImpl calculates a new J = dF/dy + gamma*dF/dy'
 * and does an LU factorization of P.
 * -----------------------------------------------------------------
 */

static int cpBandPrecSetupImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                               realtype gamma, void *bp_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CPBandPrecData pdata;
  CPodeMem cp_mem;
  int ier, retval;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  retval = cpBandPDQJacImpl(pdata, t, gamma, y, yp, r, tmp1, tmp2, tmp3);
  if (retval < 0) {
    cpProcessError(cp_mem, -1, "CPBANDPRE", "cpBandPrecSetupImpl", MSGBP_FUNC_FAILED);
    return(-1);
  }
  if (retval > 0) {
    return(1);
  }

  /* Do LU factorization of matrix. */
  ier = BandGBTRF(savedP, pivots);
 
  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandPrecSolveImpl
 * -----------------------------------------------------------------
 * cpBandPrecSolveImpl solves a linear system P x = b, where P is the
 * matrix computed by cpBandPrecSetupImpl.
 *
 * The parameters of cpBandPrecSolveImpl used here are as follows:
 *
 * b       is the right-hand side vector of the linear system.
 *
 * bp_data is a pointer to preconditioner data - the same as the
 *         bp_data parameter passed to CPSp*.
 *
 * x       is the output vector computed by cpBandPrecSolveImpl.
 *
 * The value returned by the cpBandPrecSolveImpl function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */ 

static int cpBandPrecSolveImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                               N_Vector b, N_Vector x,
                               realtype gamma, realtype delta, 
                               void *bp_data, N_Vector tmp)
{
  CPBandPrecData pdata;
  realtype *xd;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  /* Copy b to x. */
  N_VScale(ONE, b, x);

  /* Do band backsolve on the vector z. */
  xd = N_VGetArrayPointer(x);

  BandGBTRS(savedP, pivots, xd);

  return(0);
}


/*
 * -----------------------------------------------------------------
 * Function to free memory allocated by CPBandPrecInit
 * -----------------------------------------------------------------
 */

static int cpBandPrecFree(CPodeMem cp_mem)
{
  CPSpilsMem cpspils_mem;
  CPBandPrecData pdata;

  if (cp_mem->cp_lmem == NULL) return(0);
  cpspils_mem = (CPSpilsMem) cp_mem->cp_lmem;
  
  if (cpspils_mem->s_P_data == NULL) return(0);
  pdata = (CPBandPrecData) cpspils_mem->s_P_data;

  DestroyMat(savedJ);
  DestroyMat(savedP);
  DestroyArray(pivots);

  free(pdata);
  pdata = NULL;

  return(0);
}

/* 
 * =================================================================
 * DQ LOCAL JACOBIAN APROXIMATIONS
 * =================================================================
 */

#define ewt       (cp_mem->cp_ewt)
#define uround    (cp_mem->cp_uround)
#define h         (cp_mem->cp_h)
#define fe        (cp_mem->cp_fe)
#define fi        (cp_mem->cp_fi)
#define user_data (cp_mem->cp_user_data)

/*
 * -----------------------------------------------------------------
 * cpBandPDQJacExpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a band matrix of type
 * BandMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

static int cpBandPDQJacExpl(CPBandPrecData pdata, 
                            realtype t, N_Vector y, N_Vector fy, 
                            N_Vector ftemp, N_Vector ytemp)
{
  CPodeMem cp_mem;
  realtype fnorm, minInc, inc, inc_inv, srur;
  int group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int retval;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp. */
  ewt_data   = N_VGetArrayPointer(ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector. */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f. */
  srur = SUNRsqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * SUNRabs(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing. */
  width = ml + mu + 1;
  ngroups = SUNMIN(width, N);

  /* Loop over column groups. */
  for (group = 1; group <= ngroups; group++) {
    
    /* Increment all y_j in group. */
    for(j = group-1; j < N; j += width) {
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y. */

    retval = fe(t, ytemp, ftemp, user_data);
    nfeBP++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients. */
    for (j = group-1; j < N; j += width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(savedJ,j);
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mu);
      i2 = SUNMIN(j+ml, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandPDQJacImpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation
 * to the Jacobian dF/dy + gamma*dF/dy' and loads it directly into
 * savedP. It assumes that a band matrix of type BandMat is stored
 * column-wise, and that elements within each column are contiguous.
 * This makes it possible to get the address of a column of J via 
 * the macro BAND_COL and to write a simple for loop to set each of
 * the elements of a column in succession.
 * -----------------------------------------------------------------
 */

static int cpBandPDQJacImpl(CPBandPrecData pdata, 
                            realtype t, realtype gamma,
                            N_Vector y, N_Vector yp, N_Vector r, 
                            N_Vector ftemp, N_Vector ytemp, N_Vector yptemp)

{
  CPodeMem cp_mem;
  realtype inc, inc_inv, yj, ypj, srur, ewtj;
  realtype *y_data, *yp_data, *ewt_data;
  realtype *ytemp_data, *yptemp_data, *ftemp_data, *r_data, *col_j;
  int i, j, i1, i2, width, group, ngroups;
  int retval = 0;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Obtain pointers to the data for all vectors used.  */
  ewt_data    = N_VGetArrayPointer(ewt);
  r_data      = N_VGetArrayPointer(r);
  y_data      = N_VGetArrayPointer(y);
  yp_data     = N_VGetArrayPointer(yp);
  ftemp_data  = N_VGetArrayPointer(ftemp);
  ytemp_data  = N_VGetArrayPointer(ytemp);
  yptemp_data = N_VGetArrayPointer(yptemp);

  /* Initialize ytemp and yptemp. */
  N_VScale(ONE, y, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Compute miscellaneous values for the Jacobian computation. */
  srur = SUNRsqrt(uround);
  width = ml + mu + 1;
  ngroups = SUNMIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {

    /* Increment all y[j] and yp[j] for j in this group. */
    for (j=group-1; j<N; j+=width) {
        yj = y_data[j];
        ypj = yp_data[j];
        ewtj = ewt_data[j];

        /* Set increment inc to yj based on sqrt(uround)*abs(yj), with
           adjustments using ypj and ewtj if this is small, and a further
           adjustment to give it the same sign as h*ypj. */
        inc = SUNMAX( srur * SUNMAX( SUNRabs(yj), SUNRabs(h*ypj) ) , ONE/ewtj );

        if (h*ypj < ZERO) inc = -inc;
        inc = (yj + inc) - yj;

        /* Increment yj and ypj. */
        ytemp_data[j]  += gamma*inc;
        yptemp_data[j] += inc;
    }

    /* Call ODE fct. with incremented arguments. */
    retval = fi(t, ytemp, yptemp, ftemp, user_data);
    nfeBP++;
    if (retval != 0) break;

    /* Loop over the indices j in this group again. */
    for (j=group-1; j<N; j+=width) {

      /* Reset ytemp and yptemp components that were perturbed. */
      yj = ytemp_data[j]  = y_data[j];
      ypj = yptemp_data[j] = yp_data[j];
      col_j = BAND_COL(savedP, j);
      ewtj = ewt_data[j];
      
      /* Set increment inc exactly as above. */
      inc = SUNMAX( srur * SUNMAX( SUNRabs(yj), SUNRabs(h*ypj) ) , ONE/ewtj );
      if (h*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      
      /* Load the difference quotient Jacobian elements for column j. */
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mu);
      i2 = SUNMIN(j+ml,N-1);
      for (i=i1; i<=i2; i++) 
        BAND_COL_ELEM(col_j,i,j) = inc_inv*(ftemp_data[i]-r_data[i]);
    }
    
  }
  
  return(retval);
}
