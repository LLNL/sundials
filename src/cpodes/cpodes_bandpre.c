/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:07:05 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This file contains implementations of the banded difference
 * quotient Jacobian-based preconditioner and solver routines for
 * use with the CPSPILS linear solvers..
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_impl.h"
#include "cpodes_bandpre_impl.h"

#include <cpodes/cpodes_sptfqmr.h>
#include <cpodes/cpodes_spbcgs.h>
#include <cpodes/cpodes_spgmr.h>

#include <sundials/sundials_math.h>

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* Prototypes of CPBandPrecSetup and CPBandPrecSolve */
  
static int CPBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int CPBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bp_data, N_Vector tmp);

/* Prototype for difference quotient Jacobian calculation routine */

static int CPBandPDQJac(CPBandPrecData pdata, 
                        realtype t, N_Vector y, N_Vector fy, 
                        N_Vector ftemp, N_Vector ytemp);

/* Redability replacements */

#define vec_tmpl (cp_mem->cp_tempv)

/*
 * -----------------------------------------------------------------
 * Malloc, Free, and Get Functions
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPBandPrecAlloc will
 *       first test for a compatible N_Vector internal representation
 *       by checking that the function N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */

void *CPBandPrecAlloc(void *cpode_mem, long int N, 
                      long int mu, long int ml)
{
  CPodeMem cp_mem;
  CPBandPrecData pdata;
  long int mup, mlp, storagemu;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_CPMEM_NULL);
    return(NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the BAND preconditioner */
  if(vec_tmpl->ops->nvgetarraypointer == NULL) {
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_BAD_NVECTOR);
    return(NULL);
  }

  pdata = NULL;
  pdata = (CPBandPrecData) malloc(sizeof *pdata);  /* Allocate data memory */
  if (pdata == NULL) {
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_MEM_FAIL);
    return(NULL);
  }

  /* Load pointers and bandwidths into pdata block. */
  pdata->cpode_mem = cpode_mem;
  pdata->N = N;
  pdata->mu = mup = MIN(N-1, MAX(0,mu));
  pdata->ml = mlp = MIN(N-1, MAX(0,ml));

  /* Initialize nfeBP counter */
  pdata->nfeBP = 0;

  /* Allocate memory for saved banded Jacobian approximation. */
  pdata->savedJ = NULL;
  pdata->savedJ = BandAllocMat(N, mup, mlp, mup);
  if (pdata->savedJ == NULL) {
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_MEM_FAIL);
    return(NULL);
  }

  /* Allocate memory for banded preconditioner. */
  storagemu = MIN(N-1, mup+mlp);
  pdata->savedP = NULL;
  pdata->savedP = BandAllocMat(N, mup, mlp, storagemu);
  if (pdata->savedP == NULL) {
    BandFreeMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_MEM_FAIL);
    return(NULL);
  }

  /* Allocate memory for pivot array. */
  pdata->pivots = NULL;
  pdata->pivots = BandAllocPiv(N);
  if (pdata->savedJ == NULL) {
    BandFreeMat(pdata->savedP);
    BandFreeMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_MEM_FAIL);
    return(NULL);
  }

  return((void *) pdata);
}

int CPBPSptfqmr(void *cpode_mem, int pretype, int maxl, void *p_data)
{
  CPodeMem cp_mem;
  int flag;

  flag = CPSptfqmr(cpode_mem, pretype, maxl);
  if(flag != CPSPILS_SUCCESS) return(flag);
  
  cp_mem = (CPodeMem) cpode_mem;

  if ( p_data == NULL ) {
    cpProcessError(cp_mem, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBPSptfqmr", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  } 

  flag = CPSpilsSetPreconditioner(cpode_mem, CPBandPrecSetup, CPBandPrecSolve, p_data);
  if(flag != CPSPILS_SUCCESS) return(flag);

  return(CPSPILS_SUCCESS);
}

int CPBPSpbcg(void *cpode_mem, int pretype, int maxl, void *p_data)
{
  CPodeMem cp_mem;
  int flag;

  flag = CPSpbcg(cpode_mem, pretype, maxl);
  if(flag != CPSPILS_SUCCESS) return(flag);

  cp_mem = (CPodeMem) cpode_mem;

  if ( p_data == NULL ) {
    cpProcessError(cp_mem, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBPSpbcg", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  } 

  flag = CPSpilsSetPreconditioner(cpode_mem, CPBandPrecSetup, CPBandPrecSolve, p_data);
  if(flag != CPSPILS_SUCCESS) return(flag);

  return(CPSPILS_SUCCESS);
}

int CPBPSpgmr(void *cpode_mem, int pretype, int maxl, void *p_data)
{
  CPodeMem cp_mem;
  int flag;

  flag = CPSpgmr(cpode_mem, pretype, maxl);
  if(flag != CPSPILS_SUCCESS) return(flag);

  cp_mem = (CPodeMem) cpode_mem;

  if ( p_data == NULL ) {
    cpProcessError(cp_mem, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBPSpgmr", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  }

  flag = CPSpilsSetPreconditioner(cpode_mem, CPBandPrecSetup, CPBandPrecSolve, p_data);
  if(flag != CPSPILS_SUCCESS) return(flag);

  return(CPSPILS_SUCCESS);
}

void CPBandPrecFree(void **bp_data)
{
  CPBandPrecData pdata;

  if (*bp_data == NULL) return;

  pdata = (CPBandPrecData) (*bp_data);
  BandFreeMat(pdata->savedJ);
  BandFreeMat(pdata->savedP);
  BandFreePiv(pdata->pivots);

  free(*bp_data);
  *bp_data = NULL;

}

int CPBandPrecGetWorkSpace(void *bp_data, long int *lenrwBP, long int *leniwBP)
{
  CPBandPrecData pdata;
  long int N, ml, mu, smu;

  if ( bp_data == NULL ) {
    cpProcessError(NULL, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBandPrecGetWorkSpace", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  } 

  pdata = (CPBandPrecData) bp_data;

  N   = pdata->N;
  mu  = pdata->mu;
  ml  = pdata->ml;
  smu = MIN( N-1, mu + ml);

  *leniwBP = pdata->N;
  *lenrwBP = N * ( 2*ml + smu + mu + 2 );

  return(CPBANDPRE_SUCCESS);
}

int CPBandPrecGetNumRhsEvals(void *bp_data, long int *nfevalsBP)
{
  CPBandPrecData pdata;

  if (bp_data == NULL) {
    cpProcessError(NULL, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBandPrecGetNumRhsEvals", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  } 

  pdata = (CPBandPrecData) bp_data;

  *nfevalsBP = pdata->nfeBP;

  return(CPBANDPRE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPBandPrecGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *CPBandPrecGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CPBANDPRE_SUCCESS:
    sprintf(name,"CPBANDPRE_SUCCESS");
    break;   
  case CPBANDPRE_PDATA_NULL:
    sprintf(name,"CPBANDPRE_PDATA_NULL");
    break;
  case CPBANDPRE_RHSFUNC_UNRECVR:
    sprintf(name,"CPBANDPRE_RHSFUNC_UNRECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

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
 * CPBandPrecSetup
 * -----------------------------------------------------------------
 * Together CPBandPrecSetup and CPBandPrecSolve use a banded
 * difference quotient Jacobian to create a preconditioner.
 * CPBandPrecSetup calculates a new J, if necessary, then
 * calculates P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of CPBandPrecSetup are as follows:
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
 *           jok == FALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == TRUE means that Jacobian data from the
 *                  previous PrecSetup call will be reused
 *                  (with the current value of gamma).
 *         A CPBandPrecSetup call with jok == TRUE should only
 *         occur after a call with jok == FALSE.
 *
 * *jcurPtr is a pointer to an output integer flag which is
 *          set by CPBandPrecond as follows:
 *            *jcurPtr = TRUE if Jacobian data was recomputed.
 *            *jcurPtr = FALSE if Jacobian data was not recomputed,
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
 * The value to be returned by the CPBandPrecSetup function is
 *   0  if successful, or
 *   1  if the band factorization failed.
 * -----------------------------------------------------------------
 */

static int CPBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  long int ier;
  CPBandPrecData pdata;
  CPodeMem cp_mem;
  int retval;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  if (jok) {

    /* If jok = TRUE, use saved copy of J. */
    *jcurPtr = FALSE;
    BandCopy(savedJ, savedP, mu, ml);

  } else {

    /* If jok = FALSE, call CPBandPDQJac for new J value. */
    *jcurPtr = TRUE;
    BandZero(savedJ);

    retval = CPBandPDQJac(pdata, t, y, fy, tmp1, tmp2);
    if (retval < 0) {
      cpProcessError(cp_mem, CPBANDPRE_RHSFUNC_UNRECVR, "CPBANDPRE", "CPBandPrecSetup", MSGBP_RHSFUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    BandCopy(savedJ, savedP, mu, ml);

  }
  
  /* Scale and add I to get savedP = I - gamma*J. */
  BandScale(-gamma, savedP);
  BandAddI(savedP);
 
  /* Do LU factorization of matrix. */
  ier = BandGBTRF(savedP, pivots);
 
  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * CPBandPrecSolve
 * -----------------------------------------------------------------
 * CPBandPrecSolve solves a linear system P z = r, where P is the
 * matrix computed by CPBandPrecond.
 *
 * The parameters of CPBandPrecSolve used here are as follows:
 *
 * r       is the right-hand side vector of the linear system.
 *
 * bp_data is a pointer to preconditioner data - the same as the
 *         bp_data parameter passed to CPSp*.
 *
 * z       is the output vector computed by CPBandPrecSolve.
 *
 * The value returned by the CPBandPrecSolve function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */ 

static int CPBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bp_data, N_Vector tmp)
{
  CPBandPrecData pdata;
  realtype *zd;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  /* Copy r to z. */
  N_VScale(ONE, r, z);

  /* Do band backsolve on the vector z. */
  zd = N_VGetArrayPointer(z);

  BandGBTRS(savedP, pivots, zd);

  return(0);
}

#define ewt    (cp_mem->cp_ewt)
#define uround (cp_mem->cp_uround)
#define h      (cp_mem->cp_h)
#define f      (cp_mem->cp_f)
#define f_data (cp_mem->cp_f_data)

/*
 * -----------------------------------------------------------------
 * CPBandPDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a band matrix of type
 * BandMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

static int CPBandPDQJac(CPBandPrecData pdata, 
                        realtype t, N_Vector y, N_Vector fy, 
                        N_Vector ftemp, N_Vector ytemp)
{
  CPodeMem cp_mem;
  realtype fnorm, minInc, inc, inc_inv, srur;
  long int group, i, j, width, ngroups, i1, i2;
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
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing. */
  width = ml + mu + 1;
  ngroups = MIN(width, N);
  
  for (group = 1; group <= ngroups; group++) {
    
    /* Increment all y_j in group. */
    for(j = group-1; j < N; j += width) {
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y. */

    retval = f(t, ytemp, ftemp, f_data);
    nfeBP++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients. */
    for (j = group-1; j < N; j += width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(savedJ,j);
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mu);
      i2 = MIN(j+ml, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

  return(0);
}
