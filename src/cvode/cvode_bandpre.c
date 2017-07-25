/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This file contains implementations of the banded difference
 * quotient Jacobian-based preconditioner and solver routines for
 * use with the CVSPILS linear solvers..
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "cvode_bandpre_impl.h"
#include "cvode_spils_impl.h"

#include <cvode/cvode_sptfqmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_spgmr.h>

#include <sundials/sundials_math.h>

#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* Prototypes of CVBandPrecSetup and CVBandPrecSolve */
  
static int CVBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int CVBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bp_data, N_Vector tmp);

/* Prototype for CVBandPrecFree */

static int CVBandPrecFree(CVodeMem cv_mem);

/* Prototype for difference quotient Jacobian calculation routine */

static int CVBandPDQJac(CVBandPrecData pdata,
                        realtype t, N_Vector y, N_Vector fy, 
                        N_Vector ftemp, N_Vector ytemp);

/*
 * -----------------------------------------------------------------
 * Initialization, Free, and Get Functions
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVBandPrecInit will
 *       first test for a compatible N_Vector internal representation
 *       by checking that the function N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */

int CVBandPrecInit(void *cvode_mem, sunindextype N, sunindextype mu, sunindextype ml)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBandPrecData pdata;
  sunindextype mup, mlp, storagemu;
  int flag;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBANDPRE", "CVBandPrecInit", MSGBP_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Test if the NVECTOR package is compatible with the BAND preconditioner */
  if(cv_mem->cv_tempv->ops->nvgetarraypointer == NULL) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVBANDPRE", "CVBandPrecInit", MSGBP_BAD_NVECTOR);
    return(CVSPILS_ILL_INPUT);
  }

  pdata = NULL;
  pdata = (CVBandPrecData) malloc(sizeof *pdata);  /* Allocate data memory */
  if (pdata == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Load pointers and bandwidths into pdata block. */
  pdata->cvode_mem = cvode_mem;
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
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Allocate memory for banded preconditioner. */
  storagemu = SUNMIN(N-1, mup+mlp);
  pdata->savedP = NULL;
  pdata->savedP = NewBandMat(N, mup, mlp, storagemu);
  if (pdata->savedP == NULL) {
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Allocate memory for pivot array. */
  pdata->lpivots = NULL;
  pdata->lpivots = NewIndexArray(N);
  if (pdata->lpivots == NULL) {
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* make sure s_P_data is free from any previous allocations */
  if (cvspils_mem->s_pfree != NULL) {
    cvspils_mem->s_pfree(cv_mem);
  }

  /* Point to the new P_data field in the SPILS memory */
  cvspils_mem->s_P_data = pdata;

  /* Attach the pfree function */
  cvspils_mem->s_pfree = CVBandPrecFree;

  /* Attach preconditioner solve and setup functions */
  flag = CVSpilsSetPreconditioner(cvode_mem, CVBandPrecSetup, CVBandPrecSolve);

  return(flag);
}

int CVBandPrecGetWorkSpace(void *cvode_mem, sunindextype *lenrwBP, sunindextype *leniwBP)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBandPrecData pdata;
  sunindextype N, ml, mu, smu;

  
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBANDPRE", "CVBandPrecGetWorkSpace", MSGBP_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBANDPRE", "CVBandPrecGetWorkSpace", MSGBP_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  if (cvspils_mem->s_P_data == NULL) {
    cvProcessError(cv_mem, CVSPILS_PMEM_NULL, "CVBANDPRE", "CVBandPrecGetWorkSpace", MSGBP_PMEM_NULL);
    return(CVSPILS_PMEM_NULL);
  } 
  pdata = (CVBandPrecData) cvspils_mem->s_P_data;

  N   = pdata->N;
  mu  = pdata->mu;
  ml  = pdata->ml;
  smu = SUNMIN( N-1, mu + ml);

  *leniwBP = pdata->N;
  *lenrwBP = N * ( 2*ml + smu + mu + 2 );

  return(CVSPILS_SUCCESS);
}

int CVBandPrecGetNumRhsEvals(void *cvode_mem, long int *nfevalsBP)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBandPrecData pdata;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBANDPRE", "CVBandPrecGetNumRhsEvals", MSGBP_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBANDPRE", "CVBandPrecGetNumRhsEvals", MSGBP_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  if (cvspils_mem->s_P_data == NULL) {
    cvProcessError(cv_mem, CVSPILS_PMEM_NULL, "CVBANDPRE", "CVBandPrecGetNumRhsEvals", MSGBP_PMEM_NULL);
    return(CVSPILS_PMEM_NULL);
  } 
  pdata = (CVBandPrecData) cvspils_mem->s_P_data;

  *nfevalsBP = pdata->nfeBP;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVBandPrecSetup
 * -----------------------------------------------------------------
 * Together CVBandPrecSetup and CVBandPrecSolve use a banded
 * difference quotient Jacobian to create a preconditioner.
 * CVBandPrecSetup calculates a new J, if necessary, then
 * calculates P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of CVBandPrecSetup are as follows:
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
 *         A CVBandPrecSetup call with jok == TRUE should only
 *         occur after a call with jok == FALSE.
 *
 * *jcurPtr is a pointer to an output integer flag which is
 *          set by CVBandPrecond as follows:
 *            *jcurPtr = TRUE if Jacobian data was recomputed.
 *            *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                       but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * bp_data is a pointer to preconditoner data (set by CVBandPrecInit)
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *           for vectors of length N for work space. This
 *           routine uses only tmp1 and tmp2.
 *
 * The value to be returned by the CVBandPrecSetup function is
 *   0  if successful, or
 *   1  if the band factorization failed.
 * -----------------------------------------------------------------
 */

static int CVBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CVBandPrecData pdata;
  CVodeMem cv_mem;
  int retval;
  sunindextype ier;

  /* Assume matrix and lpivots have already been allocated. */
  pdata = (CVBandPrecData) bp_data;

  cv_mem = (CVodeMem) pdata->cvode_mem;

  if (jok) {

    /* If jok = TRUE, use saved copy of J. */
    *jcurPtr = FALSE;
    BandCopy(pdata->savedJ, pdata->savedP, pdata->mu, pdata->ml);

  } else {

    /* If jok = FALSE, call CVBandPDQJac for new J value. */
    *jcurPtr = TRUE;
    SetToZero(pdata->savedJ);

    retval = CVBandPDQJac(pdata, t, y, fy, tmp1, tmp2);
    if (retval < 0) {
      cvProcessError(cv_mem, -1, "CVBANDPRE", "CVBandPrecSetup", MSGBP_RHSFUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    BandCopy(pdata->savedJ, pdata->savedP, pdata->mu, pdata->ml);

  }
  
  /* Scale and add I to get savedP = I - gamma*J. */
  BandScale(-gamma, pdata->savedP);
  AddIdentity(pdata->savedP);
 
  /* Do LU factorization of matrix. */
  ier = BandGBTRF(pdata->savedP, pdata->lpivots);
 
  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * CVBandPrecSolve
 * -----------------------------------------------------------------
 * CVBandPrecSolve solves a linear system P z = r, where P is the
 * matrix computed by CVBandPrecond.
 *
 * The parameters of CVBandPrecSolve used here are as follows:
 *
 * r is the right-hand side vector of the linear system.
 *
 * bp_data is a pointer to preconditoner data (set by CVBandPrecInit)
 *
 * z is the output vector computed by CVBandPrecSolve.
 *
 * The value returned by the CVBandPrecSolve function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */ 

static int CVBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bp_data, N_Vector tmp)
{
  CVBandPrecData pdata;
  realtype *zd;

  /* Assume matrix and lpivots have already been allocated. */
  pdata = (CVBandPrecData) bp_data;

  /* Copy r to z. */
  N_VScale(ONE, r, z);

  /* Do band backsolve on the vector z. */
  zd = N_VGetArrayPointer(z);

  BandGBTRS(pdata->savedP, pdata->lpivots, zd);

  return(0);
}


static int CVBandPrecFree(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;
  CVBandPrecData pdata;

  if (cv_mem->cv_lmem == NULL) return(0);
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  
  if (cvspils_mem->s_P_data == NULL) return(0);
  pdata = (CVBandPrecData) cvspils_mem->s_P_data;

  DestroyMat(pdata->savedJ);
  DestroyMat(pdata->savedP);
  DestroyArray(pdata->lpivots);

  free(pdata);
  pdata = NULL;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * CVBandPDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a band matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

static int CVBandPDQJac(CVBandPrecData pdata, 
                        realtype t, N_Vector y, N_Vector fy, 
                        N_Vector ftemp, N_Vector ytemp)
{
  CVodeMem cv_mem;
  realtype fnorm, minInc, inc, inc_inv, srur;
  sunindextype group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int retval;

  cv_mem = (CVodeMem) pdata->cvode_mem;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp. */
  ewt_data   = N_VGetArrayPointer(cv_mem->cv_ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector. */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f. */
  srur = SUNRsqrt(cv_mem->cv_uround);
  fnorm = N_VWrmsNorm(fy, cv_mem->cv_ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * SUNRabs(cv_mem->cv_h) * cv_mem->cv_uround * pdata->N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing. */
  width = pdata->ml + pdata->mu + 1;
  ngroups = SUNMIN(width, pdata->N);
  
  for (group = 1; group <= ngroups; group++) {
    
    /* Increment all y_j in group. */
    for(j = group-1; j < pdata->N; j += width) {
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y. */

    retval = cv_mem->cv_f(t, ytemp, ftemp, cv_mem->cv_user_data);
    pdata->nfeBP++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients. */
    for (j = group-1; j < pdata->N; j += width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(pdata->savedJ,j);
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-pdata->mu);
      i2 = SUNMIN(j+pdata->ml, pdata->N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

  return(0);
}
