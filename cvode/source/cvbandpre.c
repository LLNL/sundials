/*******************************************************************
 *                                                                 *
 * File          : cvbandpre.c                                     *
 * Programmers   : Michael Wittman, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvode/LICENSE                         *
 *-----------------------------------------------------------------*
 * This file contains implementations of the banded difference     *
 * quotient Jacobian-based preconditioner and solver routines for  *
 * use with CVSpgmr.                                               *
 *                                                                 *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cvode.h"
#include "cvbandpre.h"
#include "cvspgmr.h"
#include "sundialsmath.h"
#include "nvector.h"

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* Error Messages */

#define CVBALLOC        "CVBandPreAlloc-- "
#define MSG_CVMEM_NULL  CVBALLOC "Integrator memory is NULL.\n\n"
#define MSG_WRONG_NVEC  CVBALLOC "Incompatible NVECTOR implementation.\n\n"

#define MSG_PDATA_NULL "CVBandPrecGet*-- BandPrecData is NULL. \n\n"

#define MSG_NO_PDATA   "CVBPSpgmr-- BandPrecData is NULL. \n\n"

/* Prototypes of CVBandPrecSetup and CVBandPrecSolve */
  
static int CVBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int CVBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bp_data, N_Vector tmp);

/* Prototype for difference quotient Jacobian calculation routine */

static void CVBandPDQJac(CVBandPrecData pdata, 
                         realtype t, N_Vector y, N_Vector fy, 
                         N_Vector ftemp, N_Vector ytemp);

/* Redability replacements */
#define nvspec (cv_mem->cv_nvspec)
#define errfp  (cv_mem->cv_errfp)

/*********************  Malloc, Free, and Get Functions  *******************/

void *CVBandPrecAlloc(void *cvode_mem, long int N, 
                      long int mu, long int ml)
{
  CVodeMem cv_mem;
  CVBandPrecData pdata;
  long int mup, mlp, storagemu;

  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_CVMEM_NULL);
    return(NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the BAND preconditioner */
  if ((strcmp(nvspec->tag,"serial")) || 
      nvspec->ops->nvgetdata == NULL || 
      nvspec->ops->nvsetdata == NULL) {
    fprintf(errfp, MSG_WRONG_NVEC);
    return(NULL);
  }

  pdata = (CVBandPrecData) malloc(sizeof *pdata);  /* Allocate data memory */
  if (pdata == NULL) return(NULL);

  /* Load pointers and bandwidths into pdata block. */
  pdata->cv_mem = cv_mem;
  pdata->N = N;
  pdata->mu = mup = MIN( N-1, MAX(0,mu) );
  pdata->ml = mlp = MIN( N-1, MAX(0,ml) );

  /* Initialize nfeBP counter */
  pdata->nfeBP = 0;

  /* Allocate memory for saved banded Jacobian approximation. */
  pdata->savedJ = BandAllocMat(N, mup, mlp, mup);
  if (pdata->savedJ == NULL) {
    free(pdata);
    return(NULL);
  }

  /* Allocate memory for banded preconditioner. */
  storagemu = MIN( N-1, mup + mlp);
  pdata->savedP = BandAllocMat(N, mup, mlp, storagemu);
  if (pdata->savedP == NULL) {
    BandFreeMat(pdata->savedJ);
    free(pdata);
    return(NULL);
  }

  /* Allocate memory for pivot array. */
  pdata->pivots = BandAllocPiv(N);
  if (pdata->savedJ == NULL) {
    BandFreeMat(pdata->savedP);
    BandFreeMat(pdata->savedJ);
    free(pdata);
    return(NULL);
  }

  return((void *) pdata);
}

int CVBPSpgmr(void *cvode_mem, int pretype, int maxl, void *p_data)
{
  int flag;

  if ( p_data == NULL ) {
    fprintf(stdout, MSG_NO_PDATA);
    return(BP_NO_PDATA);
  } 

  flag = CVSpgmr(cvode_mem, pretype, maxl);
  if(flag != SUCCESS) return(flag);

  flag = CVSpgmrSetPrecData(cvode_mem, p_data);
  if(flag != SUCCESS) return(flag);

  flag = CVSpgmrSetPrecSetupFn(cvode_mem, CVBandPrecSetup);
  if(flag != SUCCESS) return(flag);

  flag = CVSpgmrSetPrecSolveFn(cvode_mem, CVBandPrecSolve);
  if(flag != SUCCESS) return(flag);

  return(SUCCESS);
}

void CVBandPrecFree(void *bp_data)
{
  CVBandPrecData pdata;

  if ( bp_data != NULL ) {
    pdata = (CVBandPrecData) bp_data;
    BandFreeMat(pdata->savedJ);
    BandFreeMat(pdata->savedP);
    BandFreePiv(pdata->pivots);
    free(pdata);
  }
}

int CVBandPrecGetIntWorkSpace(void *bp_data, long int *leniwBP)
{
  CVBandPrecData pdata;

  if ( bp_data == NULL ) {
    fprintf(stdout, MSG_PDATA_NULL);
    return(BP_NO_PDATA);
  } 

  pdata = (CVBandPrecData) bp_data;

  *leniwBP = pdata->N;

  return(OKAY);
}


int CVBandPrecGetRealWorkSpace(void *bp_data, long int *lenrwBP)
{
  CVBandPrecData pdata;
  long int N, ml, mu, smu;

  if ( bp_data == NULL ) {
    fprintf(stdout, MSG_PDATA_NULL);
    return(BP_NO_PDATA);
  } 

  pdata = (CVBandPrecData) bp_data;
  N   = pdata->N;
  mu  = pdata->mu;
  ml  = pdata->ml;
  smu = MIN( N-1, mu + ml);

  *lenrwBP = N * ( 2*ml + smu + mu + 2 );

  return(OKAY);
}

int CVBandPrecGetNumRhsEvals(void *bp_data, long int *nfevalsBP)
{
  CVBandPrecData pdata;

  if ( bp_data == NULL ) {
    fprintf(stdout, MSG_PDATA_NULL);
    return(BP_NO_PDATA);
  } 

  pdata = (CVBandPrecData) bp_data;

  *nfevalsBP = pdata->nfeBP;

  return(OKAY);
}



/***************** Preconditioner setup and solve functions *******/


/* Readability Replacements */

#define N         (pdata->N)
#define mu        (pdata->mu)
#define ml        (pdata->ml)
#define pivots    (pdata->pivots)
#define savedJ    (pdata->savedJ)
#define savedP    (pdata->savedP)
#define nfeBP     (pdata->nfeBP)

/****************** CVBandPrecSetup *******************************

  Together CVBandPrecSetup and CVBandPrecSolve use a banded           
  difference quotient Jacobian to create a preconditioner.       
  CVBandPrecSetup calculates a new J, if necessary, then           
  calculates P = I - gamma*J, and does an LU factorization of P. 
                                                                 
  The parameters of CVBandPrecSetup are as follows:                
                                                                 
  t       is the current value of the independent variable.      
                                                                 
  y       is the current value of the dependent variable vector, 
            namely the predicted value of y(t).                  
                                                                 
  fy      is the vector f(t,y).                                  
                                                                 
  jok     is an input flag indicating whether Jacobian-related   
          data needs to be recomputed, as follows:               
            jok == FALSE means recompute Jacobian-related data   
                   from scratch.                                 
            jok == TRUE  means that Jacobian data from the       
                   previous PrecSetup call will be reused          
                   (with the current value of gamma).            
          A CVBandPrecSetup call with jok == TRUE should only      
          occur after a call with jok == FALSE.                  
                                                                 
  *jcurPtr is a pointer to an output integer flag which is        
           set by CVBandPrecond as follows:                       
            *jcurPtr = TRUE if Jacobian data was recomputed.     
            *jcurPtr = FALSE if Jacobian data was not recomputed,
                       but saved data was reused.                
                                                                 
  gamma   is the scalar appearing in the Newton matrix.          
                                                                 
  bp_data is a pointer to preconditoner data - the same as the   
             bp_data parameter passed to CVSpgmr.                
                                                                 
  tmp1, tmp2, and tmp3 are pointers to memory allocated    
            for vectors of length N for work space.  This        
            routine uses only tmp1 and tmp2.                 
                                                                 
                                                                 
  The value to be returned by the CVBandPrecSetup function is      
    0  if successful, or                                         
    1  if the band factorization failed.                         

 ******************************************************************/

static int CVBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  long int ier;
  CVBandPrecData pdata;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CVBandPrecData) bp_data;

  if (jok) {
    /* If jok = TRUE, use saved copy of J. */
    *jcurPtr = FALSE;
    BandCopy(savedJ, savedP, mu, ml);
  } else {
    /* If jok = FALSE, call CVBandPDQJac for new J value. */
    *jcurPtr = TRUE;
    BandZero(savedJ);
    CVBandPDQJac(pdata, t, y, fy, tmp1, tmp2);
    BandCopy(savedJ, savedP, mu, ml);
  }
  
  /* Scale and add I to get savedP = I - gamma*J. */
  BandScale(-gamma, savedP);
  BandAddI(savedP);
 
  /* Do LU factorization of matrix. */
  ier = BandFactor(savedP, pivots);
 
  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}


/******************* CVBandPrecSolve ******************************

  CVBandPrecSolve solves a linear system P z = r, where P is the
  matrix computed by CVBandPrecond.

  The parameters of CVBandPrecSolve used here are as follows:

  r       is the right-hand side vector of the linear system.

  bp_data is a pointer to preconditioner data - the same as the
           bp_data parameter passed to CVSpgmr.

  z       is the output vector computed by CVBandPrecSolve.

  The value returned by the CVBandPrecSolve function is always 0,
  indicating success.
 
 ******************************************************************/

static int CVBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bp_data, N_Vector tmp)
{
  CVBandPrecData pdata;
  realtype *zd;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CVBandPrecData) bp_data;

  /* Copy r to z. */
  N_VScale(ONE, r, z);

  /* Do band backsolve on the vector z. */
  zd = N_VGetData(z);
  BandBacksolve(savedP, pivots, zd);
  N_VSetData(zd, z);

  return(0);
}

/*************** CVBandPDQJac ****************************************

 This routine generates a banded difference quotient approximation to
 the Jacobian of f(t,y).  It assumes that a band matrix of type
 BandMat is stored column-wise, and that elements within each column
 are contiguous. This makes it possible to get the address of a column
 of J via the macro BAND_COL and to write a simple for loop to set
 each of the elements of a column in succession.

**********************************************************************/

#define ewt    (cv_mem->cv_ewt)
#define uround (cv_mem->cv_uround)
#define h      (cv_mem->cv_h)
#define f      (cv_mem->cv_f)
#define f_data (cv_mem->cv_f_data)

static void CVBandPDQJac(CVBandPrecData pdata, 
                         realtype t, N_Vector y, N_Vector fy, 
                         N_Vector ftemp, N_Vector ytemp)
{
  CVodeMem cv_mem;
  realtype    fnorm, minInc, inc, inc_inv, srur;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;

  cv_mem = pdata->cv_mem;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp. */
  ewt_data   = N_VGetData(ewt);
  fy_data    = N_VGetData(fy);
  ftemp_data = N_VGetData(ftemp);
  y_data     = N_VGetData(y);
  ytemp_data = N_VGetData(ytemp);

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
    N_VSetData(ytemp_data, ytemp);
    f(t, ytemp, ftemp, f_data);
    nfeBP++;
    ftemp_data = N_VGetData(ftemp);

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
}
