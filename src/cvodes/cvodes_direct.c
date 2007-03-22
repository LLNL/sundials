/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2007-03-22 18:05:52 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CVSDIRECT linear solvers
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_impl.h"
#include "cvodes_direct_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PRIVATE FUNCTION PROTOTYPES
 * =================================================================
 */

static int cvDlsDenseJacBWrapper(int nB, realtype t,
                                 N_Vector yB, N_Vector fyB, 
                                 DlsMat JB, void *cvode_mem,
                                 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

static int cvDlsBandJacBWrapper(int nB, int mupperB, int mlowerB, 
                                realtype t, N_Vector yB, N_Vector fyB, 
                                DlsMat Jac, void *cvode_mem, 
                                N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);


/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define f              (cv_mem->cv_f)
#define f_data         (cv_mem->cv_f_data)
#define uround         (cv_mem->cv_uround)
#define nst            (cv_mem->cv_nst)
#define tn             (cv_mem->cv_tn)
#define h              (cv_mem->cv_h)
#define gamma          (cv_mem->cv_gamma)
#define gammap         (cv_mem->cv_gammap)
#define gamrat         (cv_mem->cv_gamrat)
#define ewt            (cv_mem->cv_ewt)

#define lmem           (cv_mem->cv_lmem)

#define mtype          (cvdls_mem->d_type)
#define n              (cvdls_mem->d_n)
#define ml             (cvdls_mem->d_ml)
#define mu             (cvdls_mem->d_mu)
#define smu            (cvdls_mem->d_smu)
#define djac           (cvdls_mem->d_djac)
#define bjac           (cvdls_mem->d_bjac)
#define M              (cvdls_mem->d_M)
#define savedJ         (cvdls_mem->d_savedJ)
#define pivots         (cvdls_mem->d_pivots)
#define nstlj          (cvdls_mem->d_nstlj)
#define nje            (cvdls_mem->d_nje)
#define nfeDQ          (cvdls_mem->d_nfeDQ)
#define J_data         (cvdls_mem->d_J_data)
#define last_flag      (cvdls_mem->d_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS (FORWARD INTEGRATION)
 * =================================================================
 */
              
/*
 * CVDlsSetJacFn specifies the (dense or band) Jacobian function.
 */
int CVDlsSetJacFn(void *cvode_mem, void *jac, void *jac_data)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_MEM_NULL, "CVSDIRECT", "CVDlsSetJacFn", MSGD_CVMEM_NULL);
    return(CVDIRECT_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVDIRECT_LMEM_NULL, "CVSDIRECT", "CVDlsSetJacFn", MSGD_LMEM_NULL);
    return(CVDIRECT_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) lmem;

  if (mtype == SUNDIALS_DENSE)
    djac = (CVDlsDenseJacFn) jac;
  else if (mtype == SUNDIALS_BAND)
    bjac = (CVDlsBandJacFn) jac;

  J_data = jac_data;

  return(CVDIRECT_SUCCESS);
}

/*
 * CVDlsGetWorkSpace returns the length of workspace allocated for the
 * CVDIRECT linear solver.
 */
int CVDlsGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_MEM_NULL, "CVSDIRECT", "CVDlsGetWorkSpace", MSGD_CVMEM_NULL);
    return(CVDIRECT_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVDIRECT_LMEM_NULL, "CVSDIRECT", "CVDlsGetWorkSpace", MSGD_LMEM_NULL);
    return(CVDIRECT_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) lmem;

  if (mtype == SUNDIALS_DENSE) {
    *lenrwLS = 2*n*n;
    *leniwLS = n;
  } else if (mtype == SUNDIALS_BAND) {
    *lenrwLS = n*(smu + mu + 2*ml + 2);
    *leniwLS = n;
  }

  return(CVDIRECT_SUCCESS);
}

/*
 * CVDlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_MEM_NULL, "CVSDIRECT", "CVDlsGetNumJacEvals", MSGD_CVMEM_NULL);
    return(CVDIRECT_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVDIRECT_LMEM_NULL, "CVSDIRECT", "CVDlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(CVDIRECT_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) lmem;

  *njevals = nje;

  return(CVDIRECT_SUCCESS);
}

/*
 * CVDlsGetNumRhsEvals returns the number of calls to the ODE function
 * needed for the DQ Jacobian approximation.
 */
int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_MEM_NULL, "CVSDIRECT", "CVDlsGetNumRhsEvals", MSGD_CVMEM_NULL);
    return(CVDIRECT_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVDIRECT_LMEM_NULL, "CVSDIRECT", "CVDlsGetNumRhsEvals", MSGD_LMEM_NULL);
    return(CVDIRECT_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) lmem;

  *nfevalsLS = nfeDQ;

  return(CVDIRECT_SUCCESS);
}

/*
 * CVDlsGetReturnFlagName returns the name associated with a CVDIRECT
 * return value.
 */
char *CVDlsGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CVDIRECT_SUCCESS:
    sprintf(name,"CVDIRECT_SUCCESS");
    break;   
  case CVDIRECT_MEM_NULL:
    sprintf(name,"CVDIRECT_MEM_NULL");
    break;
  case CVDIRECT_LMEM_NULL:
    sprintf(name,"CVDIRECT_LMEM_NULL");
    break;
  case CVDIRECT_ILL_INPUT:
    sprintf(name,"CVDIRECT_ILL_INPUT");
    break;
  case CVDIRECT_MEM_FAIL:
    sprintf(name,"CVDIRECT_MEM_FAIL");
    break;
  case CVDIRECT_JACFUNC_UNRECVR:
    sprintf(name,"CVDIRECT_JACFUNC_UNRECVR");
    break;
  case CVDIRECT_JACFUNC_RECVR:
    sprintf(name,"CVDIRECT_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * CVDlsGetLastFlag returns the last flag set in a CVDIRECT function.
 */
int CVDlsGetLastFlag(void *cvode_mem, int *flag)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_MEM_NULL, "CVSDIRECT", "CVDlsGetLastFlag", MSGD_CVMEM_NULL);
    return(CVDIRECT_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVDIRECT_LMEM_NULL, "CVSDIRECT", "CVDlsGetLastFlag", MSGD_LMEM_NULL);
    return(CVDIRECT_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) lmem;

  *flag = last_flag;

  return(CVDIRECT_SUCCESS);
}

/* 
 * =================================================================
 * DQ JACOBIAN APPROXIMATIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cvDlsDenseDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a dense matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */ 

int cvDlsDenseDQJac(int N, realtype t,
                    N_Vector y, N_Vector fy, 
                    DlsMat Jac, void *jac_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype fnorm, minInc, inc, inc_inv, yjsaved, srur;
  realtype *tmp2_data, *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  int j;
  int retval = 0;

  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* jac_data points to cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvdls_mem = (CVDlsMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */

    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = MAX(srur*ABS(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = f(t, y, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;
    
    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

    DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}

/*
 * -----------------------------------------------------------------
 * cvDlsBandDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y).  It assumes that a band matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

int cvDlsBandDQJac(int N, int mupper, int mlower,
                   realtype t, N_Vector y, N_Vector fy, 
                   DlsMat Jac, void *jac_data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  N_Vector ftemp, ytemp;
  realtype fnorm, minInc, inc, inc_inv, srur;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int group, i, j, width, ngroups, i1, i2;
  int retval = 0;

  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* jac_dat points to cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvdls_mem = (CVDlsMem) lmem;

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */

    retval = f(tn, ytemp, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(Jac,j);
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }
  
  return(retval);
}

/* 
 * =================================================================
 * BACKWARD INTEGRATION SUPPORT
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Additional readability replacements 
 * -----------------------------------------------------------------
 */

#define ytmp        (ca_mem->ca_ytmp)
#define getY        (ca_mem->ca_getY)

#define mtypeB      (cvdlsB_mem->d_typeB)
#define djacB       (cvdlsB_mem->d_djacB)
#define bjacB       (cvdlsB_mem->d_bjacB)
#define jac_data_B  (cvdlsB_mem->d_jac_dataB)

/*
 * -----------------------------------------------------------------
 * EXPORTED FUNCTIONS
 * -----------------------------------------------------------------
 */

int CVDlsSetJacFnB(void *cvode_mem, int which, void *jacB, void *jac_dataB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVDlsMemB cvdlsB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDIRECT_MEM_NULL, "CVSDIRECT", "CVDlsSetJacFnB", MSGD_CVMEM_NULL);
    return(CVDIRECT_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CVDIRECT_NO_ADJ, "CVSDIRECT", "CVDlsSetJacFnB", MSGD_NO_ADJ);
    return(CVDIRECT_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    CVProcessError(cv_mem, CVDIRECT_ILL_INPUT, "CVSDIRECT", "CVDlsSetJacFnB", MSGD_BAD_WHICH);
    return(CVDIRECT_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  if (cvB_mem->cv_lmem == NULL) {
    CVProcessError(cv_mem, CVDIRECT_LMEMB_NULL, "CVSDIRECT", "CVDlsSetJacFnB", MSGD_LMEMB_NULL);
    return(CVDIRECT_LMEMB_NULL);
  }
  cvdlsB_mem = (CVDlsMemB) (cvB_mem->cv_lmem);

  switch (mtypeB) {
  case SUNDIALS_DENSE:
    djacB = (CVDlsDenseJacFnB) jacB;
    flag = CVDlsSetJacFn(cvodeB_mem, (void *)cvDlsDenseJacBWrapper, cvode_mem);
    break;
  case SUNDIALS_BAND:
    bjacB = (CVDlsBandJacFnB) jacB;
    flag = CVDlsSetJacFn(cvodeB_mem, (void *)cvDlsBandJacBWrapper, cvode_mem);
    break;
  }

  jac_data_B = jac_dataB;

  return(flag);
}

/*
 * -----------------------------------------------------------------
 * PRIVATE INTERFACE FUNCTIONS
 * -----------------------------------------------------------------
 */

/*
 * cvDlsDenseJacBWrapper
 *
 * This routine interfaces to the CVDlsDenseJacFnB routine provided 
 * by the user. cvDlsDenseJacBWrapper is of type CVDlsDenseJacFn.
 * NOTE: jac_data actually contains cvode_mem
 */


static int cvDlsDenseJacBWrapper(int nB, realtype t,
                                 N_Vector yB, N_Vector fyB, 
                                 DlsMat JB, void *cvode_mem,
                                 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVDlsMemB cvdlsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvdlsB_mem = (CVDlsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  flag = getY(cv_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    CVProcessError(cv_mem, -1, "CVSDIRECT", "cvDlsDenseJacBWrapper", MSGD_BAD_TINTERP);
    return(-1);
  }

  /* Call user's adjoint dense djacB routine (of type CVDlsDenseJacFnB) */
  retval = djacB(nB, t, ytmp, yB, fyB, JB, jac_data_B, 
                 tmp1B, tmp2B, tmp3B);

  return(retval);

}


/*
 * cvDlsBandJacBWrapper
 *
 * This routine interfaces to the CVBandJacFnB routine provided 
 * by the user. cvDlsBandJacBWrapper is of type CVDlsBandJacFn.
 * NOTE: jac_data actually contains cvode_mem
 */

static int cvDlsBandJacBWrapper(int nB, int mupperB, int mlowerB, 
                                realtype t, N_Vector yB, N_Vector fyB, 
                                DlsMat JB, void *cvode_mem, 
                                N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVDlsMemB cvdlsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvdlsB_mem = (CVDlsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  flag = getY(cv_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    CVProcessError(cv_mem, -1, "CVSDIRECT", "cvDlsBandJacBWrapper", MSGD_BAD_TINTERP);
    return(-1);
  }

  /* Call user's adjoint band bjacB routine (of type CVDlsBandJacFnB) */
  retval = bjacB(nB, mupperB, mlowerB, t, ytmp, yB, fyB, JB, jac_data_B,
                 tmp1B, tmp2B, tmp3B);

  return(retval);
}

